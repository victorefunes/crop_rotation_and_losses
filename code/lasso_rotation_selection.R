## ============================================================================
## lasso_rotation_selection.R
## Replaces the PCA feature-space analysis with a LASSO-based selection of
## which rotation sequences actually predict corn / soybean yields.
##
## Design, matching the request exactly:
##   1. Partial out field (tile_field_ID) and year fixed effects from the
##      outcome AND from every sequence dummy, via the Frisch-Waugh-Lovell
##      (FWL) demeaning step -- NOT by putting field dummies in the LASSO
##      design matrix. tile_field_ID has enormous cardinality; including it
##      as regressors would be both computationally infeasible for glmnet/
##      rlasso and statistically wrong (LASSO would shrink the FE dummies
##      too, reintroducing the omitted-fixed-effect bias the FE model exists
##      to remove). fixest::demean() runs the same alternating-projections
##      algorithm feols() uses internally, so this is numerically consistent
##      with every other regression in the paper.
##   2. LASSO of the demeaned outcome on the demeaned sequence-dummy matrix.
##      Uses hdm::rlasso (already a dependency in rotation_setup.R), which
##      sets a theoretically-justified, heteroskedasticity-robust penalty
##      (Belloni-Chernozhukov-Hansen) rather than a CV-tuned one -- more
##      defensible for selection consistency than cv.glmnet here, though a
##      cv.glmnet cross-check is included below.
##   3. Keep only the nonzero-coefficient sequences.
##   4. Refit via feols on just the selected sequences, with the paper's
##      normal FE/clustering setup, to get valid, uncorrected-for-shrinkage
##      standard errors for reporting ("post-selection OLS"). LASSO
##      coefficients themselves are shrunk and are not the right thing to
##      put in a results table.
##
## Requires: fixest, hdm, glmnet (optional cross-check), tidyverse, data.table.
## Assumes corn_jp_data / soy_jp_data and all_controls already exist, as built
## in rotation_setup.R / main_analysis.R.
## ============================================================================

library(fixest)
library(data.table)
library(hdm)
library(glmnet)   # optional cross-check only
library(arrow)
library(statar)
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

# ── Core function ──────────────────────────────────────────────────────────
# data        : corn_jp_data or soy_jp_data
# yvar        : "corn_yield" or "soy_yield"
# seqvar      : the rotation sequence factor (default "rot_crop"; reference
#               level should already be set via relevel(), as in the rest of
#               the pipeline -- monoculture for corn, soy monoculture for soy)
# fe_vars     : fixed effects to partial out before LASSO
# cluster_fml : clustering formula for the post-selection refit
# controls    : optional additional regressors to partial out alongside the
#               FEs before running LASSO (see note below -- NULL by default,
#               matching "partialling out year and field FE" literally)
#
# Returns a list with the rlasso fit, the selected sequence names, and the
# post-selection feols refit (with and without the full control vector).

lasso_select_sequences <- function(data, yvar, seqvar = "rot_crop",
                                   fe_vars = c("tile_field_ID", "year"),
                                   cluster_fml = ~COUNTY_FIPS,
                                   controls = NULL) {
 
  data <- as.data.table(data)
  stopifnot(is.factor(data[[seqvar]]))
  ref_level <- levels(data[[seqvar]])[1]
  cat(sprintf("Reference sequence (dropped): %s\n", ref_level))
 
  # --- 0. Build ONE complete-case frame up front, and derive y / X_seq / W
  # all from it. Do NOT extract y, X_seq, and W as three separate operations
  # on `data` -- model.matrix() silently drops NA rows via na.omit while plain
  # column selection does not, so three separate extractions can silently
  # disagree on row count (this is what produced the cbind error). Making the
  # NA-drop explicit here, once, before anything else is built, removes that
  # entire bug class rather than papering over one instance of it.
  needed_cols <- unique(c(yvar, seqvar, fe_vars, controls))
  analysis_dt <- data[, ..needed_cols]
 
  n_before <- nrow(analysis_dt)
  complete <- stats::complete.cases(analysis_dt)
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    cat(sprintf("Dropping %d of %d rows (%.2f%%) with missing values in %s or controls.\n",
                n_dropped, n_before, 100 * n_dropped / n_before, yvar))
  }
  analysis_dt <- analysis_dt[complete]
  analysis_dt[[seqvar]] <- droplevels(analysis_dt[[seqvar]])
  if (!(ref_level %in% levels(analysis_dt[[seqvar]]))) {
    stop("Reference level '", ref_level, "' has no complete-case observations left; ",
         "choose a different reference or investigate why it's entirely missing controls.")
  }
 
  # --- 1. Build the sequence dummy matrix (one column per non-reference level) ---
  X_seq <- model.matrix(as.formula(paste0("~", seqvar, " - 1")), data = analysis_dt)
  ref_col <- paste0(seqvar, ref_level)
  X_seq <- X_seq[, colnames(X_seq) != ref_col, drop = FALSE]
  cat(sprintf("Candidate sequences for selection: %d\n", ncol(X_seq)))
 
  # --- 2. Partial out field + year FE (and, optionally, controls) via FWL ---
  # fixest::demean() runs the same alternating-projections sweep feols() uses,
  # so this is numerically consistent with the rest of the paper's models.
  # NOTE on `controls`: as specified, only field + year FE are partialled out
  # by default. If you want the LASSO to select sequences net of weather/soil
  # effects too (recommended -- otherwise a selected sequence could partly be
  # standing in for the weather/soil conditions correlated with it), pass the
  # numeric columns from all_controls_cols here; they get demeaned alongside
  # the FEs and dropped from the design matrix before selection (they are not
  # candidates for selection themselves, only controls being partialled out).
  fe_data <- analysis_dt[, ..fe_vars]
  y <- analysis_dt[[yvar]]

  if (!is.null(controls)) {
    W <- as.matrix(analysis_dt[, ..controls])
  }

  # Guaranteed to hold by construction (all three come from analysis_dt), but
  # asserted explicitly so any future refactor that breaks this fails loudly
  # here instead of at an opaque cbind() call three lines down.
  stopifnot(nrow(X_seq) == length(y))
  if (!is.null(controls)) stopifnot(nrow(W) == length(y))

  # X_seq can easily be 300+ dummy columns; with millions of rows a single
  # cbind(y, X_seq, W) is a many-GB matrix that can fail to allocate even
  # when plenty of (fragmented) system RAM is free ("cannot allocate vector
  # of size N Gb"). Demean y + W together (small), then demean X_seq in
  # column batches, writing each batch's result directly into a
  # preallocated X_dm -- this keeps peak transient memory to one batch's
  # worth instead of a second full-size copy of the whole matrix.
  yW <- if (!is.null(controls)) cbind(y = y, W) else cbind(y = y)
  yW_dm <- fixest::demean(X = yW, f = fe_data)
  y_dm <- yW_dm[, "y"]

  n_candidate <- ncol(X_seq)
  X_dm <- matrix(NA_real_, nrow = nrow(X_seq), ncol = n_candidate,
                 dimnames = list(NULL, colnames(X_seq)))
  batch_size <- 40L
  n_batches <- ceiling(n_candidate / batch_size)
  for (b in seq_len(n_batches)) {
    cols <- ((b - 1L) * batch_size + 1L):min(b * batch_size, n_candidate)
    X_dm[, cols] <- fixest::demean(X = X_seq[, cols, drop = FALSE], f = fe_data)
  }
  rm(X_seq); gc()

  if (!is.null(controls)) {
    W_dm <- yW_dm[, colnames(W), drop = FALSE]
    # Partial the (demeaned) controls out of y and of every sequence column
    # too, via a second FWL step, so LASSO sees pure sequence variation.
    y_dm <- residuals(lm.fit(W_dm, y_dm))
    X_dm <- apply(X_dm, 2, function(col) residuals(lm.fit(W_dm, col)))
  }
 
  # --- 3. LASSO selection ---
  # hdm::rlasso: plug-in penalty robust to heteroskedasticity, no CV needed.
  # post = TRUE returns post-lasso (OLS refit on selected set internally) --
  # we still do our own post-selection feols refit below for clustered SEs
  # and to keep the reported coefficients on the paper's usual scale/table
  # format, but rlasso's post-lasso coefficients are a useful cross-check.
  fit <- hdm::rlasso(x = X_dm, y = y_dm, post = TRUE)
 
  nonzero <- names(coef(fit))[coef(fit) != 0]
  nonzero <- setdiff(nonzero, "(Intercept)")
  selected_sequences <- sub(paste0("^", seqvar), "", nonzero)
  cat(sprintf("Selected (nonzero) sequences: %d of %d\n",
              length(selected_sequences), n_candidate))
 
  # --- 3b. Optional cross-check: cv.glmnet with clustered folds ---
  # LASSO's CV-tuned lambda tends to select more variables than rlasso's
  # plug-in penalty; useful as a sensitivity check on how selection-count-
  # dependent the results are, not as the primary specification.
  cv_fit <- tryCatch({
    cv.glmnet(x = X_dm, y = y_dm, alpha = 1, standardize = TRUE)
  }, error = function(e) { message("cv.glmnet cross-check failed: ", conditionMessage(e)); NULL })
  cv_selected <- if (!is.null(cv_fit)) {
    cf <- coef(cv_fit, s = "lambda.1se")
    nm <- rownames(cf)[which(cf[, 1] != 0)]
    setdiff(sub(paste0("^", seqvar), "", nm), c("(Intercept)", ""))
  } else NULL
 
  # --- 4. Post-selection refit for valid inference ---
  # Refit on `data` (the original, full-column frame), not analysis_dt, so the
  # reported model has access to every column etable()/downstream code might
  # expect -- but apply the SAME complete-case row filter first, so the refit
  # sample matches exactly the sample the selection was estimated on. Re-level
  # rot_crop to only the selected sequences (same reference), so the reported
  # coefficients are on the paper's native scale, with the usual clustered
  # SEs -- not the shrunk LASSO estimates.
  data_sel <- data[complete][data[complete][[seqvar]] %in% c(ref_level, selected_sequences)]
  data_sel[[seqvar]] <- factor(data_sel[[seqvar]], levels = c(ref_level, selected_sequences))
 
  fml_nc <- as.formula(paste(yvar, "~", seqvar, "| tile_field_ID + year"))
  refit_nc <- feols(fml_nc, data = data_sel, cluster = cluster_fml)
 
  refit_full <- NULL
  if (!is.null(controls)) {
    fml_full <- make_jp_formula(yvar, seqvar, all_controls)  # from rotation_setup.R
    refit_full <- feols(fml_full, data = data_sel, cluster = cluster_fml)
  }
 
  list(
    rlasso_fit          = fit,
    selected_sequences  = selected_sequences,
    n_candidate         = n_candidate,
    n_selected          = length(selected_sequences),
    cv_glmnet_selected  = cv_selected,   # cross-check; compare length/overlap
    cv_glmnet_fit       = cv_fit,        # full path; plot(cv_glmnet_fit) for lambda vs CV error / nvars
    refit_no_controls   = refit_nc,
    refit_full_controls = refit_full,
    reference_level     = ref_level,
    X_dm                = X_dm,          # demeaned (and control-partialled) sequence dummies fed to LASSO
    y_dm                = y_dm           # demeaned (and control-partialled) outcome fed to LASSO
  )
}

source("rotation_setup_wa.R")

# ── Load corn data ────────────────────────────────────────────────────────────
# Load, correct, and add degree-day variables in one pipeline.
# corn_df is kept as the raw data source; corn_jp_data is the analysis sample.

cat("Loading corn data...\n")
corn_df <- read_parquet(
  "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet")

corn_df <- corn_df |>
  filter(STATE_ABBR == "IL") |>
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         corn_yield = corn_yield / 62.77)  |>
  arrange(tile_field_ID, year)

cat("Corn raw rows:", nrow(corn_df), "\n")

## Add present year crop variable
setDT(corn_df)

add_crop_year <- function(dt) {
  stopifnot(is.data.table(dt))

  # Drop any stale crop_year column from a previous run, so the regex below
  # can't accidentally re-capture it
  if ("crop_year" %in% names(dt)) {
    dt[, crop_year := NULL]
  }

  # Match crop_YYYY columns only (exactly 4 digits) -- excludes crop_year,
  # rot_seq, or anything else that happens to start with "crop_". This is a
  # pure year-lookup: it works the same whether crop_YYYY holds raw CDL names
  # ("Corn"), numeric CDL codes, or functional-group labels ("corn") from
  # recode_cdl_functional() -- it never inspects the values, only picks the
  # column matching dt$year. Forced to character so the lookup is well-defined
  # regardless of which of those three a given caller has already applied.
  crop_cols  <- grep("^crop_[0-9]{4}$", names(dt), value = TRUE)
  crop_years <- as.integer(sub("crop_", "", crop_cols))

  stopifnot(length(crop_cols) > 0, !anyNA(crop_years))

  m       <- as.matrix(dt[, lapply(.SD, as.character), .SDcols = crop_cols])
  col_idx <- match(dt$year, crop_years)

  stopifnot(length(col_idx) == nrow(dt))

  dt[, crop_year := m[cbind(seq_len(.N), col_idx)]]
  invisible(dt)
}

add_crop_year(corn_df)

corn_df <- corn_df |>
  rename(crop_0 = crop_year,
         crop_1 = prioryr_crop,
         crop_2 = prior2yr_crop,
         crop_3 = prior3yr_crop,
         crop_4 = prior4yr_crop,
         crop_5 = prior5yr_crop,
         crop_6 = prior6yr_crop)

source("rci_vectorized.R")

corn_df <- corn_df |>
    mutate(RCI = rci(crop_0, crop_1, crop_2, crop_3, crop_4, crop_5))

source("cdl_functional_recode.R")

recode_cdl_functional(corn_df, cols = paste0("crop_", 0:6))

corn_df <- corn_df |>
  #rename(RCI = annual_RCI) |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-",
          crop_2, "-", crop_1, "-", crop_0)) |>
  #rci_correction() |>
  add_degree_days()

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.
#
# Rotation universe: any COMPLETE 6-year functional-group history (crop_0..
# crop_5 all agricultural, i.e. non-NA after recode_cdl_functional()), not
# just the 29 hand-curated corn/soy/wheat sequences corn_soy_patterns still
# encodes -- lasso_rotation_selection.R now does the dozens-to-a-handful
# selection on this broader candidate set instead of a pre-filtered pattern
# list. MIN_PATTERN_FREQ drops sequences so rare they'd be perfectly
# collinear with field FE and add nothing but noise to the LASSO design
# matrix; raise/lower it there, not by hand-editing a pattern list.
# (corn_soy_patterns / rot_features are still used below, unchanged, to
# attach rot_index -- that PCA is fit only on the original 29 sequences, so
# rows outside that set simply get rot_index = NA, which is fine: only the
# fml_corn_index model actually uses rot_index, and it drops NA rows itself.)
corn_df |> 
  tab(rot_crop) |>
  data.frame() |>
  arrange(desc(Freq.)) |>
  head(20)

MIN_PATTERN_FREQ <- 250
pattern_freq     <- corn_df |> count(rot_crop, name = "N")
common_patterns  <- pattern_freq$rot_crop[pattern_freq$N >= MIN_PATTERN_FREQ]
cat(sprintf("Rotation patterns: %d distinct, %d kept at N >= %d field-years.\n",
            nrow(pattern_freq), length(common_patterns), MIN_PATTERN_FREQ))

corn_jp_data <- corn_df |>
  filter(if_all(paste0("crop_", 0:5), ~ !is.na(.))) |>
  filter(rot_crop %in% common_patterns) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")),
    # lasso_select_sequences() requires rot_crop to be a factor with the
    # reference level (continuous corn monoculture) already set via relevel()
    # -- see stopifnot(is.factor(...)) at the top of that function.
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "corn-corn-corn-corn-corn-corn"))

cat("Corn analysis sample:", nrow(corn_jp_data), "rows\n")

# Free raw data — no longer needed
#rm(corn_df);
gc()

# ── Apply to corn and soybean ────────────────────────────────────────────────
# Pass all_controls_cols (numeric control column names, no I() terms) to also
# partial out weather/soil before selection -- recommended; see note above.
# Using cluster = ~COUNTY_FIPS to match the existing corn_rot / soy_rot tables.
corn_lasso <- lasso_select_sequences(
  data        = corn_jp_data,
  yvar        = "corn_yield",
  controls    = all_controls_cols,
  cluster_fml = ~COUNTY_FIPS
)

#soy_lasso <- lasso_select_sequences(
#  data        = soy_jp_data,
#  yvar        = "soy_yield",
#  controls    = all_controls_cols,
#  cluster_fml = ~COUNTY_FIPS
#)

cat(sprintf("\nCorn: %d of %d sequences selected (rlasso); %d selected under cv.glmnet cross-check.\n",
            corn_lasso$n_selected, corn_lasso$n_candidate, length(corn_lasso$cv_glmnet_selected)))
#cat(sprintf("Soy:  %d of %d sequences selected (rlasso); %d selected under cv.glmnet cross-check.\n",
#            soy_lasso$n_selected, soy_lasso$n_candidate, length(soy_lasso$cv_glmnet_selected)))

# ── Save tables in the same etable/dict_corn style as the rest of the paper ──
# `keep` is built from functional_letters (rotation_setup_wa.R) rather than a
# hardcoded "^[CSW]-", because rot_crop's candidate universe is no longer
# corn/soy/wheat-only: recode_cdl_functional() over the full crop history
# means a selected sequence can contain any of the project's functional-group
# letters (A = ley, L = annual_legume, B = annual_broadleaf, F = fallow, on
# top of C/S/W) -- a literal "^[CSW]-" would silently drop those rows from
# the table.
seq_letter_class <- paste0("^[", paste(unname(functional_letters), collapse = ""), "]-")

etable(corn_lasso$refit_full_controls, tex = TRUE, cluster = ~COUNTY_FIPS,
       dict = dict_corn, keep = seq_letter_class,
       file = paste0(tab_dir, "corn_lasso.tex"), replace = TRUE,
       title = "LASSO-selected rotation sequence effects on corn yield")

#etable(soy_lasso$refit_full_controls, tex = TRUE, cluster = ~COUNTY_FIPS,
#       dict = dict_soy, keep = seq_letter_class,
#       file = paste0(tab_dir, "soy_lasso.tex"), replace = TRUE,
#       title = "LASSO-selected rotation sequence effects on soybean yield")
length(corn_lasso$selected_sequences)   # n_selected
corn_lasso$n_candidate                  # out of how many candidates
setdiff(corn_lasso$selected_sequences, corn_lasso$cv_glmnet_selected)  # rlasso-only picks
setdiff(corn_lasso$cv_glmnet_selected, corn_lasso$selected_sequences)  # cv.glmnet-only picks

fwrite(data.table(rot_crop = corn_lasso$selected_sequences), "corn_lasso_selected_sequences.csv")
