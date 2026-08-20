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
 
  to_demean <- if (!is.null(controls)) cbind(y = y, X_seq, W) else cbind(y = y, X_seq)
  demeaned <- fixest::demean(X = to_demean, f = fe_data)
 
  y_dm <- demeaned[, "y"]
  X_dm <- demeaned[, colnames(X_seq), drop = FALSE]
  if (!is.null(controls)) {
    W_dm <- demeaned[, colnames(W), drop = FALSE]
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
              length(selected_sequences), ncol(X_seq)))
 
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
    n_candidate         = ncol(X_seq),
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
 
# ── Apply to corn and soybean ────────────────────────────────────────────────
# Pass all_controls_cols (numeric control column names, no I() terms) to also
# partial out weather/soil before selection -- recommended; see note above.
# Using cluster = ~COUNTY_FIPS to match the existing corn_rot / soy_rot tables.
#corn_lasso <- lasso_select_sequences(
#  data        = corn_jp_data,
#  yvar        = "corn_yield",
#  controls    = all_controls_cols,
#  cluster_fml = ~COUNTY_FIPS
#)

#soy_lasso <- lasso_select_sequences(
#  data        = soy_jp_data,
#  yvar        = "soy_yield",
#  controls    = all_controls_cols,
#  cluster_fml = ~COUNTY_FIPS
#)

#cat(sprintf("\nCorn: %d of %d sequences selected (rlasso); %d selected under cv.glmnet cross-check.\n",
#            corn_lasso$n_selected, corn_lasso$n_candidate, length(corn_lasso$cv_glmnet_selected)))
#cat(sprintf("Soy:  %d of %d sequences selected (rlasso); %d selected under cv.glmnet cross-check.\n",
#            soy_lasso$n_selected, soy_lasso$n_candidate, length(soy_lasso$cv_glmnet_selected)))

# ── Save tables in the same etable/dict_corn style as the rest of the paper ──
#etable(corn_lasso$refit_full_controls, tex = TRUE, cluster = ~COUNTY_FIPS,
#       dict = dict_corn, keep = "^[CSW]-",
#       file = paste0(tab_dir, "corn_lasso.tex"), replace = TRUE,
#       title = "LASSO-selected rotation sequence effects on corn yield")

#etable(soy_lasso$refit_full_controls, tex = TRUE, cluster = ~COUNTY_FIPS,
#       dict = dict_soy, keep = "^[CSW]-",
#       file = paste0(tab_dir, "soy_lasso.tex"), replace = TRUE,
#       title = "LASSO-selected rotation sequence effects on soybean yield")
