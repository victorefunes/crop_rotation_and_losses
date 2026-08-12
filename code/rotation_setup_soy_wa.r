## ============================================================================
## rotation_setup_soy_wa.R
## Shared setup for corn and soy Just-Pope rotation analysis — SOYBEAN version.
## Source this file at the top of the soy analysis script:
##   source("rotation_setup_soy_wa.R")
##
## Rotation universe restricted to sequences ending in soybeans (5), using
## only codes 1 (corn), 5 (soybeans), 24 (winter wheat) — i.e. the soybean
## analog of the 49 corn-terminal sequences in rotation_setup_wa.r.
##
## TODO: soy_patterns below covers only ~80% of soy field-years (from the
## first 106 rows of `soy_df |> tab(rot_crop)`, filtered to 1/5/24-only
## codes). Extend with more rows from that same frequency table until
## cumulative coverage reaches ~99%, matching the corn file's standard.
## ============================================================================

library(tidyverse)
library(data.table)
library(statar)
library(fixest)
library(broom)
library(haven)
library(marginaleffects)
library(hdm)
library(dotwhisker)
library(ggfortify)
library(ggrepel)
library(patchwork)
library(knitr)
library(kableExtra)
library(sf)
library(usmap)
theme_set(theme_bw())


tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"
fig_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/figures/"

# ── Rotation patterns ─────────────────────────────────────────────────────────
# Sequences ending in soybeans (5) that account for ~80% of field-years so far
# (partial — see TODO above). Codes are CDL (1 corn, 5 soybeans, 24 winter
# wheat). Kept as numeric-code strings so rot_crop needs no recoding;
# readable labels are supplied by dict_soy below.

#soy_patterns <- tibble(pattern = c(
#  "1-5-1-5-1-5",   "5-1-1-5-1-5",   "1-1-1-5-1-5",   "5-1-5-1-5-5",
#  "5-5-1-5-1-5",   "1-1-5-1-1-5",   "5-1-5-1-1-5",   "5-1-5-5-1-5",
#  "1-5-1-1-1-5",   "5-5-5-5-5-5",   "1-5-5-5-1-5",   "1-1-1-1-1-5",
#  "1-5-1-5-5-5",   "5-1-1-1-1-5",   "1-5-5-1-5-5",   "5-1-5-5-5-5",
#  "5-5-5-5-1-5",   "5-5-5-1-5-5",   "5-5-1-5-5-5",   "1-5-5-5-5-5",
#  "5-24-1-5-1-5",  "1-1-5-5-1-5",   "1-1-5-1-5-5",   "5-1-5-24-1-5",
#  "1-5-5-1-1-5",   "1-5-1-1-5-5",   "1-1-1-1-5-5",   "5-5-1-1-1-5",
#  "24-1-5-24-1-5", "5-5-5-1-1-5",   "5-1-1-1-5-5",   "24-5-1-5-1-5",
#  "1-5-24-5-1-5",  "5-1-1-5-5-5",   "1-5-1-5-24-5",  "1-1-5-24-1-5",
#  "1-1-5-5-5-5",   "1-5-24-1-1-5",  "5-5-1-1-5-5",   "1-1-1-5-5-5",
#  "1-24-1-5-1-5"
#))

soy_patterns <- tibble(pattern = c(
  "1-1-1-1-1-5", "1-1-1-5-1-5", "1-1-5-1-1-5", "1-1-5-1-5-5",
  "1-1-5-24-1-5", "1-5-1-1-1-5", "1-5-1-1-5-5", "1-5-1-5-1-5",
  "1-5-1-5-24-5", "1-5-1-5-5-5", "1-5-24-1-5-5", "1-5-5-1-1-5",
  "1-5-5-1-5-5", "1-5-5-5-5-5", "24-1-1-1-1-5", "24-1-5-1-5-5",
  "36-36-36-5-1-5", "5-1-1-1-1-5", "5-1-1-5-1-5", "5-1-5-1-1-5",
  "5-1-5-1-5-5", "5-1-5-5-5-5", "5-5-1-5-5-5", "5-5-5-1-5-5"
))

# ── RCI correction ────────────────────────────────────────────────────────────
# Removes mismatched RCI / rot_crop combinations (data quality fix). Shared
# with the corn setup — if a wheat sequence shows an analogous RCI mismatch,
# add it here the same way.

rci_correction <- function(df) {
  df |>
    mutate(data_rm = case_when(
      rot_crop == "5-1-5-1-5-1" & RCI == 3.24 ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 3    ~ 1,
      rot_crop == "5-1-5-1-1-1" & RCI == 2.24 ~ 1,
      rot_crop == "5-1-1-5-1-5" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-5-1-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-1-5-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 2    ~ 1,
      rot_crop == "1-5-1-1-1-5" & RCI == 1.73 ~ 1,
      rot_crop == "1-1-1-5-1-5" & RCI == 0    ~ 1,
      rot_crop == "1-1-1-1-1-5" & RCI == 0    ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 2.45 ~ 1,
      rot_crop == "1-5-1-1-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-1-5-1-5" & RCI == 0    ~ 1,
      rot_crop == "1-5-1-5-1-5" & RCI == 2.45 ~ 1,
      .default = 0)) |>
    filter(data_rm == 0) |>
    select(-data_rm)
}

# ── Degree-day functions (Schlenker-Roberts 2009) ─────────────────────────────

compute_gdd <- function(tmin, tmax, base = 8, cap = 29) {
  pmax(0, pmin((tmax + tmin) / 2, cap) - base)
}
compute_edd <- function(tmin, tmax, threshold = 29) {
  pmax(0, (tmax + tmin) / 2 - threshold)
}

add_degree_days <- function(df) {
  df |>
    mutate(EDD_6 = compute_edd(tmmx_6, tmmn_6),
           EDD_7 = compute_edd(tmmx_7, tmmn_7),
           EDD_8 = compute_edd(tmmx_8, tmmn_8),
           GDD_6 = compute_gdd(tmmx_6, tmmn_6),
           GDD_7 = compute_gdd(tmmx_7, tmmn_7),
           GDD_8 = compute_gdd(tmmx_8, tmmn_8))
}

# ── Control vector ────────────────────────────────────────────────────────────
# nccpi3all_mean and soc0_100_mean excluded — collinear with field/FIPS FEs.

all_controls <- c(
  "pr_6", "pr_7", "pr_8",
  "I(pr_6^2)", "I(pr_7^2)", "I(pr_8^2)",
  "cGDD_6m", "cGDD_7m", "cGDD_8m",
  "EDD_6", "EDD_7", "EDD_8",
  "vpd_6", "vpd_7", "vpd_8",
  "soil_6", "soil_7", "soil_8",
  "rootznaws_mean"
)

# Controls for FGLS bootstrap — drop collinear soil vars
all_controls_fgls <- setdiff(all_controls, c("nccpi3all_mean", "soc0_100_mean"))

# Base column names that exist in the data — for NA filtering.
# Strips formula expressions like I(pr_6^2) which are not actual columns.
all_controls_cols <- all_controls_fgls[!grepl("^I\\(", all_controls_fgls)]

# ── Formula helper ────────────────────────────────────────────────────────────

make_jp_formula <- function(lhs, rot_var, controls,
                            fe = "tile_field_ID + year") {
  rhs <- paste(c(rot_var, controls), collapse = " + ")
  as.formula(paste(lhs, "~", rhs, "|", fe))
}

# ── Coefficient label dictionaries ────────────────────────────────────────────
# rot_crop levels stay as numeric-code strings; these map coefficient names to
# readable C/S/W labels. dict_soy is generated from the sequence list so it
# stays in sync automatically.

to_letters <- function(x) {                 # 1->C, 5->S, 24->W, 36->A
  x <- gsub("24", "W", x, fixed = TRUE)
  x <- gsub("36", "A", x, fixed = TRUE)
  x <- gsub("5",  "S", x, fixed = TRUE)
  x <- gsub("1",  "C", x, fixed = TRUE)
  x
}

make_dict <- function(patterns, ref) {
  seqs   <- setdiff(patterns, ref)
  labels <- to_letters(seqs)
  # rot_crop is factored on letter-coded levels (see to_letters() in the
  # analysis scripts), so dict keys must be "rot_crop" + letter sequence,
  # not the raw numeric-code pattern, or etable won't match the coefficients.
  setNames(labels, paste0("rot_crop", labels))
}

# Soy analysis: reference = continuous soybeans. All sequences here end in soy.
dict_soy <- make_dict(soy_patterns$pattern, ref = "5-5-5-5-5-5")

# Corn analysis: build the same way from the corn script's own 99% list
# (sequences ending in corn) — see rotation_setup_wa.r.
dict_corn <- NULL

dict_rci <- c(
  "RCI1.41" = "RCI = 1.41", "RCI1.73" = "RCI = 1.73",
  "RCI2"    = "RCI = 2",    "RCI2.24" = "RCI = 2.24",
  "RCI2.45" = "RCI = 2.45", "RCI2.65" = "RCI = 2.65",
  "RCI2.74" = "RCI = 2.74", "RCI2.83" = "RCI = 2.83",
  "RCI3"    = "RCI = 3",    "RCI3.24" = "RCI = 3.24",
  "RCI3.46" = "RCI = 3.46", "RCI3.67" = "RCI = 3.67",
  "RCI3.74" = "RCI = 3.74", "RCI4"    = "RCI = 4",
  "RCI4.24" = "RCI = 4.24", "RCI4.47" = "RCI = 4.47",
  "RCI4.74" = "RCI = 4.74", "RCI5.2"  = "RCI = 5.2"
)

dict_vpd <- c(
  "vpd_namesomewhatdry"           = "Somewhat dry season",
  "vpd_namedry"                   = "Dry season",
  "RCI x vpd_namesomewhatdry"     = "RCI x Somewhat dry season",
  "RCI x vpd_namedry"             = "RCI x Dry season"
)

# ── PCA of rotation features ──────────────────────────────────────────────────
# Features defined over the CDL codes so wheat and alfalfa contribute distinctly:
#   legume years    = soybeans (5) + WinWht/Soy dbl (26) + alfalfa (36)
#   small-grain yrs  = winter wheat (24) + WinWht/Soy dbl (26)
# (26 and 36 don't occur in the soy_patterns set, so those codes are inert
# here; they keep the code correct if the sequence list is ever widened.)
# Fit is on soy_patterns — i.e. the observed soy-terminal set — so rot_index
# reflects the rotations that actually occur. Sequences are weighted equally,
# not by frequency.

legume_codes      <- c(5L, 26L, 36L)
small_grain_codes <- c(24L, 26L)

parse_rot   <- function(x) as.integer(strsplit(x, "-", fixed = TRUE)[[1]])
shannon_div <- function(v) { p <- table(v) / length(v); -sum(p * log(p)) }
longest_run <- function(v) max(rle(v)$lengths)
count_in    <- function(v, codes) sum(v %in% codes)
gap_mean    <- function(v, codes) { p <- which(v %in% codes); if (length(p) < 2) 0 else mean(diff(p)) }
gap_min     <- function(v, codes) { p <- which(v %in% codes); if (length(p) < 2) 0 else min(diff(p)) }

build_rot_features <- function(patterns) {
  tibble(pattern = unique(patterns)) |>
    mutate(
      rot_vec    = lapply(pattern, parse_rot),
      n_legume   = sapply(rot_vec, count_in, codes = legume_codes),
      n_grain    = sapply(rot_vec, count_in, codes = small_grain_codes),
      diversity  = sapply(rot_vec, shannon_div),
      max_run    = sapply(rot_vec, longest_run),
      no_mono    = 6 - max_run,
      leg_gap    = sapply(rot_vec, gap_mean, codes = legume_codes),
      leg_mingap = sapply(rot_vec, gap_min,  codes = legume_codes),
      free_leg   = ifelse(leg_gap    == 0, 0, 1 / leg_gap),
      tight_leg  = ifelse(leg_mingap == 0, 0, 1 / leg_mingap)
    ) |>
    select(pattern, n_legume, n_grain, diversity, no_mono, free_leg, tight_leg)
}

rot_pca_features <- c("n_legume", "n_grain", "diversity",
                      "no_mono", "free_leg", "tight_leg")

build_rot_pca <- function(patterns) {
  feats <- build_rot_features(patterns)
  X <- as.data.frame(feats[, rot_pca_features])
  # Drop any zero-variance column so scale. = TRUE can't divide by zero.
  keep <- vapply(X, function(z) stats::sd(z) > 0, logical(1))
  pca  <- prcomp(X[, keep, drop = FALSE], scale. = TRUE)
  feats$rot_index <- pca$x[, 1]
  attr(feats, "pca")      <- pca
  attr(feats, "pca_cols") <- names(keep)[keep]
  feats
}

rot_features <- build_rot_pca(soy_patterns$pattern)
pca_rot      <- attr(rot_features, "pca")

cat("Rotation PCA — variance explained by first two PCs:\n")
print(summary(pca_rot)$importance[, 1:2])

# CHECK THE PC1 SIGN against pca_rot$rotation before using rot_index as a
# regressor — prcomp's sign is arbitrary. Flip with rot_index <- -rot_index if
# "more complex" should be the high end.

# Figure: rot_pca_plot — PCA biplot (rotation_plots.R has the fuller version).
autoplot(pca_rot, data = rot_features,
         loadings = TRUE, loadings.label = TRUE, loadings.label.repel = TRUE) +
  labs(title   = "PCA of rotation features (corn, soy, wheat)",
       caption = "Each point is one of the soy-terminal sequences in soy_patterns.") +
  theme_bw() +
  theme(legend.position = "none") ->
  rot_pca_plot
ggsave(paste0(fig_dir, "rot_pca_plot_soy.png"), rot_pca_plot,
       width = 9, height = 7, dpi = 300)

# ── Bootstrap function (shared) ───────────────────────────────────────────────

boot_jp_fgls <- function(dt, fml_mean, fml_var, fml_var_b,
                          B = 499, seed = 42, n_workers = NULL) {
  dt <- as.data.table(dt)

  # Remove singletons iteratively until stable
  repeat {
    n_before <- nrow(dt)
    dt <- dt[dt[, .N, by = tile_field_ID][N > 1], on = "tile_field_ID"]
    dt <- dt[dt[, .N, by = year][N > 1],          on = "year"]
    if (nrow(dt) == n_before) break
  }
  cat("Rows after singleton removal:", nrow(dt), "\n")

  # Keep only the columns the formulas actually reference. dt gets broadcast
  # to every bootstrap worker process below, so trimming unused columns here
  # (e.g. RCI, rotation lag dummies) bounds peak memory across workers.
  required_cols <- unique(c(all.vars(fml_mean), all.vars(fml_var),
                             all.vars(fml_var_b), "tile_field_ID", "year"))
  dt <- dt[, .SD, .SDcols = intersect(required_cols, names(dt))]
  cat("Columns kept for bootstrap:", ncol(dt), "\n")

  # nthreads: threads used *within* one feols call (shared memory, cheap).
  # n_workers: separate R processes for the parallel bootstrap below — each
  # holds its own full copy of dt, so this is the memory-multiplying knob.
  # Default to half the cores rather than detectCores()-1 to leave headroom;
  # lower it further (or slim dt more) if you still see out-of-memory errors.
  n_threads <- parallel::detectCores() - 1L
  if (is.null(n_workers)) n_workers <- max(1L, parallel::detectCores() %/% 2L)
  set.seed(seed)
  fields    <- unique(dt$tile_field_ID)
  n_fields  <- length(fields)
  field_idx <- match(dt$tile_field_ID, fields)

  # Point estimates on full data
  s1 <- feols(fml_mean, data = dt,
              cluster  = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, resid_sq := NA_real_]
  dt[obs(s1), resid_sq := residuals(s1)^2]

  s2a <- feols(fml_var, data = dt,
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, h_hat := NA_real_]
  dt[obs(s2a), h_hat := pmax(fitted(s2a), 1e-6)]

  s2b <- feols(fml_var, data = dt,
               weights  = ~I(1 / h_hat),
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)

  coef_hat <- coef(s2b)

  # ── Bootstrap, parallelized across iterations (one core per draw) ──────────
  boot_one <- function(b, dt, n_fields, field_idx, fml_mean, fml_var_b) {
    draws  <- tabulate(sample.int(n_fields, n_fields, replace = TRUE), n_fields)
    boot_w <- draws[field_idx]
    keep   <- boot_w > 0

    tryCatch({
      dt_b <- dt[keep]
      dt_b[, boot_w := boot_w[keep]]

      b1 <- feols(fml_mean, data = dt_b, weights = ~boot_w,
                  vcov = "iid", nthreads = 1, warn = FALSE, notes = FALSE)
      dt_b[, resid_sq_b := NA_real_]
      dt_b[obs(b1), resid_sq_b := residuals(b1)^2]

      b2a <- feols(fml_var_b, data = dt_b, weights = ~boot_w,
                   vcov = "iid", nthreads = 1, warn = FALSE, notes = FALSE)
      dt_b[, h_hat_b    := NA_real_]
      dt_b[obs(b2a), h_hat_b := pmax(fitted(b2a), 1e-6)]
      dt_b[, combined_w := boot_w / h_hat_b]

      b2b <- feols(fml_var_b, data = dt_b, weights = ~combined_w,
                   vcov = "iid", nthreads = 1, warn = FALSE, notes = FALSE)

      coef(b2b)
    }, error = function(e) {
      cat("Bootstrap iteration", b, "failed:", conditionMessage(e), "\n")
      NULL
    })
  }

  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterEvalQ(cl, { library(data.table); library(fixest) })
  parallel::clusterSetRNGStream(cl, seed)

  cat("Running", B, "bootstrap iterations on", n_workers, "workers...\n")
  boot_results <- parallel::parLapply(
    cl, seq_len(B), boot_one,
    dt = dt, n_fields = n_fields, field_idx = field_idx,
    fml_mean = fml_mean, fml_var_b = fml_var_b
  )
  cat("Bootstrap complete.\n")

  coef_boot <- matrix(NA_real_, nrow = B, ncol = length(coef_hat),
                       dimnames = list(NULL, names(coef_hat)))
  for (b in seq_len(B)) {
    res <- boot_results[[b]]
    if (!is.null(res)) coef_boot[b, names(res)] <- res
  }

  boot_se <- apply(coef_boot, 2, sd,       na.rm = TRUE)
  boot_ci <- apply(coef_boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  list(coef = coef_hat, se = boot_se, ci = boot_ci,
       fit_mean = s1, fit_var_ols = s2a, fit_var_fgls = s2b,
       boot_draws = coef_boot)
}
