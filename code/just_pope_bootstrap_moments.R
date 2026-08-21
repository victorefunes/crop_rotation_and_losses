# ==============================================================================
# Corrected moment bootstrap for the Just-Pope decomposition
# ------------------------------------------------------------------------------
# Drop-in replacement for boot_jp_fgls() in just_pope.r. Extends the existing
# field-level (cluster) pairs bootstrap from two moments (mean + variance) to
# three (mean + variance + skewness).
#
# Why a bootstrap is required for the third moment
# ------------------------------------------------------------------------------
# The stage-2 (resid^2) and stage-3 (resid^3) regressors are GENERATED: they are
# built from the stage-1 fitted mean. The analytic clustered SEs that feols
# reports for a standalone feols(resid_cube ~ ...) treat resid_cube as data and
# therefore understate uncertainty, because they ignore the sampling variability
# of the stage-1 fit that produced the residuals. The fix is to re-estimate
# stage 1 inside every bootstrap replicate and rebuild BOTH the squared and the
# cubed residual from that same replicate-specific fit, so the first-stage
# estimation error propagates into the stage-2 and stage-3 coefficient
# distributions. This is the standard Murphy-Topel / generated-regressor concern;
# resampling the whole pipeline sidesteps the need for an analytic correction.
#
# Estimator choices (kept consistent with the existing code)
# ------------------------------------------------------------------------------
#  * Resampling unit: field (tile_field_ID), via frequency weights (no row
#    duplication), exactly as in boot_jp_fgls().
#  * Stage 2 (variance): FGLS, weighting by 1/h_hat as in the original.
#  * Stage 3 (skewness): OLS moment regression of resid^3 on the same RHS
#    (Antle 1983). There is no natural FGLS weight for the third central moment,
#    so it is estimated by OLS within each replicate and its inference comes
#    entirely from the bootstrap. Do NOT report the analytic clustered SE for the
#    skewness stage -- it is the quantity this function exists to replace.
# ==============================================================================

library(fixest)
library(data.table)

# Rebuild formulas with cleaned controls (as in just_pope.r)
fml_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)
fml_var  <- make_jp_formula("resid_sq",   "rot_crop", all_controls)
fml_skew <- make_jp_formula("resid_cube", "rot_crop", all_controls)

boot_jp_moments <- function(dt, fml_mean, fml_var, fml_skew,
                            B = 999, seed = 42) {

  dt <- as.data.table(dt)

  # ── Remove singletons (drop fields/years with a single observation) ─────────
  repeat {
    n_before <- nrow(dt)
    dt <- dt[dt[, .N, by = tile_field_ID][N > 1], on = "tile_field_ID"]
    dt <- dt[dt[, .N, by = year][N > 1],          on = "year"]
    if (nrow(dt) == n_before) break
  }
  cat("Rows after singleton removal:", nrow(dt), "\n")

  n_threads <- parallel::detectCores() - 1L
  set.seed(seed)
  fields   <- unique(dt$tile_field_ID)
  n_fields <- length(fields)

  # ── Parallel formulas with bootstrap LHS names ──────────────────────────────
  fml_var_b  <- make_jp_formula("resid_sq_b",   "rot_crop", all_controls)
  fml_skew_b <- make_jp_formula("resid_cube_b", "rot_crop", all_controls)

  # ── Point estimates on the full data ────────────────────────────────────────
  # Stage 1: conditional mean
  s1 <- feols(fml_mean, data = dt,
              cluster  = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)

  # Rebuild both residual powers from the SAME stage-1 fit
  resid_full <- residuals(s1)
  dt[, `:=`(resid_sq = NA_real_, resid_cube = NA_real_)]
  dt[obs(s1), `:=`(resid_sq   = resid_full^2,
                   resid_cube = resid_full^3)]

  # Stage 2a: OLS variance -> fitted h_hat for the FGLS weights
  s2a <- feols(fml_var, data = dt,
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, h_hat := NA_real_]
  dt[obs(s2a), h_hat := pmax(fitted(s2a), 1e-6)]

  # Stage 2b: FGLS variance (point estimate reported for the variance stage)
  s2b <- feols(fml_var, data = dt,
               weights  = ~I(1 / h_hat),
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)

  # Stage 3: OLS skewness (third central moment). The clustered SE attached here
  # is the NAIVE one and is deliberately overwritten by the bootstrap below.
  s3 <- feols(fml_skew, data = dt,
              cluster  = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)

  coef_hat_mean <- coef(s1)
  coef_hat_var  <- coef(s2b)
  coef_hat_skew <- coef(s3)

  coef_boot_mean <- matrix(NA_real_, nrow = B, ncol = length(coef_hat_mean),
                           dimnames = list(NULL, names(coef_hat_mean)))
  coef_boot_var  <- matrix(NA_real_, nrow = B, ncol = length(coef_hat_var),
                           dimnames = list(NULL, names(coef_hat_var)))
  coef_boot_skew <- matrix(NA_real_, nrow = B, ncol = length(coef_hat_skew),
                           dimnames = list(NULL, names(coef_hat_skew)))

  # ── Bootstrap ───────────────────────────────────────────────────────────────
  for (b in seq_len(B)) {

    cat("\rBootstrap iteration", b, "of", B, "  ")

    # Field-level frequency weights (resample fields with replacement)
    drawn     <- table(sample(fields, n_fields, replace = TRUE))
    drawn_tbl <- data.table(
      tile_field_ID = names(drawn),
      boot_w        = as.numeric(drawn)
    )
    dt[drawn_tbl, boot_w := i.boot_w, on = "tile_field_ID"]
    dt[is.na(boot_w), boot_w := 0L]

    tryCatch({

      dt_b <- dt[boot_w > 0]

      # Stage 1 (replicate) — the single fit that feeds BOTH higher moments
      b1 <- feols(fml_mean, data = dt_b,
                  weights  = ~boot_w, vcov = "iid",
                  nthreads = n_threads, warn = FALSE, notes = FALSE)

      # Store the replicate mean coefficients. Not needed to make the mean SEs
      # valid (the mean is estimated on observed yield, so its analytic clustered
      # SEs are already correct), but capturing it here gives (i) bootstrap mean
      # SEs for table consistency and (ii) the mean's covariance with the higher
      # moments, since all three come from the same resampled draw.
      coef_boot_mean[b, names(coef(b1))] <- coef(b1)

      resid_b <- residuals(b1)
      dt_b[, `:=`(resid_sq_b = NA_real_, resid_cube_b = NA_real_)]
      dt_b[obs(b1), `:=`(resid_sq_b   = resid_b^2,
                         resid_cube_b = resid_b^3)]

      # Stage 2a (replicate): OLS variance -> h_hat_b
      b2a <- feols(fml_var_b, data = dt_b,
                   weights  = ~boot_w, vcov = "iid",
                   nthreads = n_threads, warn = FALSE, notes = FALSE)
      dt_b[, h_hat_b := NA_real_]
      dt_b[obs(b2a), h_hat_b := pmax(fitted(b2a), 1e-6)]
      dt_b[, combined_w := boot_w / h_hat_b]

      # Stage 2b (replicate): FGLS variance
      b2b <- feols(fml_var_b, data = dt_b,
                   weights  = ~combined_w, vcov = "iid",
                   nthreads = n_threads, warn = FALSE, notes = FALSE)

      # Stage 3 (replicate): OLS skewness on the cubed replicate residual.
      # Frequency-weighted by boot_w only; no h_hat weighting for the 3rd moment.
      b3 <- feols(fml_skew_b, data = dt_b,
                  weights  = ~boot_w, vcov = "iid",
                  nthreads = n_threads, warn = FALSE, notes = FALSE)

      # Align by coefficient name in case a replicate drops a collinear level
      coef_boot_var[b,  names(coef(b2b))] <- coef(b2b)
      coef_boot_skew[b, names(coef(b3))]  <- coef(b3)

    }, error = function(e) {
      cat("Bootstrap iteration", b, "failed:", conditionMessage(e), "\n")
    })

    dt[, boot_w := NULL]

    if (b %% 100 == 0) {
      cat("Bootstrap iteration", b, "of", B, "\n")
      saveRDS(list(mean = coef_boot_mean,
                   var  = coef_boot_var,
                   skew = coef_boot_skew),
              paste0("boot_progress_", b, ".rds"))
      gc()
    }
  }

  # ── Summaries ────────────────────────────────────────────────────────────────
  summarise <- function(coef_hat, draws) {
    list(
      coef = coef_hat,
      se   = apply(draws, 2, sd, na.rm = TRUE),
      ci   = apply(draws, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE),
      # bootstrap two-sided p-value: 2 * min(share below 0, share above 0)
      p    = apply(draws, 2, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(NA_real_)
        2 * min(mean(x <= 0), mean(x >= 0))
      })
    )
  }

  list(
    mean       = summarise(coef_hat_mean, coef_boot_mean),
    variance   = summarise(coef_hat_var,  coef_boot_var),
    skewness   = summarise(coef_hat_skew, coef_boot_skew),
    fit_mean   = s1,
    fit_var    = s2b,
    fit_skew   = s3,
    boot_mean  = coef_boot_mean,
    boot_var   = coef_boot_var,
    boot_skew  = coef_boot_skew
  )
}

# ── Run ─────────────────────────────────────────────────────────────────────
gc()
boot_moments <- boot_jp_moments(corn_jp_data, fml_mean, fml_var, fml_skew,
                                B = 999, seed = 42)

saveRDS(boot_moments,
        "D:/Crop data/boot_moments.rds", compress = TRUE)

# ── Assemble the skewness (stage-3) table with CORRECTED SEs ──────────────────
# Keep only rotation-sequence coefficients; attach bootstrap SE / CI / p-value.
skew_tab <- data.table(
  term      = names(boot_moments$skewness$coef),
  estimate  = boot_moments$skewness$coef,
  boot_se   = boot_moments$skewness$se,
  ci_lower  = boot_moments$skewness$ci[1, ],
  ci_upper  = boot_moments$skewness$ci[2, ],
  boot_p    = boot_moments$skewness$p
)[grepl("rot_crop", term)][
  , term := gsub("rot_crop", "", term)][order(-estimate)]

print(skew_tab)