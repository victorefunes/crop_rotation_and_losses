# ==============================================================================
# Three-stage bootstrap for the parsimonious (Z-vector) corn variance stage
# ------------------------------------------------------------------------------
# Table tab:zvector column (2) (resid_sq_z ~ late_soy + soy_gap + soy_cons + ...)
# is currently reported with analytic two-way-clustered SEs. Like the
# sequence-dummy Just-Pope variance stage, its dependent variable is a GENERATED
# regressor (squared residual from the Z-vector stage-1 mean fit), so the
# analytic SEs understate uncertainty. This script reproduces the field-level
# pairs cluster bootstrap of just_pope_bootstrap_moments.R for the 3-feature RHS.
#
# Usage
# ------------------------------------------------------------------------------
#   source("rotation_setup_wa.R")          # make_jp_formula, all_controls_fgls
#   # ... build corn_jp_data with late_soy / soy_gap / soy_cons (as in
#   #     corn_analysis_full.R lines ~572-590) ...
#   source("zvector_bootstrap_var.R")
#
# Produces:
#   z_boot_var           list(coef, se, ci, p) for the stage-2 coefficients
#   tables/zvector_boot_var.txt   coef / boot SE / boot 95% CI / boot p, one row
#                                 per Z feature (for pasting into tab:zvector notes)
# ==============================================================================

library(fixest)
library(data.table)

boot_zvar <- function(dt, feats = c("late_soy", "soy_gap", "soy_cons"),
                      controls = all_controls_fgls, B = 999, seed = 42) {

  dt <- as.data.table(dt)
  rhs <- paste(feats, collapse = " + ")
  fml_mean <- make_jp_formula("corn_yield", rhs, controls)
  fml_var  <- make_jp_formula("resid_sq_z", rhs, controls)

  # ── Singleton removal (match the moment bootstrap) ──────────────────────────
  repeat {
    n_before <- nrow(dt)
    dt <- dt[dt[, .N, by = tile_field_ID][N > 1], on = "tile_field_ID"]
    dt <- dt[dt[, .N, by = year][N > 1],          on = "year"]
    if (nrow(dt) == n_before) break
  }
  cat("Rows after singleton removal:", nrow(dt), "\n")

  n_threads <- max(1L, parallel::detectCores() - 1L)
  set.seed(seed)
  fields   <- unique(dt$tile_field_ID)
  n_fields <- length(fields)

  # ── Point estimates on the full sample ─────────────────────────────────────
  s1 <- feols(fml_mean, data = dt, cluster = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, resid_sq_z := NA_real_]
  dt[obs(s1), resid_sq_z := residuals(s1)^2]

  s2 <- feols(fml_var, data = dt, cluster = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)

  coef_hat <- coef(s2)
  coef_boot <- matrix(NA_real_, nrow = B, ncol = length(coef_hat),
                      dimnames = list(NULL, names(coef_hat)))

  # ── Bootstrap: resample fields, rebuild resid_sq_z inside every replicate ───
  for (b in seq_len(B)) {
    cat("\rZ-vector var bootstrap", b, "of", B, "  ")
    drawn <- table(sample(fields, n_fields, replace = TRUE))
    dt[data.table(tile_field_ID = names(drawn), bw = as.numeric(drawn)),
       boot_w := i.bw, on = "tile_field_ID"]
    dt[is.na(boot_w), boot_w := 0]

    tryCatch({
      dt_b <- dt[boot_w > 0]
      b1 <- feols(fml_mean, data = dt_b, weights = ~boot_w, vcov = "iid",
                  nthreads = n_threads, warn = FALSE, notes = FALSE)
      dt_b[, resid_sq_z := NA_real_]
      dt_b[obs(b1), resid_sq_z := residuals(b1)^2]
      b2 <- feols(fml_var, data = dt_b, weights = ~boot_w, vcov = "iid",
                  nthreads = n_threads, warn = FALSE, notes = FALSE)
      coef_boot[b, names(coef(b2))] <- coef(b2)
    }, error = function(e)
      cat("\niter", b, "failed:", conditionMessage(e), "\n"))

    dt[, boot_w := NULL]
    if (b %% 100 == 0) { cat("  [", b, "]\n"); gc() }
  }

  se <- apply(coef_boot, 2, sd, na.rm = TRUE)
  ci <- apply(coef_boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  p  <- apply(coef_boot, 2, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    2 * min(mean(x <= 0), mean(x >= 0))
  })

  list(coef = coef_hat, se = se, ci = ci, p = p,
       fit_mean = s1, fit_var = s2, draws = coef_boot)
}

z_boot_var <- boot_zvar(corn_jp_data)

z_feat <- c("late_soy", "soy_gap", "soy_cons")
z_tab <- data.frame(
  term      = z_feat,
  estimate  = z_boot_var$coef[z_feat],
  analytic_se = summary(z_boot_var$fit_var)$se[z_feat],
  boot_se   = z_boot_var$se[z_feat],
  boot_ci_lo = z_boot_var$ci[1, z_feat],
  boot_ci_hi = z_boot_var$ci[2, z_feat],
  boot_p    = z_boot_var$p[z_feat],
  row.names = NULL
)
print(z_tab)
write.table(z_tab, file.path(tab_dir, "zvector_boot_var.txt"),
            row.names = FALSE)
cat("\nWrote", file.path(tab_dir, "zvector_boot_var.txt"), "\n")
cat("If boot_p for soy_gap crosses 0.10, drop the marginal-significance note in\n",
    "just_pope_rotations.tex Section~\\ref{subsec:rci_effects} intro (the paragraph\n",
    "beginning 'The variance equation for corn is considerably weaker').\n")
