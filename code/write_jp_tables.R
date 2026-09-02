# ==============================================================================
# Export bootstrap Just-Pope stages to LaTeX tables matching corn_jp_var.tex
# ------------------------------------------------------------------------------
# Regenerates tables/corn_jp_var.tex from a boot_jp_moments() object, so the
# table and the manuscript text rest on the SAME (bootstrap) inference. Point
# estimates are the full-sample stage coefficients; SEs and significance stars
# come from the bootstrap draws, NOT from the analytic clustered vcov that
# fixest attaches to a single feols() fit.
#
#   format_jp_moment_tex()  pure string builder (no fixest dependency)
#   write_jp_var_tex()      variance stage  -> tables/corn_jp_var.tex
#
# The stage-3 skewness table is NOT written here. The manuscript reports the
# standardized skewness as column (2) of tables/corn_jp_moments.tex (written by
# corn_analysis_full.R, already on bootstrap SEs). The old write_jp_skew_tex()
# emitted tables/corn_jp_skew.tex on the RAW resid_cube scale (~s^3 larger) under
# a "standardized" caption, so it was deleted (2026-09-02) rather than reconciled.
#
# IMPORTANT -- variance stage point estimates change from OLS to FGLS:
#   The bootstrap stores the FGLS variance coefficients (coef(s2b)), which is the
#   efficient Just-Pope variance estimator and the quantity the bootstrap
#   actually resamples. The variance table currently in the repo is captioned
#   "(OLS)" and shows the stage-2a OLS coefficients. Regenerating with
#   write_jp_var_tex() therefore changes the coefficients, not just the SEs, and
#   the caption is updated to "(FGLS, bootstrap SE)". If you intend to keep the
#   OLS point estimates in the paper, do NOT regenerate that table from the
#   bootstrap -- report the analytic OLS table and cite the bootstrap only for
#   the stages whose regressors are generated (variance-as-generated / skewness).
# ==============================================================================

# ── 4-significant-figure formatter, matching etable's default ────────────────
.fmt_sig <- function(x, digits = 4L) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    formatC(v, digits = digits, format = "g")
  }, character(1))
}

# ── Significance stars from a p-value (thresholds match the paper: .01/.05/.10)
.stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.01, "$^{***}$",
  ifelse(p < 0.05, "$^{**}$",
  ifelse(p < 0.10, "$^{*}$", ""))))
}

# ── Pure LaTeX builder ────────────────────────────────────────────────────────
# coef, se, pval : named numeric vectors over the SAME terms, restricted to the
#                  rotation-sequence coefficients, names = sequence labels with
#                  any "rot_crop" prefix stripped.
# fe_vars        : fixed-effect variable names (raw); underscores escaped here.
format_jp_moment_tex <- function(coef, se, pval, n_obs, r2, wr2,
                                 dep_label,
                                 caption,
                                 label,
                                 fe_vars   = c("tile_field_ID", "year"),
                                 col_width = 32L) {

  pad <- function(s) formatC(s, flag = "-", width = col_width)

  est  <- .fmt_sig(coef)
  sev  <- .fmt_sig(se)
  star <- .stars(pval)
  coef_rows <- sprintf("      %s & %s%s (%s)\\\\   ",
                       pad(names(coef)), est, star, sev)

  n_fmt   <- formatC(n_obs, big.mark = ",", format = "d")
  r2_fmt  <- formatC(r2,  digits = 5, format = "f")
  wr2_fmt <- formatC(wr2, digits = 5, format = "f")

  fe_rows <- sprintf("      %s & $\\checkmark$\\\\   ",
                     pad(paste0(gsub("_", "\\\\_", fe_vars), " fixed effects")))

  c(
    "",
    "\\begin{table}[H]",
    sprintf("   \\caption{\\label{%s} %s}", label, caption),
    "   \\bigskip",
    "   \\centering",
    "   \\begin{tabular}{lc}",
    "      \\toprule",
    sprintf("      %s & %s\\\\   ", pad(""), dep_label),
    sprintf("      %s & (1)\\\\  ", pad("")),
    "      \\midrule ",
    coef_rows,
    "       \\\\",
    sprintf("      %s & %s\\\\  ", pad("Observations"), n_fmt),
    sprintf("      %s & %s\\\\  ", pad("R$^2$"),        r2_fmt),
    sprintf("      %s & %s\\\\  ", pad("Within R$^2$"), wr2_fmt),
    "       \\\\",
    fe_rows,
    "      \\bottomrule",
    "   \\end{tabular}",
    "\\end{table}",
    "", ""
  )
}

# ── Shared core: pull stats off a fixest fit, emit the file ──────────────────
.emit_jp_tex <- function(stat, fit, path, dep_label, caption, label,
                         term_prefix = "rot_crop", ...) {

  keep  <- grepl(term_prefix, names(stat$coef))
  clean <- function(v) { v <- v[keep]; names(v) <- gsub(term_prefix, "", names(v)); v }

  n_obs   <- fit$nobs
  r2_val  <- fixest::r2(fit, "r2")
  wr2_val <- fixest::r2(fit, "wr2")
  fe_vars <- if (!is.null(fit$fixef_vars)) fit$fixef_vars else c("tile_field_ID", "year")

  tex <- format_jp_moment_tex(clean(stat$coef), clean(stat$se), clean(stat$p),
                              n_obs = n_obs, r2 = r2_val, wr2 = wr2_val,
                              dep_label = dep_label, caption = caption,
                              label = label, fe_vars = fe_vars, ...)

  writeLines(tex, path)
  message("Wrote ", path, " (", sum(keep), " sequence coefficients, bootstrap SEs)")
  invisible(path)
}

# Absolute output directory -- avoids the "../tables/" vs "tables/" mismatch
# that breaks these functions depending on the caller's working directory.
jp_tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"

# ── Variance stage -> tables/corn_jp_var.tex (FGLS point est. + bootstrap SE) ──
write_jp_var_tex <- function(boot_obj,
                             fit  = boot_obj$fit_var,
                             path = paste0(jp_tab_dir, "corn_jp_var.tex"),
                             caption = paste("Stage 2 --- Corn yield conditional",
                                             "variance (FGLS, bootstrap SE)"),
                             dep_label = "resid\\_sq",
                             label = "tab:corn_jp_var", ...) {
  .emit_jp_tex(boot_obj$variance, fit, path, dep_label, caption, label, ...)
}

# ── First-stage (mean) bootstrap vcov -> tables/corn_jp_mean_vcov.rds ─────────
# The mean-stage LaTeX table (tab:corn_jp_mean) is built elsewhere from the
# analytic clustered SEs (those are already valid -- see the note at the top of
# just_pope_bootstrap_moments.R -- the mean is estimated on observed yield, not
# a generated regressor). What the analytic output does NOT give is the
# covariance BETWEEN rotation-sequence coefficients, which is needed for the
# SE of any linear combination of them (e.g. the rotation "score" in
# corn_analysis_full.R). boot_obj$boot_mean holds the B replicate draws, so
# their sample covariance is the bootstrap estimate of that full vcov.
write_jp_mean_vcov <- function(boot_obj,
                               draws = boot_obj$boot_mean,
                               path  = paste0(jp_tab_dir, "corn_jp_mean_vcov.rds"),
                               term_prefix = "rot_crop") {

  keep      <- grepl(term_prefix, colnames(draws))
  draws_kp  <- draws[, keep, drop = FALSE]
  colnames(draws_kp) <- gsub(term_prefix, "", colnames(draws_kp))

  vcov_mat <- cov(draws_kp, use = "pairwise.complete.obs")

  saveRDS(vcov_mat, path)
  message("Wrote ", path, " (", nrow(vcov_mat), " x ", ncol(vcov_mat),
          " bootstrap vcov, ", sum(complete.cases(draws_kp)), " complete replicates)")
  invisible(vcov_mat)
}

# ── Usage (after boot_moments <- boot_jp_moments(...)) ───────────────────────
# source("code/write_jp_tables.R")
write_jp_var_tex(boot_moments)                  # tables/corn_jp_var.tex  (OLS->FGLS!)
write_jp_mean_vcov(boot_moments)                # tables/corn_jp_mean_vcov.rds
#
# NOTE: with B = 499 the bootstrap p-value is granular in steps of ~2/499 ~ 0.004,
# fine for the .05/.10 stars but coarse for the .01 threshold; run B >= 999 for
# stable *** classifications.
