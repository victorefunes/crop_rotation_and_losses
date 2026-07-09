## ============================================================================
## save_models_lean.R
## Helper function to extract only what tables_combined.R needs from a
## fitted feols object and a bootstrap result, discarding everything else.
## Source this in corn_analysis.R and soy_analysis.R before saving.
## ============================================================================

lean_feols <- function(fit) {
  ## Keep only the coefficient vector, SE vector, and key scalars.
  ## Strips the model matrix, data, residuals, and fixest internals.
  list(
    coef     = coef(fit),
    se       = se(fit),
    tstat    = tstat(fit),
    pvalue   = pvalue(fit),
    ci       = confint(fit),
    nobs     = nobs(fit),
    r2       = r2(fit, type = "r2"),
    wr2      = r2(fit, type = "wr2"),
    fml      = formula(fit),
    fixef_vars = fixef_vars(fit),
    call     = fit$call
  )
}

lean_boot <- function(boot_result) {
  ## From the boot list returned by boot_jp_fgls(), keep:
  ##   - point estimates and bootstrap SEs/CIs  (small vectors)
  ##   - boot_draws matrix                       (B x n_coef — keep for CI)
  ##   - lean versions of the three fitted models
  ## Drop: full feols objects with embedded data references.
  list(
    coef        = boot_result$coef,
    se          = boot_result$se,
    ci          = boot_result$ci,
    boot_draws  = boot_result$boot_draws,       # B x n_coef matrix
    fit_mean    = lean_feols(boot_result$fit_mean),
    fit_var_ols = lean_feols(boot_result$fit_var_ols),
    fit_var_fgls = lean_feols(boot_result$fit_var_fgls)
  )
}
