# ==============================================================================
# Bootstrap vcov matrices for the three Just-Pope stages (mean, variance,
# skewness), for use as etable()'s `vcov` argument in place of the analytic
# clustered vcov that fixest attaches to a standalone feols() fit.
#
# boot_moments$boot_mean / boot_var / boot_skew (from just_pope_bootstrap_moments.R)
# are each a B x k matrix of replicate coefficients over the FULL RHS of their
# stage (rotation sequences + controls) -- not restricted to rot_crop -- so
# their sample covariance already has the coefficient coverage etable() needs:
# a vcov matrix passed to etable()/summary.fixest() must contain every
# coefficient in the corresponding model, not a subset.
# ==============================================================================

library(fixest)
library(tidyverse)

boot_stage_vcov <- function(draws) {
  cov(draws, use = "pairwise.complete.obs")
}

# Sourcing this file with `JP_VCOV_LOAD_ONLY <- TRUE` set beforehand loads
# boot_stage_vcov() only and SKIPS the read/write/print below (used by
# run_jp_bootstrap_standalone.R). Run the file directly to execute as before.
if (!exists("JP_VCOV_LOAD_ONLY") || !isTRUE(JP_VCOV_LOAD_ONLY)) {

boot_moments <- readRDS("D:/Crop data/boot_moments.rds")

jp_vcov_mean <- boot_stage_vcov(boot_moments$boot_mean)
jp_vcov_var  <- boot_stage_vcov(boot_moments$boot_var)
jp_vcov_skew <- boot_stage_vcov(boot_moments$boot_skew)

write.table(jp_vcov_mean, "D:/Crop data/jp_vcov_mean.txt")
write.table(jp_vcov_var, "D:/Crop data/jp_vcov_var.txt")
write.table(jp_vcov_skew, "D:/Crop data/jp_vcov_skew.txt")

# ── Usage ──────────────────────────────────────────────────────────────────
# etable() takes one vcov spec per model, in the same order as the models.
etable(boot_moments$fit_mean, boot_moments$fit_var, boot_moments$fit_skew,
       vcov = list(jp_vcov_mean, jp_vcov_var, jp_vcov_skew),
       keep = "^rot_crop",
       headers = c("Mean", "Variance", "Skewness"),
       title   = "Just-Pope moments (bootstrap SE)")

}  # end if (!JP_VCOV_LOAD_ONLY)

rm(boot_moments, jp_vcov_mean, jp_vcov_var, jp_vcov_skew)
gc()
