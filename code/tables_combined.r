## ============================================================================
## tables_combined.R
## Produces tables that require results from both corn_analysis.R and
## soy_analysis.R. Run this after both analysis scripts have completed
## and saved their model objects to disk.
##
## Prerequisites:
##   corn_z_models.rds  — from corn_analysis.R
##   soy_z_models.rds   — from soy_analysis.R
##
## Tables produced:
##   tab:zvector   — Effect of rotation patterns on yields (main result, Table 1)
## ============================================================================
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")
source("rotation_setup.R")

# ── Check prerequisites before loading ───────────────────────────────────────

rds_dir <- "C:/Users/vf006/Documents/"

corn_rds <- paste0(rds_dir, "corn_z_models.rds")
soy_rds  <- paste0(rds_dir, "soy_z_models.rds")

missing <- c(
  if (!file.exists(corn_rds)) corn_rds,
  if (!file.exists(soy_rds))  soy_rds
)

if (length(missing) > 0) {
  stop(
    "Missing RDS file(s). Run the following scripts first:\n",
    if (!file.exists(corn_rds)) "  -> corn_analysis.R  (saves corn_z_models.rds)\n",
    if (!file.exists(soy_rds))  "  -> soy_analysis.R   (saves soy_z_models.rds)\n",
    call. = FALSE
  )
}

# ── Load saved model objects ──────────────────────────────────────────────────

corn_z <- readRDS(corn_rds)
soy_z  <- readRDS(soy_rds)

corn_z_s1 <- corn_z$s1
corn_z_s2 <- corn_z$s2
soy_z_s1  <- soy_z$s1
soy_z_s2  <- soy_z$s2

cat("Corn Z-vector stage 1 obs:", nobs(corn_z_s1), "\n")
cat("Corn Z-vector stage 2 obs:", nobs(corn_z_s2), "\n")
cat("Soy  Z-vector stage 1 obs:", nobs(soy_z_s1),  "\n")
cat("Soy  Z-vector stage 2 obs:", nobs(soy_z_s2),  "\n")

# ── Coefficient label dictionary ──────────────────────────────────────────────

dict_z <- c(
  "late_soy" = "Latest soy rotation",
  "soy_gap"  = "Soy gap",
  "soy_cons" = "Number of consecutive soy years",
  "nsoy"     = "Number of soy harvests"
)

# ── Table: tab:zvector ────────────────────────────────────────────────────────
# Replicates Table 1 from the draft PDF.
# Columns: QDANN mean, QDANN variance, SCYM mean, SCYM variance
# corn_z_s1 = QDANN mean equation
# corn_z_s2 = QDANN variance equation
# soy_z_s1  = SCYM mean equation  (rename headers accordingly if needed)
# soy_z_s2  = SCYM variance equation

etable(corn_z_s1, corn_z_s2, soy_z_s1, soy_z_s2,
       tex       = TRUE,
       dict      = dict_z,
       headers   = c("QDANN mean", "QDANN variance",
                     "SCYM mean",  "SCYM variance"),
       keep      = c("late_soy", "soy_gap", "soy_cons", "nsoy"),
       se.below  = FALSE,
       style.tex = style.tex("aer"),
       replace   = TRUE,
       title     = "Effect of rotation patterns on yields",
       label     = "tab:zvector",
       extralines = list(
         "_Controls" = c("Yes", "Yes", "Yes", "Yes"),
         "_Year FE"  = c("\\checkmark", "\\checkmark",
                         "\\checkmark", "\\checkmark"),
         "_FIPS FE"  = c("\\checkmark", "\\checkmark",
                         "\\checkmark", "\\checkmark")
       ),
       file = paste0(tab_dir, "zvector.tex"))

cat("Z-vector table saved to:", paste0(tab_dir, "zvector.tex"), "\n")

# ── Sanity check: compare to draft PDF Table 1 ────────────────────────────────
# Expected values from draft (QDANN mean column):
#   late_soy:  2.847***
#   soy_gap:   0.311***
#   nsoy:     -0.976***
#   soy_cons: -4.671***

cat("\nZ-vector coefficients vs draft PDF Table 1 (QDANN mean):\n")
broom::tidy(corn_z_s1) |>
  filter(term %in% c("late_soy", "soy_gap", "nsoy", "soy_cons")) |>
  select(term, estimate, std.error, p.value) |>
  mutate(
    draft = c(2.847, 0.311, -0.976, -4.671)[
      match(term, c("late_soy", "soy_gap", "nsoy", "soy_cons"))],
    diff  = round(estimate - draft, 4)
  ) |>
  print()