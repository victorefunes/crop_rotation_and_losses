## ============================================================================
## tables_combined.R
## Produces tables that require results from corn_analysis.R.
## Run this after corn_analysis.R has completed and saved its model
## objects to disk.
##
## Prerequisites:
##   corn_z_models.rds  — from corn_analysis.R
##
## Tables produced:
##   tab:zvector   — Effect of rotation patterns on yields (main result, Table 1)
## ============================================================================
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")
source("rotation_setup_wa.R")

# ── Check prerequisites before loading ───────────────────────────────────────

rds_dir <- "D:/Crop data/"

corn_rds <- paste0(rds_dir, "corn_z_models.rds")

if (!file.exists(corn_rds)) {
  stop(
    "Missing RDS file. Run the following script first:\n",
    "  -> corn_analysis.R  (saves corn_z_models.rds)\n",
    call. = FALSE
  )
}

# ── Load saved model objects ──────────────────────────────────────────────────

corn_z <- readRDS(corn_rds)

corn_z_s1       <- corn_z$z_s1
corn_z_s2       <- corn_z$z_s2

cat("Corn Z-vector stage 1 obs:", nobs(corn_z_s1), "\n")
cat("Corn Z-vector stage 2 obs:", nobs(corn_z_s2), "\n")

# ── Coefficient label dictionary ──────────────────────────────────────────────

dict_z <- c(
  "late_soy" = "Latest soy rotation",
  "soy_gap"  = "Soy gap",
  "soy_cons" = "Number of consecutive soy years",
  "nsoy"     = "Number of soy harvests"
)

# ── Table: tab:zvector ────────────────────────────────────────────────────────
# Replicates Table 1 from the draft PDF.
# Columns: QDANN mean, QDANN variance
# corn_z_s1 = QDANN mean equation
# corn_z_s2 = QDANN variance equation
tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"

etable(corn_z_s1, corn_z_s2,
       tex       = TRUE,
       dict      = dict_z,
       headers   = c("Corn mean", "Corn variance"),
       keep_raw  = c("late_soy", "soy_gap", "soy_cons", "nsoy"),
       se.below  = FALSE,
       style.tex = style.tex("aer"),
       replace   = TRUE,
       title     = "Effect of rotation patterns on yields",
       label     = "tab:zvector",
       extralines = list(
         "_Controls" = c("Yes", "Yes"),
         "_Year FE"  = c("\\checkmark", "\\checkmark"),
         "_FIPS FE"  = c("\\checkmark", "\\checkmark")
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
