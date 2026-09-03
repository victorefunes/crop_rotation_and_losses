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
#
# IMPORTANT: column (2) (Corn variance) has a generated-regressor dependent
# variable (squared stage-1 residual), so the analytic two-way (county, year)
# clustered SEs that etable() prints are NOT the inference of record. The
# post-processing block below overwrites the column (2) SEs and significance
# stars in tables/zvector.tex with the B=999 field-level pairs cluster bootstrap
# results from code/zvector_bootstrap_var.R (tables/zvector_boot_var.txt).
# Point estimates are identical; only the column (2) SEs/stars change.
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
         "_Controls" = c("Yes", "Yes")
       ),
       file = paste0(tab_dir, "zvector.tex"))

cat("Z-vector table saved to:", paste0(tab_dir, "zvector.tex"), "\n")

# ── Overwrite column (2) SEs with the field-pairs bootstrap ───────────────────
# (generated-regressor inference of record; see code/zvector_bootstrap_var.R)
boot_path <- paste0(tab_dir, "zvector_boot_var.txt")
if (file.exists(boot_path)) {
  bv    <- read.table(boot_path, header = TRUE, stringsAsFactors = FALSE)
  stars <- function(p) {
    if (is.na(p))      ""
    else if (p < 0.01) "$^{***}$"
    else if (p < 0.05) "$^{**}$"
    else if (p < 0.10) "$^{*}$"
    else               ""
  }
  g4  <- function(x) formatC(signif(x, 4), format = "g")
  lbl <- c(late_soy = "Latest soy rotation",
           soy_gap  = "Soy gap",
           soy_cons = "Number of consecutive soy years")
  zt  <- readLines(paste0(tab_dir, "zvector.tex"))
  for (i in seq_len(nrow(bv))) {
    term  <- bv$term[i]
    if (!term %in% names(lbl)) next
    cell  <- sprintf("%s%s (%s)", g4(bv$estimate[i]), stars(bv$boot_p[i]),
                     g4(bv$boot_se[i]))
    rows  <- grep(paste0("^\\s*", lbl[[term]], "\\s+&"), zt)
    for (r in rows) {
      parts <- strsplit(zt[r], "&", fixed = TRUE)[[1]]
      if (length(parts) == 3) {
        parts[3] <- sprintf(" %s\\\\", cell)      # trailing LaTeX row break
        zt[r]    <- paste(parts, collapse = "&")
      }
    }
  }
  writeLines(zt, paste0(tab_dir, "zvector.tex"))
  cat("Column (2) SEs/stars replaced with field-pairs bootstrap",
      "(tables/zvector_boot_var.txt).\n")
} else {
  cat("WARNING:", boot_path, "not found -- column (2) left with analytic SEs.\n",
      "Run code/zvector_bootstrap_var.R first.\n")
}

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
