## ============================================================================
## run_jp_bootstrap_standalone.R
## ----------------------------------------------------------------------------
## Standalone driver for the three-stage Just-Pope moment bootstrap.
##
## Runs, end to end, without corn_analysis_full.R / main_analysis.R:
##   1. loads shared setup (rotation_setup_wa.R) and builds corn_jp_data
##      exactly as corn_analysis_full.R does (same filters, same recode);
##   2. runs boot_jp_moments() from just_pope_bootstrap_moments.R
##      (B replicates, field-level pairs bootstrap over all three stages);
##   3. writes boot_moments.rds;
##   4. builds the bootstrap vcov matrices with boot_stage_vcov() from
##      jp_boot_vcov.R and writes them as BOTH .txt (pipeline inputs read by
##      corn_analysis_full.R) and .csv;
##   5. writes the per-stage moment summary tables (estimate / boot SE / 95%
##      CI / bootstrap p) as .csv, and the three-stage etable as .tex.
##
## Nothing here depends on a prior run: just_pope_bootstrap_moments.R and
## jp_boot_vcov.R are sourced with load-only flags so only their functions
## (boot_jp_moments, boot_stage_vcov) are pulled in; the driver owns the run
## and every output path.
##
## Usage:  Rscript code/run_jp_bootstrap_standalone.R
##   or:   source("code/run_jp_bootstrap_standalone.R")  from an R session
## Runtime: ~hours at B = 999 on the full IL corn sample (~1.8M field-years).
## ============================================================================

# ── Config ──────────────────────────────────────────────────────────────────
CODE_DIR   <- "C:/Users/vf006/Box/crop_rotations_and_losses/code"
DATA_PARQUET <- "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet"
RDS_PATH   <- "D:/Crop data/boot_moments.rds"   # kept where jp_boot_vcov.R expects it
B          <- 999
SEED       <- 42

setwd(CODE_DIR)

library(arrow)
library(statar)

source("rotation_setup_wa.R")   # libs, make_jp_formula(), all_controls,
                                # all_controls_cols, corn_soy_patterns,
                                # rot_features, to_letters(), tab_dir, fig_dir

OUT_DIR <- tab_dir              # tables/ ; corn_analysis_full.R reads the .txt here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load + build corn_jp_data  (verbatim from corn_analysis_full.R) ──────────
cat("Loading corn data...\n")
corn_df <- read_parquet(DATA_PARQUET)

corn_df <- corn_df |>
  filter(STATE_ABBR == "IL") |>
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         corn_yield = corn_yield / 62.77) |>
  arrange(tile_field_ID, year)

cat("Corn raw rows:", nrow(corn_df), "\n")

setDT(corn_df)

add_crop_year <- function(dt) {
  stopifnot(is.data.table(dt))
  if ("crop_year" %in% names(dt)) dt[, crop_year := NULL]
  crop_cols  <- grep("^crop_[0-9]{4}$", names(dt), value = TRUE)
  crop_years <- as.integer(sub("crop_", "", crop_cols))
  stopifnot(length(crop_cols) > 0, !anyNA(crop_years))
  m       <- as.matrix(dt[, ..crop_cols])
  col_idx <- match(dt$year, crop_years)
  stopifnot(length(col_idx) == nrow(dt))
  dt[, crop_year := m[cbind(seq_len(.N), col_idx)]]
  invisible(dt)
}
add_crop_year(corn_df)

corn_df <- corn_df |>
  rename(crop_0 = crop_year,
         crop_1 = prioryr_crop,
         crop_2 = prior2yr_crop,
         crop_3 = prior3yr_crop,
         crop_4 = prior4yr_crop,
         crop_5 = prior5yr_crop,
         crop_6 = prior6yr_crop)

source("rci_vectorized.R")
corn_df <- corn_df |>
  mutate(RCI = rci(crop_0, crop_1, crop_2, crop_3, crop_4, crop_5))

source("cdl_recode.R")
recode_cdl(corn_df, cols = paste0("crop_", 0:6))

corn_df <- corn_df |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-",
                           crop_2, "-", crop_1, "-", crop_0)) |>
  add_degree_days()

corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  mutate(
    seq_num  = strsplit(rot_crop, "-", fixed = TRUE),
    soy_lag1 = sapply(seq_num, \(v) as.integer(v[5] == "5")),
    soy_lag2 = sapply(seq_num, \(v) as.integer(v[4] == "5")),
    soy_lag3 = sapply(seq_num, \(v) as.integer(v[3] == "5")),
    soy_lag4 = sapply(seq_num, \(v) as.integer(v[2] == "5")),
    soy_lag5 = sapply(seq_num, \(v) as.integer(v[1] == "5")),
    rot_crop = to_letters(rot_crop),
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")
  ) |>
  select(-seq_num) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))

cat("Corn analysis sample:", nrow(corn_jp_data), "rows\n")
rm(corn_df); gc()

# ── Load boot_jp_moments() and boot_stage_vcov() without running their tails ──
JP_BOOT_LOAD_ONLY <- TRUE
source("just_pope_bootstrap_moments.R")   # -> boot_jp_moments(), fml_mean/var/skew

JP_VCOV_LOAD_ONLY <- TRUE
source("jp_boot_vcov.R")                  # -> boot_stage_vcov()

# ── 1. Run the bootstrap ────────────────────────────────────────────────────
gc()
boot_moments <- boot_jp_moments(corn_jp_data, fml_mean, fml_var, fml_skew,
                                B = B, seed = SEED)

saveRDS(boot_moments, RDS_PATH, compress = TRUE)
cat("Wrote", RDS_PATH, "\n")

# ── 2. Bootstrap vcov matrices  (as in jp_boot_vcov.R) ──────────────────────
jp_vcov_mean <- boot_stage_vcov(boot_moments$boot_mean)
jp_vcov_var  <- boot_stage_vcov(boot_moments$boot_var)
jp_vcov_skew <- boot_stage_vcov(boot_moments$boot_skew)

# .txt  (read back by corn_analysis_full.R::read_vcov_txt from tab_dir)
write.table(jp_vcov_mean, file.path(OUT_DIR, "jp_vcov_mean.txt"))
write.table(jp_vcov_var,  file.path(OUT_DIR, "jp_vcov_var.txt"))
write.table(jp_vcov_skew, file.path(OUT_DIR, "jp_vcov_skew.txt"))

# .csv  (same matrices, comma-separated, row names in column 1)
write.csv(jp_vcov_mean, file.path(OUT_DIR, "jp_vcov_mean.csv"))
write.csv(jp_vcov_var,  file.path(OUT_DIR, "jp_vcov_var.csv"))
write.csv(jp_vcov_skew, file.path(OUT_DIR, "jp_vcov_skew.csv"))

# ── 3. Per-stage moment summary tables -> .csv ──────────────────────────────
# All RHS terms; the rotation-sequence rows are those matching ^rot_crop.
moment_csv <- function(s, path) {
  out <- data.frame(
    term      = names(s$coef),
    estimate  = as.numeric(s$coef),
    boot_se   = as.numeric(s$se),
    ci_lower  = as.numeric(s$ci[1, ]),
    ci_upper  = as.numeric(s$ci[2, ]),
    boot_p    = as.numeric(s$p),
    row.names = NULL,
    check.names = FALSE
  )
  out$is_sequence <- grepl("^rot_crop", out$term)
  out$term <- sub("^rot_crop", "", out$term)
  write.csv(out, path, row.names = FALSE)
  invisible(out)
}
moment_csv(boot_moments$mean,     file.path(OUT_DIR, "jp_boot_moments_mean.csv"))
moment_csv(boot_moments$variance, file.path(OUT_DIR, "jp_boot_moments_variance.csv"))
moment_csv(boot_moments$skewness, file.path(OUT_DIR, "jp_boot_moments_skewness.csv"))

# ── 4. Three-stage etable with bootstrap vcov -> .tex ───────────────────────
etable(boot_moments$fit_mean, boot_moments$fit_var, boot_moments$fit_skew,
       vcov    = list(jp_vcov_mean, jp_vcov_var, jp_vcov_skew),
       keep    = "^rot_crop",
       headers = c("Mean", "Variance", "Skewness"),
       title   = "Just-Pope moments (bootstrap SE)",
       tex     = TRUE,
       replace = TRUE,
       file    = file.path(OUT_DIR, "jp_boot_moments.tex"))

# console copy
print(etable(boot_moments$fit_mean, boot_moments$fit_var, boot_moments$fit_skew,
             vcov    = list(jp_vcov_mean, jp_vcov_var, jp_vcov_skew),
             keep    = "^rot_crop",
             headers = c("Mean", "Variance", "Skewness"),
             title   = "Just-Pope moments (bootstrap SE)"))

cat("\nDone. Outputs in:", OUT_DIR, "\n",
    "  jp_vcov_{mean,var,skew}.{txt,csv}\n",
    "  jp_boot_moments_{mean,variance,skewness}.csv\n",
    "  jp_boot_moments.tex\n",
    "and", RDS_PATH, "\n")
gc()