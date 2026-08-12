library(arrow)
library(tidyverse)
library(statar)
library(fixest)
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

source("rotation_setup_soy_wa.R")

# ── Load soy data ────────────────────────────────────────────────────────────
# Load, correct, and add degree-day variables in one pipeline.
# soy_df is kept as the raw data source; soy_jp_data is the analysis sample.

cat("Loading soy data...\n")
soy_df <- read_parquet(
  "D:/Crop data/d_igis13soy_11_30_2025.with_rci.parquet") 

soy_df <- soy_df |>
  filter(STATE_ABBR == "IL") |>
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         soy_yield = soy_yield / 67.25)  |>  # kg/ha -> bu/ac
  arrange(tile_field_ID, year)

cat("Soy raw rows:", nrow(soy_df), "\n")

## Add present year crop variable
setDT(soy_df)  

add_crop_year <- function(dt) {
  stopifnot(is.data.table(dt))
  
  # Drop any stale crop_year column from a previous run, so the regex below
  # can't accidentally re-capture it
  if ("crop_year" %in% names(dt)) {
    dt[, crop_year := NULL]
  }
  
  # Match crop_YYYY columns only (exactly 4 digits) -- excludes crop_year,
  # rot_seq, or anything else that happens to start with "crop_"
  crop_cols  <- grep("^crop_[0-9]{4}$", names(dt), value = TRUE)
  crop_years <- as.integer(sub("crop_", "", crop_cols))
  
  stopifnot(length(crop_cols) > 0, !anyNA(crop_years))
  
  m       <- as.matrix(dt[, ..crop_cols])
  col_idx <- match(dt$year, crop_years)
  
  stopifnot(length(col_idx) == nrow(dt))
  
  dt[, crop_year := m[cbind(seq_len(.N), col_idx)]]
  invisible(dt)
}

add_crop_year(soy_df)

soy_df <- soy_df |>
  rename(crop_0 = crop_year,
         crop_1 = prioryr_crop,
         crop_2 = prior2yr_crop,
         crop_3 = prior3yr_crop,
         crop_4 = prior4yr_crop,
         crop_5 = prior5yr_crop,
         crop_6 = prior6yr_crop) |>
  filter(!(is.na(crop_5))) |>
  filter(!(is.na(crop_4))) |>
  filter(!(is.na(crop_3)))     

source("cdl_recode.R")

recode_cdl(soy_df, cols = paste0("crop_", 0:6))

soy_df <- soy_df |>
  rename(RCI = annual_RCI) |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-", 
          crop_2, "-", crop_1, "-", crop_0)) |>
  rci_correction() |>
  add_degree_days()

  # ── Rotation patterns ─────────────────────────────────────────────────────────

expand.grid(crop_0 = c("1","5"),
            crop_1 = c("1","5","24","36"),
            crop_2 = c("1","5","24","36"),
            crop_3 = c("1","5","24","36"),
            crop_4 = c("1","5","24","36"),
            crop_5 = c("1","5","24","36")) |>
  data.frame() |>
  mutate(pattern = paste(crop_5, crop_4, crop_3, crop_2, crop_1,
                         crop_0, sep = "-")) ->
  corn_soy_wheat

soy_jp_data <- soy_df |>
  filter(rot_crop %in% corn_soy_wheat$pattern) |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  filter(!is.na(soy_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                    ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))
 
cat("Soy analysis sample:", nrow(soy_jp_data), "rows\n")
 
# Free raw data — no longer needed
rm(soy_df); gc()
   
soy_jp_data |>
  mutate(rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "5-5-5-5-5-5")) ->
  soy_jp_data

# Remove rotation sequences that are too rare to estimate effects reliably. 
# The threshold is set by `min_freq`.
min_freq <- 100
seq_counts <- table(soy_jp_data$rot_crop)
keep_levels <- names(seq_counts)[seq_counts >= min_freq]   # e.g. min_freq = 30
soy_jp_data[["rot_crop"]] <- droplevels(factor(soy_jp_data[["rot_crop"]], levels = intersect(levels(soy_jp_data[["rot_crop"]]), keep_levels)))

source("lasso_rotation_selection.R")

soy_lasso <- lasso_select_sequences(
  data        = soy_jp_data,
  yvar        = "soy_yield",
  controls    = all_controls_cols,
  cluster_fml = ~COUNTY_FIPS
)

soy_lasso$selected_sequences

seq_names   <- grep("^rot_crop", names(coef(soy_lasso$refit_full_controls)), value = TRUE)
rename_dict <- setNames(to_letters(sub("^rot_crop", "", seq_names)), seq_names)

etable(soy_lasso$refit_full_controls, 
       tex = TRUE, 
       dict = rename_dict,
       cluster = ~COUNTY_FIPS,
       file = paste0(tab_dir, "soy_lasso.tex"), replace = TRUE,
       title = "LASSO-selected rotation sequence effects on soy yield")