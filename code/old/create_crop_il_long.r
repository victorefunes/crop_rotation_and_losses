library(tidyverse)
library(data.table)
library(statar)
library(fixest)
library(broom)
library(haven)
library(marginaleffects)
library(furrr)
library(hdm)
library(dotwhisker)
library(ggfortify)
theme_set(theme_bw())

#plan(multisession, workers = max(1, parallel::detectCores() - 1))

setwd("C:/Users/vf006/Box/Economic Analysis of Soil Health Practices")

corn_df <- fread("./Data and Data Descriptions/clean/corn_rci_il_long.csv")

## Single-out corn-soy rotations
expand.grid(crop_0 = c("1","5"), 
            crop_1 = c("1","5"), 
            crop_2 = c("1","5"), 
            crop_3 = c("1","5"), 
            crop_4 = c("1","5"),
            crop_5 = c("1","5")) |>
  data.frame() |>
  mutate(pattern = paste(crop_5, crop_4, crop_3, crop_2, crop_1, 
                         crop_0, sep = "-")) ->
  corn_soy_patterns

 ## RCI correctiom
 corn_df |>
  mutate(data_rm = case_when(
    rot_crop == "5-1-5-1-5-1" & RCI == 3.24 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 3 ~ 1,
    rot_crop == "5-1-5-1-1-1" & RCI == 2.24 ~ 1,
    rot_crop == "5-1-1-5-1-5" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-5-1-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-1-5-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 2 ~ 1,
    rot_crop == "1-5-1-1-1-5" & RCI == 1.73 ~ 1,
    rot_crop == "1-1-1-5-1-5" & RCI == 0 ~ 1,
    rot_crop == "1-1-1-1-1-5" & RCI == 0 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 2.45 ~ 1,
    rot_crop == "1-5-1-1-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-1-5-1-5" & RCI == 0 ~ 1,
    rot_crop == "1-5-1-5-1-5" & RCI == 2.45 ~ 1,
    .default = 0)) |>
  filter(data_rm == 0) |>
  select(-data_rm) ->
  corn_df 

gc()   

soy_df <- fread("./Data and Data Descriptions/clean/soy_rci_il_long.csv")   

# Order columns
soy_df <- soy_df[order(tile_field_ID, year)]

 ## RCI correctiom
soy_df |>
  mutate(data_rm = case_when(
    rot_crop == "5-1-5-1-5-1" & RCI == 3.24 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 3 ~ 1,
    rot_crop == "5-1-5-1-1-1" & RCI == 2.24 ~ 1,
    rot_crop == "5-1-1-5-1-5" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-5-1-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-1-5-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 2 ~ 1,
    rot_crop == "1-5-1-1-1-5" & RCI == 1.73 ~ 1,
    rot_crop == "1-1-1-5-1-5" & RCI == 0 ~ 1,
    rot_crop == "1-1-1-1-1-5" & RCI == 0 ~ 1,
    rot_crop == "5-1-5-1-5-1" & RCI == 2.45 ~ 1,
    rot_crop == "1-5-1-1-5-1" & RCI == 2.24 ~ 1,
    rot_crop == "1-5-1-5-1-5" & RCI == 0 ~ 1,
    rot_crop == "1-5-1-5-1-5" & RCI == 2.45 ~ 1,
    .default = 0)) |>
  filter(data_rm == 0) |>
  select(-data_rm) ->
  soy_df 

gc() 

crop_df <- full_join(corn_df, soy_df, by = c("tile_field_ID", "year", "lat", "lon",
                        "NAME", "STATE_ABBR", "COUNTY_FIPS", "STATE_FIPS")) 

crop_df |>
    mutate(vpdmax_4 = ifelse(is.na(vpdmax_4.x), vpdmax_4.y, vpdmax_4.x),
           vpdmax_5 = ifelse(is.na(vpdmax_5.x), vpdmax_5.y, vpdmax_5.x),
           vpdmax_6 = ifelse(is.na(vpdmax_6.x), vpdmax_6.y, vpdmax_6.x),
           vpdmax_7 = ifelse(is.na(vpdmax_7.x), vpdmax_7.y, vpdmax_7.x),
           vpdmax_8 = ifelse(is.na(vpdmax_8.x), vpdmax_8.y, vpdmax_8.x),
           vpdmax_9 = ifelse(is.na(vpdmax_9.x), vpdmax_9.y, vpdmax_9.x),
           def_1 = ifelse(is.na(def_1.x), def_1.y, def_1.x),
           def_2 = ifelse(is.na(def_2.x), def_2.y, def_2.x),
           def_3 = ifelse(is.na(def_3.x), def_3.y, def_3.x),
           def_4 = ifelse(is.na(def_4.x), def_4.y, def_4.x),
           def_5 = ifelse(is.na(def_5.x), def_5.y, def_5.x),
           def_6 = ifelse(is.na(def_6.x), def_6.y, def_6.x),
           def_7 = ifelse(is.na(def_7.x), def_7.y, def_7.x),
           def_8 = ifelse(is.na(def_8.x), def_8.y, def_8.x),
           def_9 = ifelse(is.na(def_9.x), def_9.y, def_9.x),   
           soil_1 = ifelse(is.na(soil_1.x), soil_1.y, soil_1.x),
           soil_2 = ifelse(is.na(soil_2.x), soil_2.y, soil_2.x),
           soil_3 = ifelse(is.na(soil_3.x), soil_3.y, soil_3.x),
           soil_4 = ifelse(is.na(soil_4.x), soil_4.y, soil_4.x),
           soil_5 = ifelse(is.na(soil_5.x), soil_5.y, soil_5.x),
           soil_6 = ifelse(is.na(soil_6.x), soil_6.y, soil_6.x),
           soil_7 = ifelse(is.na(soil_7.x), soil_7.y, soil_7.x),
           soil_8 = ifelse(is.na(soil_8.x), soil_8.y, soil_8.x),
           soil_9 = ifelse(is.na(soil_9.x), soil_9.y, soil_9.x),
           tmmn_1 = ifelse(is.na(tmmn_1.x), tmmn_1.y, tmmn_1.x),
           tmmn_2 = ifelse(is.na(tmmn_2.x), tmmn_2.y, tmmn_2.x),
           tmmn_3 = ifelse(is.na(tmmn_3.x), tmmn_3.y, tmmn_3.x),
           tmmn_4 = ifelse(is.na(tmmn_4.x), tmmn_4.y, tmmn_4.x),
           tmmn_5 = ifelse(is.na(tmmn_5.x), tmmn_5.y, tmmn_5.x),
           tmmn_6 = ifelse(is.na(tmmn_6.x), tmmn_6.y, tmmn_6.x),
           tmmn_7 = ifelse(is.na(tmmn_7.x), tmmn_7.y, tmmn_7.x),
           tmmn_8 = ifelse(is.na(tmmn_8.x), tmmn_8.y, tmmn_8.x),
           tmmn_9 = ifelse(is.na(tmmn_9.x), tmmn_9.y, tmmn_9.x),
           tmmx_1 = ifelse(is.na(tmmx_1.x), tmmx_1.y, tmmx_1.x),
           tmmx_2 = ifelse(is.na(tmmx_2.x), tmmx_2.y, tmmx_2.x),
           tmmx_3 = ifelse(is.na(tmmx_3.x), tmmx_3.y, tmmx_3.x),
           tmmx_4 = ifelse(is.na(tmmx_4.x), tmmx_4.y, tmmx_4.x),
           tmmx_5 = ifelse(is.na(tmmx_5.x), tmmx_5.y, tmmx_5.x),
           tmmx_6 = ifelse(is.na(tmmx_6.x), tmmx_6.y, tmmx_6.x),
           tmmx_7 = ifelse(is.na(tmmx_7.x), tmmx_7.y, tmmx_7.x),
           tmmx_8 = ifelse(is.na(tmmx_8.x), tmmx_8.y, tmmx_8.x),
           tmmx_9 = ifelse(is.na(tmmx_9.x), tmmx_9.y, tmmx_9.x),
           vpd_1 = ifelse(is.na(vpd_1.x), vpd_1.y, vpd_1.x),
           vpd_2 = ifelse(is.na(vpd_2.x), vpd_2.y, vpd_2.x),
           vpd_3 = ifelse(is.na(vpd_3.x), vpd_3.y, vpd_3.x),
           vpd_4 = ifelse(is.na(vpd_4.x), vpd_4.y, vpd_4.x),
           vpd_5 = ifelse(is.na(vpd_5.x), vpd_5.y, vpd_5.x),
           vpd_6 = ifelse(is.na(vpd_6.x), vpd_6.y, vpd_6.x),
           vpd_7 = ifelse(is.na(vpd_7.x), vpd_7.y, vpd_7.x),
           vpd_8 = ifelse(is.na(vpd_8.x), vpd_8.y, vpd_8.x),
           vpd_9 = ifelse(is.na(vpd_9.x), vpd_9.y, vpd_9.x),
           pr_1 = ifelse(is.na(pr_1.x), pr_1.y, pr_1.x),
           pr_2 = ifelse(is.na(pr_2.x), pr_2.y, pr_2.x),
           pr_3 = ifelse(is.na(pr_3.x), pr_3.y, pr_3.x),
           pr_4 = ifelse(is.na(pr_4.x), pr_4.y, pr_4.x),
           pr_5 = ifelse(is.na(pr_5.x), pr_5.y, pr_5.x),
           pr_6 = ifelse(is.na(pr_6.x), pr_6.y, pr_6.x),
           pr_7 = ifelse(is.na(pr_7.x), pr_7.y, pr_7.x),
           pr_8 = ifelse(is.na(pr_8.x), pr_8.y, pr_8.x),
           pr_9 = ifelse(is.na(pr_9.x), pr_9.y, pr_9.x),
           cGDD_4m = ifelse(is.na(cGDD_4m.x), cGDD_4m.y, cGDD_4m.x),
           cGDD_5m = ifelse(is.na(cGDD_5m.x), cGDD_5m.y, cGDD_5m.x),
           cGDD_6m = ifelse(is.na(cGDD_6m.x), cGDD_6m.y, cGDD_6m.x),
           cGDD_7m = ifelse(is.na(cGDD_7m.x), cGDD_7m.y, cGDD_7m.x),
           cGDD_8m = ifelse(is.na(cGDD_8m.x), cGDD_8m.y, cGDD_8m.x),
           cGDD_9m = ifelse(is.na(cGDD_9m.x), cGDD_9m.y, cGDD_9m.x),
           CC_probability = ifelse(is.na(CC_probability.x), CC_probability.y, CC_probability.x),
           ccyear = ifelse(is.na(ccyear.x), ccyear.y, ccyear.x),
           RCI = ifelse(is.na(RCI.x), RCI.y, RCI.x),
           rot_crop = ifelse(is.na(rot_crop.x), rot_crop.y, rot_crop.x),
           rot_fg = ifelse(is.na(rot_fg.x), rot_fg.y, rot_fg.x),
           ncrops = ifelse(is.na(ncrops.x), ncrops.y, ncrops.x),
           crop = ifelse(is.na(crop.x), crop.y, crop.x),
           fg = ifelse(is.na(fg.x), fg.y, fg.x)) |>
           rename(rootznaws_mean = rootznaws_mean.x,
                  soc0_100_mean = soc0_100_mean.x,
                  nccpi3all_mean = nccpi3all_mean.x) |>
    select(-ends_with(".x"), -ends_with(".y")) ->
    crop_df                        

gc()
rm(corn_df, soy_df)

library(arrow)

crop_df |>
    write_parquet("C:/Users/vf006/Documents/crop_il_long.parquet")

## ── Summary statistics by rotation type ──────────────────────────────────────
# Corn monoculture, perfect corn-soy rotation, and transitioning sequences,
# restricted to fields following one of the 64 pure corn/soybean 6-year
# sequences already enumerated in corn_soy_patterns above. crop_df already
# pools corn-year and soy-year observations (full_join on field + year), so
# corn_yield and soy_yield are each populated only on the years they apply.

library(gt)

compute_gdd <- function(tmin, tmax, base = 8, cap = 29) {
  pmax(0, pmin((tmax + tmin) / 2, cap) - base)
}

compute_edd <- function(tmin, tmax, threshold = 29) {
  pmax(0, (tmax + tmin) / 2 - threshold)
}

count_soy_years <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  sum(rotation %in% c(2L, 5L))
}

consecutive_soy <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_runs <- rle(rotation)
  lengths  <- soy_runs$lengths[soy_runs$values %in% c(2L, 5L)]
  if (is_empty(lengths)) 0L else ifelse(max(lengths) > 1, 1L, 0L)
}

soy_gap <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_years <- which(rotation %in% c(2L, 5L))
  if (length(soy_years) < 2) return(0)
  mean(diff(soy_years))
}

smallest_soy_gap <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_years <- which(rotation %in% c(2L, 5L))
  if (length(soy_years) < 2) return(0)
  min(diff(soy_years))
}

parse_rot <- function(x) {
  v <- as.integer(strsplit(x, "-")[[1]])
  ifelse(v == 5L, 2L, v)
}

# Rotation-structure features (PC1 rotation-intensity index), same
# construction as code/just_pope.r, precomputed on the 64 unique patterns.
rot_features <- corn_soy_patterns |>
  mutate(
    rot_vec     = lapply(pattern, parse_rot),
    n_soy       = sapply(rot_vec, count_soy_years),
    consec_soy  = sapply(rot_vec, consecutive_soy),
    gap_soy     = sapply(rot_vec, soy_gap),
    min_gap_soy = sapply(rot_vec, smallest_soy_gap)
  ) |>
  select(pattern, n_soy, consec_soy, gap_soy, min_gap_soy)

rot_pca_in <- rot_features |>
  transmute(
    n_soy     = n_soy,
    no_consec = 1 - consec_soy,
    freq_soy  = ifelse(gap_soy     == 0, 0, 1 / gap_soy),
    tight_soy = ifelse(min_gap_soy == 0, 0, 1 / min_gap_soy)
  )

pca_rot <- prcomp(rot_pca_in, scale. = TRUE)
rot_features <- rot_features |> mutate(rot_index = pca_rot$x[, 1])

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  left_join(rot_features, by = c("rot_crop" = "pattern")) |>
  mutate(
    EDD_6    = compute_edd(tmmx_6, tmmn_6),
    EDD_7    = compute_edd(tmmx_7, tmmn_7),
    EDD_8    = compute_edd(tmmx_8, tmmn_8),
    soy_lag1 = as.integer(substr(rot_crop, 9, 9) == "5"),   # t-1
    soy_lag2 = as.integer(substr(rot_crop, 7, 7) == "5"),   # t-2
    soy_lag3 = as.integer(substr(rot_crop, 5, 5) == "5"),   # t-3
    soy_lag4 = as.integer(substr(rot_crop, 3, 3) == "5"),   # t-4
    soy_lag5 = as.integer(substr(rot_crop, 1, 1) == "5"),   # t-5
    rotation_type = case_when(
      rot_crop == "1-1-1-1-1-1" ~ "Corn monoculture",
      rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5") ~ "Perfect corn-soy rotation",
      rot_crop == "5-5-5-5-5-5" ~ NA_character_,           # soy monoculture: out of scope
      .default = "Transitioning sequences"
    )
  ) |>
  filter(!is.na(rotation_type)) |>
  mutate(rotation_type = factor(rotation_type,
           levels = c("Corn monoculture", "Perfect corn-soy rotation", "Transitioning sequences"))) ->
  rotation_summary_df

yield_vars    <- c("corn_yield", "soy_yield")
rotation_vars <- c("RCI", "n_soy", "consec_soy", "gap_soy", "min_gap_soy", "rot_index",
                    "soy_lag1", "soy_lag2", "soy_lag3", "soy_lag4", "soy_lag5")
weather_vars  <- c("pr_6", "pr_7", "pr_8",
                    "cGDD_6m", "cGDD_7m", "cGDD_8m",
                    "EDD_6", "EDD_7", "EDD_8",
                    "vpd_6", "vpd_7", "vpd_8",
                    "soil_6", "soil_7", "soil_8")
awc_vars      <- "rootznaws_mean"
summary_vars  <- c(yield_vars, rotation_vars, weather_vars, awc_vars)

var_labels <- c(
  corn_yield  = "Corn yield (bu/acre)",
  soy_yield   = "Soybean yield (bu/acre)",
  RCI         = "Rotational Complexity Index (RCI)",
  n_soy       = "Soybean years (of 6)",
  consec_soy  = "Consecutive soy years (0/1)",
  gap_soy     = "Mean gap between soy years",
  min_gap_soy = "Minimum gap between soy years",
  rot_index   = "Rotation intensity index (PC1)",
  soy_lag1    = "Soy in t-1",
  soy_lag2    = "Soy in t-2",
  soy_lag3    = "Soy in t-3",
  soy_lag4    = "Soy in t-4",
  soy_lag5    = "Soy in t-5",
  pr_6        = "Precipitation, June",
  pr_7        = "Precipitation, July",
  pr_8        = "Precipitation, August",
  cGDD_6m     = "Cumulative GDD, June",
  cGDD_7m     = "Cumulative GDD, July",
  cGDD_8m     = "Cumulative GDD, August",
  EDD_6       = "Extreme degree days, June",
  EDD_7       = "Extreme degree days, July",
  EDD_8       = "Extreme degree days, August",
  vpd_6       = "VPD, June",
  vpd_7       = "VPD, July",
  vpd_8       = "VPD, August",
  soil_6      = "Soil moisture, June",
  soil_7      = "Soil moisture, July",
  soil_8      = "Soil moisture, August",
  rootznaws_mean = "Available water capacity (AWC), root zone"
)

var_category <- c(
  corn_yield = "Yields", soy_yield = "Yields",
  RCI = "Rotation variables", n_soy = "Rotation variables", consec_soy = "Rotation variables",
  gap_soy = "Rotation variables", min_gap_soy = "Rotation variables", rot_index = "Rotation variables",
  soy_lag1 = "Rotation variables", soy_lag2 = "Rotation variables", soy_lag3 = "Rotation variables",
  soy_lag4 = "Rotation variables", soy_lag5 = "Rotation variables",
  pr_6 = "Weather controls", pr_7 = "Weather controls", pr_8 = "Weather controls",
  cGDD_6m = "Weather controls", cGDD_7m = "Weather controls", cGDD_8m = "Weather controls",
  EDD_6 = "Weather controls", EDD_7 = "Weather controls", EDD_8 = "Weather controls",
  vpd_6 = "Weather controls", vpd_7 = "Weather controls", vpd_8 = "Weather controls",
  soil_6 = "Weather controls", soil_7 = "Weather controls", soil_8 = "Weather controls",
  rootznaws_mean = "Soil characteristics"
)

rotation_summary_df |>
  select(rotation_type, all_of(summary_vars)) |>
  pivot_longer(-rotation_type, names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(rotation_type, variable) |>
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") |>
  mutate(
    category = var_category[variable],
    label    = var_labels[variable],
    cell     = sprintf("%.2f (%.2f)", mean, sd)
  ) |>
  select(category, variable, label, rotation_type, cell) |>
  pivot_wider(names_from = rotation_type, values_from = cell) ->
  summary_wide

rotation_summary_df |>
  count(rotation_type) |>
  pivot_wider(names_from = rotation_type, values_from = n) |>
  mutate(across(everything(), ~ formatC(as.numeric(.x), big.mark = ",", format = "d")),
         category = "Sample size",
         label    = "Observations (field-years)") ->
  n_row

bind_rows(summary_wide, n_row) |>
  mutate(category = factor(category, levels = c("Yields", "Rotation variables",
                                                  "Weather controls", "Soil characteristics",
                                                  "Sample size"))) |>
  arrange(category) |>
  select(category, label,
         `Corn monoculture`, `Perfect corn-soy rotation`, `Transitioning sequences`) ->
  summary_stats_by_rotation

summary_stats_by_rotation |>
  gt(groupname_col = "category", rowname_col = "label") |>
  tab_header(
    title    = "Summary statistics by rotation type",
    subtitle = "Mean (SD); corn- and soybean-year observations, restricted to 6-year corn/soybean sequences"
  ) |>
  tab_source_note(
    "Corn monoculture: 1-1-1-1-1-1. Perfect corn-soy rotation: 5-1-5-1-5-1 or 1-5-1-5-1-5.
     Transitioning sequences: all other corn/soybean 6-year sequences (excludes soybean monoculture)."
  ) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) ->
  summary_stats_gt

gtsave(summary_stats_gt, "C:/Users/vf006/Box/crop_rotations_and_losses/tables/summary_stats_by_rotation.html")
gtsave(summary_stats_gt, "C:/Users/vf006/Box/crop_rotations_and_losses/tables/summary_stats_by_rotation.png")