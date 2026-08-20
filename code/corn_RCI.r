## ============================================================================
## corn_RCI.R
## Just-Pope production risk analysis — CORN
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Tables produced:
##   tab:corn_rot       — Rotation patterns and corn yields
##   tab:corn_rci       — Rotational Complexity and corn yields
##   tab:corn_rot_vpd   — Rotation x VPD interaction, corn
##   tab:corn_rci_vpd   — RCI x VPD interaction, corn
##   tab:corn_jp_mean   — Just-Pope stage 1: corn yield mean
##   tab:corn_jp_var    — Just-Pope stage 2: corn yield variance (OLS)
##   tab:corn_rci_jp    — Just-Pope factor RCI: corn yield moments
##
## Figures produced:
##   rot_pca_plot       — PCA biplot of rotation features (shared)
##   corn_rot_plot      — Response of corn yields to rotation sequences
##   corn_rci_plot      — Response of corn yields to RCI
##   corn_var_plot      — Response of std dev of corn yields to rotation sequences
##   corn_coeff_plot    — Mean vs variance coefficients scatter (corn)
##   corn_jp_plot       — Just-Pope mean-variance decomposition (corn)
##   rci_plot           — Nonlinear effect of RCI on corn yield mean and variance
##   corn_yield_map     — Spatial map of corn yields (2016)
##   rci_map            — Spatial map of RCI values (2016)
##   nccpi_corn_map     — Spatial map of NCCPI corn (2016)
## ============================================================================
library(arrow)
library(tidyverse)
library(statar)
library(fixest)
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

source("rotation_setup_wa.R")

# ── Load corn data ────────────────────────────────────────────────────────────
# Load, correct, and add degree-day variables in one pipeline.
# corn_df is kept as the raw data source; corn_jp_data is the analysis sample.

cat("Loading corn data...\n")
corn_df <- read_parquet(
  "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet") 

corn_df <- corn_df |>
  filter(STATE_ABBR == "IL") |> 
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         corn_yield = corn_yield / 62.77)  |> 
  arrange(tile_field_ID, year) 

cat("Corn raw rows:", nrow(corn_df), "\n")

## Add present year crop variable
setDT(corn_df)  

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

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.
 
corn_df <- corn_df |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))
 
cat("Corn analysis sample:", nrow(corn_df), "rows\n")

xcols <- paste0("crop_", 0:5)

corn_code     <- 1L
soy_code      <- 5L
wheat_codes   <- c(22L, 23L, 24L, 26L, 225L, 236L, 238L)   # MATCH your has_wheat definition
alfalfa_codes <- 36L
core_codes    <- c(corn_code, soy_code, wheat_codes, alfalfa_codes)

# non-crop land covers that exist in corn_df but NOT in corn_jp_data — must be excluded
nonag_codes <- c(0L, 63L, 64L, 65L, 81L, 82L, 83L, 87L, 88L, 111L, 112L,
                 121L, 122L, 123L, 124L, 131L, 141L, 142L, 143L, 152L, 190L, 195L)

all_codes        <- sort(unique(unlist(corn_df[, ..xcols], use.names = FALSE)))
other_crop_codes <- setdiff(all_codes, c(core_codes, nonag_codes))

corn_df[, has_other := Reduce(`|`, lapply(.SD, \(v) v %in% other_crop_codes)),
        .SDcols = xcols]

corn_df[, has_wheat := Reduce(`|`, lapply(.SD, \(v) v %in% wheat_codes)),
        .SDcols = xcols]

corn_df[, has_alfalfa := Reduce(`|`, lapply(.SD, \(v) v %in% alfalfa_codes)),
        .SDcols = xcols]

corn_df[, regime := fcase(
  has_other,                  "other",                   # checked first — wins
  !has_wheat & !has_alfalfa,  "corn_soy",
   has_wheat & !has_alfalfa,  "corn_soy_wheat",
  !has_wheat &  has_alfalfa,  "corn_soy_alfalfa",
   has_wheat &  has_alfalfa,  "corn_soy_wheat_alfalfa"
)]
corn_df[, regime := factor(regime, levels = c(
  "corn_soy","corn_soy_wheat","corn_soy_alfalfa","corn_soy_wheat_alfalfa","other"))]

corn_df[, .N, by = regime][order(-N)]
 
corn_df <- corn_df |>
  filter(!is.na(regime))

# Free raw data — no longer needed
gc() 


# Switching rgression with factorial interactions of RCI, has_wheat, and has_alfalfa
switch_reg_factorial <- feols(
  corn_yield ~ RCI * has_wheat * has_alfalfa +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~COUNTY_FIPS
)
etable(switch_reg_factorial, keep = "RCI")

avg_slopes(switch_reg_factorial, variables = "RCI",
          newdata = datagrid(has_wheat = c(FALSE, TRUE), 
            has_alfalfa = c(FALSE, TRUE)),
          by = c("has_wheat", "has_alfalfa"))


reg_factorial <- feols(
  corn_yield ~ RCI * regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~tile_field_ID + year
)
etable(reg_factorial, keep = "RCI")

library(splines)
reg_curve <- feols(
  corn_yield ~ ns(RCI, df = 3) * regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~COUNTY_FIPS
)
summary(reg_curve)

## Problem: some regimes have very few observations
corn_df |> 
  filter(!is.na(RCI) | !is.na(regime)) |>
  tab(regime, RCI) |>
  data.frame() |>
  pivot_wider(id_cols = "RCI", names_from = "regime", values_from = "Freq.")

# Regime-based bins
corn_df[, RCI_bin := fcase(
  RCI == 0,                 "monoculture",
  RCI > 0 & RCI < 2.24,     "low",
  RCI == 2.24,              "perfect_rotation",
  RCI > 2.24 & RCI <= 3.24, "moderate",
  RCI > 3.24,               "high"
)]
corn_df[, RCI_bin := factor(RCI_bin,
  levels = c("monoculture", "low", "perfect_rotation", "moderate", "high"))]

table(corn_df$regime, corn_df$RCI_bin)

switch_reg_factor <- feols(
  corn_yield ~ RCI_bin +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8
    | tile_field_ID + year,
  data = corn_df,
  split = ~regime,
  cluster = ~COUNTY_FIPS
)

etable(switch_reg_factor, keep = "RCI")
gc()


## Stigler's model
setorder(corn_df, tile_field_ID, year)

# 1) per-observation trigger: define the ONE event you want to anchor on.
#    simplest = "field is no longer plain corn_soy this year"
corn_df[, treated_now := regime != "corn_soy"]          # or: has_wheat | has_alfalfa | has_other

# 2) first year the trigger fires, per field; Inf if it never does
corn_df[, first_treat := {yrs <- year[treated_now]
  if (length(yrs)) min(yrs) else Inf
}, by = tile_field_ID]

est_es <- feols(
  corn_yield ~ sunab(first_treat, year) +        # Sun-Abraham: clean leads/lags
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~COUNTY_FIPS)

iplot(est_es)                                     # pre-event coefficients ARE the placebo
wald(est_es, "year::-[2-9]")                      # joint test that leads = 0