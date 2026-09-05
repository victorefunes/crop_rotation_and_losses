## ============================================================================
## corn_data_prep.R
## Shared data processing for the corn Just-Pope analysis.
## Loads raw field-level data, builds crop-lag / rotation / RCI / weather
## variables, assembles the analysis sample (corn_jp_data), and produces the
## summary-statistics table. Source this file at the top of:
##   corn_rotation_analysis.R  — rotation-sequence analyses
##   corn_rci_analysis.R       — RCI-based analyses
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Produces:
##   corn_jp_data   — analysis sample (shared)
##   vpd_controls   — weather control vector used by VPD-interaction models
##   load_corn_sf_2016() — helper to (re)build 2016 spatial data for maps
## Table produced here:
##   corn_summary_stats.tex — Summary statistics by rotation type
## ============================================================================

library(arrow)
library(statar)
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
  # rot_seq, or anything else that happens to start with "crop_". This is a
  # pure year-lookup: it works the same whether crop_YYYY holds raw CDL names
  # ("Corn"), numeric CDL codes, or functional-group labels ("corn") from
  # recode_cdl_functional() -- it never inspects the values, only picks the
  # column matching dt$year. Forced to character so the lookup is well-defined
  # regardless of which of those three a given caller has already applied.
  crop_cols  <- grep("^crop_[0-9]{4}$", names(dt), value = TRUE)
  crop_years <- as.integer(sub("crop_", "", crop_cols))

  stopifnot(length(crop_cols) > 0, !anyNA(crop_years))

  m       <- as.matrix(dt[, lapply(.SD, as.character), .SDcols = crop_cols])
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

source("cdl_functional_recode.R")

recode_cdl_functional(corn_df, cols = paste0("crop_", 0:6))

corn_df <- corn_df |>
  #rename(RCI = annual_RCI) |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-",
          crop_2, "-", crop_1, "-", crop_0)) |>
  #rci_correction() |>
  add_degree_days()

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.
#
# Rotation universe: any COMPLETE 6-year functional-group history (crop_0..
# crop_5 all agricultural, i.e. non-NA after recode_cdl_functional()), not
# just the 29 hand-curated corn/soy/wheat sequences corn_soy_patterns still
# encodes -- lasso_rotation_selection.R now does the dozens-to-a-handful
# selection on this broader candidate set instead of a pre-filtered pattern
# list. MIN_PATTERN_FREQ drops sequences so rare they'd be perfectly
# collinear with field FE and add nothing but noise to the LASSO design
# matrix; raise/lower it there, not by hand-editing a pattern list.
# (corn_soy_patterns / rot_features are still used below, unchanged, to
# attach rot_index -- that PCA is fit only on the original 29 sequences, so
# rows outside that set simply get rot_index = NA, which is fine: only the
# fml_corn_index model actually uses rot_index, and it drops NA rows itself.)
MIN_PATTERN_FREQ <- 30
pattern_freq     <- corn_df |> count(rot_crop, name = "N")
common_patterns  <- pattern_freq$rot_crop[pattern_freq$N >= MIN_PATTERN_FREQ]
cat(sprintf("Rotation patterns: %d distinct, %d kept at N >= %d field-years.\n",
            nrow(pattern_freq), length(common_patterns), MIN_PATTERN_FREQ))

corn_jp_data <- corn_df |>
  filter(if_all(paste0("crop_", 0:5), ~ !is.na(.))) |>
  filter(rot_crop %in% common_patterns) |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  mutate(
    # Lag dummies via delimiter split — robust to multi-word groups
    # (e.g. "annual_cereal"). Tokens run oldest -> newest; t[6] is the current
    # year (always corn), so lag k is token 6 - k. Soy = "soybeans" on the
    # functional-group string, after recode_cdl_functional().
    seq_num  = strsplit(rot_crop, "-", fixed = TRUE),
    soy_lag1 = sapply(seq_num, \(v) as.integer(v[5] == "soybeans")),
    soy_lag2 = sapply(seq_num, \(v) as.integer(v[4] == "soybeans")),
    soy_lag3 = sapply(seq_num, \(v) as.integer(v[3] == "soybeans")),
    soy_lag4 = sapply(seq_num, \(v) as.integer(v[2] == "soybeans")),
    soy_lag5 = sapply(seq_num, \(v) as.integer(v[1] == "soybeans")),
    rot_crop = to_letters(rot_crop),   # corn->C, soybeans->S, annual_cereal->W (helper from rotation_setup.R)
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

# Free raw data — no longer needed
rm(corn_df);
gc()

# ── Summary statistics table ──────────────────────────────────────────────────
# Table: tab:summary — means and SDs by rotation type
# Produced here because corn_jp_data contains all needed variables, and this
# table blends both rotation-type and RCI groupings.
# Written to tables/summary_stats.tex for \input{} in the manuscript.

corn_jp_data |>
  mutate(rotation_type = case_when(
    as.character(rot_crop) == "C-C-C-C-C-C"               ~ "Corn monoculture",
    as.character(rot_crop) %in% c("S-C-S-C-S-C",
                                   "C-S-C-S-C-S")         ~ "Perfect rotation",
    RCI >= 1.41 & RCI <= 2.00                              ~ "Transitioning",
    TRUE                                                   ~ "Other"
  )) |>
  filter(rotation_type != "Other") |>
  group_by(rotation_type) |>
  summarise(
    N               = n(),
    `Corn yield`    = sprintf("%.1f (%.1f)", mean(corn_yield, na.rm=TRUE),
                                              sd(corn_yield,   na.rm=TRUE)),
    `RCI`           = sprintf("%.2f (%.2f)", mean(RCI,         na.rm=TRUE),
                                              sd(RCI,          na.rm=TRUE)),
    `Precip (Jun)`  = sprintf("%.1f (%.1f)", mean(pr_6,        na.rm=TRUE),
                                              sd(pr_6,         na.rm=TRUE)),
    `GDD (Jul)`     = sprintf("%.0f (%.0f)", mean(cGDD_7m,     na.rm=TRUE),
                                              sd(cGDD_7m,      na.rm=TRUE)),
    `EDD (Jul)`     = sprintf("%.1f (%.1f)", mean(EDD_7,       na.rm=TRUE),
                                              sd(EDD_7,        na.rm=TRUE)),
    `VPD (Jul)`     = sprintf("%.2f (%.2f)", mean(vpd_7,       na.rm=TRUE),
                                              sd(vpd_7,        na.rm=TRUE)),
    `AWC (mm)`      = sprintf("%.0f (%.0f)", mean(rootznaws_mean, na.rm=TRUE),
                                              sd(rootznaws_mean,  na.rm=TRUE))
  ) |>
  mutate(rotation_type = factor(rotation_type,
                                 levels = c("Corn monoculture",
                                            "Perfect rotation",
                                            "Transitioning"))) |>
  arrange(rotation_type) |>
  rename(`Rotation type` = rotation_type) |>
  kable(format  = "latex",
        booktabs = TRUE,
        caption = "Summary statistics by rotation type",
        label   = "corn_summary",
        align   = c("l","r","r","r","r","r","r","r","r")) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  save_kable(file = paste0(tab_dir, "corn_summary_stats.tex"))

cat("Summary statistics table saved.\n")

# ── Final sample used by both downstream analyses ─────────────────────────────
# Drop the few rows with missing/out-of-range VPD (excluded from every
# VPD-interaction and Just-Pope model, so it's cleanest to drop once, here).

corn_jp_data <- corn_jp_data |>
  filter(!is.na(vpd_name))

# ── Shared constants ──────────────────────────────────────────────────────────

# VPD clashes with vpd_name, so it is removed from controls.
vpd_controls <- c("pr_6", "pr_7", "pr_8",
  "I(pr_6^2)", "I(pr_7^2)", "I(pr_8^2)",
  "cGDD_6m", "cGDD_7m", "cGDD_8m",
  "EDD_6", "EDD_7", "EDD_8",
  "soil_6", "soil_7", "soil_8")

# ── Spatial data helper (2016 cross-section) ──────────────────────────────────
# Reloads a minimal version of corn_df for spatial plotting only. Not run
# automatically — call from whichever analysis script needs a map.
# NOTE (pre-existing, carried over unchanged from corn_analysis_full.R):
# this reload path calls rci_correction() and expects an `annual_RCI` column,
# neither of which is defined/produced by the main pipeline above — this
# function will error unless rci_correction() has been defined elsewhere in
# the session (e.g. by sourcing main_analysis.R or rotation_setup.R).

#load_corn_sf_2016 <- function() {
#  cat("Reloading corn data for spatial maps...\n")
#  corn_df <- read_parquet(
#    "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet")

#  corn_df <- corn_df |>
#    filter(STATE_ABBR == "IL" & year == 2016) |>
#    mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
#           corn_yield = corn_yield / 62.77)  |>
#    arrange(tile_field_ID, year)

  ## Add present year crop variable
#  setDT(corn_df)

#  add_crop_year(corn_df)

#  corn_df <- corn_df |>
#    rename(crop_0 = crop_year,
#           crop_1 = prioryr_crop,
#           crop_2 = prior2yr_crop,
#           crop_3 = prior3yr_crop,
#           crop_4 = prior4yr_crop,
#           crop_5 = prior5yr_crop,
#           crop_6 = prior6yr_crop)

#  source("cdl_functional_recode.R")

#  recode_cdl_functional(corn_df, cols = paste0("crop_", 0:6))

#  corn_df <- corn_df |>
#    rename(RCI = annual_RCI) |>
#    mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-",
#            crop_2, "-", crop_1, "-", crop_0)) |>
#    rci_correction()

#  corn_sf <- corn_df |>
#    filter(rot_crop %in% corn_soy_patterns$pattern) |>
#    arrange(tile_field_ID, year) |>
#    st_as_sf(coords = c("lon", "lat"), crs = st_crs("EPSG:4326"))
#  rm(corn_df); gc()

#  il_map <- us_map(regions = "counties") |>
#    filter(abbr == "IL") |>
#    st_as_sf() |>
#    st_transform(st_crs("EPSG:4326"))

#  list(corn_sf = corn_sf, il_map = il_map)
#}

cat("corn_data_prep.R: shared setup complete.\n")
