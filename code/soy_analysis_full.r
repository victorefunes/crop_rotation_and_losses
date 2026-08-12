## ============================================================================
## soy_analysis.R
## Just-Pope production risk analysis — SOYBEANS
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Tables produced:
##   tab:soy_rot       — Rotation patterns and soy yields
##   tab:soy_rci       — Rotational Complexity and soy yields
##   tab:soy_rot_vpd   — Rotation x VPD interaction, soy
##   tab:soy_rci_vpd   — RCI x VPD interaction, soy
##   tab:soy_jp_mean   — Just-Pope stage 1: soy yield mean
##   tab:soy_jp_var    — Just-Pope stage 2: soy yield variance (OLS)
##   tab:soy_rci_jp    — Just-Pope factor RCI: soy yield moments
##
## Figures produced:
##   rot_pca_plot       — PCA biplot of rotation features (shared)
##   soy_rot_plot      — Response of soy yields to rotation sequences
##   soy_rci_plot      — Response of soy yields to RCI
##   soy_var_plot      — Response of std dev of soy yields to rotation sequences
##   soy_coeff_plot    — Mean vs variance coefficients scatter (soy)
##   soy_jp_plot       — Just-Pope mean-variance decomposition (soy)
##   rci_plot           — Nonlinear effect of RCI on soy yield mean and variance
##   soy_yield_map     — Spatial map of soy yields (2016)
##   rci_map            — Spatial map of RCI values (2016)
##   nccpi_soy_map     — Spatial map of NCCPI soy (2016)
## ============================================================================
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
         crop_6 = prior6yr_crop) 

source("cdl_recode.R")

recode_cdl(soy_df, cols = paste0("crop_", 0:6))

soy_df <- soy_df |>
  rename(RCI = annual_RCI) |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-", 
          crop_2, "-", crop_1, "-", crop_0)) |>
  rci_correction() |>
  add_degree_days()

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.
 
soy_jp_data <- soy_df |>
  filter(rot_crop %in% soy_patterns$pattern) |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  mutate(
    # Lag dummies via delimiter split — robust to multi-digit codes (24 = wheat).
    # Tokens run oldest -> newest; t[6] is the current year (always corn), so
    # lag k is token 6 - k. Soy = "5" on the numeric string, before recoding.
    seq_num  = strsplit(rot_crop, "-", fixed = TRUE),
    soy_lag1 = sapply(seq_num, \(v) as.integer(v[5] == "5")),
    soy_lag2 = sapply(seq_num, \(v) as.integer(v[4] == "5")),
    soy_lag3 = sapply(seq_num, \(v) as.integer(v[3] == "5")),
    soy_lag4 = sapply(seq_num, \(v) as.integer(v[2] == "5")),
    soy_lag5 = sapply(seq_num, \(v) as.integer(v[1] == "5")),
    rot_crop = to_letters(rot_crop),   # 1->C, 5->S, 24->W (helper from rotation_setup.R)
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")
  ) |>
  select(-seq_num) |>
  filter(!is.na(soy_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))
 
cat("Soy analysis sample:", nrow(soy_jp_data), "rows\n")
 
# Free raw data — no longer needed
rm(soy_df); gc()
 
# ── Summary statistics table ──────────────────────────────────────────────────
# Table: tab:summary — means and SDs by rotation type
# Produced here because soy_jp_data contains all needed variables.
# Written to tables/summary_stats.tex for \input{} in the manuscript.
 
soy_jp_data |>
  mutate(rotation_type = case_when(
    as.character(rot_crop) == "S-S-S-S-S-S"               ~ "Soy monoculture",
    as.character(rot_crop) %in% c("S-C-S-C-S-C",
                                   "C-S-C-S-C-S")         ~ "Perfect rotation",
    RCI >= 1.41 & RCI <= 2.00                              ~ "Transitioning",
    TRUE                                                   ~ "Other"
  )) |>
  filter(rotation_type != "Other") |>
  group_by(rotation_type) |>
  summarise(
    N               = n(),
    `Soy yield`     = sprintf("%.1f (%.1f)", mean(soy_yield, na.rm=TRUE),
                                              sd(soy_yield,   na.rm=TRUE)),
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
                                 levels = c("Soy monoculture",
                                            "Perfect rotation",
                                            "Transitioning"))) |>
  arrange(rotation_type) |>
  rename(`Rotation type` = rotation_type) |>
  kable(format  = "latex",
        booktabs = TRUE,
        caption = "Summary statistics by rotation type",
        label   = "summary",
        align   = c("l","r","r","r","r","r","r","r","r")) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  #footnote(general = paste0(
  #  "Means with standard deviations in parentheses. ",
  #  "Corn yield in bushels per acre (QDANN). ",
  #  "Precipitation in mm; GDD and EDD in degree-days; ",
  #  "VPD in kPa; AWC in mm. ",
  #  "Perfect rotation = alternating C-S sequence (S-C-S-C-S-C or C-S-C-S-C-S). ",
  #  "Transitioning = sequences with RCI between 1.41 and 2.00. ",
  #  "Sample: Illinois field-year observations, 2009--2022."
  #),
  #general_title = "Notes:") |>
  save_kable(file = paste0(tab_dir, "summary_stats.tex"))
 
cat("Summary statistics table saved.\n")
 
# ── 1. OLS mean model ─────────────────────────────────────────────────────────
# Table: tab:soy_rot | Figure: soy_rot_plot
 
soy_yield_formula <- make_jp_formula("soy_yield", "rot_crop", all_controls)
 
feols(soy_yield ~ rot_crop | tile_field_ID + year,
      data = soy_jp_data, cluster = ~COUNTY_FIPS) -> soy_rot_nc
 
feols(soy_yield_formula,
      data = soy_jp_data, cluster = ~COUNTY_FIPS) -> soy_rot
 
etable(soy_rot_nc, soy_rot,
       tex      = TRUE,
       dict     = dict_soy,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       fontsize = "scriptsize",
       arraystretch = 0.8,
       title    = "Rotation patterns and soy yields",
       label    = "tab:soy_rot",
       extralines = list("_Controls" = c("No", "Yes")),
       file     = paste0(tab_dir, "soy_rot.tex"))
 
# BH multiple comparisons correction
pvals_soy <- broom::tidy(soy_rot) |>
  filter(grepl("rot_crop", term)) |>
  arrange(p.value) |>
  mutate(p_adj_bh = p.adjust(p.value, method = "BH"),
         sig_raw  = p.value  < 0.05,
         sig_bh   = p_adj_bh < 0.05)
 
cat("Soy sequences significant at 5% (unadjusted):", sum(pvals_soy$sig_raw), "\n")
cat("Soy sequences significant at 5% FDR (BH):    ", sum(pvals_soy$sig_bh),  "\n")
 
# Figure: soy_rot_plot — Response of soy yields to rotation sequences
soy_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    .default              = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of soy yields to rotation sequences",
       caption = "Reference: S-S-S-S-S-S. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rot_plot
ggsave(paste0(fig_dir, "soy_rot_plot.png"), soy_rot_plot,
       width = 10, height = 7.5, dpi = 300)

library(tidytext)   # for reorder_within / scale_y_reordered

soy_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(
    system = if_else(grepl("W", term), "Corn-soy-wheat", "Corn-soy"),
    group  = case_when(
      term == "C-S-C-S-C-S" ~ "Perfect rotation",
      .default              = "Other")) |>
  ggplot(aes(x = reorder_within(term, -prms.y, system), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_x_reordered() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of soy yields to rotation sequences",
       caption = "Reference: S-S-S-S-S-S. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rot_plot
ggsave(paste0(fig_dir, "soy_rot_plot.png"), soy_rot_plot,
       width = 10, height = 7.5, dpi = 300)

 
# ── 2. RCI models — soy ──────────────────────────────────────────────────────
# Table: tab:soy_rci | Figure: soy_rci_plot
 
soy_rci_nc <- make_jp_formula("soy_yield", "RCI", NULL)
soy_rci_cs <- make_jp_formula("soy_yield", "RCI", all_controls)
 
# Remove infrequent RCI levels (fewer than 100 observations)
rci_keep <- soy_jp_data |>
  count(RCI) |>
  filter(n >= 100) |>
  pull(RCI)

soy_jp_data |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_nc, data = _, cluster = ~COUNTY_FIPS) -> soy_rci_all

soy_jp_data |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_cs, data = _, cluster = ~COUNTY_FIPS) -> soy_rci_cs
 
etable(soy_rci_all, soy_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("No controls", "Weather and soil controls"),
       keep     = "RCI",
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       title    = "Rotation Complexity and soy yields",
       label    = "tab:soy_rci",
       file     = paste0(tab_dir, "soy_rci.tex"))
 
# Figure: soy_rci_plot — Changes in RCI and soy yields
soy_rci_cs |>
  coefplot() |>
  data.frame() |>
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  mutate(term  = as.numeric(term),
         group = case_when(
           term == 2.24 ~ "Perfect rotation",
           .default     = "Other")) |>
  filter(term < 5.2) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = "Rotational Complexity Index (RCI)", y = "Coefficient Estimate",
       title   = "Changes in RCI and soy yields",
       caption = "Soy fields. Reference: RCI = 0. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rci_plot
ggsave(paste0(fig_dir, "soy_rci_plot.png"), soy_rci_plot,
       width = 10, height = 7.5, dpi = 300)
 
rm(soy_rci_all, soy_rci_cs); gc()
 
# ── 3. VPD interaction models — soy ─────────────────────────────────────────
# Tables: tab:soy_rot_vpd, tab:soy_rci_vpd
 
soy_vpd_formula     <- make_jp_formula("soy_yield", "rot_crop + vpd_name",
                                        all_controls)
soy_rci_vpd_formula <- make_jp_formula("soy_yield", "RCI * vpd_name",
                                        all_controls)
 
feols(soy_yield ~ rot_crop + vpd_name | tile_field_ID + year,
      data = soy_jp_data, cluster = ~COUNTY_FIPS) -> soy_rot_vpd_nc                                         
 
soy_jp_data |>
  feols(soy_vpd_formula, data = _, cluster = ~COUNTY_FIPS) -> soy_rot_vpd
 
soy_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) -> soy_rci_vpd
 
etable(soy_rot_vpd,
       tex      = TRUE,
       dict     = c(dict_soy, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Effect of weather and rotation sequences on soy yields",
       label    = "tab:soy_rot_vpd",
       file     = paste0(tab_dir, "soy_rot_vpd.tex"))
 
etable(soy_rci_vpd,
       tex      = TRUE,
       dict     = dict_vpd,
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "RCI x drought interaction effects on soy yields",
       label    = "tab:soy_rci_vpd",
       file     = paste0(tab_dir, "soy_rci_vpd.tex"))     
 
rm(soy_rot_vpd, soy_rci_vpd); 
gc()
 
# ── 4. Just-Pope stage 1 — soy ───────────────────────────────────────────────
# Table: tab:soy_jp_mean | Figures: soy_rot_plot (mean), soy_jp_plot
 
# Lag comparison models (confirmatory Z-vector pre-registration)
fml_soy_lag1      <- make_jp_formula("soy_yield", "soy_lag1", all_controls_fgls)
fml_soy_lag1_lag2 <- make_jp_formula("soy_yield", "soy_lag1 + soy_lag2",
                                       all_controls_fgls)
fml_soy_index     <- make_jp_formula("soy_yield", "rot_index", all_controls_fgls)
fml_soy_mean      <- make_jp_formula("soy_yield", "rot_crop", all_controls_fgls)
 
feols(fml_soy_lag1,      data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s1_lag1
feols(fml_soy_lag1_lag2, data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s1_lag2
feols(fml_soy_index,     data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s1_idx
feols(fml_soy_mean,      data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s1
 
bh_note_soy <- paste0(
  "Benjamini-Hochberg FDR correction across ", nrow(pvals_soy),
  " rotation-sequence coefficients: ",
  sum(pvals_soy$sig_raw), " significant at 5\\% (unadjusted); ",
  sum(pvals_soy$sig_bh),  " significant at 5\\% FDR."
)
 
rot_dict <- setNames(
  chartr("15", "CS", sub("^rot_crop", "", grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE))),
  grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE)
)
 
# Table: tab:soy_jp_mean
etable(soy_jp_s1, soy_jp_s1_lag1, soy_jp_s1_lag2, soy_jp_s1_idx,
       tex      = TRUE,
       keep     = c("^[CSW]-", "soy_lag", "rot_index"),
       dict     = rot_dict,
       notes    = bh_note_soy,
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Soy yield: full sequences vs lag summary variables",
       label    = "tab:soy_jp_mean",
       file     = paste0(tab_dir, "soy_jp_mean.tex"))
 
rm(soy_jp_s1_lag1, soy_jp_s1_lag2, soy_jp_s1_idx, soy_rot_nc, soy_rot)
gc()
 
# ── Z-vector model (confirmatory spec) ───────────────────────────────────────
# Table: tab:zvector — Effect of rotation patterns on yields
# This is the paper's main result table (Table 1 in the draft PDF).
# Four structural features: late_soy, soy_gap, soy_cons, nsoy
# Run on both corn and soy data, both mean and variance stages.
# NOTE: requires late_soy, soy_gap, soy_cons, nsoy to be in soy_jp_data.
# These must be constructed before this chunk runs.
 
# Construct Z-vector variables if not already present
# Current crop is always soy (S) here — every sequence in soy_patterns ends
# in S — so a "late_soy" defined the corn-script way would be constant (-1)
# for every row and get dropped for collinearity. Track CORN recency
# instead, since corn is the crop that actually varies across these
# soy-terminal sequences.
# late_corn: negative integer = how many periods ago was the last corn year
soy_jp_data <- soy_jp_data |>
  mutate(
    # Parse rotation sequence to find last corn year
    seq_vec   = strsplit(as.character(rot_crop), "-"),
    late_corn = sapply(seq_vec, function(v) {
      corn_pos <- which(rev(v) == "C")   # positions from most recent (1=t-1)
      if (length(corn_pos) == 0) 0L else -min(corn_pos)
    }),
    corn_cons = sapply(seq_vec, function(v) {
      runs <- rle(v)
      as.integer(any(runs$lengths[runs$values == "C"] >= 2))
    }),
    corn_gap  = sapply(seq_vec, function(v) {
      pos <- which(v == "C")
      if (length(pos) < 2) 0L else min(diff(pos))
    }),
    ncorn     = sapply(seq_vec, function(v) sum(v == "C"))
  ) |>
  select(-seq_vec)
 
# Z-vector stage 1 — soy mean
fml_z_soy_mean <- make_jp_formula("soy_yield",
                                   "late_corn + corn_gap + corn_cons + ncorn",
                                   all_controls_fgls)

feols(fml_z_soy_mean, data = soy_jp_data,
      cluster = ~tile_field_ID + year) -> soy_z_s1

# Z-vector stage 2 — soy variance
soy_jp_data <- soy_z_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid_sq_z = (soy_yield - .fitted)^2) |>
  select(-starts_with("."))

fml_z_soy_var <- make_jp_formula("resid_sq_z",
                                   "late_corn + corn_gap + corn_cons + ncorn",
                                   all_controls_fgls)
 
feols(fml_z_soy_var, data = soy_jp_data,
      cluster = ~tile_field_ID + year) -> soy_z_s2
 
# ── Save Z-vector models for tables_combined.R ────────────────────────────────
#saveRDS(
#  list(z_s1        = soy_z_s1,
#       z_s2        = soy_z_s2,
#       rot_vpd_nc  = soy_rot_vpd_nc,
#       rot_vpd     = soy_rot_vpd),
#  file     = "D:/Crop data/soy_z_models.rds",
#  compress = "bzip2"
#)  
 
# ── 5. Just-Pope stage 2 — soy ───────────────────────────────────────────────
# Table: tab:soy_jp_var | Figures: soy_var_plot, soy_coeff_plot, soy_jp_plot
 
soy_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_jp_data
 
fml_soy_var <- make_jp_formula("resid_sq", "rot_crop", all_controls_fgls)
feols(fml_soy_var, data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s2
 
# Table: tab:soy_jp_var
etable(soy_jp_s2,
       tex      = TRUE,
       keep     = "^[CSW]-",
       dict     = rot_dict,
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2 — Soy yield conditional variance (OLS)",
       label    = "tab:soy_jp_var",
       file     = paste0(tab_dir, "soy_jp_var.tex"))
 
# ── Figure: score_yield — Response of yields to rotation score values ─────────
# soy_jp_s1 / soy_jp_s2 contain the sequence-level coefficients.
# soy_z_s1 supplies the Z-vector coefficients used to compute the score.
 
# Step 1: compute score for each unique sequence from Z-vector coefficients
z_coefs <- coef(soy_z_s1)

# feols silently drops any term that's collinear with the others (common
# here since late_corn/corn_gap/corn_cons/ncorn are derived from the same
# small set of ~41 sequences). z_coefs["name"] returns NA for a dropped
# term, and NA * anything poisons the whole score sum -- so fall back to 0
# for any missing term instead, and warn which one(s) got dropped.
get_z_coef <- function(name) {
  if (name %in% names(z_coefs)) return(unname(z_coefs[name]))
  warning("Z-vector term '", name, "' was dropped from soy_z_s1 ",
          "(collinear with the others) -- treating its score contribution as 0.")
  0
}

score_df <- soy_jp_data |>
  group_by(rot_crop) |>
  summarise(
    late_corn = mean(late_corn, na.rm = TRUE),
    corn_gap  = mean(corn_gap,  na.rm = TRUE),
    corn_cons = mean(corn_cons, na.rm = TRUE),
    ncorn     = mean(ncorn,     na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    score = get_z_coef("late_corn") * late_corn +
            get_z_coef("corn_gap")  * corn_gap  +
            get_z_coef("corn_cons") * corn_cons +
            get_z_coef("ncorn")     * ncorn
  )
 
# Step 2: extract sequence-level coefficients from the FULL sequence models
soy_s1_coef <- broom::tidy(soy_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    mean_est = estimate,
    mean_se  = std.error
  )
 
soy_s2_coef <- broom::tidy(soy_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    var_est  = estimate,
    var_se   = std.error
  )
 
# Step 3: join score onto sequence-level coefficients
score_plot_df <- score_df |>
  left_join(soy_s1_coef, by = "rot_crop") |>
  left_join(soy_s2_coef, by = "rot_crop") |>
  filter(!is.na(mean_est)) |>
  arrange(score) |>
  mutate(score_label = factor(round(score, 2),
                               levels = unique(round(score, 2))))
 
# Step 4: two-panel plot sorted by score
p_mean <- ggplot(score_plot_df,
                 aes(x = mean_est,
                     y = reorder(score_label, score))) +
  geom_point(colour = "#e05c5c", size = 2) +
  geom_errorbar(aes(xmin = mean_est - 1.96 * mean_se,
                    xmax = mean_est + 1.96 * mean_se),
                orientation = "y", width = 0, colour = "#e05c5c") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "Coefficient Estimate", y = "Rotation score", title = "Mean") +
  theme_bw()
 
p_var <- ggplot(score_plot_df,
                aes(x = var_est,
                    y = reorder(score_label, score))) +
  geom_point(colour = "#5c9ee0", size = 2) +
  geom_errorbar(aes(xmin = var_est - 1.96 * var_se,
                    xmax = var_est + 1.96 * var_se),
                orientation = "y", width = 0, colour = "#5c9ee0") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "Coefficient Estimate", y = NULL, title = "Variance") +
  theme_bw()
 
score_yield <- p_mean + p_var +
  patchwork::plot_annotation(
    title   = "Effect of rotation score on yields",
    caption = paste0(
      "Reference: soy monoculture (score = 0). Two-way clustering: field + year.\n",
      "Score = Z-vector coefficients applied to sequence-level feature averages.\n",
      "Left: stage-1 mean coefficients. Right: stage-2 variance coefficients."
    )
  )
 
ggsave(paste0(fig_dir, "soy_score_yield.png"), score_yield,
       width = 12, height = 8, dpi = 300)
cat("score_yield figure saved.\n")
 
 
# Summary data frame
soy_s1_coef <- broom::tidy(soy_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term),
            mean_est = estimate, mean_se = std.error, mean_p = p.value)
 
soy_s2_coef <- broom::tidy(soy_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term),
            var_est  = estimate, var_se  = std.error, var_p  = p.value)
 
soy_jp_summary <- soy_s1_coef |>
  left_join(soy_s2_coef, by = "rot_crop") |>
  mutate(mean_sig  = mean_p < 0.05,
         var_sig   = var_p  < 0.05,
         dominates = mean_est > 0 & var_est < 0,
         tradeoff  = mean_est > 0 & var_est > 0) |>
  arrange(desc(dominates), desc(mean_est))
 
cat("Soy sequences dominating monoculture:\n")
soy_jp_summary |> 
  filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |> 
  print(n = Inf)
 
# Figure: soy_var_plot — Response of std dev of soy yields to rotation sequences
soy_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(
    system = if_else(grepl("W", term), "Corn-soy-wheat", "Corn-soy"),
    group  = case_when(
      term == "C-S-C-S-C-S" ~ "Perfect rotation",
      .default              = "Other")) |>
  arrange(system, prms.y) |>
  mutate(term = factor(term, levels = unique(term))) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of standard deviation of soy yields to rotation sequences",
       caption = "Reference: S-S-S-S-S-S. Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_var_plot
ggsave(paste0(fig_dir, "soy_var_plot.png"), soy_var_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: soy_coeff_plot — Mean vs variance coefficients
soy_jp_summary |>
  mutate(pr = factor(as.integer(rot_crop %in% c("C-S-C-S-C-S")))) |>
  ggplot(aes(x = var_est, y = mean_est, label = rot_crop)) +
  geom_jitter(aes(color = pr), size = 2) +
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  labs(x = "Soy standard deviation model coefficients",
       y = "Soy mean yield model coefficients") ->
  soy_coeff_plot
ggsave(paste0(fig_dir, "soy_coeff_plot.png"), soy_coeff_plot,
       width = 10, height = 10, dpi = 300)

# Figure: soy_jp_plot — Just-Pope mean-variance decomposition
soy_jp_summary |>
  mutate(
    type = case_when(
      dominates           ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0        ~ "Worse than monoculture",
      TRUE                ~ "No significant difference"),
    perfect = rot_crop %in% c("C-S-C-S-C-S")
  ) |>
  ggplot(aes(x = mean_est, y = var_est, colour = type, shape = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8) +
  geom_errorbar(aes(xmin = mean_est - 1.96*mean_se,
                    xmax = mean_est + 1.96*mean_se),
                orientation = "y", width = 0, alpha = 0.4) +
  geom_errorbar(aes(ymin = var_est - 1.96*var_se,
                    ymax = var_est + 1.96*var_se),
                width = 0, alpha = 0.4) +
  ggrepel::geom_label_repel(data = ~filter(., perfect),
                             aes(label = rot_crop), size = 3,
                             show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Dominates monoculture"     = "#2ca02c",
    "Mean-variance trade-off"   = "#ff7f0e",
    "Worse than monoculture"    = "#d62728",
    "No significant difference" = "grey60")) +
  labs(x      = "Stage 1: effect on mean soy yield (bu/acre)",
       y      = "Stage 2: effect on yield variance",
       colour = NULL, shape = NULL,
       title  = "Just-Pope decomposition: soy rotation effects on mean and variance",
       caption = "Reference: S-S-S-S-S-S. Two-way clustering: field + year.\nQuadrant IV: higher mean, lower variance — unambiguously better.") +
  theme_bw() + theme(legend.position = "bottom") ->
  soy_jp_plot
ggsave(paste0(fig_dir, "soy_jp_plot.png"), soy_jp_plot,
       width = 9, height = 7, dpi = 300)

rm(soy_jp_s2, soy_s1_coef, soy_s2_coef, soy_jp_summary, soy_jp_plot,
   soy_coeff_plot, soy_var_plot); gc()
 
# ── 6. Just-Pope factor RCI — soy ───────────────────────────────────────────
# Table: tab:soy_rci_jp | Figure: rci_plot
 
soy_rci_jp_data <- soy_jp_data |>
  select(-any_of(c("resid_sq", "h_hat"))) |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI))
 
fml_soy_rci_mean <- make_jp_formula("soy_yield", "RCI", all_controls_fgls)
feols(fml_soy_rci_mean, data = soy_rci_jp_data,
      cluster = ~tile_field_ID+year) -> soy_rci_jp_s1
 
soy_rci_jp_s1 |>
  augment(newdata = soy_rci_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_rci_jp_data
 
fml_soy_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls_fgls)
feols(fml_soy_rci_var, data = soy_rci_jp_data,
      cluster = ~tile_field_ID+year) -> soy_rci_jp_s2
 
# Table: tab:soy_rci_jp
etable(soy_rci_jp_s1, soy_rci_jp_s2,
       tex      = TRUE,
       keep     = "^RCI",
       dict     = dict_rci,
       headers  = c("Mean", "Variance"),
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Just-Pope: factor RCI and soy yield moments",
       label    = "tab:soy_rci_jp",
       file     = paste0(tab_dir, "soy_rci_jp.tex"))
 
# Figure: rci_plot — Nonlinear effect of RCI on soy yield mean and variance
rci_plot_df <- bind_rows(
  broom::tidy(soy_rci_jp_s1) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI","",term)),
              est = estimate, se = std.error, moment = "Mean"),
  broom::tidy(soy_rci_jp_s2) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI","",term)),
              est = estimate, se = std.error, moment = "Variance")
)
 
ggplot(rci_plot_df, aes(x = rci, y = est, colour = moment, fill = moment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = est - 1.96*se, ymax = est + 1.96*se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~moment, scales = "free_y") +
  scale_colour_manual(values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_fill_manual(  values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_x_continuous(breaks = sort(unique(rci_plot_df$rci))) +
  labs(x      = "Rotational Complexity Index (RCI)",
       y      = "Coefficient relative to RCI = 0",
       title  = "Nonlinear effect of RCI on soy yield mean and variance",
       caption = "Soy, corn, and wheat fields. Reference: RCI = 0. Two-way clustering.") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) ->
  rci_plot
ggsave(paste0(fig_dir, "soy_rci_plot.png"), rci_plot,
       width = 9, height = 7, dpi = 300)
 
rm(soy_rci_jp_data, soy_rci_jp_s1, soy_rci_jp_s2, rci_plot_df, rci_plot); gc()
 
# ── 7. FGLS + bootstrap — soy ────────────────────────────────────────────────
 
fml_mean  <- make_jp_formula("soy_yield", "rot_crop", all_controls_fgls)
fml_var   <- make_jp_formula("resid_sq",   "rot_crop", all_controls_fgls)
fml_var_b <- make_jp_formula("resid_sq_b", "rot_crop", all_controls_fgls)
 
# Stage 2a to get h_hat for FGLS weights
soy_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_jp_data
 
feols(fml_var, data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s2a
 
soy_jp_s2a |>
  augment(newdata = soy_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |>
  select(-starts_with(".")) ->
  soy_jp_data
 
cat("Obs hitting h_hat floor:", sum(soy_jp_data$h_hat == 1e-6), "\n")
 
feols(fml_var, data = soy_jp_data,
      weights = ~I(1/h_hat),
      cluster = ~tile_field_ID+year) -> soy_jp_s2b
 
etable(soy_jp_s1, soy_jp_s2a, soy_jp_s2b,
       keep  = "rot_crop",
       title = "Just-Pope FGLS: soy rotation effects on mean and variance")
 
rm(soy_jp_s2a); gc()
 
source("save_models_lean.R")
boot_soy <- boot_jp_fgls(soy_jp_data, fml_mean, fml_var, fml_var_b,
                          B = 499, seed = 42, n_workers = 2)
 
saveRDS(
  list(
    # Z-vector models — lean extracts only
    z_s1         = lean_feols(soy_z_s1),
    z_s2         = lean_feols(soy_z_s2),
    # Full sequence models needed for score_yield figure
    jp_s1        = lean_feols(soy_jp_s1),
    jp_s2        = lean_feols(soy_jp_s2),
    # VPD models
    rot_vpd_nc   = lean_feols(soy_rot_vpd_nc),
    rot_vpd      = lean_feols(soy_rot_vpd),
    # Bootstrap — lean strips embedded data from feols objects inside
    boot         = lean_boot(boot_soy),
    # Score construction inputs — just the coefficient vector
    z_coefs      = coef(soy_z_s1),
    # Score data frame — sequence-level feature averages (small)
    score_df     = soy_jp_data |>
                     group_by(rot_crop) |>
                     summarise(late_soy = mean(late_soy, na.rm = TRUE),
                               soy_gap  = mean(soy_gap,  na.rm = TRUE),
                               soy_cons = mean(soy_cons, na.rm = TRUE),
                               nsoy     = mean(nsoy,     na.rm = TRUE),
                               .groups  = "drop")
  ),
  file     = "D:/Crop data/soy_z_models.rds",
  compress = "bzip2"
)
 
cat("Soy lean models saved. File size:",
    round(file.size("D:/Crop data/soy_z_models.rds") / 1e6, 1),
    "MB\n")
 
rm(soy_jp_data, soy_jp_s1, soy_jp_s2b, boot_soy); gc()
 
# ── 8. Spatial maps — soy ────────────────────────────────────────────────────
# Figures: soy_yield_map, rci_map, nccpi_soy_map
# Reload a minimal version of soy_df for spatial plotting only
 
cat("Reloading soy data for spatial maps...\n")
soy_df <- read_parquet(
  "D:/Crop data/d_igis13soy_11_30_2025.with_rci.parquet") 

soy_sf <- soy_df |>
  filter(STATE_ABBR == "IL" & year == 2016) |>
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         soy_yield = soy_yield / 67.25)  |>  # kg/ha -> bu/ac
  arrange(tile_field_ID, year) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs("EPSG:4326"))
  rm(soy_df); gc()
 
il_map <- us_map(regions = "counties") |>
  filter(abbr == "IL") |>
  st_as_sf() |>
  st_transform(st_crs("EPSG:4326"))
 
# Figure: soy_yield_map — Soy yields across Illinois (2016)
soy_sf |>
  filter(!is.na(soy_yield)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = soy_yield), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "Soy yields — QDANN (2016)",
       caption = "Source: Ma et al. (2024)") ->
  soy_yield_map
ggsave(paste0(fig_dir, "soy_yield_map.png"), soy_yield_map,
       width = 10, height = 10, dpi = 300)
 
# Figure: rci_map — RCI values across Illinois (2016)
soy_sf |>
  rename(RCI = annual_RCI) |>
  filter(!is.na(RCI)) |>
  mutate(RCI = factor(RCI)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = RCI), size = 0.2) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "bottom") +
  labs(title   = "Rotational Complexity Index (2016)",
       caption = "Source: CDL / Socolar et al. (2021)") ->
  rci_map
ggsave(paste0(fig_dir, "rci_map.png"), rci_map,
       width = 10, height = 10, dpi = 300)
 
# Figure: nccpi_soy_map — NCCPI soy productivity index (2016)
soy_sf |>
  filter(!is.na(nccpi3soy_mean)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = nccpi3soy_mean), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "NCCPI — Soy (2016)",
       caption = "Source: gSSURGO") ->
  nccpi_soy_map
ggsave(paste0(fig_dir, "nccpi_soy_map.png"), nccpi_soy_map,
       width = 10, height = 10, dpi = 300)
 
rm(soy_sf, il_map, soy_yield_map, rci_map, nccpi_soy_map); gc()
 
# ── 9. County-level yields and maximum entropy distribution ───────────────────
 
cat("Reloading soy data for county-level analysis...\n")
soy_df <- read_parquet(
  "D:/Crop data/d_igis13soy_11_30_2025.with_rci.parquet") 

soy_county <- soy_df |>
  filter(STATE_ABBR == "IL") |>
  select(COUNTY_FIPS, STATE_FIPS, year, soy_yield) |>
  mutate(soy_yield = soy_yield / 67.25) |>  # kg/ha -> bu/ac
  group_by(COUNTY_FIPS, STATE_FIPS, year) |>
  summarise(yield = mean(soy_yield, na.rm = TRUE), .groups = "drop") |>
  rename(state = STATE_FIPS, county = COUNTY_FIPS)
rm(soy_df); gc() 
 
source("./max_entropy/maxent_tack2013.r")
soy_det <- detrend_yields_proportional(soy_county)
dist_m   <- tack_table2(soy_det$yield_norm)
 
dist_m$table 
rm(soy_county, soy_det); gc()
 
cat("\n=== Soy analysis complete ===\n")
cat("All tables saved to:", tab_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
 