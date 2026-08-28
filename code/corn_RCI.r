## ============================================================================
## corn_RCI.R
## Just-Pope production risk analysis — CORN
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

corn_df[, has_soy := Reduce(`|`, lapply(.SD, \(v) v %in% soy_code)),
        .SDcols = xcols]

# Perfect corn/soy rotation = strict yearly alternation over the 6-yr window
# (crop_0 is always corn, so this means crop_1,3,5 = soy and crop_2,4 = corn)
corn_df[, perfect_cs := crop_0 == corn_code & crop_1 == soy_code &
                         crop_2 == corn_code & crop_3 == soy_code &
                         crop_4 == corn_code & crop_5 == soy_code]

# corn monoculture, corn/soy rotations (perfect vs other), corn/soy/wheat, other
corn_df[, regime := fcase(
  has_other | has_alfalfa,    "other",                   # checked first — wins
  has_wheat,                  "corn_soy_wheat",
  has_soy & perfect_cs,       "corn_soy_perfect",
  has_soy,                    "corn_soy_other",
  default =                   "corn_monoculture"
)]
corn_df[, regime := factor(regime, levels = c(
  "corn_monoculture","corn_soy_perfect","corn_soy_other","corn_soy_wheat","other"))]

corn_df[, .N, by = regime][order(-N)]
 
corn_df <- corn_df |>
  filter(!is.na(regime))

# Free raw data — no longer needed
gc() 

source("rci_shapley_decomp.R")

out <- decompose_rci(corn_df)                       # per field-year contributions
out[!is.na(dRCI), .(max_abs_check = max(abs(check)))]        # should be ~1e-12

# # sample-level attribution: of the average RCI change, how much is each margin?
out[dRCI > 0, lapply(.SD, mean), .SDcols = patterns("^shap_")]
#
# # share of total ΔRCI variation carried by each component:
out[!is.na(dRCI), lapply(.SD, function(s) sum(s)/sum(dRCI)),
     .SDcols = patterns("^shap_")]

corn_df <- corn_df |>
  left_join(out, by = c("tile_field_ID", "year")) 

rm(out); gc()

corn_df[!is.na(dRCI), lapply(.SD, function(s) sum(s)/sum(dRCI)),
     .SDcols = patterns("^shap_"), by = regime]

corn_df[!is.na(dRCI), lapply(.SD, function(s) mean(s)),
     .SDcols = patterns("^shap_"), by = regime]


# Categories:.default# 1. Corn monoculture
# 2. Corn/soy rotations
# 3. Corn/soy and wheat rotations

#Q1 what is the correct specification?

reg_pure <- feols(corn_yield~RCI:regime + regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + 
    rootznaws_mean| tile_field_ID + year, corn_df, 
    cluster=~COUNTY_FIPS)
reg_full <- feols(corn_yield~RCI*regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + 
    rootznaws_mean| tile_field_ID + year, corn_df, 
    cluster=~COUNTY_FIPS)
wald(reg_full, "^regime") 

# the regime main effects are not redundant with tile_field_ID fixed effects — 
# meaning regime does vary meaningfully within fields over time (or at least isn't fully 
# absorbed by the FE), and dropping it (as RCI:regime-only does) would misattribute genuine 
# regime-level yield differences into the RCI interaction slopes

reg_factorial <- feols(
  corn_yield ~ RCI*regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~COUNTY_FIPS
)
etable(reg_factorial, keep = c("RCI", "regime"))

 # This confirms the full-interaction spec was the right call — both the regime main effects and the 
 # RCI×regime slopes are substantively different from each other and mostly significant:

# Baseline regime (implicit reference category, presumably continuous corn) has RCI slope = 0.740.
# corn_soy_perfect: RCI effect = 0.740 + 0.896 = 1.64 — much steeper RCI-yield sensitivity than baseline.
# corn_soy_other: RCI effect = 0.740 + 1.582 = 2.32 — the steepest.
# corn_soy_wheat: RCI effect = 0.740 − 0.694 ≈ 0.046 — essentially flat, RCI barely matters under this regime.
# other: RCI effect = 0.740 + 0.600 = 1.34.

# The interaction terms being significant (and of varying sign/magnitude relative to baseline) is 
# exactly the pattern the earlier Wald test predicted — regime materially moderates the RCI-yield 
# relationship, not just the intercept. This is good evidence you'd have gotten a distorted, blended 
# RCI slope had you used RCI:regime alone without letting the intercepts (main effects) shift freely.


library(marginaleffects)

# Average marginal effect of RCI, by regime, this computes computes ∂corn_yield/∂RCI for each regime
mfx <- slopes(reg_factorial, variables = "RCI", by = "regime",
              newdata = datagrid(regime = unique(corn_df$regime)),
              vcov = ~COUNTY_FIPS)
mfx

ggplot(mfx, aes(x = regime, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  labs(y = "Marginal effect of RCI on corn yield", x = "Regime") +
  theme_minimal()


# Cluster selection
reg_factorial_alt <- feols(
  corn_yield ~ RCI*regime +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8 + rootznaws_mean
    | tile_field_ID + year,
  data = corn_df, cluster = ~COUNTY_FIPS + year
)
etable(reg_factorial, reg_factorial_alt, keep = c("RCI", "regime"))

V <- vcov(reg_factorial_alt, cluster = ~COUNTY_FIPS + year)

# This is a known artifact of the Cameron-Gelbach-Miller (CGM) multiway clustering 
# formula — it's not additive across dimensions the way one-way clustering is, and 
# with a large number of fixed effects (tile_field_ID has presumably thousands of levels) 
# relative to the cluster counts (102 counties × 15 years), the resulting covariance estimate 
# can end up not positive semi-definite. fixest automatically applies an eigenvalue correction 
# (clips negative eigenvalues to zero) to make it usable — this is standard practice, the same 
# fix Stata's reghdfe/cgmreg and other implementations use, so it's not a sign your model is broken.

mfx <- slopes(reg_factorial, variables = "RCI", by = "regime",
              newdata = datagrid(regime = unique(corn_df$regime)),
              vcov = V)
mfx

ggplot(mfx, aes(x = regime, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  labs(y = "Marginal effect of RCI on corn yield", x = "Regime") +
  theme_minimal()

# Why not FEs at the tile_field_ID+year level?
# switching to tile_field_ID + year would likely understate your SEs relative to COUNTY_FIPS + year,
# it drops the meaningful spatial clustering in favor of a within-field-only dimension that's largely 
# redundant with your FE structure.



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

# what crops drive the results in switch-reg_factor?
# Decompse RCi by changes in iots elements

## Stigler's model
setorder(corn_df, tile_field_ID, year)

# 1) per-observation trigger: define the ONE event you want to anchor on.
#    simplest = "field is no longer plain corn_soy this year"
corn_df[, treated_now := regime != "corn_soy"]          # or: has_wheat | has_alfalfa | has_other

# 2) first year the trigger fires, per field; Inf if it never does
corn_df[, first_treat := {yrs <- year[treated_now]
  if (length(yrs)) min(yrs) else Inf
}, by = tile_field_ID]

# Pre-compute the relative period so binning `first_treat` (below) can't
# interfere with sunab's automatic period-from-cohort calculation.
corn_df[, rel_period := year - first_treat]

# The full dynamic Sun-Abraham decomposition (cohort x relative-period
# interactions) is unidentified with tile_field_ID FE: ~300K field levels at
# ~7.5 obs/field leaves too little within-field variation to separately pin
# down ~20-30 interaction parameters (confirmed by comparing to a single
# static treated_now dummy, which IS precisely estimated under the same FE).
# So: (1) headline effect from a static ATT, (2) a coarse 3-bin event study
# (pre / event-year / post) as a lightweight pre-trends check, both cheap
# enough in parameters to be identified at this FE cardinality.

# (1) Static ATT — headline effect
est_att <- feols(
  corn_yield ~ treated_now +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8
    | tile_field_ID + year,
  data = corn_df, vcov = "hetero", mem.clean = TRUE)
etable(est_att)

# (2) Coarse pre/event/post event study — pre-trends validity check.
# rootznaws_mean is dropped: it's time-invariant per field, so it's always
# fully absorbed by tile_field_ID FE (confirmed collinear in every field-FE
# spec tested).
corn_df[, event_bin := fcase(
  is.infinite(first_treat), "control",
  rel_period <= -1,         "pre",
  rel_period == 0,          "event",
  rel_period >= 1,          "post"
)]
corn_df[, event_bin := factor(event_bin,
  levels = c("control", "pre", "event", "post"))]

# "control" fields never leave the reference category, so within tile_field_ID
# FE they carry zero identifying variation for pre/event/post — they only
# inflate the sample without helping. Restrict to switcher fields, where the
# comparison is a genuine within-field transition (pre -> event -> post).
switchers_df <- corn_df[event_bin != "control"]
switchers_df[, event_bin := droplevels(event_bin)]   # drops the unused "control"
                                                       # level so pre/event/post
                                                       # don't sum to a constant

est_es <- feols(
  corn_yield ~ event_bin +
    pr_6 + pr_7 + pr_8 + I(pr_6^2) + I(pr_7^2) + I(pr_8^2) +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8
    | tile_field_ID + year,
  data = switchers_df, vcov = "hetero", mem.clean = TRUE)

etable(est_es)

## Callaway & Sant'Anna group-time ATT — alternative estimator that computes
## each (cohort, calendar-year) effect from its own 2x2 DiD against a
## not-yet-treated comparison group, instead of one giant two-way-FE
## regression. Sidesteps the tile_field_ID cardinality problem entirely.
library(did)

corn_df[, first_treat_cs := fifelse(is.infinite(first_treat), 0, first_treat)]
corn_df[, field_id_num := .GRP, by = tile_field_ID]

cs_out <- att_gt(
  yname = "corn_yield",
  tname = "year",
  idname = "field_id_num",
  gname = "first_treat_cs",
  xformla = ~ pr_6 + pr_7 + pr_8 +
    cGDD_6m + cGDD_7m + cGDD_8m + EDD_6 + EDD_7 + EDD_8 +
    vpd_6 + vpd_7 + vpd_8 + soil_6 + soil_7 + soil_8,
  data = corn_df,
  control_group = "notyettreated",
  allow_unbalanced_panel = TRUE,
  bstrap = FALSE,
  panel = TRUE
)
summary(cs_out)

cs_dynamic <- aggte(cs_out, type = "dynamic", na.rm = TRUE, min_e = -4, max_e = 4)
summary(cs_dynamic)

cs_simple <- aggte(cs_out, type = "simple", na.rm = TRUE)
summary(cs_simple)
