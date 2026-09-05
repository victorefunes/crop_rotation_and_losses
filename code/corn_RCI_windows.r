## ============================================================================
## corn_RCI_windows.R
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

corn_df[, origin := min(year), by = tile_field_ID]   # or your actual per-field origin
corn_df[, `:=`(
  window_id  = (year - origin) %/% 6L,
  window_pos = (year - origin) %%  6L
)]

hist <- melt(corn_df,
             id.vars       = c("tile_field_ID", "year"),
             measure.vars  = paste0("crop_", 0:6),
             variable.name = "lag", value.name = "crop")
hist[, k := as.integer(sub("crop_", "", lag))]
hist[, obs_year := year - k]

# if two corn-obs rows disagree about the same field-year, keep the nearer (smaller k) one
setorder(hist, tile_field_ID, obs_year, k)
long <- unique(hist[!is.na(crop)], by = c("tile_field_ID", "obs_year"))[
          , .(tile_field_ID, year = obs_year, crop)]

stopifnot(long[, .N, by = .(tile_field_ID, year)][N > 1L, .N] == 0L)  # no contradictions

source("rci_shapley_decomp.R")
source("decompose_rci_eras.R")

setnames(long, "year", "obs_year", skip_absent=TRUE)

# Windows anchored at origin = 2003 so that:
#   window 0 (t0) = 2003-2008, window 1 (t1) = 2009-2014, window 2 (t2) = 2015-2020
RCI_ORIGIN <- 2003L

# one county per field (guard against a field straddling a county line)
xwalk <- unique(corn_df[, .(tile_field_ID, COUNTY_FIPS)])
stopifnot(xwalk[, .N, by = tile_field_ID][N > 1L, .N] == 0L)

# ── Helper: decompose + aggregate by county + plot for one era-to-era comparison ──
decompose_and_plot <- function(anchor_A, anchor_B, period_label) {

  eras <- decompose_rci_eras(long, id = "tile_field_ID", year = "obs_year", crop = "crop",
                              origin = RCI_ORIGIN, anchor_A = anchor_A, anchor_B = anchor_B,
                              merge_perennial = TRUE)

  cat("\n==", period_label, "==\n")
  print(eras[, .(.N, ok = sum(!is.na(check)), na = sum(is.na(check)))])   # how many usable
  print(eras[!is.na(check), .(max_abs_check = max(abs(check)))])         # ~1e-12
  print(eras[!is.na(dRCI),
             lapply(.SD, function(s) sum(s) / sum(dRCI)),
             .SDcols = patterns("^shap_")])

  ec <- merge(eras, xwalk, by = "tile_field_ID", all.x = TRUE)

  by_county <- ec[!is.na(dRCI), .(
    n_fields    = .N,
    dRCI        = mean(dRCI),
    shap_number = mean(shap_number),
    shap_t1c    = mean(shap_t1c),
    shap_t2     = mean(shap_t2)
  ), by = COUNTY_FIPS][order(COUNTY_FIPS)]

  # exact by construction: per-field dRCI = sum of its shap_*, so the means add up too
  by_county[, check := shap_number + shap_t1c + shap_t2 - dRCI]     # ~1e-15

  # each component's share of the county's mean ΔRCI
  by_county[, `:=`(sh_number = shap_number / dRCI,
                   sh_t1c    = shap_t1c    / dRCI,
                   sh_t2     = shap_t2     / dRCI)]

  print(by_county[])

  print(
    by_county |>
      ggplot(aes(x = dRCI)) +
      geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
      theme_minimal() +
      labs(title = paste("Distribution of Mean ΔRCI by County —", period_label),
           x = "Mean ΔRCI",
           y = "Number of Counties")
  )

  print(
    by_county |>
      pivot_longer(cols = starts_with("sh_"), names_to = "component", values_to = "share") |>
      ggplot(aes(x = share, fill = component)) +
      geom_histogram(binwidth = 10, color = "black") +
      theme_bw() +
      facet_wrap(. ~ component, scales = "free") +
      labs(title = paste("Distribution of ΔRCI Shares by County —", period_label),
           x = "Share of Mean ΔRCI",
           y = "Number of Counties") +
      theme(legend.position = "bottom")
  )

  print(
    by_county |>
      pivot_longer(cols = dRCI:shap_t2, names_to = "component", values_to = "share") |>
      ggplot(aes(x = share, fill = component)) +
      geom_density(color = "black") +
      theme_bw() +
      facet_wrap(. ~ component, scales = "free") +
      geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
      labs(title = paste("Distribution of ΔRCI Shares by County —", period_label),
           x = "Share of Mean ΔRCI",
           y = "Number of Counties") +
      theme(legend.position = "bottom")
  )

  list(eras = eras, by_county = by_county)
}

# t0 (2003-2008) -> t1 (2009-2014)
res_t0_t1 <- decompose_and_plot(anchor_A = 0L, anchor_B = 1L, period_label = "t0 (2003-2008) -> t1 (2009-2014)")
eras3_t0_t1     <- res_t0_t1$eras
by_county_t0_t1 <- res_t0_t1$by_county

# t1 (2009-2014) -> t2 (2015-2020) -- the next period
res_t1_t2 <- decompose_and_plot(anchor_A = 1L, anchor_B = 2L, period_label = "t1 (2009-2014) -> t2 (2015-2020)")
eras3_t1_t2     <- res_t1_t2$eras
by_county_t1_t2 <- res_t1_t2$by_county

res_t0_t2 <- decompose_and_plot(anchor_A = 0L, anchor_B = 2L, period_label = "t0 (2003-2008) -> t2 (2015-2020)")
eras3_t0_t2     <- res_t0_t2$eras
by_county_t0_t2 <- res_t0_t2$by_county

# ── Stack weather/soil covariates into a long panel matching the RCI decompositions ──
# Granularity = field x transition (t0->t1, t1->t2, t0->t2), matching dRCI's granularity:
# each covariate is summarized as its window-A / window-B mean plus the within-field
# change dX = X_B - X_A, so it lines up with dRCI as a "before vs after" comparison.

weather_soil_cols <- all_controls_cols   # pr_*, cGDD_*, EDD_*, vpd_*, soil_* (see rotation_setup_wa.R)
panel_vars        <- c(weather_soil_cols, "corn_yield")

corn_df[, rci_window := (year - RCI_ORIGIN) %/% 6L]

window_means <- corn_df[, lapply(.SD, mean, na.rm = TRUE),
                         by = .(tile_field_ID, rci_window),
                         .SDcols = panel_vars]

# attach window-A / window-B covariate means + deltas to one era's decomposition
build_covariate_panel <- function(eras, anchor_A, anchor_B, period_label) {
  wa <- window_means[rci_window == anchor_A]; setnames(wa, panel_vars, paste0(panel_vars, "_A"))
  wb <- window_means[rci_window == anchor_B]; setnames(wb, panel_vars, paste0(panel_vars, "_B"))
  wa[, rci_window := NULL]; wb[, rci_window := NULL]

  panel <- merge(eras, wa, by = "tile_field_ID", all.x = TRUE)
  panel <- merge(panel, wb, by = "tile_field_ID", all.x = TRUE)

  for (v in panel_vars)
    panel[, (paste0("d_", v)) := get(paste0(v, "_B")) - get(paste0(v, "_A"))]

  # log-yield growth, kept alongside the level difference d_corn_yield
  panel[, d_log_corn_yield := log(corn_yield_B) - log(corn_yield_A)]

  panel[, period := period_label]
  panel[]
}

panel_t0_t1 <- build_covariate_panel(eras3_t0_t1, 0L, 1L, "t0_t1")
panel_t1_t2 <- build_covariate_panel(eras3_t1_t2, 1L, 2L, "t1_t2")
panel_t0_t2 <- build_covariate_panel(eras3_t0_t2, 0L, 2L, "t0_t2")

panel_long <- rbind(panel_t0_t1, panel_t1_t2, panel_t0_t2, fill = TRUE)
panel_long <- merge(panel_long, xwalk, by = "tile_field_ID", all.x = TRUE)   # attach county for FE/clustering

cat("Long panel (field x transition):", nrow(panel_long), "rows,",
    uniqueN(panel_long$tile_field_ID), "fields,",
    uniqueN(panel_long$period), "periods\n")

rm(panel_t0_t1, panel_t1_t2, panel_t0_t2, eras3_t0_t1, eras3_t1_t2, eras3_t0_t2); gc()

corn_df |>
  select(tile_field_ID, STATE_FIPS, lat, lon, nccpi3corn_mean,
          soc0_100_mean, rootznaws_mean) |>
  distinct(.keep_all = TRUE) ->
  id_df 

panel_long <- left_join(panel_long, id_df, by = "tile_field_ID")           

rm(corn_df, id_df);gc() 

feols(d_corn_yield ~ dRCI + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period, 
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d1

feols(d_corn_yield ~ shap_number + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d2

feols(d_corn_yield ~ shap_t1c + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d3


feols(d_corn_yield ~ shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d4

feols(d_corn_yield ~ shap_number+shap_t1c+shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d5
etable(corn_d1, corn_d2, corn_d3, corn_d4, corn_d5, 
        keep = c("dRCI", "shap_number", "shap_t1c", "shap_t2"), 
        digits = 3)


feols(d_corn_yield ~ dRCI + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        tile_field_ID+period, 
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d1

feols(d_corn_yield ~ shap_number + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        tile_field_ID+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d2

feols(d_corn_yield ~ shap_t1c + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        tile_field_ID+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d3


feols(d_corn_yield ~ shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        tile_field_ID+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d4

feols(d_corn_yield ~ shap_number+shap_t1c+shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        tile_field_ID+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_d5
etable(corn_d1, corn_d2, corn_d3, corn_d4, corn_d5, 
        keep = c("dRCI", "shap_number", "shap_t1c", "shap_t2"), 
        digits = 3)        

## In logs
feols(d_log_corn_yield ~ dRCI + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period, 
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_log_d1

feols(d_log_corn_yield ~ shap_number + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_log_d2

feols(d_log_corn_yield ~ shap_t1c + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_log_d3


feols(d_log_corn_yield ~ shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_log_d4

feols(d_log_corn_yield ~ shap_number+shap_t1c+shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + 
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 + 
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean| 
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~COUNTY_FIPS) -> 
        corn_log_d5
etable(corn_log_d1, corn_log_d2, corn_log_d3, corn_log_d4, corn_log_d5, 
        keep = c("dRCI", "shap_number", "shap_t1c", "shap_t2"), 
        digits = 3)        

# ── Robustness: cluster by field instead of county ──────────────────────────
# County clustering assumes independence across ~100 counties; field clustering
# is more conservative (>1 obs per field-period across the 3 stacked transitions)
# and checks whether significance survives a finer clustering level.

feols(d_corn_yield ~ dRCI + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m +
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean|
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~tile_field_ID) ->
        corn_d1_fe

feols(d_corn_yield ~ shap_number + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m +
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean|
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~tile_field_ID) ->
        corn_d2_fe

feols(d_corn_yield ~ shap_t1c + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m +
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean|
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~tile_field_ID) ->
        corn_d3_fe

feols(d_corn_yield ~ shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m +
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean|
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~tile_field_ID) ->
        corn_d4_fe

feols(d_corn_yield ~ shap_number+shap_t1c+shap_t2 + d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m +
        d_cGDD_8m + d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
        d_soil_6 + d_soil_7 + d_soil_8+nccpi3corn_mean+soc0_100_mean+rootznaws_mean|
        COUNTY_FIPS+period,
        data = panel_long, cluster = ~tile_field_ID) ->
        corn_d5_fe

etable(corn_d1_fe, corn_d2_fe, corn_d3_fe, corn_d4_fe, corn_d5_fe, 
        keep = c("dRCI", "shap_number", "shap_t1c", "shap_t2"), 
        digits = 3)

cor(panel_long[, .(shap_number, shap_t1c, shap_t2)], use = "complete.obs")
lm_d5 <- lm(d_corn_yield ~ shap_number + shap_t1c + shap_t2 +
              d_pr_6 + d_pr_7 + d_pr_8 + d_cGDD_6m + d_cGDD_7m + d_cGDD_8m +
              d_EDD_6 + d_EDD_7 + d_EDD_8 + d_vpd_6 + d_vpd_7 + d_vpd_8 +
              d_soil_6 + d_soil_7 + d_soil_8 +
              nccpi3corn_mean + soc0_100_mean + rootznaws_mean +
              factor(COUNTY_FIPS) + factor(period),
            data = panel_long)
car::vif(lm_d5)

## Try Lewbel IV in this contxt, using the heteroskedasticity-based 
# instruments for the three RCI components