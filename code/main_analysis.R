## ============================================================================
## Crop Rotations and Yield Risk: Just-Pope Production Risk Analysis
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas — Department of Agricultural Economics and Agribusiness
## ============================================================================
## File structure:
##   0. Setup and data loading
##   1. Rotation sequence patterns and RCI correction
##   2. PCA of rotation features (Figure: rot_pca_plot)
##   3. OLS mean models — corn (Table: tab:corn_rot)
##   4. OLS mean models — soy (Table: tab:soy_rot)
##   5. RCI mean models — corn and soy (Tables: tab:corn_rci, tab:soy_rci)
##   6. VPD interaction models (Tables: tab:corn_rot_vpd, tab:soy_rot_vpd,
##                                      tab:corn_rci_vpd, tab:soy_rci_vpd)
##   7. Just-Pope stage 1 + stage 2 — corn (Tables: tab:corn_jp_mean, tab:corn_jp_var)
##      Mean-variance scatter (Figure: corn_jp_plot)
##      Coefficient plots (Figures: corn_rot_plot, corn_var_plot, corn_coeff_plot)
##   8. Just-Pope stage 1 + stage 2 — soy (Tables: tab:soy_jp_mean, tab:soy_jp_var)
##      Mean-variance scatter (Figure: soy_jp_plot)
##      Coefficient plots (Figures: soy_rot_plot, soy_var_plot, soy_coeff_plot)
##   9. Just-Pope factor RCI (Tables: tab:corn_rci_jp, tab:soy_rci_jp)
##      RCI nonlinearity figure (Figure: rci_plot)
##  10. FGLS estimation + bootstrap inference
##  11. Spatial maps (Figures: rci_map, corn_yield_map, soy_yield_map,
##                             nccpi_corn_map, nccpi_soy_map)
##  12. County-level yields and maximum entropy distribution
## ============================================================================

# ── 0. Setup ─────────────────────────────────────────────────────────────────

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
library(ggrepel)
library(sf)
library(usmap)
theme_set(theme_bw())

# plan(multisession, workers = max(1, parallel::detectCores() - 1))

setwd("C:/Users/vf006/Box/Economic Analysis of Soil Health Practices")

tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"
fig_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/figures/"

# Load separate corn and soy files (as per current pipeline)
corn_df <- fread("./Data and Data Descriptions/clean/corn_rci_il_long.csv")
soy_df  <- fread("./Data and Data Descriptions/clean/soy_rci_il_long.csv")

corn_df <- corn_df[order(tile_field_ID, year)]
soy_df  <- soy_df[order(tile_field_ID, year)]

# ── 1. Rotation patterns and RCI correction ───────────────────────────────────

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

# RCI correction — remove mismatched RCI/rot_crop combinations
# (same logic applied to both corn_df and soy_df)
rci_correction <- function(df) {
  df |>
    mutate(data_rm = case_when(
      rot_crop == "5-1-5-1-5-1" & RCI == 3.24 ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 3    ~ 1,
      rot_crop == "5-1-5-1-1-1" & RCI == 2.24 ~ 1,
      rot_crop == "5-1-1-5-1-5" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-5-1-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-1-5-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 2    ~ 1,
      rot_crop == "1-5-1-1-1-5" & RCI == 1.73 ~ 1,
      rot_crop == "1-1-1-5-1-5" & RCI == 0    ~ 1,
      rot_crop == "1-1-1-1-1-5" & RCI == 0    ~ 1,
      rot_crop == "5-1-5-1-5-1" & RCI == 2.45 ~ 1,
      rot_crop == "1-5-1-1-5-1" & RCI == 2.24 ~ 1,
      rot_crop == "1-5-1-5-1-5" & RCI == 0    ~ 1,
      rot_crop == "1-5-1-5-1-5" & RCI == 2.45 ~ 1,
      .default = 0)) |>
    filter(data_rm == 0) |>
    select(-data_rm)
}

corn_df <- rci_correction(corn_df)
soy_df  <- rci_correction(soy_df)

# Transition sequences for reference (RCI values)
# C-C-C-C-C-C: 0.00  |  S-S-S-S-S-S: 0.00
# C-C-C-C-C-S: 1.41  |  S-S-S-S-S-C: 1.41
# C-C-C-C-S-C: 1.73  |  S-S-S-S-C-S: 1.73
# C-C-C-S-C-S: 2.00  |  S-S-S-C-S-C: 2.00
# C-C-S-C-S-C: 2.24  |  S-S-C-S-C-S: 2.24
# C-S-C-S-C-S: 2.24  |  S-C-S-C-S-C: 2.24

# ── Controls ──────────────────────────────────────────────────────────────────

# Degree-day variables (Schlenker-Roberts 2009)
compute_gdd <- function(tmin, tmax, base = 8, cap = 29) {
  pmax(0, pmin((tmax + tmin) / 2, cap) - base)
}
compute_edd <- function(tmin, tmax, threshold = 29) {
  pmax(0, (tmax + tmin) / 2 - threshold)
}

add_degree_days <- function(df) {
  df |>
    mutate(EDD_6 = compute_edd(tmmx_6, tmmn_6),
           EDD_7 = compute_edd(tmmx_7, tmmn_7),
           EDD_8 = compute_edd(tmmx_8, tmmn_8),
           GDD_6 = compute_gdd(tmmx_6, tmmn_6),
           GDD_7 = compute_gdd(tmmx_7, tmmn_7),
           GDD_8 = compute_gdd(tmmx_8, tmmn_8))
}

corn_df <- add_degree_days(corn_df)
soy_df  <- add_degree_days(soy_df)

# Revised control vector: quadratic precipitation, EDD, drop tmmx/tmmn
# nccpi3all_mean and soc0_100_mean dropped — collinear with field FEs
all_controls <- c(
  "pr_6", "pr_7", "pr_8",
  "I(pr_6^2)", "I(pr_7^2)", "I(pr_8^2)",
  "cGDD_6m", "cGDD_7m", "cGDD_8m",
  "EDD_6", "EDD_7", "EDD_8",
  "vpd_6", "vpd_7", "vpd_8",
  "soil_6", "soil_7", "soil_8",
  "rootznaws_mean"
)

# Formula helper — used throughout
make_jp_formula <- function(lhs, rot_var, controls,
                            fe = "tile_field_ID + year") {
  rhs <- paste(c(rot_var, controls), collapse = " + ")
  as.formula(paste(lhs, "~", rhs, "|", fe))
}

# Coefficient label dictionaries
dict_corn <- c(
  "rot_cropC-C-C-C-S-C" = "C-C-C-C-S-C", "rot_cropC-C-C-S-C-C" = "C-C-C-S-C-C",
  "rot_cropC-C-C-S-S-C" = "C-C-C-S-S-C", "rot_cropC-C-S-C-C-C" = "C-C-S-C-C-C",
  "rot_cropC-C-S-C-S-C" = "C-C-S-C-S-C", "rot_cropC-C-S-S-C-C" = "C-C-S-S-C-C",
  "rot_cropC-C-S-S-S-C" = "C-C-S-S-S-C", "rot_cropC-S-C-C-C-C" = "C-S-C-C-C-C",
  "rot_cropC-S-C-C-S-C" = "C-S-C-C-S-C", "rot_cropC-S-C-S-C-C" = "C-S-C-S-C-C",
  "rot_cropC-S-C-S-S-C" = "C-S-C-S-S-C", "rot_cropC-S-S-C-C-C" = "C-S-S-C-C-C",
  "rot_cropC-S-S-C-S-C" = "C-S-S-C-S-C", "rot_cropC-S-S-S-C-C" = "C-S-S-S-C-C",
  "rot_cropC-S-S-S-S-C" = "C-S-S-S-S-C", "rot_cropS-C-C-C-C-C" = "S-C-C-C-C-C",
  "rot_cropS-C-C-C-S-C" = "S-C-C-C-S-C", "rot_cropS-C-C-S-C-C" = "S-C-C-S-C-C",
  "rot_cropS-C-C-S-S-C" = "S-C-C-S-S-C", "rot_cropS-C-S-C-C-C" = "S-C-S-C-C-C",
  "rot_cropS-C-S-C-S-C" = "S-C-S-C-S-C", "rot_cropS-C-S-S-C-C" = "S-C-S-S-C-C",
  "rot_cropS-C-S-S-S-C" = "S-C-S-S-S-C", "rot_cropS-S-C-C-C-C" = "S-S-C-C-C-C",
  "rot_cropS-S-C-C-S-C" = "S-S-C-C-S-C", "rot_cropS-S-C-S-C-C" = "S-S-C-S-C-C",
  "rot_cropS-S-C-S-S-C" = "S-S-C-S-S-C", "rot_cropS-S-S-C-C-C" = "S-S-S-C-C-C",
  "rot_cropS-S-S-C-S-C" = "S-S-S-C-S-C", "rot_cropS-S-S-S-C-C" = "S-S-S-S-C-C",
  "rot_cropS-S-S-S-S-C" = "S-S-S-S-S-C"
)

dict_soy <- c(
  "rot_cropC-C-C-C-C-S" = "C-C-C-C-C-S", "rot_cropC-C-C-C-S-S" = "C-C-C-C-S-S",
  "rot_cropC-C-C-S-C-S" = "C-C-C-S-C-S", "rot_cropC-C-C-S-S-S" = "C-C-C-S-S-S",
  "rot_cropC-C-S-C-C-S" = "C-C-S-C-C-S", "rot_cropC-C-S-C-S-S" = "C-C-S-C-S-S",
  "rot_cropC-C-S-S-C-S" = "C-C-S-S-C-S", "rot_cropC-C-S-S-S-S" = "C-C-S-S-S-S",
  "rot_cropC-S-C-C-C-S" = "C-S-C-C-C-S", "rot_cropC-S-C-C-S-S" = "C-S-C-C-S-S",
  "rot_cropC-S-C-S-C-S" = "C-S-C-S-C-S", "rot_cropC-S-C-S-S-S" = "C-S-C-S-S-S",
  "rot_cropC-S-S-C-C-S" = "C-S-S-C-C-S", "rot_cropC-S-S-C-S-S" = "C-S-S-C-S-S",
  "rot_cropC-S-S-S-C-S" = "C-S-S-S-C-S", "rot_cropC-S-S-S-S-S" = "C-S-S-S-S-S",
  "rot_cropS-C-C-C-C-S" = "S-C-C-C-C-S", "rot_cropS-C-C-C-S-S" = "S-C-C-C-S-S",
  "rot_cropS-C-C-S-C-S" = "S-C-C-S-C-S", "rot_cropS-C-C-S-S-S" = "S-C-C-S-S-S",
  "rot_cropS-C-S-C-C-S" = "S-C-S-C-C-S", "rot_cropS-C-S-C-S-S" = "S-C-S-C-S-S",
  "rot_cropS-C-S-S-C-S" = "S-C-S-S-C-S", "rot_cropS-C-S-S-S-S" = "S-C-S-S-S-S",
  "rot_cropS-S-C-C-C-S" = "S-S-C-C-C-S", "rot_cropS-S-C-C-S-S" = "S-S-C-C-S-S",
  "rot_cropS-S-C-S-C-S" = "S-S-C-S-C-S", "rot_cropS-S-C-S-S-S" = "S-S-C-S-S-S",
  "rot_cropS-S-S-C-C-S" = "S-S-S-C-C-S", "rot_cropS-S-S-C-S-S" = "S-S-S-C-S-S",
  "rot_cropS-S-S-S-C-S" = "S-S-S-S-C-S"
)

dict_rci <- c(
  "RCI1.41" = "RCI = 1.41", "RCI1.73" = "RCI = 1.73",
  "RCI2"    = "RCI = 2",    "RCI2.24" = "RCI = 2.24",
  "RCI2.45" = "RCI = 2.45", "RCI2.65" = "RCI = 2.65",
  "RCI2.74" = "RCI = 2.74", "RCI2.83" = "RCI = 2.83",
  "RCI3"    = "RCI = 3",    "RCI3.24" = "RCI = 3.24",
  "RCI3.46" = "RCI = 3.46", "RCI3.67" = "RCI = 3.67",
  "RCI3.74" = "RCI = 3.74", "RCI4"    = "RCI = 4",
  "RCI4.24" = "RCI = 4.24", "RCI4.47" = "RCI = 4.47",
  "RCI4.74" = "RCI = 4.74", "RCI5.2"  = "RCI = 5.2"
)

# ── 2. PCA of rotation features ───────────────────────────────────────────────
# Motivates the Z-vector parameterization (Section 4.2 of paper)
# Figure: rot_pca_plot

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
    free_soy  = ifelse(gap_soy     == 0, 0, 1 / gap_soy),
    tight_soy = ifelse(min_gap_soy == 0, 0, 1 / min_gap_soy)
  )

pca_rot <- prcomp(rot_pca_in, scale. = TRUE)
cat("Rotation PCA — variance explained by first two PCs:\n")
print(summary(pca_rot)$importance[, 1:2])
cat("PC1 loadings:\n")
print(pca_rot$rotation[, 1])

rot_features <- rot_features |>
  mutate(rot_index = pca_rot$x[, 1])

# Figure: rot_pca_plot — PCA biplot (Section 4.2 / Figure PCA in paper)
autoplot(pca_rot, data = corn_soy_patterns,
         loadings = TRUE, loadings.label = TRUE) +
  labs(title = "PCA of corn-soy rotation features") +
  theme_bw() +
  theme(legend.position = "none") ->
  rot_pca_plot
ggsave(paste0(fig_dir, "rot_pca_plot.png"), rot_pca_plot,
       width = 9, height = 7, dpi = 300)

# ── 3. OLS mean models — corn ─────────────────────────────────────────────────
# Tables: tab:corn_rot (no controls / with controls)
# Figures: corn_rot_plot

corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  left_join(
    rot_features |> select(pattern, n_soy, consec_soy, gap_soy, min_gap_soy, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  mutate(
    soy_lag1 = as.integer(substr(rot_crop, 9, 9) == "5"),
    soy_lag2 = as.integer(substr(rot_crop, 7, 7) == "5"),
    soy_lag3 = as.integer(substr(rot_crop, 5, 5) == "5"),
    soy_lag4 = as.integer(substr(rot_crop, 3, 3) == "5"),
    soy_lag5 = as.integer(substr(rot_crop, 1, 1) == "5"),
    rot_crop  = gsub("1", "C", gsub("5", "S", rot_crop)),
    rot_crop  = factor(rot_crop),
    rot_crop  = relevel(rot_crop, ref = "C-C-C-C-C-C")
  ) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls), ~ !is.na(.)))

cat("Corn analysis sample:", nrow(corn_jp_data), "rows\n")

corn_yield_formula <- make_jp_formula("corn_yield", "rot_crop", all_controls)

# No-controls model
feols(corn_yield ~ rot_crop | tile_field_ID + year,
      data    = corn_jp_data,
      cluster = ~COUNTY_FIPS) ->
  corn_rot_nc

# With-controls model
feols(corn_yield_formula,
      data    = corn_jp_data,
      cluster = ~COUNTY_FIPS) ->
  corn_rot

# Table: tab:corn_rot — Rotation patterns and corn yields
etable(corn_rot_nc, corn_rot,
       dict     = dict_corn,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       se.below = FALSE,
       title    = "Rotation patterns and corn yields",
       label    = "tab:corn_rot",
       extralines = list("_Controls" = c("No", "Yes")))
etable(corn_rot_nc, corn_rot,
       tex      = TRUE,
       dict     = dict_corn,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       se.below = FALSE,
       title    = "Rotation patterns and corn yields",
       label    = "tab:corn_rot",
       extralines = list("_Controls" = c("No", "Yes")),
       file     = paste0(tab_dir, "corn_rot.tex"))

# Figure: corn_rot_plot — Response of corn yields to rotation sequences
corn_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "S-C-S-C-S-C" ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-C-S-C") ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  labs(title = "Response of corn yields to rotation sequences",
       caption = "Reference: corn monoculture (C-C-C-C-C-C). Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rot_plot
ggsave(paste0(fig_dir, "corn_rot_plot.png"), corn_rot_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 4. OLS mean models — soy ──────────────────────────────────────────────────
# Tables: tab:soy_rot
# Figures: soy_rot_plot

soy_jp_data <- soy_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(
    rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")
  ) |>
  filter(!is.na(soy_yield)) |>
  filter(if_all(all_of(all_controls), ~ !is.na(.)))

cat("Soy analysis sample:", nrow(soy_jp_data), "rows\n")

soy_yield_formula <- make_jp_formula("soy_yield", "rot_crop", all_controls)

feols(soy_yield ~ rot_crop | tile_field_ID + year,
      data    = soy_jp_data,
      cluster = ~COUNTY_FIPS) ->
  soy_rot_nc

feols(soy_yield_formula,
      data    = soy_jp_data,
      cluster = ~COUNTY_FIPS) ->
  soy_rot

# Table: tab:soy_rot — Rotation patterns and soybean yields
etable(soy_rot_nc, soy_rot,
       tex      = TRUE,
       dict     = dict_soy,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       se.below = FALSE,
       title    = "Rotation patterns and soybean yields",
       label    = "tab:soy_rot",
       extralines = list("_Controls" = c("No", "Yes")),
       file     = paste0(tab_dir, "soy_rot.tex"))

# Figure: soy_rot_plot — Response of soybean yields to rotation sequences
soy_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S",
                "S-S-S-S-C-S", "S-S-C-S-C-S") ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  labs(title = "Response of soybean yields to rotation sequences",
       caption = "Reference: soy monoculture (S-S-S-S-S-S). Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rot_plot
ggsave(paste0(fig_dir, "soy_rot_plot.png"), soy_rot_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 5. RCI models — corn and soy ──────────────────────────────────────────────
# Tables: tab:corn_rci, tab:soy_rci
# Figures: corn_rci_plot, soy_rci_plot

corn_rci_formula <- make_jp_formula("corn_yield", "RCI", all_controls)
soy_rci_formula  <- make_jp_formula("soy_yield",  "RCI", all_controls)

# All crops
corn_df |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_all

soy_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_all

# Corn-soy only
corn_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_cs

soy_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_cs

# Table: tab:corn_rci — Rotational Complexity and corn yields
etable(corn_rci_all, corn_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("All crops", "Corn and soybeans only"),
       keep     = "RCI",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       se.below = FALSE,
       title    = "Rotation Complexity and corn yields",
       label    = "tab:corn_rci",
       file     = paste0(tab_dir, "corn_rci.tex"))

# Table: tab:soy_rci — Rotational Complexity and soybean yields
etable(soy_rci_all, soy_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("All crops", "Corn and soybeans only"),
       keep     = "RCI",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       se.below = FALSE,
       title    = "Rotation Complexity and soybean yields",
       label    = "tab:soy_rci",
       file     = paste0(tab_dir, "soy_rci.tex"))

# Figure: corn_rci_plot — RCI effects on corn yields
corn_rci_cs |>
  coefplot() |>
  data.frame() |>
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  mutate(term = as.numeric(term),
         group = case_when(
           term == 2.24 ~ "Perfect rotation",
           term %in% c(1.41, 1.73, 2) ~ "Transitioning",
           .default = "Other")) |>
  filter(term < 5.2) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Rotational Complexity Index (RCI)") +
  ylab("Coefficient Estimate") +
  labs(title = "Changes in RCI and corn yields",
       caption = "Corn-soy fields only. Reference: RCI = 1.41. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rci_plot
ggsave(paste0(fig_dir, "corn_rci_plot.png"), corn_rci_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: soy_rci_plot — RCI effects on soybean yields
soy_rci_cs |>
  coefplot() |>
  data.frame() |>
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  mutate(term = as.numeric(term),
         group = case_when(
           term == 2.24 ~ "Perfect rotation",
           term %in% c(1.41, 1.73, 2) ~ "Transitioning",
           .default = "Other")) |>
  filter(term < 5.2) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Rotational Complexity Index (RCI)") +
  ylab("Coefficient Estimate") +
  labs(title = "Changes in RCI and soybean yields",
       caption = "Corn-soy fields only. Reference: RCI = 1.41. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rci_plot
ggsave(paste0(fig_dir, "soy_rci_plot.png"), soy_rci_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 6. VPD interaction models ─────────────────────────────────────────────────
# Tables: tab:corn_rot_vpd, tab:soy_rot_vpd, tab:corn_rci_vpd, tab:soy_rci_vpd

corn_df |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ 1,
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ 2,
    vpdmax_7 > 2.1                     ~ 3,
    .default = NA),
    vpd_name = factor(vpd_name, labels = c("normal", "somewhat dry", "dry"))) ->
  corn_df

soy_df |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ 1,
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ 2,
    vpdmax_7 > 2.1                     ~ 3,
    .default = NA),
    vpd_name = factor(vpd_name, labels = c("normal", "somewhat dry", "dry"))) ->
  soy_df

dict_vpd <- c("vpd_namesomewhatdry" = "Somewhat dry season",
              "vpd_namedry"         = "Dry season")

corn_vpd_formula    <- make_jp_formula("corn_yield", "rot_crop + vpd_name", all_controls)
soy_vpd_formula     <- make_jp_formula("soy_yield",  "rot_crop + vpd_name", all_controls)
corn_rci_vpd_formula <- make_jp_formula("corn_yield", "RCI * vpd_name",     all_controls)
soy_rci_vpd_formula  <- make_jp_formula("soy_yield",  "RCI * vpd_name",     all_controls)

corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot_vpd

soy_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")) |>
  feols(soy_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rot_vpd

corn_df |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_vpd

soy_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_vpd

# Table: tab:corn_rot_vpd
etable(corn_rot_vpd,
       tex      = TRUE,
       dict     = c(dict_corn, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Effect of weather and rotation sequences on corn yields",
       label    = "tab:corn_rot_vpd",
       file     = paste0(tab_dir, "corn_rot_vpd.tex"))

# Table: tab:soy_rot_vpd
etable(soy_rot_vpd,
       tex      = TRUE,
       dict     = c(dict_soy, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Effect of weather and rotation sequences on soybean yields",
       label    = "tab:soy_rot_vpd",
       file     = paste0(tab_dir, "soy_rot_vpd.tex"))

# Table: tab:corn_rci_vpd
dict_rci_vpd <- c(dict_vpd,
                  "RCI x vpd_namesomewhatdry" = "RCI x Somewhat dry season",
                  "RCI x vpd_namedry"          = "RCI x Dry season")
etable(corn_rci_vpd,
       tex      = TRUE,
       dict     = dict_rci_vpd,
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "RCI x drought interaction effects on corn yields",
       label    = "tab:corn_rci_vpd",
       file     = paste0(tab_dir, "corn_rci_vpd.tex"))

# Table: tab:soy_rci_vpd
etable(soy_rci_vpd,
       tex      = TRUE,
       dict     = dict_rci_vpd,
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_", "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "RCI x drought interaction effects on soybean yields",
       label    = "tab:soy_rci_vpd",
       file     = paste0(tab_dir, "soy_rci_vpd.tex"))

# ── 7. Just-Pope — corn ───────────────────────────────────────────────────────
# Tables: tab:corn_jp_mean (stage 1), tab:corn_jp_var (stage 2 OLS)
# Figures: corn_jp_plot (mean-variance scatter), corn_var_plot, corn_coeff_plot

# BH multiple comparisons correction
pvals_corn <- broom::tidy(corn_rot) |>
  filter(grepl("rot_crop", term)) |>
  arrange(p.value) |>
  mutate(
    p_adj_bh = p.adjust(p.value, method = "BH"),
    sig_raw  = p.value  < 0.05,
    sig_bh   = p_adj_bh < 0.05
  )
cat("Corn sequences significant at 5% (unadjusted):", sum(pvals_corn$sig_raw), "\n")
cat("Corn sequences significant at 5% FDR (BH):    ", sum(pvals_corn$sig_bh),  "\n")

bh_note_corn <- paste0(
  "Benjamini-Hochberg FDR correction across ", nrow(pvals_corn),
  " rotation-sequence coefficients: ",
  sum(pvals_corn$sig_raw), " significant at 5\\% (unadjusted); ",
  sum(pvals_corn$sig_bh),  " significant at 5\\% FDR."
)

# Stage 1 lag comparison (confirmatory Z-vector spec pre-registration)
fml_corn_lag1      <- make_jp_formula("corn_yield", "soy_lag1", all_controls)
fml_corn_lag1_lag2 <- make_jp_formula("corn_yield", "soy_lag1 + soy_lag2", all_controls)
fml_corn_index     <- make_jp_formula("corn_yield", "rot_index", all_controls)

feols(fml_corn_lag1,      data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s1_lag1
feols(fml_corn_lag1_lag2, data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s1_lag1_lag2
feols(fml_corn_index,     data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s1_index

fml_corn_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)
feols(fml_corn_mean, data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s1

rot_names  <- grep("^rot_crop", names(coef(corn_jp_s1)), value = TRUE)
rot_labels <- chartr("15", "CS", sub("^rot_crop", "", rot_names))
rot_dict   <- setNames(rot_labels, rot_names)

# Table: tab:corn_jp_mean — Stage 1 corn yield mean
etable(corn_jp_s1, corn_jp_s1_lag1, corn_jp_s1_lag1_lag2, corn_jp_s1_index,
       tex      = TRUE,
       keep     = c("^[CS]-", "soy_lag", "rot_index"),
       dict     = rot_dict,
       notes    = bh_note_corn,
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Corn yield: full sequences vs lag summary variables",
       label    = "tab:corn_jp_mean",
       file     = paste0(tab_dir, "corn_jp_mean.tex"))

# Stage 2 OLS variance
corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_jp_data

fml_corn_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)
feols(fml_corn_var, data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s2

# Table: tab:corn_jp_var — Stage 2 corn yield variance
etable(corn_jp_s2,
       tex      = TRUE,
       keep     = "rot_crop",
       dict     = rot_dict,
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2 — Corn yield conditional variance (OLS)",
       label    = "tab:corn_jp_var",
       file     = paste0(tab_dir, "corn_jp_var.tex"))

# Summary table
corn_s1_coef <- broom::tidy(corn_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop", "", term),
            mean_est = estimate, mean_se = std.error, mean_p = p.value)

corn_s2_coef <- broom::tidy(corn_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop", "", term),
            var_est  = estimate, var_se  = std.error, var_p  = p.value)

corn_jp_summary <- corn_s1_coef |>
  left_join(corn_s2_coef, by = "rot_crop") |>
  mutate(
    mean_sig  = mean_p < 0.05,
    var_sig   = var_p  < 0.05,
    dominates = mean_est > 0 & var_est < 0,
    tradeoff  = mean_est > 0 & var_est > 0
  ) |>
  arrange(desc(dominates), desc(mean_est))

cat("Corn sequences dominating monoculture (higher mean + lower variance):\n")
corn_jp_summary |> filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |> print(n = Inf)

# Figure: corn_rot_plot (variance) — Response of std dev of corn yields to rotation sequences
corn_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "S-C-S-C-S-C" ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-C-S-C") ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  labs(title = "Response of standard deviation of corn yields to rotation sequences",
       caption = "Reference: corn monoculture (C-C-C-C-C-C). Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_var_plot
ggsave(paste0(fig_dir, "corn_var_plot.png"), corn_var_plot,  # fixed: was corn_rot_plot
       width = 10, height = 7.5, dpi = 300)

# Figure: corn_coeff_plot — Mean vs variance coefficients scatter
corn_jp_summary |>
  mutate(pr = ifelse(rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5"), 1, 0),
         pr = factor(pr)) |>
  ggplot(aes(x = var_est, y = mean_est, label = rot_crop)) +
  geom_jitter(aes(color = pr), size = 2) +
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  xlab("Corn standard deviation model coefficients") +
  ylab("Corn mean yield model coefficients") +
  labs(title = "Just-Pope decomposition: corn rotation effects on mean and variance",
       caption = "Reference: corn monoculture. Quadrant IV (right, above): dominates monoculture.") ->
  corn_coeff_plot
ggsave(paste0(fig_dir, "corn_coeff_plot.png"), corn_coeff_plot,
       width = 10, height = 10, dpi = 300)

# Figure: corn_jp_plot — mean vs variance with error bars (JP section figure)
corn_jp_summary |>
  mutate(
    type = case_when(
      dominates           ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0        ~ "Worse than monoculture",
      TRUE                ~ "No significant difference"),
    perfect = rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5")
  ) |>
  ggplot(aes(x = mean_est, y = var_est, colour = type, shape = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8) +
  geom_errorbar(aes(xmin = mean_est - 1.96 * mean_se,
                    xmax = mean_est + 1.96 * mean_se),
                orientation = "y", width = 0, alpha = 0.4) +
  geom_errorbar(aes(ymin = var_est - 1.96 * var_se,
                    ymax = var_est + 1.96 * var_se),
                width = 0, alpha = 0.4) +
  ggrepel::geom_label_repel(data = ~ filter(., perfect),
                             aes(label = rot_crop), size = 3,
                             show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Dominates monoculture"     = "#2ca02c",
    "Mean-variance trade-off"   = "#ff7f0e",
    "Worse than monoculture"    = "#d62728",
    "No significant difference" = "grey60")) +
  labs(x       = "Stage 1 coefficient: effect on mean corn yield (bu/acre)",
       y       = "Stage 2 coefficient: effect on yield variance",
       colour  = NULL, shape = NULL,
       title   = "Just-Pope decomposition: corn rotation effects on mean and variance",
       caption = "Reference: corn monoculture (C-C-C-C-C-C). Two-way clustering: field + year.\nQuadrant IV (bottom-right): higher mean, lower variance — unambiguously better.") +
  theme_bw() +
  theme(legend.position = "bottom") ->
  corn_jp_plot
ggsave(paste0(fig_dir, "corn_jp_plot.png"), corn_jp_plot,
       width = 9, height = 7, dpi = 300)

# ── 8. Just-Pope — soy ───────────────────────────────────────────────────────
# Tables: tab:soy_jp_mean, tab:soy_jp_var
# Figures: soy_jp_plot, soy_var_plot, soy_coeff_plot

fml_soy_mean <- make_jp_formula("soy_yield", "rot_crop", all_controls)
feols(fml_soy_mean, data = soy_jp_data, cluster = ~tile_field_ID + year) -> soy_jp_s1

# Table: tab:soy_jp_mean
etable(soy_jp_s1,
       tex      = TRUE,
       keep     = "rot_crop",
       dict     = setNames(chartr("15","CS", sub("^rot_crop","", grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE))),
                           grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE)),
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Soybean yield conditional mean",
       label    = "tab:soy_jp_mean",
       file     = paste0(tab_dir, "soy_jp_mean.tex"))

soy_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_jp_data

fml_soy_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)
feols(fml_soy_var, data = soy_jp_data, cluster = ~tile_field_ID + year) -> soy_jp_s2

# Table: tab:soy_jp_var
etable(soy_jp_s2,
       tex      = TRUE,
       keep     = "rot_crop",
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2 — Soybean yield conditional variance (OLS)",
       label    = "tab:soy_jp_var",
       file     = paste0(tab_dir, "soy_jp_var.tex"))

soy_s1_coef <- broom::tidy(soy_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term), mean_est=estimate, mean_se=std.error, mean_p=p.value)

soy_s2_coef <- broom::tidy(soy_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term), var_est=estimate, var_se=std.error, var_p=p.value)

soy_jp_summary <- soy_s1_coef |>
  left_join(soy_s2_coef, by = "rot_crop") |>
  mutate(dominates = mean_est > 0 & var_est < 0,
         tradeoff  = mean_est > 0 & var_est > 0,
         mean_sig  = mean_p < 0.05, var_sig = var_p < 0.05) |>
  arrange(desc(dominates), desc(mean_est))

# Figure: soy_var_plot — Response of std dev of soybean yields to rotation sequences
soy_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S",
                "S-S-S-S-C-S", "S-S-C-S-C-S") ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  labs(title = "Response of standard deviation of soybean yields to rotation sequences",
       caption = "Reference: soy monoculture (S-S-S-S-S-S). Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_var_plot
ggsave(paste0(fig_dir, "soy_var_plot.png"), soy_var_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: soy_coeff_plot
soy_jp_summary |>
  mutate(pr = ifelse(rot_crop %in% c("5-1-5-1-5-1","1-5-1-5-1-5"), 1, 0),
         pr = factor(pr)) |>
  ggplot(aes(x = var_est, y = mean_est, label = rot_crop)) +
  geom_jitter(aes(color = pr), size = 2) +
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  xlab("Soybeans standard deviation model coefficients") +
  ylab("Soybeans mean yield model coefficients") ->
  soy_coeff_plot
ggsave(paste0(fig_dir, "soy_coeff_plot.png"), soy_coeff_plot,
       width = 10, height = 10, dpi = 300)

# Figure: soy_jp_plot
soy_jp_summary |>
  mutate(
    type = case_when(
      dominates           ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0        ~ "Worse than monoculture",
      TRUE                ~ "No significant difference"),
    perfect = rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5")
  ) |>
  ggplot(aes(x = mean_est, y = var_est, colour = type, shape = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8) +
  geom_errorbar(aes(xmin = mean_est - 1.96 * mean_se,
                    xmax = mean_est + 1.96 * mean_se),
                orientation = "y", width = 0, alpha = 0.4) +
  geom_errorbar(aes(ymin = var_est - 1.96 * var_se,
                    ymax = var_est + 1.96 * var_se),
                width = 0, alpha = 0.4) +
  ggrepel::geom_label_repel(data = ~ filter(., perfect),
                             aes(label = rot_crop), size = 3,
                             show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Dominates monoculture"     = "#2ca02c",
    "Mean-variance trade-off"   = "#ff7f0e",
    "Worse than monoculture"    = "#d62728",
    "No significant difference" = "grey60")) +
  labs(x      = "Stage 1 coefficient: effect on mean soy yield (bu/acre)",
       y      = "Stage 2 coefficient: effect on yield variance",
       colour = NULL, shape = NULL,
       title  = "Just-Pope decomposition: soy rotation effects on mean and variance",
       caption = "Reference: soy monoculture (S-S-S-S-S-S). Two-way clustering: field + year.") +
  theme_bw() +
  theme(legend.position = "bottom") ->
  soy_jp_plot
ggsave(paste0(fig_dir, "soy_jp_plot.png"), soy_jp_plot,
       width = 9, height = 7, dpi = 300)

# ── 9. Just-Pope factor RCI ───────────────────────────────────────────────────
# Tables: tab:corn_rci_jp, tab:soy_rci_jp
# Figure: rci_plot

corn_rci_jp_data <- corn_jp_data |>
  select(-any_of(c("resid", "resid_sq", "h_hat"))) |>
  mutate(RCI = factor(RCI))

fml_corn_rci_mean <- make_jp_formula("corn_yield", "RCI", all_controls)
feols(fml_corn_rci_mean, data = corn_rci_jp_data, cluster = ~tile_field_ID + year) -> corn_rci_jp_s1

corn_rci_jp_s1 |>
  augment(newdata = corn_rci_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_rci_jp_data

fml_corn_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls)
feols(fml_corn_rci_var, data = corn_rci_jp_data, cluster = ~tile_field_ID + year) -> corn_rci_jp_s2

soy_rci_jp_data <- soy_jp_data |>
  select(-any_of(c("resid", "resid_sq", "h_hat"))) |>
  mutate(RCI = factor(RCI))

fml_soy_rci_mean <- make_jp_formula("soy_yield", "RCI", all_controls)
feols(fml_soy_rci_mean, data = soy_rci_jp_data, cluster = ~tile_field_ID + year) -> soy_rci_jp_s1

soy_rci_jp_s1 |>
  augment(newdata = soy_rci_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_rci_jp_data

fml_soy_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls)
feols(fml_soy_rci_var, data = soy_rci_jp_data, cluster = ~tile_field_ID + year) -> soy_rci_jp_s2

# Table: tab:corn_rci_jp — Just-Pope factor RCI, corn
etable(corn_rci_jp_s1, corn_rci_jp_s2,
       tex      = TRUE,
       keep     = "^RCI",
       dict     = dict_rci,
       headers  = c("Mean", "Variance"),
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Just-Pope: factor RCI and corn yield moments (corn-soy fields)",
       label    = "tab:corn_rci_jp",
       file     = paste0(tab_dir, "corn_rci_jp.tex"))

# Table: tab:soy_rci_jp — Just-Pope factor RCI, soy
etable(soy_rci_jp_s1, soy_rci_jp_s2,
       tex      = TRUE,
       keep     = "^RCI",
       dict     = dict_rci,
       headers  = c("Mean", "Variance"),
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Just-Pope: factor RCI and soy yield moments (corn-soy fields)",
       label    = "tab:soy_rci_jp",
       file     = paste0(tab_dir, "soy_rci_jp.tex"))

# Figure: rci_plot — Nonlinear effect of RCI on corn yield mean and variance
rci_plot_df <- bind_rows(
  broom::tidy(corn_rci_jp_s1) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI","",term)),
              est = estimate, se = std.error, moment = "Mean"),
  broom::tidy(corn_rci_jp_s2) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI","",term)),
              est = estimate, se = std.error, moment = "Variance")
)

ggplot(rci_plot_df, aes(x = rci, y = est, colour = moment, fill = moment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~moment, scales = "free_y") +
  scale_colour_manual(values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_fill_manual(  values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_x_continuous(breaks = sort(unique(rci_plot_df$rci))) +
  labs(x       = "Rotational Complexity Index (RCI)",
       y       = "Coefficient relative to RCI = 1.41",
       colour  = NULL, fill = NULL,
       title   = "Nonlinear effect of RCI on corn yield mean and variance",
       caption = "Corn-soy fields only. Two-way clustering: field + year.\nReference: RCI = 1.41 (near-monoculture).") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) ->
  rci_plot
ggsave(paste0(fig_dir, "rci_plot.png"), rci_plot,
       width = 9, height = 7, dpi = 300)

# ── 10. FGLS estimation + bootstrap ──────────────────────────────────────────
# Rebuild formulas after dropping collinear soil variables
all_controls <- setdiff(all_controls, c("nccpi3all_mean", "soc0_100_mean"))
fml_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)
fml_var  <- make_jp_formula("resid_sq",   "rot_crop", all_controls)

# Complete-case filter — ensures identical sample across all three stages
corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1","C",gsub("5","S",rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(setdiff(all_controls, c("I(pr_6^2)","I(pr_7^2)","I(pr_8^2)"))),
                ~ !is.na(.)))

cat("FGLS analysis sample:", nrow(corn_jp_data), "rows\n")

# Stage 1
feols(fml_mean, data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s1

corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_jp_data

# Stage 2a OLS
feols(fml_var, data = corn_jp_data, cluster = ~tile_field_ID + year) -> corn_jp_s2a

corn_jp_s2a |>
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |>
  select(-starts_with(".")) ->
  corn_jp_data

cat("Obs hitting h_hat floor:", sum(corn_jp_data$h_hat == 1e-6), "\n")

# Stage 2b FGLS
feols(fml_var, data = corn_jp_data,
      weights = ~I(1 / h_hat),
      cluster = ~tile_field_ID + year) -> corn_jp_s2b

etable(corn_jp_s1, corn_jp_s2a, corn_jp_s2b,
       keep  = "rot_crop",
       title = "Just-Pope FGLS: corn rotation effects on mean and variance")

# Bootstrap function
boot_jp_fgls <- function(dt, fml_mean, fml_var, B = 999, seed = 42) {

  dt <- as.data.table(dt)

  repeat {
    n_before <- nrow(dt)
    dt <- dt[dt[, .N, by = tile_field_ID][N > 1], on = "tile_field_ID"]
    dt <- dt[dt[, .N, by = year][N > 1],          on = "year"]
    if (nrow(dt) == n_before) break
  }
  cat("Rows after singleton removal:", nrow(dt), "\n")

  n_threads <- parallel::detectCores() - 1L
  set.seed(seed)
  fields   <- unique(dt$tile_field_ID)
  n_fields <- length(fields)

  fml_var_b <- make_jp_formula("resid_sq_b", "rot_crop", all_controls)

  s1 <- feols(fml_mean, data = dt, cluster = c("tile_field_ID","year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, resid_sq := NA_real_]
  dt[obs(s1), resid_sq := residuals(s1)^2]

  s2a <- feols(fml_var, data = dt, cluster = c("tile_field_ID","year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, h_hat := NA_real_]
  dt[obs(s2a), h_hat := pmax(fitted(s2a), 1e-6)]

  s2b <- feols(fml_var, data = dt, weights = ~I(1/h_hat),
               cluster = c("tile_field_ID","year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)

  coef_hat  <- coef(s2b)
  coef_boot <- matrix(NA, nrow = B, ncol = length(coef_hat),
                      dimnames = list(NULL, names(coef_hat)))

  for (b in seq_len(B)) {
    cat("\rBootstrap iteration", b, "of", B, "  ")

    drawn     <- table(sample(fields, n_fields, replace = TRUE))
    drawn_tbl <- data.table(tile_field_ID = names(drawn),
                            boot_w        = as.numeric(drawn))
    dt[drawn_tbl, boot_w := i.boot_w, on = "tile_field_ID"]
    dt[is.na(boot_w), boot_w := 0L]

    tryCatch({
      dt_b <- dt[boot_w > 0]

      b1 <- feols(fml_mean, data = dt_b, weights = ~boot_w,
                  vcov = "iid", nthreads = n_threads,
                  warn = FALSE, notes = FALSE)
      dt_b[, resid_sq_b := NA_real_]
      dt_b[obs(b1), resid_sq_b := residuals(b1)^2]

      b2a <- feols(fml_var_b, data = dt_b, weights = ~boot_w,
                   vcov = "iid", nthreads = n_threads,
                   warn = FALSE, notes = FALSE)
      dt_b[, h_hat_b   := NA_real_]
      dt_b[obs(b2a), h_hat_b := pmax(fitted(b2a), 1e-6)]
      dt_b[, combined_w := boot_w / h_hat_b]

      b2b <- feols(fml_var_b, data = dt_b, weights = ~combined_w,
                   vcov = "iid", nthreads = n_threads,
                   warn = FALSE, notes = FALSE)

      coef_boot[b, ] <- coef(b2b)

    }, error = function(e) {
      cat("Bootstrap iteration", b, "failed:", conditionMessage(e), "\n")
    })

    dt[, boot_w := NULL]

    if (b %% 100 == 0) {
      saveRDS(coef_boot, paste0("boot_progress_", b, ".rds"))
      gc()
    }
  }

  boot_se <- apply(coef_boot, 2, sd,       na.rm = TRUE)
  boot_ci <- apply(coef_boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  list(coef = coef_hat, se = boot_se, ci = boot_ci,
       fit_mean = s1, fit_var_ols = s2a, fit_var_fgls = s2b,
       boot_draws = coef_boot)
}

gc()
boot_fgls <- boot_jp_fgls(corn_jp_data, fml_mean, fml_var, B = 499, seed = 42)
saveRDS(boot_fgls, "C:/Users/vf006/Documents/boot_fgls.rds", compress = "xz")

# ── 11. Spatial maps ──────────────────────────────────────────────────────────
# Figures: rci_map, corn_yield_map, soy_yield_map, nccpi_corn_map, nccpi_soy_map

crs <- st_crs("EPSG:4326")
us_map(regions = "counties") |>
  filter(abbr == "IL") |>
  st_as_sf() |>
  st_transform(crs) ->
  il_map

# Use corn_df for spatial data (both files share same non-yield columns)
corn_df |>
  st_as_sf(coords = c("lon", "lat"), crs = crs) ->
  crop_sf

# Figure: rci_map — RCI values across Illinois (2016)
crop_sf |>
  filter(year == 2016 & !is.na(RCI)) |>
  mutate(RCI = factor(RCI)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = RCI), size = 0.2) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "bottom") +
  labs(title = "Rotational Complexity Index (2016)", caption = "Source: CDL / Socolar et al. (2021)") ->
  rci_map
ggsave(paste0(fig_dir, "rci_map.png"), rci_map,
       width = 10, height = 10, dpi = 300)

# Figure: corn_yield_map — Corn yields across Illinois (2016)
crop_sf |>
  filter(year == 2016 & !is.na(corn_yield)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = corn_yield), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title = "Corn yields — QDANN (2016)", caption = "Source: Ma et al. (2024)") ->
  corn_yield_map
ggsave(paste0(fig_dir, "corn_yield_map.png"), corn_yield_map,
       width = 10, height = 10, dpi = 300)

# Figure: soy_yield_map — Soybean yields across Illinois (2016)
soy_df |>
  st_as_sf(coords = c("lon", "lat"), crs = crs) |>
  filter(year == 2016 & !is.na(soy_yield)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = soy_yield), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title = "Soybean yields — QDANN (2016)", caption = "Source: Ma et al. (2024)") ->
  soy_yield_map
ggsave(paste0(fig_dir, "soy_yield_map.png"), soy_yield_map,
       width = 10, height = 10, dpi = 300)

# Figure: nccpi_corn_map — NCCPI corn productivity index (2016)
crop_sf |>
  filter(year == 2016 & !is.na(nccpi3corn_mean)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = nccpi3corn_mean), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title = "NCCPI — Corn (2016)", caption = "Source: gSSURGO") ->
  nccpi_corn_map
ggsave(paste0(fig_dir, "nccpi_corn_map.png"), nccpi_corn_map,
       width = 10, height = 10, dpi = 300)

# Figure: nccpi_soy_map — NCCPI soybean productivity index (2016)
crop_sf |>
  filter(year == 2016 & !is.na(nccpi3soy_mean)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = nccpi3soy_mean), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title = "NCCPI — Soybeans (2016)", caption = "Source: gSSURGO") ->
  nccpi_soy_map
ggsave(paste0(fig_dir, "nccpi_soy_map.png"), nccpi_soy_map,
       width = 10, height = 10, dpi = 300)

# ── 12. County-level yields and maximum entropy distribution ──────────────────
source("C:/Users/vf006/Box/crop_rotations_and_losses/code/maxent_tack2013.r")

corn_m <- corn_jp_data |>
  group_by(COUNTY_FIPS, STATE_FIPS, year) |>
  summarise(yield = mean(corn_yield, na.rm = TRUE)) |>
  ungroup() |>
  rename(state = STATE_FIPS, county = COUNTY_FIPS)

corn_det  <- detrend_yields_proportional(corn_m)
dist_m    <- tack_table2(corn_det$yield_norm)
