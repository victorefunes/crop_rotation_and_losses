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

# Order columns
corn_df <- corn_df[order(tile_field_ID, year)]

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

# ── Formula helper ────────────────────────────────────────────────────────────
# Centralise the control vector so stage 1 and stage 2 always use
# identical RHS variables (avoids silent mismatches when copy-pasting).

weather_controls <- c(
  paste0("pr_",   6:8),
  paste0("cGDD_", 6:8, "m"),
  paste0("tmmx_", 6:8),
  paste0("tmmn_", 6:8),
  paste0("soil_", 6:8),
  paste0("vpd_",  6:8)
)

compute_gdd <- function(tmin, tmax, base = 8, cap = 29) {
  # Use sinusoidal approximation within day
  # Mean temperature within [base, cap] range
  pmax(0, pmin((tmax + tmin) / 2, cap) - base)
}

compute_edd <- function(tmin, tmax, threshold = 29) {
  # Degree days above threshold
  pmax(0, (tmax + tmin) / 2 - threshold)
}

corn_df |>
  mutate(EDD_6 = compute_edd(tmmx_6, tmmn_6),
         EDD_7 = compute_edd(tmmx_7, tmmn_7), 
         EDD_8 = compute_edd(tmmx_8, tmmn_8),
         GDD_6 = compute_gdd(tmmx_6, tmmn_6),
         GDD_7 = compute_gdd(tmmx_7, tmmn_7), 
         GDD_8 = compute_gdd(tmmx_8, tmmn_8)) ->
  corn_df

# Replace tmmx and tmmn with degree day variables (compute from GRIDMET daily)
weather_controls_revised <- c(
  # Precipitation: linear + quadratic
  "pr_6", "pr_7", "pr_8",
  "pr_6^2", "pr_7^2", "pr_8^2",

  # Beneficial heat: GDD 8-29C by month
  "cGDD_6m", "cGDD_7m", "cGDD_8m",

  # Harmful heat: EDD above 29C by month
  "EDD_6", "EDD_7", "EDD_8",

  # VPD: linear (quadratic for July optional)
  "vpd_6", "vpd_7", "vpd_8",

  # Soil moisture: keep monthly, drop PDSI
  "soil_6", "soil_7", "soil_8",

  # Cumulative GDD already covered by GDD_6/7/8 — drop cGDD to avoid collinearity
  # "cGDD_6m", "cGDD_7m", "cGDD_8m"  # redundant with GDD_6/7/8

  # Soil characteristics (identified with FIPS FEs)
  "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"
)

soil_controls <- c("nccpi3all_mean", "rootznaws_mean", "soc0_100_mean")
#all_controls  <- c(weather_controls, soil_controls)
all_controls  <- c(weather_controls_revised, "rootznaws_mean")

make_jp_formula <- function(lhs, rot_var, controls,
                            fe = "tile_field_ID + year") {
  rhs <- paste(c(rot_var, controls), collapse = " + ")
  as.formula(paste(lhs, "~", rhs, "|", fe))
}  

# ── Rotation-structure features ───────────────────────────────────────────────
# Precompute on the 64 unique patterns; join to data instead of row-wise sapply.
# The helper functions expect 1 = corn, 2 = soy; recode 5 → 2 at parse time.
count_soy_years <- function(rotation) {
  if (is.character(rotation))
    rotation <- as.integer(strsplit(rotation, "-")[[1]])
  sum(rotation %in% c(2L, 5L))
}

consecutive_soy <- function(rotation) {
  if (is.character(rotation))
    rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_runs <- rle(rotation)
  lengths  <- soy_runs$lengths[soy_runs$values %in% c(2L, 5L)]
  if (is_empty(lengths)) 0L else ifelse(max(lengths) > 1, 1L, 0L)
}

soy_gap <- function(rotation) {
  if (is.character(rotation))
    rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_years <- which(rotation %in% c(2L, 5L))
  if (length(soy_years) < 2) return(0)
  mean(diff(soy_years))
}

smallest_soy_gap <- function(rotation) {
  if (is.character(rotation))
    rotation <- as.integer(strsplit(rotation, "-")[[1]])
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

# PCA composite: re-orient all inputs so higher value = more rotation-intensive
# gap_soy / min_gap_soy return 0 when there is <=1 soy year (gap undefined);
# treat those as 0 frequency rather than NA so continuous corn scores lowest.
rot_pca_in <- rot_features |>
  transmute(
    n_soy       = n_soy,
    no_consec   = 1 - consec_soy,                                # blocking is bad
    freq_soy    = ifelse(gap_soy     == 0, 0, 1 / gap_soy),     # higher freq = better
    tight_soy   = ifelse(min_gap_soy == 0, 0, 1 / min_gap_soy)  # tighter min gap = better
  )

pca_rot <- prcomp(rot_pca_in, scale. = TRUE)
cat("Rotation PCA — variance explained by first two PCs:\n")
print(summary(pca_rot)$importance[, 1:2])
cat("PC1 loadings:\n")
print(pca_rot$rotation[, 1])

rot_features <- rot_features |>
  mutate(rot_index = pca_rot$x[, 1])

autoplot(pca_rot, data = corn_soy_patterns, loadings = TRUE, loadings.label = TRUE) +
  labs(title = "PCA of corn-soy rotation features") +
  theme_bw() +
  theme(legend.position = "none") ->
  rot_pca_plot  

ggsave("C:/Users/vf006/Box/crop_rotations_and_losses/figures/rot_pca_plot.png", rot_pca_plot,
       width = 9, height = 7, dpi = 300)

# ── Stage 1: conditional mean (corn) ─────────────────────────────────────────
# Identical specification to corn_rot above, but two-way clustered.
# rot_crop factor with corn monoculture (1-1-1-1-1-1) as reference.

corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  left_join(
    rot_features |> select(pattern, n_soy, consec_soy, gap_soy, min_gap_soy, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  mutate(
    # Binary soy-lag indicators parsed from "t5-t4-t3-t2-t1-t0" string
    soy_lag1 = as.integer(substr(rot_crop, 9, 9) == "5"),   # t-1
    soy_lag2 = as.integer(substr(rot_crop, 7, 7) == "5"),   # t-2
    soy_lag3 = as.integer(substr(rot_crop, 5, 5) == "5"),   # t-3
    soy_lag4 = as.integer(substr(rot_crop, 3, 3) == "5"),   # t-4
    soy_lag5 = as.integer(substr(rot_crop, 1, 1) == "5"),   # t-5
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "1-1-1-1-1-1")
  )

fml_corn_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)

corn_jp_s1 <- feols(
  fml_corn_mean,
  data    = corn_jp_data,
  cluster = ~tile_field_ID + year   # two-way: field + year
)

etable(corn_jp_s1, keep = "rot_crop",
       title = "Stage 1 — Corn yield conditional mean (two-way clustered)")

# ── Stage 1 (lag specifications): isolate Pattern 1 ──────────────────────────
# soy_lag1 alone captures Pattern 1 (immediate predecessor).
# Compare against soy_lag1 + soy_lag2 and the full 31-dummy spec:
# if within-R2 barely changes, Pattern 1 is the dominant source of variation.

fml_corn_lag1      <- make_jp_formula("corn_yield", "soy_lag1", all_controls)
fml_corn_lag1_lag2 <- make_jp_formula("corn_yield", "soy_lag1 + soy_lag2", all_controls)
fml_corn_index     <- make_jp_formula("corn_yield", "rot_index", all_controls)

corn_jp_s1_lag1 <- feols(
  fml_corn_lag1,
  data    = corn_jp_data,
  cluster = ~tile_field_ID + year
)

corn_jp_s1_lag1_lag2 <- feols(
  fml_corn_lag1_lag2,
  data    = corn_jp_data,
  cluster = ~tile_field_ID + year
)

corn_jp_s1_index <- feols(
  fml_corn_index,
  data    = corn_jp_data,
  cluster = ~tile_field_ID + year
)

rot_names  <- grep("^rot_crop", names(coef(corn_jp_s1)), value = TRUE)
rot_labels <- chartr("15", "CS", sub("^rot_crop", "", rot_names))
rot_dict   <- setNames(rot_labels, rot_names)

etable(
  corn_jp_s1, corn_jp_s1_lag1, corn_jp_s1_lag1_lag2, corn_jp_s1_index,
  keep    = c("rot_crop", "soy_lag", "rot_index"),
  dict    = rot_dict,
  title   = "Stage 1 — Corn yield: full sequences vs lag summary variables",
  file    = "tables/corn_stage1.tex",
  replace = TRUE
)


gc()

# ── Stage 2: conditional variance (corn) ──────────────────────────────────────
# Add squared residuals to the data, then regress on same RHS.
# Note: we attach residuals by row index — safe because feols drops NAs
# and singletons, so we match on the original row positions fixest returns.
corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid = corn_yield - .fitted) |> 
  select(tile_field_ID, year, resid) |>
  mutate(resid_sq = resid^2) ->
  corn_jp_resid

corn_jp_data <- corn_jp_data |>
  left_join(corn_jp_resid, by = c("tile_field_ID", "year"));rm(corn_jp_resid)

fml_corn_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)

corn_jp_s2 <- feols(
  fml_corn_var,
  data    = corn_jp_data,
  cluster = ~tile_field_ID + year
)

etable(corn_jp_s2, keep = "rot_crop",
       title = "Stage 2 — Corn yield conditional variance (two-way clustered)")       

# ── Summary: sequences that simultaneously raise mean AND lower variance ───────
# A rotation is unambiguously preferred to monoculture if its stage-1
# coefficient is positive (higher mean) AND its stage-2 coefficient is
# negative (lower variance). Flag these sequences for the rotation score.

corn_s1_coef <- broom::tidy(corn_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop  = gsub("rot_crop", "", term),
    mean_est  = estimate,
    mean_se   = std.error,
    mean_p    = p.value
  )

corn_s2_coef <- broom::tidy(corn_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop  = gsub("rot_crop", "", term),
    var_est   = estimate,
    var_se    = std.error,
    var_p     = p.value
  )

corn_jp_summary <- corn_s1_coef |>
  left_join(corn_s2_coef, by = "rot_crop") |>
  mutate(
    mean_sig  = mean_p < 0.05,
    var_sig   = var_p  < 0.05,
    # Unambiguously better than monoculture: higher mean AND lower variance
    dominates = mean_est > 0 & var_est < 0,
    # Risk-increasing: higher mean but also higher variance (trade-off)
    tradeoff  = mean_est > 0 & var_est > 0
  ) |>
  arrange(desc(dominates), desc(mean_est))

cat("Sequences that dominate corn monoculture (higher mean + lower variance):\n")
corn_jp_summary |>
  filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |>
  print(n = Inf)

cat("\nSequences with mean-variance trade-off (higher mean, higher variance):\n")
corn_jp_summary |>
  filter(tradeoff & mean_sig) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |>
  print(n = Inf)     

# ── Plot: mean vs variance coefficients, coloured by outcome type ─────────────
corn_jp_summary |>
  mutate(
    type = case_when(
      dominates          ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0       ~ "Worse than monoculture",
      TRUE               ~ "No significant difference"
    ),
    # Flag perfect rotation for labelling
    perfect = rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5")
  ) |>
  ggplot(aes(x = mean_est, y = var_est,
             colour = type, shape = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8) +
  geom_errorbar(aes(xmin = mean_est - 1.96 * mean_se,
                     xmax = mean_est + 1.96 * mean_se),
                     orientation = "y",
                 width = 0, alpha = 0.4) +
  geom_errorbar(aes(ymin = var_est - 1.96 * var_se,
                    ymax = var_est + 1.96 * var_se),
                width = 0, alpha = 0.4) +
  ggrepel::geom_label_repel(
    data = ~ filter(., perfect),
    aes(label = rot_crop),
    size = 3, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Dominates monoculture"       = "#2ca02c",
    "Mean-variance trade-off"     = "#ff7f0e",
    "Worse than monoculture"      = "#d62728",
    "No significant difference"   = "grey60"
  )) +
  labs(
    x      = "Stage 1 coefficient: effect on mean corn yield (bu/acre)",
    y      = "Stage 2 coefficient: effect on yield variance",
    colour = NULL, shape = NULL,
    title  = "Just-Pope decomposition: corn rotation effects on mean and variance",
    caption = "Reference: corn monoculture (1-1-1-1-1-1). Two-way clustering: field + year.\nQuadrant II (top-left): lower mean, higher variance — unambiguously worse.\nQuadrant IV (bottom-right): higher mean, lower variance — unambiguously better."
  ) +
  theme_bw() +
  theme(legend.position = "bottom") ->
  corn_jp_plot
ggsave("C:/Users/vf006/Box/crop_rotations_and_losses/figures/corn_jp_plot.png", corn_jp_plot,
       width = 9, height = 7, dpi = 300)    

## Soy data
rm(corn_df)

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

# ── Stage 1: conditional mean (soy) ──────────────────────────────────────────
# Reference: soy monoculture (5-5-5-5-5-5).

soy_jp_data <- soy_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "5-5-5-5-5-5")
  )

fml_soy_mean <- make_jp_formula("soy_yield", "rot_crop", all_controls)

soy_jp_s1 <- feols(
  fml_soy_mean,
  data    = soy_jp_data,
  cluster = ~tile_field_ID + year
)

etable(soy_jp_s1, keep = "rot_crop",
       title = "Stage 1 — Soy yield conditional mean (two-way clustered)")  

# ── Stage 2: conditional variance (soy) ──────────────────────────────────────
soy_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid = soy_yield - .fitted) |> 
  select(tile_field_ID, year, resid) |>
  mutate(resid_sq = resid^2) ->
  soy_jp_resid

soy_jp_data <- soy_jp_data |>
  left_join(soy_jp_resid, by = c("tile_field_ID", "year"));rm(soy_jp_resid)

fml_soy_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)

soy_jp_s2 <- feols(
  fml_soy_var,
  data    = soy_jp_data,
  cluster = ~tile_field_ID + year
)

etable(soy_jp_s2, keep = "rot_crop",
       title = "Stage 2 — Soy yield conditional variance (two-way clustered)")

soy_s1_coef <- broom::tidy(soy_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    mean_est = estimate,
    mean_se  = std.error,
    mean_p   = p.value
  )

soy_s2_coef <- broom::tidy(soy_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    var_est  = estimate,
    var_se   = std.error,
    var_p    = p.value
  )

soy_jp_summary <- soy_s1_coef |>
  left_join(soy_s2_coef, by = "rot_crop") |>
  mutate(
    dominates = mean_est > 0 & var_est < 0,
    tradeoff  = mean_est > 0 & var_est > 0,
    mean_sig  = mean_p < 0.05,
    var_sig   = var_p  < 0.05
  ) |>
  arrange(desc(dominates), desc(mean_est))

cat("Sequences that dominate soy monoculture (higher mean + lower variance):\n")
soy_jp_summary |>
  filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |>
  print(n = Inf)     

soy_jp_summary |>
  mutate(
    type = case_when(
      dominates          ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0       ~ "Worse than monoculture",
      TRUE               ~ "No significant difference"
    ),
    # Flag perfect rotation for labelling
    perfect = rot_crop %in% c("5-1-5-1-5-1", "1-5-1-5-1-5")
  ) |>
  ggplot(aes(x = mean_est, y = var_est,
             colour = type, shape = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3, alpha = 0.8) +
  geom_errorbar(aes(xmin = mean_est - 1.96 * mean_se,
                     xmax = mean_est + 1.96 * mean_se),
                     orientation = "y",
                 width = 0, alpha = 0.4) +
  geom_errorbar(aes(ymin = var_est - 1.96 * var_se,
                    ymax = var_est + 1.96 * var_se),
                width = 0, alpha = 0.4) +
  ggrepel::geom_label_repel(
    data = ~ filter(., perfect),
    aes(label = rot_crop),
    size = 3, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Dominates monoculture"       = "#2ca02c",
    "Mean-variance trade-off"     = "#ff7f0e",
    "Worse than monoculture"      = "#d62728",
    "No significant difference"   = "grey60"
  )) +
  labs(
    x      = "Stage 1 coefficient: effect on mean soy yield (bu/acre)",
    y      = "Stage 2 coefficient: effect on yield variance",
    colour = NULL, shape = NULL,
    title  = "Just-Pope decomposition: soy rotation effects on mean and variance",
    caption = "Reference: soy monoculture (5-5-5-5-5-5). Two-way clustering: field + year.\nQuadrant II (top-left): lower mean, higher variance — unambiguously worse.\nQuadrant IV (bottom-right): higher mean, lower variance — unambiguously better."
  ) +
  theme_bw() +
  theme(legend.position = "bottom") ->
  soy_jp_plot
ggsave("C:/Users/vf006/Box/crop_rotations_and_losses/figures/soy_jp_plot.png", soy_jp_plot,
       width = 9, height = 7, dpi = 300)    

# ── Motivation ────────────────────────────────────────────────────────────────
# The continuous RCI model above showed near-zero soy coefficients (0.004,
# 0.071) despite significant effects in the factor RCI model. This confirms
# the nonlinearity noted in the paper: a linear RCI rating approach would
# be misspecified. We therefore use factor RCI in both JP stages, which
# allows each RCI level to have its own mean and variance effect.

rm(soy_df)
# ── Corn: factor RCI ─────────────────────────────────────────────────────────
# Corn-soy fields only (consistent with corn_rci_cs above).
# Reference: corn monoculture maps to RCI = 0; or use lowest observed RCI.

corn_rci_jp_data <- corn_jp_data |>
  select(-c(resid, resid_sq)) |>
  mutate(RCI = factor(RCI))          # reference = lowest RCI level (RCI=1.41)

fml_corn_rci_mean <- make_jp_formula("corn_yield", "RCI", all_controls)

corn_rci_jp_s1 <- feols(
  fml_corn_rci_mean,
  data    = corn_rci_jp_data,
  cluster = ~tile_field_ID + year
 )

etable(corn_rci_jp_s1, keep = "^RCI",
       title = "Stage 1 — Corn yield conditional mean (two-way clustered)")   

corn_rci_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid = corn_yield - .fitted) |> 
  select(tile_field_ID, year, resid) |>
  mutate(resid_sq = resid^2) ->
  corn_rci_jp_resid

corn_rci_jp_data <- corn_rci_jp_data |>
  left_join(corn_rci_jp_resid, by = c("tile_field_ID", "year"));rm(corn_rci_jp_resid)

fml_corn_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls)

corn_rci_jp_s2 <- feols(
  fml_corn_rci_var,
  data    = corn_rci_jp_data,
  cluster = ~tile_field_ID + year
)

etable(
  corn_rci_jp_s1, corn_rci_jp_s2,
  keep  = "RCI",
  title = "Just-Pope: factor RCI and corn yield moments (corn-soy fields)"
)  

# ── Soy: factor RCI ──────────────────────────────────────────────────────────

soy_rci_jp_data <- soy_jp_data |>
  select(-c(resid, resid_sq)) |>
  mutate(RCI = factor(RCI))          # reference = lowest RCI level (RCI=1.41)

fml_soy_rci_mean <- make_jp_formula("soy_yield", "RCI", all_controls)

soy_rci_jp_s1 <- feols(
  fml_soy_rci_mean,
  data    = soy_rci_jp_data,
  cluster = ~tile_field_ID + year
)

soy_rci_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid = soy_yield - .fitted) |> 
  select(tile_field_ID, year, resid) |>
  mutate(resid_sq = resid^2) ->
  soy_rci_jp_resid

soy_rci_jp_data <- soy_rci_jp_data |>
  left_join(soy_rci_jp_resid, by = c("tile_field_ID", "year"));rm(soy_rci_jp_resid)

fml_soy_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls)

soy_rci_jp_s2 <- feols(
  fml_soy_rci_var,
  data    = soy_rci_jp_data,
  cluster = ~tile_field_ID + year
)

etable(
  soy_rci_jp_s1, soy_rci_jp_s2,
  keep  = "RCI",
  title = "Just-Pope: factor RCI and soy yield moments (corn-soy fields)"
)

# ── Figure: RCI nonlinearity in mean and variance ─────────────────────────────
# This is the key plot for the paper's claim that a linear RCI rating
# approach is misspecified. Show mean and variance coefficients as a
# function of RCI level for corn.

rci_levels <- levels(corn_rci_jp_data$RCI)[-1]   # drop reference (RCI=1.41)

rci_plot_df <- bind_rows(
  broom::tidy(corn_rci_jp_s1) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI", "", term)),
              est = estimate, se = std.error, moment = "Mean"),
  broom::tidy(corn_rci_jp_s2) |>
    filter(grepl("^RCI", term)) |>
    transmute(rci = as.numeric(gsub("RCI", "", term)),
              est = estimate, se = std.error, moment = "Variance")
)

ggplot(rci_plot_df, aes(x = rci, y = est, colour = moment, fill = moment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = est - 1.96 * se,
                  ymax = est + 1.96 * se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~moment, scales = "free_y") +
  scale_colour_manual(values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_fill_manual(  values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_x_continuous(breaks = sort(unique(rci_plot_df$rci))) +
  labs(
    x       = "Rotational Complexity Index (RCI)",
    y       = "Coefficient relative to RCI = 1.41",
    colour  = NULL, fill = NULL,
    title   = "Nonlinear effect of RCI on corn yield mean and variance",
    caption = "Corn-soy fields only. Two-way clustering: field + year.\nReference category: RCI = 1.41 (near-monoculture)."
  ) +
  theme_bw() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1)
  ) ->
  rci_plot
ggsave("C:/Users/vf006/Box/crop_rotations_and_losses/figures/rci_plot.png", rci_plot,
       width = 9, height = 7, dpi = 300)   

## Robustness: observation count reconciliation

cat("Observation counts after final crop_df cleaning:\n")
cat("  corn_jp_s1 :", formatC(nobs(corn_jp_s1),     big.mark = ","), "\n")
cat("  soy_jp_s1  :", formatC(nobs(soy_jp_s1),      big.mark = ","), "\n")



## FGLS

# ── Stage 1: conditional mean (identical to JP section) ──────────────────────
corn_jp_data <- corn_jp_data |>
  filter(!is.na(corn_yield))

fml_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)

corn_jp_s1 <- feols(fml_mean, data = corn_jp_data,
                    cluster = ~tile_field_ID + year)

# Attach residuals — safe row-aligned assignment via named index
corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid = corn_yield - .fitted) |> 
  select(tile_field_ID, year, resid) |>
  mutate(resid_sq = resid^2) ->
  corn_jp_resid

corn_jp_data <- corn_jp_data |>
  select(-c(resid, resid_sq)) |>
  left_join(corn_jp_resid, by = c("tile_field_ID", "year"));rm(corn_jp_resid)

# ── Stage 2a: OLS on squared residuals → variance function estimate ───────────
fml_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)

corn_jp_s2a <- feols(fml_var, data = corn_jp_data,
                     cluster = ~tile_field_ID + year)

# Fitted values from 2a are the estimated conditional variance h(x)
# Must be strictly positive for valid WLS weights — floor at a small positive
# value to handle any negative fitted values from the linear approximation
corn_jp_s2 |>
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |> # floor: weights must be finite
  select(tile_field_ID, year, h_hat) ->
  corn_jp_hat

corn_jp_data <- corn_jp_data |>
  left_join(corn_jp_hat, by = c("tile_field_ID", "year"));rm(corn_jp_hat) 

# ── Stage 2b: FGLS — WLS of û² on X, weighting by 1/h_hat ───────────────────
# feols accepts weights via the `weights` argument
corn_jp_s2b <- feols(fml_var, data = corn_jp_data,
                     weights = ~I(1 / h_hat),
                     cluster = ~tile_field_ID + year)       


# ── Filter to complete cases on all controls BEFORE stage 1 ──────────────────
vars_needed <- c("corn_yield", "rot_crop", "tile_field_ID", "year", all_controls)

corn_jp_data <- corn_jp_data |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls), ~ !is.na(.)))

cat("Rows after complete-case filter:", nrow(corn_jp_data), "\n")

corn_jp_s2a |>                               # was corn_jp_s2
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |>
  select(tile_field_ID, year, h_hat) ->
  corn_jp_hat

# ── Complete-case filter — identical sample across all three stages ────────────
corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls), ~ !is.na(.)))

cat("Analysis sample:", nrow(corn_jp_data), "rows\n")

# ── Stage 1 ───────────────────────────────────────────────────────────────────
fml_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)

corn_jp_s1 <- feols(fml_mean, data = corn_jp_data,
                    cluster = ~tile_field_ID + year)

corn_jp_data <- corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with("."))                    # drop .fitted, .resid etc.

# ── Stage 2a: OLS variance ────────────────────────────────────────────────────
fml_var <- make_jp_formula("resid_sq", "rot_crop", all_controls)

corn_jp_s2a <- feols(fml_var, data = corn_jp_data,
                     cluster = ~tile_field_ID + year)

corn_jp_data <- corn_jp_s2a |>
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6, na.rm = TRUE)) |>
  select(-starts_with("."))

cat("Obs hitting h_hat floor:", sum(corn_jp_data$h_hat == 1e-6), "\n")

# ── Stage 2b: FGLS ───────────────────────────────────────────────────────────
corn_jp_s2b <- feols(fml_var, data = corn_jp_data,
                     weights = ~I(1 / h_hat),
                     cluster = ~tile_field_ID + year)  
                     
etable(corn_jp_s1, corn_jp_s2a, corn_jp_s2b,
       keep = "rot_crop",
       title = "Just-Pope decomposition: corn rotation effects on mean and variance (FGLS)")

# ── Remove time-invariant soil vars collinear with field FEs ─────────────────
all_controls <- setdiff(all_controls, c("nccpi3all_mean", "soc0_100_mean"))

# Rebuild formulas with cleaned controls
fml_mean <- make_jp_formula("corn_yield", "rot_crop", all_controls)
fml_var  <- make_jp_formula("resid_sq",   "rot_crop", all_controls)       

boot_jp_fgls <- function(dt, fml_mean, fml_var, B = 999, seed = 42) {

  dt <- as.data.table(dt)

  # ── Remove singletons ────────────────────────────────────────────────────
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

  # ── Bootstrap variance formula uses resid_sq_b as LHS ───────────────────
  # fml_var has "resid_sq" as LHS — create a parallel formula for the boot
  fml_var_b <- make_jp_formula("resid_sq_b", "rot_crop", all_controls)

  # ── Point estimates on full data ──────────────────────────────────────────
  s1 <- feols(fml_mean, data = dt,
              cluster  = c("tile_field_ID", "year"),
              nthreads = n_threads,
              warn = FALSE, notes = FALSE)
  dt[, resid_sq := NA_real_]
  dt[obs(s1), resid_sq := residuals(s1)^2]

  s2a <- feols(fml_var, data = dt,
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads,
               warn = FALSE, notes = FALSE)
  dt[, h_hat := NA_real_]
  dt[obs(s2a), h_hat := pmax(fitted(s2a), 1e-6)]

  s2b <- feols(fml_var, data = dt,
               weights  = ~I(1 / h_hat),
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads,
               warn = FALSE, notes = FALSE)

  coef_hat  <- coef(s2b)
  coef_boot <- matrix(NA, nrow = B, ncol = length(coef_hat),
                      dimnames = list(NULL, names(coef_hat)))

  # ── Bootstrap ─────────────────────────────────────────────────────────────
  for (b in seq_len(B)) {

    cat("\rBootstrap iteration", b, "of", B, "  ")

    # Field-level frequency weights — no row duplication
    drawn     <- table(sample(fields, n_fields, replace = TRUE))
    drawn_tbl <- data.table(
      tile_field_ID = names(drawn),
      boot_w        = as.numeric(drawn)
    )
    dt[drawn_tbl, boot_w := i.boot_w, on = "tile_field_ID"]
    dt[is.na(boot_w), boot_w := 0L]

    tryCatch({

      # Work on the active subset only — avoids repeated subsetting in feols
      dt_b <- dt[boot_w > 0]

      # Stage 1
      b1 <- feols(fml_mean, data = dt_b,
                  weights  = ~boot_w,
                  vcov     = "iid",
                  nthreads = n_threads,
                  warn = FALSE, notes = FALSE)
      dt_b[, resid_sq_b := NA_real_]
      dt_b[obs(b1), resid_sq_b := residuals(b1)^2]

      # Stage 2a
      b2a <- feols(fml_var_b, data = dt_b,
                   weights  = ~boot_w,
                   vcov     = "iid",
                   nthreads = n_threads,
                   warn = FALSE, notes = FALSE)
      dt_b[, h_hat_b    := NA_real_]
      dt_b[obs(b2a), h_hat_b := pmax(fitted(b2a), 1e-6)]
      dt_b[, combined_w := boot_w / h_hat_b]

      # Stage 2b — only coef() is used; VCOV not needed inside the loop
      b2b <- feols(fml_var_b, data = dt_b,
                   weights  = ~combined_w,
                   vcov     = "iid",
                   nthreads = n_threads,
                   warn = FALSE, notes = FALSE)

      coef_boot[b, ] <- coef(b2b)

    }, error = function(e) {
      cat("Bootstrap iteration", b, "failed:", conditionMessage(e), "\n")
    })

    # Reset boot_w — dt_b is a copy so no cleanup needed there
    dt[, boot_w := NULL]

    if (b %% 100 == 0) {
      cat("Bootstrap iteration", b, "of", B, "\n")
      saveRDS(coef_boot, paste0("boot_progress_", b, ".rds"))
      gc()
    }
  }

  boot_se <- apply(coef_boot, 2, sd,       na.rm = TRUE)
  boot_ci <- apply(coef_boot, 2, quantile,
                   probs = c(0.025, 0.975), na.rm = TRUE)

  list(coef         = coef_hat,
       se           = boot_se,
       ci           = boot_ci,
       fit_mean     = s1,
       fit_var_ols  = s2a,
       fit_var_fgls = s2b,
       boot_draws   = coef_boot)
}

gc()
boot_fgls <- boot_jp_fgls(corn_jp_data, fml_mean, fml_var, B = 499, seed = 42)
boot_fgls |>
  saveRDS("C:/Users/vf006/Documents//boot_fgls.rds", compress = "xz")

## County-level yields
source("C:/Users/vf006/Box/crop_rotations_and_losses/code/maxent_tack2013.r")

corn_m <- corn_jp_data |>
  group_by(COUNTY_FIPS, STATE_FIPS, year) |>
  summarise(yield = mean(corn_yield, na.rm = TRUE)) |>
  ungroup() |>
  rename(state = STATE_FIPS, 
        county = COUNTY_FIPS)

corn_det <- detrend_yields_proportional(corn_m)
dist_m <- tack_table2(corn_det$yield_norm)