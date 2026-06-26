library(tidyverse)
library(data.table)
library(statar)
library(fixest)
library(broom)
library(haven)
library(marginaleffects)
library(furrr)
library(dotwhisker)
theme_set(theme_bw())

plan(multisession, workers = max(1, parallel::detectCores() - 1))

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
soil_controls <- c("nccpi3all_mean", "rootznaws_mean", "soc0_100_mean")
all_controls  <- c(weather_controls, soil_controls)

make_jp_formula <- function(lhs, rot_var, controls,
                            fe = "tile_field_ID + year") {
  rhs <- paste(c(rot_var, controls), collapse = " + ")
  as.formula(paste(lhs, "~", rhs, "|", fe))
}  

# ── Stage 1: conditional mean (corn) ─────────────────────────────────────────
# Identical specification to corn_rot above, but two-way clustered.
# rot_crop factor with corn monoculture (1-1-1-1-1-1) as reference.

corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(
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
  theme(legend.position = "bottom")   

## SOy data
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
  theme(legend.position = "bottom")      

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
  )

## Robustness: observation count reconciliation

cat("Observation counts after final crop_df cleaning:\n")
cat("  corn_jp_s1 :", formatC(nobs(corn_jp_s1),     big.mark = ","), "\n")
cat("  corn_rot   :", formatC(nobs(corn_rot),        big.mark = ","), "\n")
cat("  soy_jp_s1  :", formatC(nobs(soy_jp_s1),      big.mark = ","), "\n")
cat("  soy_rot    :", formatC(nobs(soy_rot),         big.mark = ","), "\n")

# Flag if corn or soy JP obs counts deviate from the corresponding mean models
# (they should match since the only difference is two-way vs county clustering)
if (nobs(corn_jp_s1) != nobs(corn_rot)) {
  warning("Corn JP and corn_rot observation counts differ — check crop_df state.")
}
if (nobs(soy_jp_s1) != nobs(soy_rot)) {
  warning("Soy JP and soy_rot observation counts differ — check crop_df state.")
}