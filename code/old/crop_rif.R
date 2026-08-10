library(tidyverse)
library(data.table)
library(statar)
library(fixest)
library(broom)
library(haven)
library(marginaleffects)
library(furrr)
library(dotwhisker)
library(dineq)
library(rifreg)
library(TTR)
theme_set(theme_bw())

plan(multisession, workers = max(1, parallel::detectCores() - 1))

setwd("C:/Users/vf006/Box/Economic Analysis of Soil Health Practices")

tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"
fig_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/figures/"

crop_df <- fread("./Data and Data Descriptions/clean/data_il_long.csv")

# Order columns
crop_df <- crop_df[order(tile_field_ID, year)]

# List of all corn and soy patterns
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

# RCI values correction
crop_df |>
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
  crop_df

prob <- seq(0.1, 0.9, 0.1)

crop_df |>
  filter(!is.na(corn_yield)) |>
  group_by(year) |>
  select(tile_field_ID, year, corn_yield) |>
  mutate(rifs = get_rif(dep_var = corn_yield, statistic = "quantiles", 
                        probs = prob) |> as_tibble()) |> 
  unnest(cols = rifs, names_sep = "_") |>
  select(-c(corn_yield, rifs_weights)) |>
  ungroup() ->
  corn_df

crop_df |>
  filter(!is.na(soy_yield)) |>
  group_by(year) |>
  select(tile_field_ID, year, soy_yield) |>
  mutate(rifs = get_rif(dep_var = soy_yield, statistic = "quantiles", 
                        probs = prob) |> as_tibble()) |> 
  unnest(cols = rifs, names_sep = "_") |>
  select(-c(soy_yield, rifs_weights)) |>
  ungroup() ->
  soy_df


names(corn_df)[3:11] <- paste0("corn_rif_", prob)
names(soy_df)[3:11] <- paste0("soy_rif_", prob)

crop_df |>
  left_join(corn_df, by = c("tile_field_ID", "year")) |>
  left_join(soy_df, by = c("tile_field_ID", "year")) ->
  crop_df

rm(corn_df, soy_df)

## RCI regressions
corn_RCI_rif <- paste("c(", paste("corn_rif_", prob, collapse = ",", sep = ""), 
                      ")~RCI+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

dict <- c("RCI1.41" = "RCI = 1.41", "RCI1.73" = "RCI = 1.73",
          "RCI2" = "RCI = 2", "RCI2.24" = "RCI = 2.24",
          "RCI2.45" = "RCI = 2.45", "RCI2.65" = "RCI = 2.65",
          "RCI2.74" = "RCI = 2.74", "RCI2.83" = "RCI = 2.83",
          "RCI3" = "RCI = 3","RCI3.24" = "RCI = 3.24", 
          "RCI3.46" = "RCI = 3.46", "RCI3.67" = "RCI = 3.67", 
          "RCI3.74" = "RCI = 3.74", "RCI4" = "RCI = 4",
          "RCI4.24" = "RCI = 4.24", "RCI4.47" = "RCI = 4.47",
          "RCI4.74" = "RCI = 4.74", "RCI5.2" = "RCI = 5.2")

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_RCI_rif, data = _, cluster = ~COUNTY_FIPS) ->
  corn_RCI_rif
etable(corn_RCI_rif, keep = "RCI",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

soy_RCI_rif <- paste("c(", paste("soy_rif_", prob, collapse = ",", sep = ""), 
                      ")~RCI+", 
                      paste("pr_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                      "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                      "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                      "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_RCI_rif, data = _, cluster = ~COUNTY_FIPS) ->
  soy_RCI_rif
etable(soy_RCI_rif, keep = "RCI",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

## Mean regressions
corn_RCI_formula <- paste("corn_yield~RCI+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_mean

soy_RCI_formula <- paste("soy_yield~RCI+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_mean
etable(corn_rci_mean, soy_rci_mean, keep = "RCI",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

coef(corn_RCI_rif, long = TRUE) |>
  filter(str_detect(coefficient, "RCI")) |>
  separate(lhs, into = c("temp", "quantile"), sep = 9) |>
  select(-temp) |>
  mutate(quantile = as.numeric(quantile)) |>
  separate(coefficient, into = c("temp", "RCI"), sep = 3) |>
  select(-c(id, temp)) ->
  coefs

coefs |>
  ggplot(aes(x = RCI, y = estimate)) +
  geom_point(aes(color = quantile)) +
  geom_line(aes(color = quantile), group = 1) + 
  facet_wrap(quantile~., scales = "fixed")

crop_df |>
  filter(!is.na(corn_yield)) |>
  group_by(year) |>
  select(tile_field_ID, year, corn_yield) |>
  mutate(rifs = get_rif(dep_var = corn_yield, statistic = "variance") |> 
           as_tibble()) |> 
  unnest(cols = rifs, names_sep = "_") |>
  select(-c(corn_yield, rifs_weights)) |>
  rename(rif_var_corn = rifs_rif_variance) |>
  ungroup() ->
  corn_df

crop_df |>
  filter(!is.na(soy_yield)) |>
  group_by(year) |>
  select(tile_field_ID, year, soy_yield) |>
  mutate(rifs = get_rif(dep_var = soy_yield, statistic = "variance") |> 
           as_tibble()) |> 
  unnest(cols = rifs, names_sep = "_") |>
  select(-c(soy_yield, rifs_weights)) |>
  rename(rif_var_soy = rifs_rif_variance) |>
  ungroup() ->
  soy_df

crop_df |>
  left_join(corn_df, by = c("tile_field_ID", "year")) |>
  left_join(soy_df, by = c("tile_field_ID", "year")) ->
  crop_df

rm(corn_df, soy_df)

corn_RCI_var <- paste("rif_var_corn~RCI+", 
                      paste("pr_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                      "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                      "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                      "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_RCI_var, data = _, cluster = ~COUNTY_FIPS) ->
  corn_RCI_var


soy_RCI_var <- paste("rif_var_soy~RCI+", 
                      paste("pr_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                      "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                      "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                      "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_RCI_var, data = _, cluster = ~COUNTY_FIPS) ->
  soy_RCI_var
etable(corn_RCI_var, soy_RCI_var, keep = "RCI",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

coefplot(corn_RCI_var, keep = "RCI")
coefplot(soy_RCI_var, keep = "RCI")

# Perfect rotation vs. monoculture
crop_df |>
  mutate(perfect_rot = case_when(
    rot_crop == "5-1-5-1-5-1" ~ 1,
    rot_crop == "1-5-1-5-1-5" ~ 1,
    rot_crop == "1-1-1-1-1-1" ~ 0,
    rot_crop == "5-5-5-5-5-5" ~ 0,
    .default = NA)) ->
  crop_df

corn_rot <- paste("corn_yield~perfect_rot+", 
                      paste("pr_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                      "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                      "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                      "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(!is.na(perfect_rot)) |>
  feols(corn_rot, data = _, cluster = ~COUNTY_FIPS) ->
  corn_prot

soy_rot <- paste("soy_yield~perfect_rot+", 
                  paste("pr_", 6:8, collapse = "+", sep = ""), 
                  "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                  "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                  "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                  "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                  "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                  "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(!is.na(perfect_rot)) |>
  feols(soy_rot, data = _, cluster = ~COUNTY_FIPS) ->
  soy_prot

etable(corn_prot, soy_prot, keep = "perfect_rot",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

corn_prot_rif_fml <- paste("c(", paste("corn_rif_", prob, collapse = ",", sep = ""), 
                       ")~perfect_rot+", 
                       paste("pr_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                       "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                       "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                       "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(!is.na(perfect_rot)) |>
  feols(corn_prot_rif_fml, data = _, cluster = ~COUNTY_FIPS) ->
  corn_prot_rif

soy_prot_rif_fml <- paste("c(", paste("soy_rif_", prob, collapse = ",", sep = ""), 
                       ")~perfect_rot+", 
                       paste("pr_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                       "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                       "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                       "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(!is.na(perfect_rot)) |>
  feols(soy_prot_rif_fml, data = _, cluster = ~COUNTY_FIPS) ->
  soy_prot_rif

etable(corn_prot_rif, soy_prot_rif, keep = "perfect_rot",
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)

confint(corn_prot_rif, long = TRUE, level = 0.99) |>
  filter(coefficient == "perfect_rot") |>
  separate(lhs, into = c("temp", "quantile"), sep = 9) |>
  select(-c(temp, coefficient)) |>
  data.frame() |>
  rename(ci_low = `X0.5..`,
         ci_high = `X99.5..`) ->
  corn_prot_rif_ci

coef(corn_prot_rif, long = TRUE) |>
  filter(coefficient == "perfect_rot") |>
  separate(lhs, into = c("temp", "quantile"), sep = 9) |>
  mutate(quantile = factor(quantile)) |>
  select(-c(id, temp)) |>
  left_join(corn_prot_rif_ci, by = "quantile") |>
  ggplot(aes(x = quantile, y = estimate)) +
  geom_hline(yintercept = 3.825, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_ribbon(aes(x = quantile, ymin = ci_low, ymax = ci_high), 
              fill = "lightblue", alpha = 0.5, group = 1) +
  geom_line(group = 1) 

confint(soy_prot_rif, long = TRUE, level = 0.99) |>
  filter(coefficient == "perfect_rot") |>
  separate(lhs, into = c("temp", "quantile"), sep = 9) |>
  select(-c(temp, coefficient)) |>
  data.frame() |>
  rename(ci_low = `X0.5..`,
         ci_high = `X99.5..`) ->
  soy_prot_rif_ci

coef(soy_prot_rif, long = TRUE) |>
  filter(coefficient == "perfect_rot") |>
  separate(lhs, into = c("temp", "quantile"), sep = 9) |>
  mutate(quantile = factor(quantile)) |>
  select(-c(id, temp)) |>
  left_join(soy_prot_rif_ci, by = "quantile") |>
  ggplot(aes(x = quantile, y = estimate)) +
  geom_hline(yintercept = 0.969, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_ribbon(aes(x = quantile, ymin = ci_low, ymax = ci_high), 
              fill = "lightblue", alpha = 0.5, group = 1) +
  geom_line(group = 1) 