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

# Transition out of corn/soy monoculture
# C-C-C-C-C-C: 0.00
# C-C-C-C-C-S: 1.41
# C-C-C-C-S-C: 1.73
# C-C-C-S-C-S: 2.00
# C-C-S-C-S-C: 2.24
# C-S-C-S-C-S: 2.24

# S-S-S-S-S-S: 0.00
# S-S-S-S-S-C: 1.41
# S-S-S-S-C-S: 1.73
# S-S-S-C-S-C: 2.00
# S-S-C-S-C-S: 2.24
# S-C-S-C-S-C: 2.24

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

crop_df |> 
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  tab(rot_crop) |> 
  data.frame() |>
  arrange(desc(Freq.)) |>
  select(-`Cum.`) |>
  mutate(cum = cumsum(Percent))

## Corn rotations using corn monoculture as the base category
corn_yield_formula <- paste("corn_yield ~ rot_crop+", 
                            paste("pr_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                            "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                            "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                            "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_yield~rot_crop|tile_field_ID+year, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot_nc

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot

dict_corn <- c("rot_cropC-C-C-C-S-C" = "C-C-C-C-S-C", "rot_cropC-C-C-S-C-C" = "C-C-C-S-C-C",
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
          "rot_cropS-S-S-S-S-C" = "S-S-S-S-S-C")

etable(corn_rot_nc, corn_rot, dict = dict_corn,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       extralines = list("_Controls" = c("No", "Yes")))
etable(corn_rot_nc, corn_rot, 
       tex = TRUE,
       dict = dict_corn,
       headers = c("No controls", "With controls"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       se.below = FALSE,
       title = "Rotation patterns and corn yields",
       label = "tab:table1",
       extralines = list("_Controls" = c("No", "Yes")),
       file = paste0(tab_dir, "corn_rot.tex"))


corn_rot |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |> 
  separate(term, into = c("temp", "term"), sep = 8) |>
  arrange(estimate) |>
  dwplot(style = "dotwhisker",
         ci = 0.99, 
         dodge_size = 0.4,
         dist_args = list(alpha = 0.75),
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  aes(color = group) +
  scale_color_discrete() +
  theme_bw() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  corn_rot_plot
corn_rot_plot
#ggsave(corn_rot_plot, 
#        filename = paste0(fig_dir, "corn_rot_plot.png"), 
#        width = 10, height = 7.5)

corn_rot |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-C-S-C-S-C" ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C", "S-S-S-S-S-C",
                "S-S-S-C-S-C") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  corn_rot_plot
corn_rot_plot
ggsave(corn_rot_plot, 
       filename = paste0(fig_dir, "corn_rot_plot.png"), 
       width = 10, height = 7.5)


## Soybean rotations 
soy_yield_formula <- paste("soy_yield ~ rot_crop+", 
                           paste("pr_", 6:8, collapse = "+", sep = ""), 
                           "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                           "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                           "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                           "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                           "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                           "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()


dict_soy <- c("rot_cropC-C-C-C-C-S" = "C-C-C-C-C-S", "rot_cropC-C-C-C-S-S" = "C-C-C-C-S-S",
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
          "rot_cropS-S-S-S-C-S" = "S-S-S-S-C-S")

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")) |>
  feols(soy_yield~rot_crop|tile_field_ID+year, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rot_nc

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")) |>
  feols(soy_yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rot
etable(soy_rot_nc, soy_rot, dict = dict_soy,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       extralines = list("_Controls" = c("No", "Yes")))
etable(soy_rot_nc, soy_rot,
       tex = TRUE,
       dict = dict_soy,
       headers = c("No controls", "With controls"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       se.below = FALSE,
       title = "Rotation patterns and soybean yields",
       label = "tab:table2",
       extralines = list("_Controls" = c("No", "Yes")),
       file = paste0(tab_dir, "soy_rot.tex"))

soy_rot |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |> 
  separate(term, into = c("temp", "term"), sep = 8) |>
  arrange(estimate) |>
  dwplot(style = "dotwhisker",
         dist_args = list(color = "black", alpha = 0.75),
         ci = 0.99,
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  theme_bw() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  soy_rot_plot
soy_rot_plot
#ggsave(soy_rot_plot, 
#       filename = paste0(fig_dir, "soy_rot_plot.png"), 
#       width = 10, height = 7.5)

soy_rot |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S", "S-S-S-S-C-S",
                "S-S-C-S-C-S") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  soy_rot_plot
soy_rot_plot
ggsave(soy_rot_plot, 
       filename = paste0(fig_dir, "soy_rot_plot.png"), 
       width = 10, height = 7.5)


## RCI regressions
corn_RCI_formula <- paste("corn_yield~RCI+", 
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
  feols(corn_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_all

crop_df |>
  mutate(RCI = factor(RCI)) |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  feols(corn_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_cs

etable(corn_rci_all, corn_rci_cs, keep = "RCI",
       headers = c("All crops", "Corn and soybeans only"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)
etable(corn_rci_all, corn_rci_cs, 
       tex = TRUE,
       dict = dict,
       headers = c("All crops", "Corn and soybeans only"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       se.below = FALSE,
       title = "Rotation Complexity and corn yields",
       label = "tab:table3",
       extralines = list("_Controls" = c("No", "Yes")),
       file = paste0(tab_dir, "corn_rci.tex"))


corn_rci_all |> 
  tidy() |> 
  filter(grepl("RCI", term)) |> 
  separate(term, into = c("temp", "term"), sep = 3) |>
  mutate(term = as.numeric(term))|>
  filter(term < 5) |>
  arrange(desc(term)) |>
  dwplot(style = "dotwhisker",
         dist_args = list(color = "black", alpha = 0.75),
         ci = 0.99,
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("Rotation Complexity Index") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  corn_rci_plot
corn_rci_plot
#ggsave(corn_rci_plot, 
#       filename = paste0(fig_dir, "corn_rci_plot.png"), 
#       width = 10, height = 7.5)

corn_rci_all |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  filter(term < 5.2) |>
  mutate(group = case_when(
    term == 2.24 ~ "Perfect rotation",
    term %in% c(1.41, 1.73, 2) ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  corn_rci_plot
corn_rci_plot
ggsave(corn_rci_plot, 
       filename = paste0(fig_dir, "corn_rci_plot.png"), 
       width = 10, height = 7.5)

soy_RCI_formula <- paste("soy_yield~RCI+", 
                         paste("pr_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                         "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                         "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_all

crop_df |>
  mutate(RCI = factor(RCI)) |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  feols(soy_RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_cs

etable(soy_rci_all, soy_rci_cs, keep = "RCI",
       headers = c("All crops", "Corn and soybeans only"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       dict = dict)
etable(soy_rci_all, soy_rci_cs, 
       tex = TRUE,
       dict = dict,
       headers = c("All crops", "Corn and soybeans only"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       se.below = FALSE,
       title = "Rotation Complexity and soybean yields",
       label = "tab:table4",
       extralines = list("_Controls" = c("No", "Yes")),
       file = paste0(tab_dir, "soy_rci.tex"))

soy_rci_all |> 
  tidy() |> 
  filter(grepl("RCI", term)) |> 
  separate(term, into = c("temp", "term"), sep = 3) |>
  mutate(term = as.numeric(term))|>
  arrange(term) |>
  filter(term < 5) |>
  arrange(desc(term)) |>
  dwplot(style = "dotwhisker",
         dist_args = list(color = "black", alpha = 0.75),
         ci = 0.99,
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("Rotation Complexity Index") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  soy_rci_plot
soy_rci_plot
#ggsave(soy_rci_plot, 
#       filename = paste0(fig_dir, "soy_rci_plot.png"), 
#       width = 10, height = 7.5)

soy_rci_all |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  filter(term < 5.2) |>
  mutate(group = case_when(
    term == 2.24 ~ "Perfect rotation",
    term %in% c(1.41, 1.73, 2) ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  soy_rci_plot
soy_rci_plot
ggsave(soy_rci_plot, 
       filename = paste0(fig_dir, "soy_rci_plot.png"), 
       width = 10, height = 7.5)


## VPDmax_7 as a category:
# Normal       0<=vpdmax_7<1.9
# somewhat dry 1.9<=vpdmax_7<=2.1
# dry          2.1<vpdmax_7<= 4.6 

crop_df |>
  mutate(vpd_name = case_when(
    vpdmax_7 >=0 & vpdmax_7 < 1.9 ~ 1,
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ 2,
    vpdmax_7 > 2.1 ~ 3,
    .default = NA),
    vpd_name = factor(vpd_name, labels = c("normal", "somewhat dry", "dry"))) ->
  crop_df

corn_vpd_formula <- paste("corn_yield ~ rot_crop+vpd_name+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot_vpd
etable(corn_rot_vpd,
       dict = c(dict_corn, "vpd_namesomewhatdry" = "Somewhat dry season",
                "vpd_namedry" = "Dry season"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))
etable(corn_rot_vpd, 
       tex = TRUE,
       dict = c(dict_corn, "vpd_namesomewhatdry" = "Somewhat dry season",
                "vpd_namedry" = "Dry season"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       title = "Effect of weather and rotation sequences on corn yields",
       label = "tab:table5",
       file = paste0(tab_dir, "corn_rot_vpd.tex"))

soy_vpd_formula <- paste("soy_yield ~ rot_crop+vpd_name+", 
                         paste("pr_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                         "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                         "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                         "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")) |>
  feols(soy_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rot_vpd
etable(soy_rot_vpd,
       dict = c(dict_soy, "vpd_namesomewhatdry" = "Somewhat dry season",
                "vpd_namedry" = "Dry season"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))
etable(soy_rot_vpd, 
       tex = TRUE,
       dict = c(dict_soy, "vpd_namesomewhatdry" = "Somewhat dry season",
                "vpd_namedry" = "Dry season"),
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       title = "Effect of weather and rotation sequences on corn yields",
       label = "tab:table6",
       file = paste0(tab_dir, "soy_rot_vpd.tex"))


corn_rci_vpd_formula <- paste("corn_yield ~ RCI*vpd_name+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

dict_vpd <- c("vpd_namesomewhatdry" = "Somewhat dry season",
              "vpd_namedry" = "Dry season", 
              "RCI x vpd_namesomewhatdry" = "RCI x Somewhat dry season",
              "RCI x vpd_namedry" = "RCI x Dry season")

crop_df |>
  feols(corn_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rci_vpd
etable(corn_rci_vpd, 
       dict = dict_vpd,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))
etable(corn_rci_vpd, 
       tex = TRUE,
       dict = dict_vpd,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       title = "Effect of weather and rotation sequences on corn yields",
       label = "tab:table7",
       file = paste0(tab_dir, "corn_rci_vpd.tex"))

soy_rci_vpd_formula <- paste("soy_yield ~ RCI*vpd_name+", 
                              paste("pr_", 6:8, collapse = "+", sep = ""), 
                              "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                              "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                              "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                              "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                              "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                              "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  feols(soy_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rci_vpd
etable(soy_rci_vpd, 
       dict = dict_vpd,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))
etable(soy_rci_vpd, 
       tex = TRUE,
       dict = dict_vpd,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"),
       style.tex = style.tex("aer"),
       replace = TRUE,
       title = "Effect of weather and rotation sequences on corn yields",
       label = "tab:table8",
       file = paste0(tab_dir, "soy_rci_vpd.tex"))

## Conditional standard error regressions
corn_rot |>
  augment(newdata = crop_df |> 
            filter(rot_crop %in% corn_soy_patterns$pattern) |>
            mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)))) |>
  rename(corn_yield_pred = .fitted) |>
  mutate(corn_yield_res = corn_yield-corn_yield_pred,
         scaled_res = scale(corn_yield_res)) ->
  corn_res

soy_rot |>
  augment(newdata = crop_df |> 
            filter(rot_crop %in% corn_soy_patterns$pattern) |>
            mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)))) |>
  rename(soy_yield_pred = .fitted) |>
  mutate(soy_yield_res = soy_yield-soy_yield_pred,
         scaled_res = scale(soy_yield_res)) ->
  soy_res

corn_var_formula <- paste("res_var ~ rot_crop+", 
                          paste("pr_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                          "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                          "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                          "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

corn_res |>
  mutate(rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C"),
         res_var = sqrt(corn_yield_res^2)) |>
  feols(corn_var_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_var
etable(corn_var, keep = "rot_crop")

corn_var |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |> 
  separate(term, into = c("temp", "term"), sep = 8) |>
  arrange(estimate) |>
  dwplot(style = "dotwhisker",
         dist_args = list(color = "black", alpha = 0.75),
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  theme_bw() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  corn_var_plot
corn_var_plot
#ggsave(corn_var_plot, 
#       filename = paste0(fig_dir, "corn_var_plot.png"), 
#       width = 10, height = 7.5)

corn_var |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-C-S-C-S-C" ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C", "S-S-S-S-S-C",
                "S-S-S-C-S-C") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  corn_var_plot
corn_var_plot
ggsave(corn_rot_plot, 
       filename = paste0(fig_dir, "corn_var_plot.png"), 
       width = 10, height = 7.5)


soy_var_formula <- paste("res_var ~ rot_crop+", 
                         paste("pr_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                         "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                         "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

soy_res |>
  mutate(rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S"),
         res_var = sqrt(soy_yield_res^2)) |>
  feols(soy_var_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_var
etable(soy_var, keep = "rot_crop")

soy_var |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |> 
  separate(term, into = c("temp", "term"), sep = 8) |>
  arrange(estimate) |>
  dwplot(style = "dotwhisker",
         dist_args = list(color = "black", alpha = 0.75),
         vline = geom_vline(
           xintercept = 0, 
           colour = "grey60", 
           linetype = 2)) +
  theme_bw() +
  theme(legend.position = "none") + 
  xlab("Coefficient Estimate") + ylab("") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) ->
  soy_var_plot
soy_var_plot
#ggsave(soy_var_plot, 
#       filename = paste0(fig_dir, "soy_var_plot.png"), 
#       width = 10, height = 7.5)

soy_var |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S", "S-S-S-S-C-S",
                "S-S-C-S-C-S") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  soy_var_plot
soy_var_plot
ggsave(soy_var_plot, 
       filename = paste0(fig_dir, "soy_var_plot.png"), 
       width = 10, height = 7.5)

## Coefficient plots
corn_rot |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |>  
  separate(term, into = c("temp", "term"), sep = 8) |>
  select(c(term, estimate)) |>
  rename(mean = estimate) ->
  corn_m

corn_var |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |>  
  separate(term, into = c("temp", "term"), sep = 8) |>
  select(c(term, estimate)) |>
  rename(sd = estimate) ->
  corn_sd
corn_coeff <- left_join(corn_m, corn_sd, by = "term")

corn_coeff |>
  mutate(pr = ifelse(term == "S-C-S-C-S-C", 1, 0),
         pr = factor(pr)) |>
  ggplot(aes(x = sd, y = mean, label = term)) + 
  geom_jitter(aes(color = pr), size = 2) + 
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") + 
  xlab("Corn standard deviation model coefficients") +
  ylab("Corn mean yield model coefficients") ->
  corn_coeff_plot
corn_coeff_plot
ggsave(corn_coeff_plot, 
       filename = paste0(fig_dir, "corn_coeff_plot.png"), 
       width = 10, height = 10)

soy_rot |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |>  
  separate(term, into = c("temp", "term"), sep = 8) |>
  select(c(term, estimate)) |>
  rename(mean = estimate) ->
  soy_m

soy_var |> 
  tidy() |> 
  filter(grepl("rot_crop", term)) |>  
  separate(term, into = c("temp", "term"), sep = 8) |>
  select(c(term, estimate)) |>
  rename(sd = estimate) ->
  soy_sd
soy_coeff <- left_join(soy_m, soy_sd, by = "term")

soy_coeff |>
  mutate(pr = ifelse(term == "C-S-C-S-C-S", 1, 0),
         pr = factor(pr)) |>
  ggplot(aes(x = sd, y = mean, label = term)) + 
  geom_jitter(aes(color = pr), size = 2) + 
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none")  + 
  xlab("Soybeans standard deviation model coefficients") +
  ylab("Soybeans mean yield model coefficients") ->
  soy_coeff_plot
soy_coeff_plot
ggsave(soy_coeff_plot, 
       filename = paste0(fig_dir, "soy_coeff_plot.png"), 
       width = 10, height = 10)

## Corn with beginning-of-season precipitation
corn_yield_formula_pr <- paste("corn_yield ~ rot_crop+", 
                            paste("pr_", 3:8, collapse = "+", sep = ""), 
                            "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                            "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                            "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                            "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_yield_formula_pr, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot_pr
etable(corn_rot, corn_rot_pr)

### Cover crop dummy
corn_yield_cc <- paste("corn_yield ~ cov_crop+", 
                               paste("pr_", 3:8, collapse = "+", sep = ""), 
                               "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                               "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                               "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                               "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                               "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                               "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

corn_yield_cc_rci <- paste("corn_yield ~ cov_crop*RCI+", 
                       paste("pr_", 3:8, collapse = "+", sep = ""), 
                       "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                       "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                       "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                       "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(cov_crop = ifelse(CC_probability_corn > 58, 1, 0),
         RCI = factor(RCI)) |>
  feols(corn_yield_cc, data = _, cluster = ~COUNTY_FIPS) ->
  corn_cc

crop_df |>
  mutate(cov_crop = ifelse(CC_probability_corn > 58, 1, 0),
         RCI = factor(RCI)) |>
  feols(corn_yield_cc_rci, data = _, cluster = ~COUNTY_FIPS) ->
  corn_cc_rci
etable(corn_cc, corn_cc_rci, keep = c("RCI", "cov_crop"))

corn_cc_rci |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("RCI", prms.estimate_names)) |> 
  mutate(treat = ifelse(prms.x < 40, 0, 1),
         treat = factor(treat, labels = c("No cover crop", "Cover crop")),
         RCI = case_when(
           treat == "No cover crop" ~ gsub("RCI", "", prms.estimate_names),
           treat == "Cover crop" ~ gsub("cov_crop:RCI", "", prms.estimate_names))) |>
  filter(RCI < 5) |>
  ggplot(aes(x = RCI, y = prms.y, groups = treat)) +
  geom_point(size = 2, aes(color = treat)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = treat), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Coefficient Estimate") + ylab("Rotation Complexity Index x Cover Crop") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  corn_cc_rci_plot

soy_yield_cc <- paste("soy_yield ~ cov_crop+", 
                       paste("pr_", 3:8, collapse = "+", sep = ""), 
                       "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                       "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                       "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                       "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

soy_yield_cc_rci <- paste("soy_yield ~ cov_crop*RCI+", 
                      paste("pr_", 3:8, collapse = "+", sep = ""), 
                      "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                      "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                      "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                      "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                      "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  mutate(cov_crop = ifelse(CC_probability_soy > 58, 1, 0)) |>
  feols(soy_yield_cc, data = _, cluster = ~COUNTY_FIPS) ->
  soy_cc
etable(soy_cc)

crop_df |>
  mutate(cov_crop = ifelse(CC_probability_corn > 58, 1, 0),
         RCI = factor(RCI)) |>
  feols(soy_yield_cc_rci, data = _, cluster = ~COUNTY_FIPS) ->
  soy_cc_rci
etable(soy_cc, soy_cc_rci, keep = c("RCI", "cov_crop"))

soy_cc_rci |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("RCI", prms.estimate_names)) |> 
  mutate(treat = ifelse(prms.x < 40, 0, 1),
         treat = factor(treat, labels = c("No cover crop", "Cover crop")),
         RCI = case_when(
           treat == "No cover crop" ~ gsub("RCI", "", prms.estimate_names),
           treat == "Cover crop" ~ gsub("cov_crop:RCI", "", prms.estimate_names))) |>
  filter(RCI < 5) |>
  ggplot(aes(x = RCI, y = prms.y, groups = treat)) +
  geom_point(size = 2, aes(color = treat)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = treat), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  xlab("Coefficient Estimate") + ylab("Rotation Complexity Index x Cover Crop") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom") ->
  soy_cc_rci_plot

## Cover crops interacted with VPD


## Rotation models with scaled outcomes
crop_df |> 
  mutate(pr_6 = scale(pr_6),
         pr_7 = scale(pr_7),
         pr_8 = scale(pr_8),
         cGDD_6m = scale(cGDD_6m),
         cGDD_7m = scale(cGDD_7m),
         cGDD_8m = scale(cGDD_8m),
         tmmx_6 = scale(tmmx_6),
         tmmx_7 = scale(tmmx_7),
         tmmx_8 = scale(tmmx_8),
         tmmn_6 = scale(tmmn_6),
         tmmn_7 = scale(tmmn_7),
         tmmn_8 = scale(tmmn_8),
         soil_6 = scale(soil_6),
         soil_7 = scale(soil_7),
         soil_8 = scale(soil_8),
         vpd_6 = scale(vpd_6),
         vpd_7 = scale(vpd_7),
         vpd_8 = scale(vpd_8),
         nccpi3all_mean = scale(nccpi3all_mean),
         rootznaws_mean = scale(rootznaws_mean),
         soc0_100_mean = scale(soc0_100_mean),
         corn_yield = scale(corn_yield)) |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_rot_scaled
etable(corn_rot_scaled)

corn_rot_scaled |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-C-S-C-S-C" ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C", 
                "S-S-S-S-S-C", "S-S-S-C-S-C") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() + ylim(c(-0.1,0.2)) +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")

crop_df |>
  mutate(pr_6 = scale(pr_6),
         pr_7 = scale(pr_7),
         pr_8 = scale(pr_8),
         cGDD_6m = scale(cGDD_6m),
         cGDD_7m = scale(cGDD_7m),
         cGDD_8m = scale(cGDD_8m),
         tmmx_6 = scale(tmmx_6),
         tmmx_7 = scale(tmmx_7),
         tmmx_8 = scale(tmmx_8),
         tmmn_6 = scale(tmmn_6),
         tmmn_7 = scale(tmmn_7),
         tmmn_8 = scale(tmmn_8),
         soil_6 = scale(soil_6),
         soil_7 = scale(soil_7),
         soil_8 = scale(soil_8),
         vpd_6 = scale(vpd_6),
         vpd_7 = scale(vpd_7),
         vpd_8 = scale(vpd_8),
         nccpi3all_mean = scale(nccpi3all_mean),
         rootznaws_mean = scale(rootznaws_mean),
         soc0_100_mean = scale(soc0_100_mean),
         soy_yield = scale(soy_yield)) |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")) |>
  feols(soy_yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_rot_scaled
etable(soy_rot_scaled)

soy_rot_scaled |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "C-S-C-S-C-S" ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S", 
                "S-S-S-S-C-S", "S-S-C-S-C-S") ~ "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")


crop_df |>
  feols(corn_yield~nccpi3corn_mean|tile_field_ID+year, data = _) ->
  soil_corn

crop_df |>
  feols(soy_yield~nccpi3soy_mean|tile_field_ID+year, data = _) ->
  soil_soy
etable(soil_corn, soil_soy)

library(sf)
library(usmap)

crs <- st_crs("EPSG:4326")
us_map(regions = "counties") |>
  filter(abbr == "IL") |>
  st_as_sf() |>
  st_transform(crs) ->
  map

crop_df |>
  st_as_sf(coords = c("lon", "lat"), crs = crs) ->
  crop_sf

crop_sf |>
  filter(year == 2016 & !is.na(RCI)) |>
  mutate(RCI = factor(RCI)) |>
  ggplot() + 
  geom_sf(data = map) + 
  geom_sf(aes(color = RCI), size = 0.2, fill = "black") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position = "bottom") + 
  ggtitle("RCI values (2016)") ->
  rci_map
rci_map
ggsave(rci_map, 
       filename = paste0(fig_dir, "rci_map.png"), 
       width = 10, height = 10)

# Corn yields
crop_sf |>
  filter(year == 2016 & !is.na(corn_yield)) |>
  ggplot() + 
  geom_sf(data = map) + 
  geom_sf(aes(color = corn_yield), size = 0.2, fill = "black") +
  scale_color_viridis_c()  +
  theme(legend.position = "bottom") + 
  ggtitle("Corn yields (2016)") ->
  corn_yield_map
corn_yield_map
ggsave(corn_yield_map, 
       filename = paste0(fig_dir, "corn_yield_map.png"), 
       width = 10, height = 10)

# Soybean yields
crop_sf |>
  filter(year == 2016 & !is.na(soy_yield)) |>
  ggplot() + 
  geom_sf(data = map) + 
  geom_sf(aes(color = soy_yield), size = 0.2, fill = "black") +
  scale_color_viridis_c()  +
  theme(legend.position = "bottom") + 
  ggtitle("Soybeans yields (2016)") ->
  soy_yield_map
soy_yield_map
ggsave(soy_yield_map, 
       filename = paste0(fig_dir, "soy_yield_map.png"), 
       width = 10, height = 10)

## NCCPI
crop_sf |>
  filter(year == 2016 & !is.na(nccpi3corn_mean)) |>
  ggplot() + 
  geom_sf(data = map) + 
  geom_sf(aes(color = nccpi3corn_mean), size = 0.2, fill = "black") +
  scale_color_viridis_c()  +
  theme(legend.position = "bottom") + 
  ggtitle("NCCPI Corn (2016)") ->
  nccpi_corn_map
nccpi_corn_map
ggsave(nccpi_corn_map, 
       filename = paste0(fig_dir, "nccpi_corn_map.png"), 
       width = 10, height = 10)

crop_sf |>
  filter(year == 2016 & !is.na(nccpi3soy_mean)) |>
  ggplot() + 
  geom_sf(data = map) + 
  geom_sf(aes(color = nccpi3soy_mean), size = 0.2, fill = "black") +
  scale_color_viridis_c()  +
  theme(legend.position = "bottom") + 
  ggtitle("NCCPI Soybeans (2016)") ->
  nccpi_soy_map
nccpi_soy_map
ggsave(nccpi_soy_map, 
       filename = paste0(fig_dir, "nccpi_soy_map.png"), 
       width = 10, height = 10)

### Both crops scaled
crop_df |>
  mutate(corn_yield = scale(corn_yield),
         soy_yield = scale(soy_yield)) |>
  pivot_longer(c(corn_yield, soy_yield), names_to = "crop_name", values_to = "yield",
               names_repair = "minimal") ->
  crop_long

yield_formula <- paste("yield ~ rot_crop+", 
                       paste("pr_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                       "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                       "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                       "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

dict <- c("rot_cropC-C-C-C-S-C" = "C-C-C-C-S-C", "rot_cropC-C-C-S-C-C" = "C-C-C-S-C-C",
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
          "rot_cropS-S-S-S-S-C" = "S-S-S-S-S-C", "rot_cropC-C-C-C-C-S" = "C-C-C-C-C-S", 
          "rot_cropC-C-C-C-S-S" = "C-C-C-C-S-S", "rot_cropS-S-S-S-S-S" = "S-S-S-S-S-S", 
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
          "rot_cropS-S-S-S-C-S" = "S-S-S-S-C-S", 
          "rot_cropCS_perfect_rotation" = "CS_perfect_rotation")

crop_long |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
         rot_crop = ifelse(rot_crop %in% c("C-S-C-S-C-S", "S-C-S-C-S-C"), 
                           "CS_perfect_rotation", rot_crop),
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  crop_rot

etable(crop_rot, 
       dict = dict,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))

# Transition out of corn/soy monoculture
# C-C-C-C-C-C: 0.00
# C-C-C-C-C-S: 1.41
# C-C-C-C-S-C: 1.73
# C-C-C-S-C-S: 2.00
# C-C-S-C-S-C: 2.24
# C-S-C-S-C-S: 2.24

# S-S-S-S-S-S: 0.00
# S-S-S-S-S-C: 1.41
# S-S-S-S-C-S: 1.73
# S-S-S-C-S-C: 2.00
# S-S-C-S-C-S: 2.24
# S-C-S-C-S-C: 2.24

crop_rot |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-S-S-S-S-S" ~ "Soy Monoculture",
    term == "CS_perfect_rotation" ~ "Perfect Rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-C-S-C", "C-C-C-S-C-S", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-S-C-S", "S-S-S-C-S-C", "S-S-C-S-C-S") ~ 
      "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")

crop_rot |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-S-S-S-S-S" ~ "Soy Monoculture",
    term == "CS_perfect_rotation" ~ "Perfect Rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-C-S-C", "C-C-C-S-C-S", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-S-C-S", "S-S-S-C-S-C", "S-S-C-S-C-S") ~ 
      "Transitioning",
    .default = 'Other')) |>
  filter(group != "Other") |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")

crop_rot |>
  augment(newdata = crop_long |> 
            filter(rot_crop %in% corn_soy_patterns$pattern) |>
            mutate(rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
                   rot_crop = ifelse(rot_crop %in% c("C-S-C-S-C-S", "S-C-S-C-S-C"), 
                                     "CS_perfect_rotation", rot_crop),
                   rot_crop = factor(rot_crop),
                   rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C"))) |>
  rename(yield_pred = .fitted) |>
  mutate(yield_res = yield-yield_pred,
         var_yield = yield_res^2) ->
  crop_res

var_formula <- paste("var_yield ~ rot_crop+", 
                     paste("pr_", 6:8, collapse = "+", sep = ""), 
                     "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                     "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                     "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                     "+", paste("soil_", 6:8, collapse = "+", sep = ""), 
                     "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                     "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_res |>
  feols(var_formula, data = _, cluster = ~COUNTY_FIPS) ->
  crop_var

etable(crop_var, 
       dict = dict,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))

crop_var |>
  coefplot() |>
  data.frame() |> 
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group =  case_when(
    term == "S-S-S-S-S-S" ~ "Soy Monoculture",
    term == "CS_perfect_rotation" ~ "Perfect Rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-C-S-C", "C-C-C-S-C-S", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-S-C-S", "S-S-S-C-S-C", "S-S-C-S-C-S") ~ 
      "Transitioning",
    .default = 'Other')) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("Coefficient Estimate") + ylab("Crop sequence") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")


RCI_formula <- paste("yield~RCI+", 
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


crop_long |>
  mutate(RCI = factor(RCI)) |>
  feols(RCI_formula, data = _, cluster = ~COUNTY_FIPS) ->
  rci_cs

etable(rci_cs, 
       dict = dict,
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))


## Add third crop sequences
# Ley
expand.grid(crop_0 = c("1","5","6"), 
            crop_1 = c("1","5","6"), 
            crop_2 = c("1","5","6"), 
            crop_3 = c("1","5","6"), 
            crop_4 = c("1","5","6"),
            crop_5 = c("1","5","6")) |>
  data.frame() |>
  mutate(pattern = paste(crop_5, crop_4, crop_3, crop_2, crop_1, 
                         crop_0, sep = "-")) ->
  cs_ley_patterns

corn_fg_formula <- paste("corn_yield ~ rot_fg+", 
                            paste("pr_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                            "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                            "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                            "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                            "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_fg %in% cs_ley_patterns$pattern) |>
  mutate(rot_fg = gsub("1", "C", gsub("2", "S",gsub("6", "L", rot_fg))),  
         rot_fg = factor(rot_fg),
         rot_fg = relevel(rot_fg, ref = "C-C-C-C-C-C")) |>
  feols(corn_fg_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_ley_rot

etable(corn_ley_rot, 
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))

corn_ley_rot |> 
  tidy() |> 
  filter(grepl("rot_fg", term)) |> 
  separate(term, into = c("temp", "term"), sep = 6) |>
  select(-temp) |>
  arrange(desc(estimate)) |>
  data.frame()

soy_fg_formula <- paste("soy_yield ~ rot_fg+", 
                         paste("pr_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("cGDD_", 6:8, "m", collapse = "+", sep = ""), 
                         "+", paste("tmmx_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("tmmn_", 6:8, collapse = "+", sep = ""), 
                         "+", paste("soil_", 6:8, collapse = "+", sep = ""),
                         "+", paste("vpd_", 6:8, collapse = "+", sep = ""),
                         "+nccpi3all_mean+rootznaws_mean+soc0_100_mean|tile_field_ID+year") |>
  as.formula()

crop_df |>
  filter(rot_fg %in% cs_ley_patterns$pattern) |>
  mutate(rot_fg = gsub("1", "C", gsub("2", "S",gsub("6", "L", rot_fg))),  
         rot_fg = factor(rot_fg),
         rot_fg = relevel(rot_fg, ref = "S-S-S-S-S-S")) |>
  feols(soy_fg_formula, data = _, cluster = ~COUNTY_FIPS) ->
  soy_ley_rot

etable(soy_ley_rot, 
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))

soy_ley_rot |> 
  tidy() |> 
  filter(grepl("rot_fg", term)) |> 
  separate(term, into = c("temp", "term"), sep = 6) |>
  select(-temp) |>
  arrange(desc(estimate)) |>
  data.frame()

# Wheat: repalce all categories with 99
## Durum Wheat: 22 
## Spring Wheat: 23 
## Winter Wheat: 24
## Dbl Crop WinWht/Soybeans: 26

# This specification crashes due to dimensionality!
expand.grid(crop_0 = c("1","5","99"), 
            crop_1 = c("1","5","99"), 
            crop_2 = c("1","5","99"), 
            crop_3 = c("1","5","99"), 
            crop_4 = c("1","5","99"),
            crop_5 = c("1","5","99")) |>
  data.frame() |>
  mutate(pattern = paste(crop_5, crop_4, crop_3, crop_2, crop_1, 
                         crop_0, sep = "-")) ->
  csw_patterns

intersect(csw_patterns$pattern, corn_soy_patterns$pattern)

crop_df |>
  mutate(rot_crop = gsub("22", "99", gsub("23", "99", gsub("24", "99", gsub("26", "99", rot_crop))))) |>
  filter(rot_crop %in% csw_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", gsub("99", "W", rot_crop))),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |> 
  tab(rot_crop) |> 
  arrange(desc(Freq.)) -> 
  rot_crop_wheat

crop_df |>
  mutate(rot_crop = gsub("22", "99", gsub("23", "99", gsub("24", "99", gsub("26", "99", rot_crop))))) |>
  filter(rot_crop %in% csw_patterns$pattern) |>
  mutate(rot_crop = gsub("1", "C", gsub("5", "S", gsub("99", "W", rot_crop))),  
         rot_crop = factor(rot_crop),
         rot_crop = relevel(rot_crop, ref = "C-C-C-C-C-C")) |>
  feols(corn_yield_formula, data = _, cluster = ~COUNTY_FIPS) ->
  corn_w_rot

etable(corn_w_rot, 
       drop = c("pr_", "cGDD_", "tmmx_", "tmmn_", "soil_", "vpd_", "Constant",
                "nccpi3all_mean", "rootznaws_mean", "soc0_100_mean"))


## Add frequency tables or rotations
## Add extra RCI column
## Order coefficient effects of crop sequences from monoculture to perfect rotation