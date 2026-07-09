## ============================================================================
## corn_analysis.R
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

setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

source("rotation_setup.R")

# ── Load corn data ────────────────────────────────────────────────────────────
# Load, correct, and add degree-day variables in one pipeline.
# corn_df is kept as the raw data source; corn_jp_data is the analysis sample.

cat("Loading corn data...\n")
corn_df <- fread(
  "./Data and Data Descriptions/clean/corn_rci_il_long.csv"
) |>
  (\(df) df[order(df$tile_field_ID, df$year)])() |>
  rci_correction() |>
  add_degree_days()

cat("Corn raw rows:", nrow(corn_df), "\n")

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.

corn_jp_data <- corn_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  left_join(
    rot_features |> select(pattern, n_soy, consec_soy, gap_soy,
                           min_gap_soy, rot_index),
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
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))

cat("Corn analysis sample:", nrow(corn_jp_data), "rows\n")

# Free raw data — no longer needed
rm(corn_df); gc()

# ── Summary statistics table ──────────────────────────────────────────────────
# Table: tab:summary — means and SDs by rotation type
# Produced here because corn_jp_data contains all needed variables.
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
# Table: tab:corn_rot | Figure: corn_rot_plot

corn_yield_formula <- make_jp_formula("corn_yield", "rot_crop", all_controls)

feols(corn_yield ~ rot_crop | tile_field_ID + year,
      data = corn_jp_data, cluster = ~COUNTY_FIPS) -> corn_rot_nc

feols(corn_yield_formula,
      data = corn_jp_data, cluster = ~COUNTY_FIPS) -> corn_rot

etable(corn_rot_nc, corn_rot,
       tex      = TRUE,
       dict     = dict_corn,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       title    = "Rotation patterns and corn yields",
       label    = "tab:corn_rot",
       extralines = list("_Controls" = c("No", "Yes")),
       file     = paste0(tab_dir, "corn_rot.tex"))

# BH multiple comparisons correction
pvals_corn <- broom::tidy(corn_rot) |>
  filter(grepl("rot_crop", term)) |>
  arrange(p.value) |>
  mutate(p_adj_bh = p.adjust(p.value, method = "BH"),
         sig_raw  = p.value  < 0.05,
         sig_bh   = p_adj_bh < 0.05)

cat("Corn sequences significant at 5% (unadjusted):", sum(pvals_corn$sig_raw), "\n")
cat("Corn sequences significant at 5% FDR (BH):    ", sum(pvals_corn$sig_bh),  "\n")

# Figure: corn_rot_plot — Response of corn yields to rotation sequences
corn_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "S-C-S-C-S-C"                         ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C", "C-C-S-C-S-C",
                "S-S-S-S-S-C", "S-S-S-C-S-C")    ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of corn yields to rotation sequences",
       caption = "Reference: C-C-C-C-C-C. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rot_plot
ggsave(paste0(fig_dir, "corn_rot_plot.png"), corn_rot_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 2. RCI models — corn ──────────────────────────────────────────────────────
# Table: tab:corn_rci | Figure: corn_rci_plot

corn_rci_nc <- make_jp_formula("corn_yield", "RCI", NULL)
corn_rci_cs <- make_jp_formula("corn_yield", "RCI", all_controls)

corn_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_nc, data = _, cluster = ~COUNTY_FIPS) -> corn_rci_all

corn_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_cs, data = _, cluster = ~COUNTY_FIPS) -> corn_rci_cs

etable(corn_rci_all, corn_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("All crops", "Corn and soybeans only"),
       keep     = "RCI",
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       title    = "Rotation Complexity and corn yields",
       label    = "tab:corn_rci",
       file     = paste0(tab_dir, "corn_rci.tex"))

# Figure: corn_rci_plot — Changes in RCI and corn yields
corn_rci_cs |>
  coefplot() |>
  data.frame() |>
  filter(grepl("RCI", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 3) |>
  select(-temp) |>
  mutate(term  = as.numeric(term),
         group = case_when(
           term == 2.24              ~ "Perfect rotation",
           term %in% c(1.41,1.73,2) ~ "Transitioning",
           .default = "Other")) |>
  filter(term < 5.2) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = "Rotational Complexity Index (RCI)", y = "Coefficient Estimate",
       title   = "Changes in RCI and corn yields",
       caption = "Corn-soy fields only. Reference: RCI = 1.41. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rci_plot
ggsave(paste0(fig_dir, "corn_rci_plot.png"), corn_rci_plot,
       width = 10, height = 7.5, dpi = 300)

rm(corn_rci_all, corn_rci_cs); gc()

# ── 3. VPD interaction models — corn ─────────────────────────────────────────
# Tables: tab:corn_rot_vpd, tab:corn_rci_vpd

corn_vpd_formula     <- make_jp_formula("corn_yield", "rot_crop + vpd_name",
                                         all_controls)
corn_rci_vpd_formula <- make_jp_formula("corn_yield", "RCI * vpd_name",
                                         all_controls)

feols(corn_yield ~ rot_crop + vpd_name | tile_field_ID + year,
      data = corn_jp_data, cluster = ~COUNTY_FIPS) -> corn_rot_vpd_nc                                         

corn_jp_data |>
  feols(corn_vpd_formula, data = _, cluster = ~COUNTY_FIPS) -> corn_rot_vpd

corn_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS) -> corn_rci_vpd

etable(corn_rot_vpd,
       tex      = TRUE,
       dict     = c(dict_corn, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Effect of weather and rotation sequences on corn yields",
       label    = "tab:corn_rot_vpd",
       file     = paste0(tab_dir, "corn_rot_vpd.tex"))

etable(corn_rci_vpd,
       tex      = TRUE,
       dict     = dict_vpd,
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "RCI x drought interaction effects on corn yields",
       label    = "tab:corn_rci_vpd",
       file     = paste0(tab_dir, "corn_rci_vpd.tex"))     

#rm(corn_rot_vpd, corn_rci_vpd); 
gc()

# ── 4. Just-Pope stage 1 — corn ───────────────────────────────────────────────
# Table: tab:corn_jp_mean | Figures: corn_rot_plot (mean), corn_jp_plot

# Lag comparison models (confirmatory Z-vector pre-registration)
fml_corn_lag1      <- make_jp_formula("corn_yield", "soy_lag1", all_controls_fgls)
fml_corn_lag1_lag2 <- make_jp_formula("corn_yield", "soy_lag1 + soy_lag2",
                                       all_controls_fgls)
fml_corn_index     <- make_jp_formula("corn_yield", "rot_index", all_controls_fgls)
fml_corn_mean      <- make_jp_formula("corn_yield", "rot_crop", all_controls_fgls)

feols(fml_corn_lag1,      data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s1_lag1
feols(fml_corn_lag1_lag2, data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s1_lag2
feols(fml_corn_index,     data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s1_idx
feols(fml_corn_mean,      data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s1

bh_note_corn <- paste0(
  "Benjamini-Hochberg FDR correction across ", nrow(pvals_corn),
  " rotation-sequence coefficients: ",
  sum(pvals_corn$sig_raw), " significant at 5\\% (unadjusted); ",
  sum(pvals_corn$sig_bh),  " significant at 5\\% FDR."
)

rot_dict <- setNames(
  chartr("15", "CS", sub("^rot_crop", "", grep("^rot_crop", names(coef(corn_jp_s1)), value=TRUE))),
  grep("^rot_crop", names(coef(corn_jp_s1)), value=TRUE)
)

# Table: tab:corn_jp_mean
etable(corn_jp_s1, corn_jp_s1_lag1, corn_jp_s1_lag2, corn_jp_s1_idx,
       tex      = TRUE,
       keep     = c("^[CS]-", "soy_lag", "rot_index"),
       dict     = rot_dict,
       notes    = bh_note_corn,
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Corn yield: full sequences vs lag summary variables",
       label    = "tab:corn_jp_mean",
       file     = paste0(tab_dir, "corn_jp_mean.tex"))

rm(corn_jp_s1_lag1, corn_jp_s1_lag2, corn_jp_s1_idx, corn_rot_nc, corn_rot)
gc()

# ── Z-vector model (confirmatory spec) ───────────────────────────────────────
# Table: tab:zvector — Effect of rotation patterns on yields
# This is the paper's main result table (Table 1 in the draft PDF).
# Four structural features: late_soy, soy_gap, soy_cons, nsoy
# Run on both corn and soy data, both mean and variance stages.
# NOTE: requires late_soy, soy_gap, soy_cons, nsoy to be in corn_jp_data.
# These must be constructed before this chunk runs.

# Construct Z-vector variables if not already present
# late_soy: negative integer = how many periods ago was the last soy harvest
corn_jp_data <- corn_jp_data |>
  mutate(
    # Parse rotation sequence to find last soy year
    seq_vec   = strsplit(as.character(rot_crop), "-"),
    late_soy  = sapply(seq_vec, function(v) {
      soy_pos <- which(rev(v) == "S")   # positions from most recent (1=t-1)
      if (length(soy_pos) == 0) 0L else -min(soy_pos)
    }),
    soy_cons  = sapply(seq_vec, function(v) {
      runs <- rle(v)
      as.integer(any(runs$lengths[runs$values == "S"] >= 2))
    }),
    soy_gap   = sapply(seq_vec, function(v) {
      pos <- which(v == "S")
      if (length(pos) < 2) 0L else min(diff(pos))
    }),
    nsoy      = sapply(seq_vec, function(v) sum(v == "S"))
  ) |>
  select(-seq_vec)

# Z-vector stage 1 — corn mean
fml_z_corn_mean <- make_jp_formula("corn_yield",
                                    "late_soy + soy_gap + soy_cons + nsoy",
                                    all_controls_fgls)

feols(fml_z_corn_mean, data = corn_jp_data,
      cluster = ~tile_field_ID + year) -> corn_z_s1

# Z-vector stage 2 — corn variance
corn_jp_data <- corn_z_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq_z = (corn_yield - .fitted)^2) |>
  select(-starts_with("."))

fml_z_corn_var <- make_jp_formula("resid_sq_z",
                                   "late_soy + soy_gap + soy_cons + nsoy",
                                   all_controls_fgls)

feols(fml_z_corn_var, data = corn_jp_data,
      cluster = ~tile_field_ID + year) -> corn_z_s2

# ── Save Z-vector models for tables_combined.R ────────────────────────────────
#saveRDS(
#  list(z_s1        = corn_z_s1,
#       z_s2        = corn_z_s2,
#       rot_vpd_nc  = corn_rot_vpd_nc,
#       rot_vpd     = corn_rot_vpd),
#  file     = "C:/Users/vf006/Documents/corn_z_models.rds",
#  compress = "bzip2"
#)  

# ── 5. Just-Pope stage 2 — corn ───────────────────────────────────────────────
# Table: tab:corn_jp_var | Figures: corn_var_plot, corn_coeff_plot, corn_jp_plot

corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_jp_data

fml_corn_var <- make_jp_formula("resid_sq", "rot_crop", all_controls_fgls)
feols(fml_corn_var, data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s2

# Table: tab:corn_jp_var
etable(corn_jp_s2,
       tex      = TRUE,
       keep     = "^[CS]-",
       dict     = rot_dict,
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2 — Corn yield conditional variance (OLS)",
       label    = "tab:corn_jp_var",
       file     = paste0(tab_dir, "corn_jp_var.tex"))

# ── Figure: score_yield — Response of yields to rotation score values ─────────
# corn_jp_s1 / corn_jp_s2 contain the sequence-level coefficients.
# corn_z_s1 supplies the Z-vector coefficients used to compute the score.

# Step 1: compute score for each unique sequence from Z-vector coefficients
z_coefs <- coef(corn_z_s1)

score_df <- corn_jp_data |>
  group_by(rot_crop) |>
  summarise(
    late_soy = mean(late_soy, na.rm = TRUE),
    soy_gap  = mean(soy_gap,  na.rm = TRUE),
    soy_cons = mean(soy_cons, na.rm = TRUE),
    nsoy     = mean(nsoy,     na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    score = z_coefs["late_soy"] * late_soy +
            z_coefs["soy_gap"]  * soy_gap  +
            z_coefs["soy_cons"] * soy_cons +
            z_coefs["nsoy"]     * nsoy
  )

# Step 2: extract sequence-level coefficients from the FULL sequence models
corn_s1_coef <- broom::tidy(corn_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    mean_est = estimate,
    mean_se  = std.error
  )

corn_s2_coef <- broom::tidy(corn_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(
    rot_crop = gsub("rot_crop", "", term),
    var_est  = estimate,
    var_se   = std.error
  )

# Step 3: join score onto sequence-level coefficients
score_plot_df <- score_df |>
  left_join(corn_s1_coef, by = "rot_crop") |>
  left_join(corn_s2_coef, by = "rot_crop") |>
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
      "Reference: corn monoculture (score = 0). Two-way clustering: field + year.\n",
      "Score = Z-vector coefficients applied to sequence-level feature averages.\n",
      "Left: stage-1 mean coefficients. Right: stage-2 variance coefficients."
    )
  )

ggsave(paste0(fig_dir, "score_yield.png"), score_yield,
       width = 12, height = 8, dpi = 300)
cat("score_yield figure saved.\n")


# Summary data frame
corn_s1_coef <- broom::tidy(corn_jp_s1) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term),
            mean_est = estimate, mean_se = std.error, mean_p = p.value)

corn_s2_coef <- broom::tidy(corn_jp_s2) |>
  filter(grepl("rot_crop", term)) |>
  transmute(rot_crop = gsub("rot_crop","",term),
            var_est  = estimate, var_se  = std.error, var_p  = p.value)

corn_jp_summary <- corn_s1_coef |>
  left_join(corn_s2_coef, by = "rot_crop") |>
  mutate(mean_sig  = mean_p < 0.05,
         var_sig   = var_p  < 0.05,
         dominates = mean_est > 0 & var_est < 0,
         tradeoff  = mean_est > 0 & var_est > 0) |>
  arrange(desc(dominates), desc(mean_est))

cat("Corn sequences dominating monoculture:\n")
corn_jp_summary |> filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |> print(n = Inf)

# Figure: corn_var_plot — Response of std dev of corn yields to rotation sequences
corn_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "S-C-S-C-S-C"                      ~ "Perfect rotation",
    term %in% c("C-C-C-C-S-C","C-C-S-C-S-C",
                "S-S-S-S-S-C","S-S-S-C-S-C")  ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of standard deviation of corn yields to rotation sequences",
       caption = "Reference: C-C-C-C-C-C. Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_var_plot
ggsave(paste0(fig_dir, "corn_var_plot.png"), corn_var_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: corn_coeff_plot — Mean vs variance coefficients
corn_jp_summary |>
  mutate(pr = factor(as.integer(rot_crop %in% c("5-1-5-1-5-1","1-5-1-5-1-5")))) |>
  ggplot(aes(x = var_est, y = mean_est, label = rot_crop)) +
  geom_jitter(aes(color = pr), size = 2) +
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  labs(x = "Corn standard deviation model coefficients",
       y = "Corn mean yield model coefficients") ->
  corn_coeff_plot
ggsave(paste0(fig_dir, "corn_coeff_plot.png"), corn_coeff_plot,
       width = 10, height = 10, dpi = 300)

# Figure: corn_jp_plot — Just-Pope mean-variance decomposition
corn_jp_summary |>
  mutate(
    type = case_when(
      dominates           ~ "Dominates monoculture",
      tradeoff & mean_sig ~ "Mean-variance trade-off",
      mean_est < 0        ~ "Worse than monoculture",
      TRUE                ~ "No significant difference"),
    perfect = rot_crop %in% c("5-1-5-1-5-1","1-5-1-5-1-5")
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
  labs(x      = "Stage 1: effect on mean corn yield (bu/acre)",
       y      = "Stage 2: effect on yield variance",
       colour = NULL, shape = NULL,
       title  = "Just-Pope decomposition: corn rotation effects on mean and variance",
       caption = "Reference: C-C-C-C-C-C. Two-way clustering: field + year.\nQuadrant IV: higher mean, lower variance — unambiguously better.") +
  theme_bw() + theme(legend.position = "bottom") ->
  corn_jp_plot
ggsave(paste0(fig_dir, "corn_jp_plot.png"), corn_jp_plot,
       width = 9, height = 7, dpi = 300)

rm(corn_jp_s2, corn_s1_coef, corn_s2_coef, corn_jp_summary, corn_jp_plot,
   corn_coeff_plot, corn_var_plot); gc()

# ── 6. Just-Pope factor RCI — corn ───────────────────────────────────────────
# Table: tab:corn_rci_jp | Figure: rci_plot

corn_rci_jp_data <- corn_jp_data |>
  select(-any_of(c("resid_sq", "h_hat"))) |>
  mutate(RCI = factor(RCI))

fml_corn_rci_mean <- make_jp_formula("corn_yield", "RCI", all_controls_fgls)
feols(fml_corn_rci_mean, data = corn_rci_jp_data,
      cluster = ~tile_field_ID+year) -> corn_rci_jp_s1

corn_rci_jp_s1 |>
  augment(newdata = corn_rci_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_rci_jp_data

fml_corn_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls_fgls)
feols(fml_corn_rci_var, data = corn_rci_jp_data,
      cluster = ~tile_field_ID+year) -> corn_rci_jp_s2

# Table: tab:corn_rci_jp
etable(corn_rci_jp_s1, corn_rci_jp_s2,
       tex      = TRUE,
       keep     = "^RCI",
       dict     = dict_rci,
       headers  = c("Mean", "Variance"),
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Just-Pope: factor RCI and corn yield moments",
       label    = "tab:corn_rci_jp",
       file     = paste0(tab_dir, "corn_rci_jp.tex"))

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
  geom_ribbon(aes(ymin = est - 1.96*se, ymax = est + 1.96*se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~moment, scales = "free_y") +
  scale_colour_manual(values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_fill_manual(  values = c("Mean" = "#1f77b4", "Variance" = "#d62728")) +
  scale_x_continuous(breaks = sort(unique(rci_plot_df$rci))) +
  labs(x      = "Rotational Complexity Index (RCI)",
       y      = "Coefficient relative to RCI = 1.41",
       title  = "Nonlinear effect of RCI on corn yield mean and variance",
       caption = "Corn-soy fields only. Reference: RCI = 1.41. Two-way clustering.") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) ->
  rci_plot
ggsave(paste0(fig_dir, "rci_plot.png"), rci_plot,
       width = 9, height = 7, dpi = 300)

rm(corn_rci_jp_data, corn_rci_jp_s1, corn_rci_jp_s2, rci_plot_df, rci_plot); gc()

# ── 7. FGLS + bootstrap — corn ────────────────────────────────────────────────

fml_mean  <- make_jp_formula("corn_yield", "rot_crop", all_controls_fgls)
fml_var   <- make_jp_formula("resid_sq",   "rot_crop", all_controls_fgls)
fml_var_b <- make_jp_formula("resid_sq_b", "rot_crop", all_controls_fgls)

# Stage 2a to get h_hat for FGLS weights
corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_jp_data

feols(fml_var, data = corn_jp_data, cluster = ~tile_field_ID+year) -> corn_jp_s2a

corn_jp_s2a |>
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |>
  select(-starts_with(".")) ->
  corn_jp_data

cat("Obs hitting h_hat floor:", sum(corn_jp_data$h_hat == 1e-6), "\n")

feols(fml_var, data = corn_jp_data,
      weights = ~I(1/h_hat),
      cluster = ~tile_field_ID+year) -> corn_jp_s2b

etable(corn_jp_s1, corn_jp_s2a, corn_jp_s2b,
       keep  = "rot_crop",
       title = "Just-Pope FGLS: corn rotation effects on mean and variance")

rm(corn_jp_s2a); gc()

source("save_models_lean.R")
boot_corn <- boot_jp_fgls(corn_jp_data, fml_mean, fml_var, fml_var_b,
                           B = 499, seed = 42, n_workers = 4)

saveRDS(
  list(
    # Z-vector models — lean extracts only
    z_s1         = lean_feols(corn_z_s1),
    z_s2         = lean_feols(corn_z_s2),
    # Full sequence models needed for score_yield figure
    jp_s1        = lean_feols(corn_jp_s1),
    jp_s2        = lean_feols(corn_jp_s2),
    # VPD models
    rot_vpd_nc   = lean_feols(corn_rot_vpd_nc),
    rot_vpd      = lean_feols(corn_rot_vpd),
    # Bootstrap — lean strips embedded data from feols objects inside
    boot         = lean_boot(boot_corn),
    # Score construction inputs — just the coefficient vector
    z_coefs      = coef(corn_z_s1),
    # Score data frame — sequence-level feature averages (small)
    score_df     = corn_jp_data |>
                     group_by(rot_crop) |>
                     summarise(late_soy = mean(late_soy, na.rm = TRUE),
                               soy_gap  = mean(soy_gap,  na.rm = TRUE),
                               soy_cons = mean(soy_cons, na.rm = TRUE),
                               nsoy     = mean(nsoy,     na.rm = TRUE),
                               .groups  = "drop")
  ),
  file     = "C:/Users/vf006/Documents/corn_z_models.rds",
  compress = "bzip2"
)

cat("Corn lean models saved. File size:",
    round(file.size("C:/Users/vf006/Documents/corn_z_models.rds") / 1e6, 1),
    "MB\n")

rm(corn_jp_data, corn_jp_s1, corn_jp_s2b, boot_corn); gc()

# ── 8. Spatial maps — corn ────────────────────────────────────────────────────
# Figures: corn_yield_map, rci_map, nccpi_corn_map
# Reload a minimal version of corn_df for spatial plotting only

cat("Reloading corn data for spatial maps...\n")
corn_sf <- fread(
  "./Data and Data Descriptions/clean/corn_rci_il_long.csv",
  select = c("tile_field_ID", "year", "lon", "lat",
             "corn_yield", "RCI", "nccpi3corn_mean")
) |>
  filter(year == 2016) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs("EPSG:4326"))

il_map <- us_map(regions = "counties") |>
  filter(abbr == "IL") |>
  st_as_sf() |>
  st_transform(st_crs("EPSG:4326"))

# Figure: corn_yield_map — Corn yields across Illinois (2016)
corn_sf |>
  filter(!is.na(corn_yield)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = corn_yield), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "Corn yields — QDANN (2016)",
       caption = "Source: Ma et al. (2024)") ->
  corn_yield_map
ggsave(paste0(fig_dir, "corn_yield_map.png"), corn_yield_map,
       width = 10, height = 10, dpi = 300)

# Figure: rci_map — RCI values across Illinois (2016)
corn_sf |>
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

# Figure: nccpi_corn_map — NCCPI corn productivity index (2016)
corn_sf |>
  filter(!is.na(nccpi3corn_mean)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = nccpi3corn_mean), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "NCCPI — Corn (2016)",
       caption = "Source: gSSURGO") ->
  nccpi_corn_map
ggsave(paste0(fig_dir, "nccpi_corn_map.png"), nccpi_corn_map,
       width = 10, height = 10, dpi = 300)

rm(corn_sf, il_map, corn_yield_map, rci_map, nccpi_corn_map); gc()

# ── 9. County-level yields and maximum entropy distribution ───────────────────

cat("Reloading corn data for county-level analysis...\n")
corn_county <- fread(
  "./Data and Data Descriptions/clean/corn_rci_il_long.csv",
  select = c("tile_field_ID", "COUNTY_FIPS", "STATE_FIPS", "year", "corn_yield")
) |>
  group_by(COUNTY_FIPS, STATE_FIPS, year) |>
  summarise(yield = mean(corn_yield, na.rm = TRUE), .groups = "drop") |>
  rename(state = STATE_FIPS, county = COUNTY_FIPS)

source("C:/Users/vf006/Box/crop_rotations_and_losses/code/maxent_tack2013.r")
corn_det <- detrend_yields_proportional(corn_county)
dist_m   <- tack_table2(corn_det$yield_norm)

rm(corn_county, corn_det); gc()

cat("\n=== Corn analysis complete ===\n")
cat("All tables saved to:", tab_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
