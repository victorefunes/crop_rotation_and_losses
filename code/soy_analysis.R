## ============================================================================
## soy_analysis.R
## Just-Pope production risk analysis — SOYBEANS
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Tables produced:
##   tab:soy_rot        — Rotation patterns and soybean yields
##   tab:soy_rci        — Rotational Complexity and soybean yields
##   tab:soy_rot_vpd    — Rotation x VPD interaction, soy
##   tab:soy_rci_vpd    — RCI x VPD interaction, soy
##   tab:soy_jp_mean    — Just-Pope stage 1: soybean yield mean
##   tab:soy_jp_var     — Just-Pope stage 2: soybean yield variance (OLS)
##   tab:soy_rci_jp     — Just-Pope factor RCI: soybean yield moments
##
## Figures produced:
##   soy_rot_plot       — Response of soybean yields to rotation sequences
##   soy_rci_plot       — Response of soybean yields to RCI
##   soy_var_plot       — Response of std dev of soybean yields to rotation sequences
##   soy_coeff_plot     — Mean vs variance coefficients scatter (soy)
##   soy_jp_plot        — Just-Pope mean-variance decomposition (soy)
##   soy_yield_map      — Spatial map of soybean yields (2016)
##   nccpi_soy_map      — Spatial map of NCCPI soybean (2016)
## ============================================================================

setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

source("rotation_setup.R")

# ── Load soy data ─────────────────────────────────────────────────────────────

cat("Loading soy data...\n")
soy_df <- fread(
  "./Data and Data Descriptions/clean/soy_rci_il_long.csv"
) |>
  (\(df) df[order(df$tile_field_ID, df$year)])() |>
  rci_correction() |>
  add_degree_days()

cat("Soy raw rows:", nrow(soy_df), "\n")

# ── Build analysis sample ─────────────────────────────────────────────────────

soy_jp_data <- soy_df |>
  filter(rot_crop %in% corn_soy_patterns$pattern) |>
  mutate(
    rot_crop = gsub("1", "C", gsub("5", "S", rot_crop)),
    rot_crop = factor(rot_crop),
    rot_crop = relevel(rot_crop, ref = "S-S-S-S-S-S")
  ) |>
  filter(!is.na(soy_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))

cat("Soy analysis sample:", nrow(soy_jp_data), "rows\n")

rm(soy_df); gc()

# ── 1. OLS mean model — soy ───────────────────────────────────────────────────
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
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       title    = "Rotation patterns and soybean yields",
       label    = "tab:soy_rot",
       extralines = list("_Controls" = c("No", "Yes")),
       file     = paste0(tab_dir, "soy_rot.tex"))

# BH correction for soy
pvals_soy <- broom::tidy(soy_rot) |>
  filter(grepl("rot_crop", term)) |>
  arrange(p.value) |>
  mutate(p_adj_bh = p.adjust(p.value, method = "BH"),
         sig_raw  = p.value  < 0.05,
         sig_bh   = p_adj_bh < 0.05)

cat("Soy sequences significant at 5% (unadjusted):", sum(pvals_soy$sig_raw), "\n")
cat("Soy sequences significant at 5% FDR (BH):    ", sum(pvals_soy$sig_bh),  "\n")

# Figure: soy_rot_plot — Response of soybean yields to rotation sequences
soy_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "C-S-C-S-C-S"                         ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S", "C-C-C-S-C-S",
                "S-S-S-S-C-S", "S-S-C-S-C-S")    ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of soybean yields to rotation sequences",
       caption = "Reference: S-S-S-S-S-S. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rot_plot
ggsave(paste0(fig_dir, "soy_rot_plot.png"), soy_rot_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 2. RCI models — soy ───────────────────────────────────────────────────────
# Table: tab:soy_rci | Figure: soy_rci_plot

soy_rci_formula <- make_jp_formula("soy_yield", "RCI", all_controls)

soy_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_formula, data = _, cluster = ~COUNTY_FIPS) -> soy_rci_all

soy_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(soy_rci_formula, data = _, cluster = ~COUNTY_FIPS) -> soy_rci_cs

etable(soy_rci_all, soy_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("All crops", "Corn and soybeans only"),
       keep     = "RCI",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       title    = "Rotation Complexity and soybean yields",
       label    = "tab:soy_rci",
       file     = paste0(tab_dir, "soy_rci.tex"))

# Figure: soy_rci_plot — Changes in RCI and soybean yields
soy_rci_cs |>
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
       title   = "Changes in RCI and soybean yields",
       caption = "Corn-soy fields only. Reference: RCI = 1.41. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_rci_plot
ggsave(paste0(fig_dir, "soy_rci_plot.png"), soy_rci_plot,
       width = 10, height = 7.5, dpi = 300)

rm(soy_rci_all, soy_rci_cs); gc()

# ── 3. VPD interaction models — soy ──────────────────────────────────────────
# Tables: tab:soy_rot_vpd, tab:soy_rci_vpd

soy_vpd_formula     <- make_jp_formula("soy_yield", "rot_crop + vpd_name",
                                        all_controls)
soy_rci_vpd_formula <- make_jp_formula("soy_yield",  "RCI * vpd_name",
                                        all_controls)

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
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Effect of weather and rotation sequences on soybean yields",
       label    = "tab:soy_rot_vpd",
       file     = paste0(tab_dir, "soy_rot_vpd.tex"))

etable(soy_rci_vpd,
       tex      = TRUE,
       dict     = dict_vpd,
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "RCI x drought interaction effects on soybean yields",
       label    = "tab:soy_rci_vpd",
       file     = paste0(tab_dir, "soy_rci_vpd.tex"))

rm(soy_rot_vpd, soy_rci_vpd); gc()

# ── 4. Just-Pope stage 1 — soy ────────────────────────────────────────────────
# Table: tab:soy_jp_mean

fml_soy_mean <- make_jp_formula("soy_yield", "rot_crop", all_controls_fgls)
feols(fml_soy_mean, data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s1

soy_rot_dict <- setNames(
  chartr("15","CS", sub("^rot_crop","", grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE))),
  grep("^rot_crop", names(coef(soy_jp_s1)), value=TRUE)
)

bh_note_soy <- paste0(
  "Benjamini-Hochberg FDR correction across ", nrow(pvals_soy),
  " rotation-sequence coefficients: ",
  sum(pvals_soy$sig_raw), " significant at 5\\% (unadjusted); ",
  sum(pvals_soy$sig_bh),  " significant at 5\\% FDR."
)

etable(soy_jp_s1,
       tex      = TRUE,
       keep     = "rot_crop",
       dict     = soy_rot_dict,
       notes    = bh_note_soy,
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Soybean yield conditional mean",
       label    = "tab:soy_jp_mean",
       file     = paste0(tab_dir, "soy_jp_mean.tex"))

rm(soy_rot_nc, soy_rot); gc()

# ── 5. Just-Pope stage 2 — soy ────────────────────────────────────────────────
# Table: tab:soy_jp_var | Figures: soy_var_plot, soy_coeff_plot, soy_jp_plot

soy_jp_s1 |>
  augment(newdata = soy_jp_data) |>
  mutate(resid_sq = (soy_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  soy_jp_data

fml_soy_var <- make_jp_formula("resid_sq", "rot_crop", all_controls_fgls)
feols(fml_soy_var, data = soy_jp_data, cluster = ~tile_field_ID+year) -> soy_jp_s2

etable(soy_jp_s2,
       tex      = TRUE,
       keep     = "rot_crop",
       dict     = soy_rot_dict,
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2 — Soybean yield conditional variance (OLS)",
       label    = "tab:soy_jp_var",
       file     = paste0(tab_dir, "soy_jp_var.tex"))

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
soy_jp_summary |> filter(dominates) |>
  select(rot_crop, mean_est, mean_se, var_est, var_se) |> print(n = Inf)

# Figure: soy_var_plot — Response of std dev of soybean yields to rotation sequences
soy_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(group = case_when(
    term == "C-S-C-S-C-S"                         ~ "Perfect rotation",
    term %in% c("C-C-C-C-C-S","C-C-C-S-C-S",
                "S-S-S-S-C-S","S-S-C-S-C-S")     ~ "Transitioning",
    .default = "Other")) |>
  ggplot(aes(x = reorder(term, -prms.y), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of standard deviation of soybean yields to rotation sequences",
       caption = "Reference: S-S-S-S-S-S. Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  soy_var_plot
ggsave(paste0(fig_dir, "soy_var_plot.png"), soy_var_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: soy_coeff_plot
soy_jp_summary |>
  mutate(pr = factor(as.integer(rot_crop %in% c("5-1-5-1-5-1","1-5-1-5-1-5")))) |>
  ggplot(aes(x = var_est, y = mean_est, label = rot_crop)) +
  geom_jitter(aes(color = pr), size = 2) +
  geom_text(check_overlap = TRUE, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  labs(x = "Soybeans standard deviation model coefficients",
       y = "Soybeans mean yield model coefficients") ->
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
  labs(x      = "Stage 1: effect on mean soy yield (bu/acre)",
       y      = "Stage 2: effect on yield variance",
       colour = NULL, shape = NULL,
       title  = "Just-Pope decomposition: soy rotation effects on mean and variance",
       caption = "Reference: S-S-S-S-S-S. Two-way clustering: field + year.\nQuadrant IV: higher mean, lower variance — unambiguously better.") +
  theme_bw() + theme(legend.position = "bottom") ->
  soy_jp_plot
ggsave(paste0(fig_dir, "soy_jp_plot.png"), soy_jp_plot,
       width = 9, height = 7, dpi = 300)

rm(soy_jp_s2, soy_s1_coef, soy_s2_coef, soy_jp_summary,
   soy_jp_plot, soy_coeff_plot, soy_var_plot); gc()

# ── 6. Just-Pope factor RCI — soy ─────────────────────────────────────────────
# Table: tab:soy_rci_jp

soy_rci_jp_data <- soy_jp_data |>
  select(-any_of(c("resid_sq","h_hat"))) |>
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

etable(soy_rci_jp_s1, soy_rci_jp_s2,
       tex      = TRUE,
       keep     = "^RCI",
       dict     = dict_rci,
       headers  = c("Mean", "Variance"),
       se.below = FALSE,
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Just-Pope: factor RCI and soybean yield moments",
       label    = "tab:soy_rci_jp",
       file     = paste0(tab_dir, "soy_rci_jp.tex"))

rm(soy_rci_jp_data, soy_rci_jp_s1, soy_rci_jp_s2); gc()

# ── 7. FGLS + bootstrap — soy ─────────────────────────────────────────────────

fml_mean  <- make_jp_formula("soy_yield",  "rot_crop", all_controls_fgls)
fml_var   <- make_jp_formula("resid_sq",   "rot_crop", all_controls_fgls)
fml_var_b <- make_jp_formula("resid_sq_b", "rot_crop", all_controls_fgls)

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

cat("Soy obs hitting h_hat floor:", sum(soy_jp_data$h_hat == 1e-6), "\n")

feols(fml_var, data = soy_jp_data,
      weights = ~I(1/h_hat),
      cluster = ~tile_field_ID+year) -> soy_jp_s2b

etable(soy_jp_s1, soy_jp_s2a, soy_jp_s2b,
       keep  = "rot_crop",
       title = "Just-Pope FGLS: soy rotation effects on mean and variance")

rm(soy_jp_s2a); gc()

gc()
boot_soy <- boot_jp_fgls(soy_jp_data, fml_mean, fml_var, fml_var_b,
                          B = 499, seed = 42)
saveRDS(boot_soy,
        "C:/Users/vf006/Documents/boot_soy.rds",
        compress = "xz")

rm(soy_jp_data, soy_jp_s1, soy_jp_s2b, boot_soy); gc()

# ── 8. Spatial maps — soy ────────────────────────────────────────────────────
# Figures: soy_yield_map, nccpi_soy_map

cat("Reloading soy data for spatial maps...\n")
soy_sf <- fread(
  "./Data and Data Descriptions/clean/soy_rci_il_long.csv",
  select = c("tile_field_ID", "year", "lon", "lat",
             "soy_yield", "nccpi3soy_mean")
) |>
  filter(year == 2016) |>
  st_as_sf(coords = c("lon", "lat"), crs = st_crs("EPSG:4326"))

il_map <- us_map(regions = "counties") |>
  filter(abbr == "IL") |>
  st_as_sf() |>
  st_transform(st_crs("EPSG:4326"))

# Figure: soy_yield_map — Soybean yields across Illinois (2016)
soy_sf |>
  filter(!is.na(soy_yield)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = soy_yield), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "Soybean yields — QDANN (2016)",
       caption = "Source: Ma et al. (2024)") ->
  soy_yield_map
ggsave(paste0(fig_dir, "soy_yield_map.png"), soy_yield_map,
       width = 10, height = 10, dpi = 300)

# Figure: nccpi_soy_map — NCCPI soybean productivity index (2016)
soy_sf |>
  filter(!is.na(nccpi3soy_mean)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = nccpi3soy_mean), size = 0.2) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "NCCPI — Soybeans (2016)",
       caption = "Source: gSSURGO") ->
  nccpi_soy_map
ggsave(paste0(fig_dir, "nccpi_soy_map.png"), nccpi_soy_map,
       width = 10, height = 10, dpi = 300)

rm(soy_sf, il_map, soy_yield_map, nccpi_soy_map); gc()

cat("\n=== Soy analysis complete ===\n")
cat("All tables saved to:", tab_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
