## ============================================================================
## corn_rci_analysis.R
## Just-Pope production risk analysis — CORN, Rotational Complexity Index (RCI)
## models.
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Requires: corn_data_prep.R (sourced below) — provides corn_jp_data,
## vpd_controls, load_corn_sf_2016(), and everything from rotation_setup_wa.R.
##
## Tables produced:
##   tab:corn_rci     — Rotational Complexity and corn yields
##   tab:corn_rci_vpd  — RCI x VPD interaction, corn
##   tab:corn_rci_jp   — Just-Pope factor RCI: corn yield moments
##
## Figures produced:
##   corn_rci_plot   — Response of corn yields to RCI (linear-in-levels model)
##   rci_plot        — Nonlinear effect of RCI on corn yield mean and variance
##                      (NOTE: both figures are saved to the same file,
##                      "corn_rci_plot.png" — rci_plot overwrites corn_rci_plot,
##                      carried over unchanged from corn_analysis_full.R)
##   rci_map         — Spatial map of RCI values (2016)
##   nccpi_corn_map  — Spatial map of NCCPI corn (2016)
## ============================================================================

source("corn_data_prep.R")

# ── 1. RCI models — corn ──────────────────────────────────────────────────────
# Table: tab:corn_rci | Figure: corn_rci_plot

corn_rci_nc <- make_jp_formula("corn_yield", "RCI", NULL)
corn_rci_cs <- make_jp_formula("corn_yield", "RCI", all_controls)

# Remove infrequent RCI levels (fewer than 100 observations)
rci_keep <- corn_jp_data |>
  count(RCI) |>
  filter(n >= 100) |>
  pull(RCI)

corn_jp_data |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_nc, data = _, cluster = ~COUNTY_FIPS+year) -> corn_rci_nc

corn_jp_data |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_cs, data = _, cluster = ~COUNTY_FIPS+year) -> corn_rci_cs

etable(corn_rci_nc, corn_rci_cs,
       tex      = TRUE,
       dict     = dict_rci,
       headers  = c("No controls", "Weather and soil controls"),
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
           term == 2.24 ~ "Perfect rotation",
           .default     = "Other")) |>
  filter(term < 5.2) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = "Rotational Complexity Index (RCI)", y = "Coefficient Estimate",
       title   = "Changes in RCI and corn yields",
       caption = "Corn, soy, and wheat fields. Reference: RCI = 0. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rci_plot
ggsave(paste0(fig_dir, "corn_rci_plot.png"), corn_rci_plot,
       width = 10, height = 7.5, dpi = 300)

rm(corn_rci_nc, corn_rci_cs); gc()

# ── 2. VPD interaction — RCI, corn ────────────────────────────────────────────
# Table: tab:corn_rci_vpd

corn_rci_vpd_formula <- make_jp_formula("corn_yield", "RCI * vpd_name",
                                         vpd_controls)

corn_jp_data |>
  mutate(RCI = factor(RCI)) |>
  feols(corn_rci_vpd_formula, data = _, cluster = ~COUNTY_FIPS+year) -> corn_rci_vpd

etable(corn_rci_vpd,
       tex      = TRUE,
       dict     = c(dict_rci, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       fontsize = "scriptsize",
       title    = "RCI x drought interaction effects on corn yields",
       label    = "tab:corn_rci_vpd",
       file     = paste0(tab_dir, "corn_rci_vpd.tex"))

gc()

# ── 3. Just-Pope factor RCI — corn ───────────────────────────────────────────
# Table: tab:corn_rci_jp | Figure: rci_plot

corn_rci_jp_data <- corn_jp_data |>
  select(-any_of(c("resid_sq", "h_hat"))) |>
  filter(RCI %in% rci_keep) |>
  mutate(RCI = factor(RCI))

fml_corn_rci_mean <- make_jp_formula("corn_yield", "RCI", all_controls_fgls)
feols(fml_corn_rci_mean, data = corn_rci_jp_data,
      cluster = ~COUNTY_FIPS+year) -> corn_rci_jp_s1

corn_rci_jp_s1 |>
  augment(newdata = corn_rci_jp_data) |>
  mutate(resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_rci_jp_data

fml_corn_rci_var <- make_jp_formula("resid_sq", "RCI", all_controls_fgls)
feols(fml_corn_rci_var, data = corn_rci_jp_data,
      cluster = ~COUNTY_FIPS+year) -> corn_rci_jp_s2

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
       y      = "Coefficient relative to RCI = 0",
       title  = "Nonlinear effect of RCI on corn yield mean and variance",
       caption = "Corn, soy, and wheat fields. Reference: RCI = 0. Two-way clustering.") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) ->
  rci_plot
ggsave(paste0(fig_dir, "corn_rci_plot.png"), rci_plot,
       width = 9, height = 7, dpi = 300)

rm(corn_rci_jp_data, corn_rci_jp_s1, corn_rci_jp_s2, rci_plot_df, rci_plot); gc()

# ── 4. Spatial maps — RCI, corn ───────────────────────────────────────────────
# Figures: rci_map, nccpi_corn_map
# See load_corn_sf_2016() in corn_data_prep.R for caveats about this reload
# path (rci_correction()/annual_RCI dependency, pre-existing).

corn_spatial <- load_corn_sf_2016()
corn_sf <- corn_spatial$corn_sf
il_map  <- corn_spatial$il_map

# Figure: rci_map — RCI values across Illinois (2016)
corn_sf |>
  filter(!is.na(RCI)) |>
  ggplot() +
  geom_sf(data = il_map) +
  geom_sf(aes(color = RCI), size = 0.2, alpha = 0.3) +
  scale_color_viridis_c() +
  theme(legend.position = "bottom") +
  labs(title   = "Rotational Complexity Index (2016)",
       caption = "Source: CDL / Socolar et al. (2021)") ->
  rci_map
ggsave(paste0(fig_dir, "rci_map.png"), rci_map,
       width = 10, height = 12.5, dpi = 600)

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
       width = 10, height = 12.5, dpi = 600)

rm(corn_spatial, corn_sf, il_map, rci_map, nccpi_corn_map); gc()

cat("\n=== Corn RCI analysis complete ===\n")
cat("All tables saved to:", tab_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
