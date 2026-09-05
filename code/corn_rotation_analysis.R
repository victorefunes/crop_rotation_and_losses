## ============================================================================
## corn_rotation_analysis.R
## Just-Pope production risk analysis — CORN, rotation-sequence models.
## Authors: Lawson Connor, Victor Funes-Leal, Eunchun Park
## University of Arkansas
## ----------------------------------------------------------------------------
## Requires: corn_data_prep.R (sourced below) — provides corn_jp_data,
## vpd_controls, load_corn_sf_2016(), and everything from rotation_setup_wa.R.
##
## Tables produced:
##   tab:corn_rot        — Rotation patterns and corn yields
##   tab:corn_rot_vpd     — Rotation x VPD interaction, corn
##   tab:corn_jp_mean     — Just-Pope stage 1: corn yield mean
##   tab:corn_jp_moments  — Just-Pope stage 2/3: corn yield variance (FGLS) and
##                          standardized skewness
##
## Figures produced:
##   corn_rot_plot_nc  — Response of corn yields to rotation sequences (no controls)
##   corn_rot_plot     — Response of corn yields to rotation sequences
##   score_yield       — Effect of rotation score (Z-vector) on yields
##   corn_var_plot     — Response of std dev of corn yields to rotation sequences
##   corn_coeff_plot   — Mean vs variance coefficients scatter (corn)
##   corn_jp_plot      — Just-Pope mean-variance decomposition (corn)
##   corn_yield_map    — Spatial map of corn yields (2016)
## ============================================================================

source("corn_data_prep.R")

# ── Shared p-value / vcov splicing helpers (rotation tables only) ────────────

fmt_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return(character(0))
  if (p < 0.001) return("[p<0.001]")
  sprintf("[p=%.3f]", p)
}

insert_pval <- function(cell, p) {
  p_str <- fmt_p(p)
  if (length(p_str) == 0) return(cell)
  sub("(\\([0-9.]+\\))", paste0("\\1 ", p_str), cell)
}

get_rot_pvals <- function(mod) {
  broom::tidy(mod) |>
    filter(grepl("^rot_crop", term)) |>
    transmute(term = gsub("^rot_crop", "", term), p.value)
}

# read_vcov_txt() -- write.table() writes coefficient names verbatim, but a
# double-wrapped "I(I(pr_6^2))" was saved instead of "I(pr_6^2)" -- fix that on
# load so the names line up exactly with the fitted models' coefficient names,
# which etable() requires for a vcov matrix to be accepted.
read_vcov_txt <- function(path) {
  m <- as.matrix(read.table(path, header = TRUE, row.names = 1,
                             check.names = FALSE))
  fix_names <- function(nm) gsub("^I\\(I\\((.*)\\)\\)$", "I(\\1)", nm)
  colnames(m) <- fix_names(colnames(m))
  rownames(m) <- fix_names(rownames(m))
  m
}

# ── 1. OLS mean model ─────────────────────────────────────────────────────────
# Table: tab:corn_rot | Figure: corn_rot_plot

corn_yield_formula <- make_jp_formula("corn_yield", "rot_crop", all_controls)

feols(corn_yield ~ rot_crop | tile_field_ID + year,
      data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_rot_nc

feols(corn_yield_formula,
      data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_rot

# Order rot_crop coefficient rows by estimate (with-controls model), decreasing
rot_order <- broom::tidy(corn_rot) |>
  filter(grepl("^rot_crop", term)) |>
  arrange(desc(estimate)) |>
  pull(term) |>
  gsub("^rot_crop", "", x = _)
rot_order_regex <- paste0("^", rot_order, "$")

etable(corn_rot_nc, corn_rot,
       tex      = TRUE,
       dict     = dict_corn,
       headers  = c("No controls", "With controls"),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_", "vpd_",
                    "rootznaws", "Constant"),
       order    = rot_order_regex,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       fontsize = "scriptsize",
       arraystretch = 0.8,
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

# Splice p-values into corn_rot.tex, next to each coefficient's SE
pvals_nc   <- get_rot_pvals(corn_rot_nc)
pvals_full <- get_rot_pvals(corn_rot)

corn_rot_tex_path <- paste0(tab_dir, "corn_rot.tex")
tex_lines <- readLines(corn_rot_tex_path)

tex_lines <- vapply(tex_lines, function(line) {
  label_match <- regmatches(line, regexpr("^\\s*[A-Za-z0-9\\-]+(?=\\s{2,})", line, perl = TRUE))
  if (length(label_match) == 0) return(line)
  term_label <- trimws(label_match)

  p_nc   <- pvals_nc$p.value[pvals_nc$term     == term_label]
  p_full <- pvals_full$p.value[pvals_full$term == term_label]
  if (length(p_nc) == 0 && length(p_full) == 0) return(line)

  parts <- strsplit(line, "&", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(line)

  parts[2] <- insert_pval(parts[2], p_nc)
  parts[3] <- insert_pval(parts[3], p_full)

  paste(parts, collapse = "&")
}, character(1), USE.NAMES = FALSE)

writeLines(tex_lines, corn_rot_tex_path)

# Figure: corn_rot_plot — Response of corn yields to rotation sequences
library(tidytext)   # for reorder_within / scale_x_reordered

corn_rot_nc |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(
    system = if_else(grepl("W", term), "Corn-soy-wheat", "Corn-soy"),
    group  = case_when(
      term == "S-C-S-C-S-C" ~ "Perfect rotation",
      .default              = "Other")) |>
  ggplot(aes(x = reorder_within(term, -prms.y, system), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_x_reordered() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of corn yields to rotation sequences (no controls)",
       caption = "Reference: C-C-C-C-C-C. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom")  ->
  corn_rot_plot_nc
ggsave(paste0(fig_dir, "corn_rot_plot_nc.png"), corn_rot_plot_nc,
       width = 10, height = 7.5, dpi = 300)

corn_rot |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(
    system = if_else(grepl("W", term), "Corn-soy-wheat", "Corn-soy"),
    group  = case_when(
      term == "S-C-S-C-S-C" ~ "Perfect rotation",
      .default              = "Other")) |>
  ggplot(aes(x = reorder_within(term, -prms.y, system), y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_x_reordered() +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of corn yields to rotation sequences",
       caption = "Reference: C-C-C-C-C-C. Clustered at COUNTY_FIPS.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_rot_plot
ggsave(paste0(fig_dir, "corn_rot_plot.png"), corn_rot_plot,
       width = 10, height = 7.5, dpi = 300)

# ── 2. VPD interaction — rotation, corn ───────────────────────────────────────
# Table: tab:corn_rot_vpd

corn_vpd_formula <- make_jp_formula("corn_yield", "rot_crop + vpd_name",
                                     vpd_controls)

feols(corn_yield ~ rot_crop + vpd_name | tile_field_ID + year,
      data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_rot_vpd_nc

corn_jp_data |>
  feols(corn_vpd_formula, data = _, cluster = ~COUNTY_FIPS+year) -> corn_rot_vpd

# Order rot_crop coefficient rows by estimate, decreasing
rot_vpd_order <- broom::tidy(corn_rot_vpd) |>
  filter(grepl("^rot_crop", term)) |>
  arrange(desc(estimate)) |>
  pull(term) |>
  gsub("^rot_crop", "", x = _)
rot_vpd_order_regex <- paste0("^", rot_vpd_order, "$")

etable(corn_rot_vpd,
       tex      = TRUE,
       dict     = c(dict_corn, dict_vpd),
       drop     = c("pr_", "cGDD_", "EDD_", "soil_",
                    "rootznaws", "Constant"),
       order    = rot_vpd_order_regex,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE, se.below = FALSE,
       fontsize = "scriptsize",
       title    = "Effect of weather and rotation sequences on corn yields",
       label    = "tab:corn_rot_vpd",
       file     = paste0(tab_dir, "corn_rot_vpd.tex"))

gc()

# ── 3. Just-Pope stage 1 — corn ───────────────────────────────────────────────
# Table: tab:corn_jp_mean | Figures: corn_rot_plot (mean), score_yield

# Lag comparison models (confirmatory Z-vector pre-registration)
fml_corn_lag1      <- make_jp_formula("corn_yield", "soy_lag1", all_controls_fgls)
fml_corn_lag1_lag2 <- make_jp_formula("corn_yield", "soy_lag1 + soy_lag2",
                                       all_controls_fgls)
fml_corn_index     <- make_jp_formula("corn_yield", "rot_index", all_controls_fgls)
fml_corn_mean      <- make_jp_formula("corn_yield", "rot_crop", all_controls_fgls)

feols(fml_corn_lag1,      data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_jp_s1_lag1
feols(fml_corn_lag1_lag2, data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_jp_s1_lag2
feols(fml_corn_index,     data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_jp_s1_idx
feols(fml_corn_mean,      data = corn_jp_data, cluster = ~COUNTY_FIPS+year) -> corn_jp_s1

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

# Bootstrap vcov matrices, from jp_boot_vcov.R -> tables/*.txt (as tab_dir).
jp_vcov_mean <- read_vcov_txt(paste0(tab_dir, "jp_vcov_mean.txt"))

# Table: tab:corn_jp_mean
etable(corn_jp_s1, corn_jp_s1_lag1, corn_jp_s1_lag2, corn_jp_s1_idx,
       tex      = TRUE,
       keep     = c("^[CSW]-", "soy_lag", "rot_index"),
       dict     = rot_dict,
       vcov     = list(jp_vcov_mean, NULL, NULL, NULL),
       notes    = bh_note_corn,
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 1 — Corn yield: full sequences vs lag summary variables",
       label    = "tab:corn_jp_mean",
       file     = paste0(tab_dir, "corn_jp_mean.tex"))

# Splice BH-adjusted p-values into corn_jp_mean.tex, next to each rot_crop coefficient's SE
pvals_jp_mean <- broom::tidy(corn_jp_s1) |>
  filter(term %in% names(rot_dict)) |>
  transmute(term_label = rot_dict[term],
            p.value = p.adjust(p.value, method = "BH"))

corn_jp_mean_tex_path <- paste0(tab_dir, "corn_jp_mean.tex")
tex_lines <- readLines(corn_jp_mean_tex_path)

tex_lines <- vapply(tex_lines, function(line) {
  label_match <- regmatches(line, regexpr("^\\s*[A-Za-z0-9\\-]+(?=\\s{2,})", line, perl = TRUE))
  if (length(label_match) == 0) return(line)
  term_label <- trimws(label_match)

  p <- pvals_jp_mean$p.value[pvals_jp_mean$term_label == term_label]
  if (length(p) == 0) return(line)

  parts <- strsplit(line, "&", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(line)

  parts[2] <- insert_pval(parts[2], p)

  paste(parts, collapse = "&")
}, character(1), USE.NAMES = FALSE)

writeLines(tex_lines, corn_jp_mean_tex_path)

rm(corn_jp_s1_lag1, corn_jp_s1_lag2, corn_jp_s1_idx, corn_rot_nc, corn_rot)
gc()

# ── Z-vector model (confirmatory spec) ───────────────────────────────────────
# Four structural features: late_soy, soy_gap, soy_cons (nsoy dropped —
# collinear with soy_gap, r=0.65). Run on corn data, both mean and variance
# stages. Feeds the rotation-score figure (score_yield) below.
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
  ) |>
  select(-seq_vec)

# Z-vector stage 1 — corn mean
fml_z_corn_mean <- make_jp_formula("corn_yield",
                                    "late_soy + soy_gap + soy_cons",
                                    all_controls_fgls)

feols(fml_z_corn_mean, data = corn_jp_data,
      cluster = ~COUNTY_FIPS+year) -> corn_z_s1

# Z-vector stage 2 — corn variance
corn_jp_data <- corn_z_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid_sq_z = (corn_yield - .fitted)^2) |>
  select(-starts_with("."))

fml_z_corn_var <- make_jp_formula("resid_sq_z",
                                   "late_soy + soy_gap + soy_cons",
                                   all_controls_fgls)

feols(fml_z_corn_var, data = corn_jp_data,
      cluster = ~COUNTY_FIPS+year) -> corn_z_s2

# ── 4. Just-Pope stage 2 — corn ───────────────────────────────────────────────
# Table: tab:corn_jp_moments | Figures: corn_var_plot, corn_coeff_plot, corn_jp_plot

corn_jp_s1 |>
  augment(newdata = corn_jp_data) |>
  mutate(resid = corn_yield - .fitted,
         resid_sq = (corn_yield - .fitted)^2) |>
  select(-starts_with(".")) ->
  corn_jp_data

fml_corn_var <- make_jp_formula("resid_sq", "rot_crop", all_controls_fgls)

# Stage 2a: OLS variance regression, used only to build the FGLS weights (h_hat)
feols(fml_corn_var, data = corn_jp_data) -> corn_jp_s2a

corn_jp_s2a |>
  augment(newdata = corn_jp_data) |>
  mutate(h_hat = pmax(.fitted, 1e-6)) |>
  select(-starts_with(".")) ->
  corn_jp_data

# Stage 2b: FGLS variance regression, weighted by 1/h_hat -- matches the
# estimator whose bootstrap vcov (jp_vcov_var, computed the same way in
# just_pope_bootstrap_moments.R) is attached below, so the SEs in
# tab:corn_jp_moments correspond to the same estimator as the point estimates.
feols(fml_corn_var, data = corn_jp_data, weights = ~I(1/h_hat)) -> corn_jp_s2

rm(corn_jp_s2a); gc()

# Table tab:corn_jp_moments is produced jointly with stage-3 skewness below
# (combined into a single table) — see the "Stage 3" section.

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
    .groups  = "drop"
  ) |>
  mutate(
    score = z_coefs["late_soy"] * late_soy +
            z_coefs["soy_gap"]  * soy_gap  +
            z_coefs["soy_cons"] * soy_cons
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
  select(rot_crop, mean_est, mean_se, var_est, var_se) |>
  print(n = Inf)

# Figure: corn_var_plot — Response of std dev of corn yields to rotation sequences
corn_jp_s2 |>
  coefplot() |>
  data.frame() |>
  filter(grepl("rot_crop", prms.estimate_names)) |>
  separate(prms.estimate_names, into = c("temp", "term"), sep = 8) |>
  select(-temp) |>
  mutate(
    system = if_else(grepl("W", term), "Corn-soy-wheat", "Corn-soy"),
    group  = case_when(
      term == "S-C-S-C-S-C" ~ "Perfect rotation",
      .default              = "Other")) |>
  arrange(system, prms.y) |>
  mutate(term = factor(term, levels = unique(term))) |>
  ggplot(aes(x = term, y = prms.y)) +
  geom_point(size = 2, aes(color = group)) +
  geom_errorbar(aes(ymin = prms.ci_low, ymax = prms.ci_high, color = group),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient Estimate", y = "Crop sequence",
       title   = "Response of standard deviation of corn yields to rotation sequences",
       caption = "Reference: C-C-C-C-C-C. Two-way clustering: field + year.") +
  theme(legend.title = element_blank(), legend.position = "bottom") ->
  corn_var_plot
ggsave(paste0(fig_dir, "corn_var_plot.png"), corn_var_plot,
       width = 10, height = 7.5, dpi = 300)

# Figure: corn_coeff_plot — Mean vs variance coefficients
corn_jp_summary |>
  mutate(pr = factor(as.integer(rot_crop %in% c("S-C-S-C-S-C")))) |>
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
    perfect = rot_crop %in% c("S-C-S-C-S-C")
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

# corn_jp_s2 is kept alive (not rm()'d here) — its variance estimate is used
# to standardize stage-3 skewness below, and both models are reported jointly
# in a single table (tab:corn_jp_moments).
corn_stage2_var <- mean(fitted(corn_jp_s2))

rm(corn_s1_coef, corn_s2_coef, corn_jp_summary, corn_jp_plot,
   corn_coeff_plot, corn_var_plot); gc()

# ── Stage 3: conditional skewness / downside risk (corn) ──────────────────────
# Antle (1983) moment-based extension. resid is already on corn_jp_data from the
# stage-2 join; the third central moment is the natural downside-risk statistic:
#   E[(y - mu)^3 | X].  Coefficient > 0  => sequence shifts mass toward the RIGHT
#   tail relative to monoculture (LESS downside risk, insurer-favorable);
#   Coefficient < 0  => heavier LEFT tail (MORE downside risk, loss-cost relevant).
corn_jp_data <- corn_jp_data |>
  mutate(resid_cube = resid^3)

# Standardize by stage-2 variance^{3/2} so coefficients/SEs are on a
# standardized-skewness scale: rescaling the LHS by a constant scales OLS
# coefficients and SEs by that same constant.
skew_scale <- corn_stage2_var^1.5

corn_jp_data <- corn_jp_data |>
  mutate(resid_cube_std = resid_cube / skew_scale)

fml_corn_skew <- make_jp_formula("resid_cube_std", "rot_crop", all_controls)

corn_jp_s3 <- feols(
  fml_corn_skew,
  data    = corn_jp_data
)

# jp_vcov_var.txt was bootstrapped off the FGLS (1/h_hat weighted) variance
# stage in just_pope_bootstrap_moments.R; corn_jp_s2 above is now fit the same
# way, so its point estimates and this vcov correspond to the same estimator.
# jp_vcov_skew.txt is on the RAW resid_cube scale; corn_jp_s3 regresses the
# standardized resid_cube_std = resid_cube / skew_scale, so Var(std coef) =
# Var(raw coef) / skew_scale^2 -- rescaled below before use.
jp_vcov_var  <- read_vcov_txt(paste0(tab_dir, "jp_vcov_var.txt"))
jp_vcov_skew <- read_vcov_txt(paste0(tab_dir, "jp_vcov_skew.txt")) / skew_scale^2

# Table: tab:corn_jp_moments — stage-2 variance and stage-3 (standardized) skewness,
# reported jointly as a single table.
etable(corn_jp_s2, corn_jp_s3,
       tex      = TRUE,
       keep     = "^[CSW]-",
       dict     = rot_dict,
       vcov     = list(jp_vcov_var, jp_vcov_skew),
       headers  = c("Variance", "Skewness (standardized)"),
       se.below = FALSE,
       placement = "H",
       style.tex = style.tex("aer"),
       replace  = TRUE,
       title    = "Stage 2/3 — Corn yield conditional variance (FGLS) and skewness",
       label    = "tab:corn_jp_moments",
       file     = paste0(tab_dir, "corn_jp_moments.tex"))

rm(corn_jp_s2, corn_jp_s3); gc()

# Collinearity check on the retained Z-vector features (documents why nsoy was
# dropped: in the 4-feature version nsoy correlated with soy_gap at r=0.65,
# making its partial sign unstable across samples).
cor(dplyr::select(corn_jp_data, late_soy, soy_cons, soy_gap))
# feature ranges
corn_jp_data |>
  dplyr::summarise(dplyr::across(c(late_soy, soy_cons, soy_gap),
                                 list(min = min, max = max, mean = mean)))

# ── 5. FGLS + bootstrap — corn (rotation-sequence stages) ────────────────────

source("jp_boot_vcov.R")

# ── 6. Spatial map — corn yields ──────────────────────────────────────────────
# Figure: corn_yield_map
# See load_corn_sf_2016() in corn_data_prep.R for caveats about this reload
# path (rci_correction()/annual_RCI dependency, pre-existing).

corn_spatial <- load_corn_sf_2016()
corn_sf <- corn_spatial$corn_sf
il_map  <- corn_spatial$il_map

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
       width = 10, height = 12.5, dpi = 600)

rm(corn_spatial, corn_sf, il_map, corn_yield_map); gc()

# ── 7. County-level yields and maximum entropy distribution ───────────────────

cat("Reloading corn data for county-level analysis...\n")
corn_df <- read_parquet(
  "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet")

corn_county <- corn_df |>
  filter(STATE_ABBR == "IL") |>
  select(COUNTY_FIPS, STATE_FIPS, year, corn_yield) |>
  mutate(corn_yield = corn_yield / 62.77) |>  # kg/ha -> bu/ac
  group_by(COUNTY_FIPS, STATE_FIPS, year) |>
  summarise(yield = mean(corn_yield, na.rm = TRUE), .groups = "drop") |>
  rename(state = STATE_FIPS, county = COUNTY_FIPS)
rm(corn_df); gc()

source("./max_entropy/maxent_tack2013.r")
corn_det <- detrend_yields_proportional(corn_county)
dist_m   <- tack_table2(corn_det$yield_norm)

dist_m$table
rm(corn_county, corn_det); gc()

cat("\n=== Corn rotation analysis complete ===\n")
cat("All tables saved to:", tab_dir, "\n")
cat("All figures saved to:", fig_dir, "\n")
