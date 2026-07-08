## ============================================================================
## rotation_setup.R
## Shared setup for corn and soy Just-Pope rotation analysis.
## Source this file at the top of each analysis script:
##   source("rotation_setup.R")
## ============================================================================

library(tidyverse)
library(data.table)
library(statar)
library(fixest)
library(broom)
library(haven)
library(marginaleffects)
library(hdm)
library(dotwhisker)
library(ggfortify)
library(ggrepel)
library(patchwork)
library(knitr)
library(kableExtra)
library(sf)
library(usmap)
theme_set(theme_bw())

setwd("C:/Users/vf006/Box/Economic Analysis of Soil Health Practices")

tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"
fig_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/figures/"

# ── Rotation patterns ─────────────────────────────────────────────────────────

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

# ── RCI correction ────────────────────────────────────────────────────────────
# Removes mismatched RCI / rot_crop combinations (data quality fix).
# Apply to any loaded crop data frame before analysis.

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

# ── Degree-day functions (Schlenker-Roberts 2009) ─────────────────────────────

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

# ── Control vector ────────────────────────────────────────────────────────────
# nccpi3all_mean and soc0_100_mean excluded — collinear with field/FIPS FEs.

all_controls <- c(
  "pr_6", "pr_7", "pr_8",
  "I(pr_6^2)", "I(pr_7^2)", "I(pr_8^2)",
  "cGDD_6m", "cGDD_7m", "cGDD_8m",
  "EDD_6", "EDD_7", "EDD_8",
  "vpd_6", "vpd_7", "vpd_8",
  "soil_6", "soil_7", "soil_8",
  "rootznaws_mean"
)

# Controls for FGLS bootstrap — drop collinear soil vars
all_controls_fgls <- setdiff(all_controls, c("nccpi3all_mean", "soc0_100_mean"))

# Base column names that exist in the data — for NA filtering.
# Strips formula expressions like I(pr_6^2) which are not actual columns.
all_controls_cols <- all_controls_fgls[!grepl("^I\\(", all_controls_fgls)]

# ── Formula helper ────────────────────────────────────────────────────────────

make_jp_formula <- function(lhs, rot_var, controls,
                            fe = "tile_field_ID + year") {
  rhs <- paste(c(rot_var, controls), collapse = " + ")
  as.formula(paste(lhs, "~", rhs, "|", fe))
}

# ── Coefficient label dictionaries ────────────────────────────────────────────

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

dict_vpd <- c(
  "vpd_namesomewhatdry"           = "Somewhat dry season",
  "vpd_namedry"                   = "Dry season",
  "RCI x vpd_namesomewhatdry"     = "RCI x Somewhat dry season",
  "RCI x vpd_namedry"             = "RCI x Dry season"
)

# ── PCA of rotation features ──────────────────────────────────────────────────
# Run once; rot_features and pca_rot are shared by both analysis files.

count_soy_years <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  sum(rotation %in% c(2L, 5L))
}
consecutive_soy <- function(rotation) {
  if (is.character(rotation)) rotation <- as.integer(strsplit(rotation, "-")[[1]])
  soy_runs <- rle(rotation)
  lengths  <- soy_runs$lengths[soy_runs$values %in% c(2L, 5L)]
  if (rlang::is_empty(lengths)) 0L else ifelse(max(lengths) > 1, 1L, 0L)
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

rot_features <- rot_features |>
  mutate(rot_index = pca_rot$x[, 1])

# Figure: rot_pca_plot — PCA biplot of rotation features (Section 4.2)
autoplot(pca_rot, data = corn_soy_patterns,
         loadings = TRUE, loadings.label = TRUE) +
  labs(title   = "PCA of corn-soy rotation features",
       caption = "Each point represents one of the 64 binary six-year corn-soy sequences.") +
  theme_bw() +
  theme(legend.position = "none") ->
  rot_pca_plot
ggsave(paste0(fig_dir, "rot_pca_plot.png"), rot_pca_plot,
       width = 9, height = 7, dpi = 300)

# ── Bootstrap function (shared) ───────────────────────────────────────────────

boot_jp_fgls <- function(dt, fml_mean, fml_var, fml_var_b,
                          B = 499, seed = 42) {
  dt <- as.data.table(dt)

  # Remove singletons iteratively until stable
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

  # Point estimates on full data
  s1 <- feols(fml_mean, data = dt,
              cluster  = c("tile_field_ID", "year"),
              nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, resid_sq := NA_real_]
  dt[obs(s1), resid_sq := residuals(s1)^2]

  s2a <- feols(fml_var, data = dt,
               cluster  = c("tile_field_ID", "year"),
               nthreads = n_threads, warn = FALSE, notes = FALSE)
  dt[, h_hat := NA_real_]
  dt[obs(s2a), h_hat := pmax(fitted(s2a), 1e-6)]

  s2b <- feols(fml_var, data = dt,
               weights  = ~I(1 / h_hat),
               cluster  = c("tile_field_ID", "year"),
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
      dt_b[, h_hat_b    := NA_real_]
      dt_b[obs(b2a), h_hat_b := pmax(fitted(b2a), 1e-6)]
      dt_b[, combined_w := boot_w / h_hat_b]

      b2b <- feols(fml_var_b, data = dt_b, weights = ~combined_w,
                   vcov = "iid", nthreads = n_threads,
                   warn = FALSE, notes = FALSE)

      coef_boot[b, ] <- coef(b2b)

    }, error = function(e) {
      cat("\nBootstrap iteration", b, "failed:", conditionMessage(e), "\n")
    })

    dt[, boot_w := NULL]

    if (b %% 100 == 0) {
      saveRDS(coef_boot, paste0("boot_progress_", b, ".rds"))
      gc()
    }
  }
  cat("\n")

  boot_se <- apply(coef_boot, 2, sd,       na.rm = TRUE)
  boot_ci <- apply(coef_boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  list(coef = coef_hat, se = boot_se, ci = boot_ci,
       fit_mean = s1, fit_var_ols = s2a, fit_var_fgls = s2b,
       boot_draws = coef_boot)
}
