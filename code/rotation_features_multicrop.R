## ============================================================================
## rotation_features_multicrop.R
## Generalized rotation-feature construction for corn / soy / alfalfa / wheat.
## Drop-in replacement for the corn-soy-specific block in rotation_setup.R
## (the expand.grid at lines 33-42 and the PCA feature block at lines 179-229).
##
## Design change vs. the binary version: the old features were all "soy vs.
## not-soy" (n_soy, consec_soy, gap_soy, min_gap_soy). With four crops there is
## no single target crop, so the feature set is rebuilt around measures that are
## defined for any number of crops -- richness, diversity, transitions, longest
## monoculture run -- plus legume-spacing terms that GENERALIZE the old soy-gap
## features (soybean and alfalfa are both N-fixers, which is the mechanism that
## matters for the corn-yield mean/variance story).
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
library(arrow)
theme_set(theme_bw())

setwd("C:/Users/vf006/Box/Economic Analysis of Soil Health Practices")

tab_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/tables/"
fig_dir <- "C:/Users/vf006/Box/crop_rotations_and_losses/figures/"

corn_data <- read_parquet("D:/Crop data/d_igis13_12_1_2025.with_rci.parquet")
corn_data |>
  filter(STATE_ABBR == "IL") ->
  corn_data

setDT(corn_data)  

add_crop_year <- function(dt) {
  stopifnot(is.data.table(dt))
  
  # Drop any stale crop_year column from a previous run, so the regex below
  # can't accidentally re-capture it
  if ("crop_year" %in% names(dt)) {
    dt[, crop_year := NULL]
  }
  
  # Match crop_YYYY columns only (exactly 4 digits) -- excludes crop_year,
  # rot_seq, or anything else that happens to start with "crop_"
  crop_cols  <- grep("^crop_[0-9]{4}$", names(dt), value = TRUE)
  crop_years <- as.integer(sub("crop_", "", crop_cols))
  
  stopifnot(length(crop_cols) > 0, !anyNA(crop_years))
  
  m       <- as.matrix(dt[, ..crop_cols])
  col_idx <- match(dt$year, crop_years)
  
  stopifnot(length(col_idx) == nrow(dt))
  
  dt[, crop_year := m[cbind(seq_len(.N), col_idx)]]
  invisible(dt)
}

add_crop_year(corn_data)

corn_data <- corn_data |>
  rename(crop_0 = crop_year,
         crop_1 = prioryr_crop,
         crop_2 = prior2yr_crop,
         crop_3 = prior3yr_crop,
         crop_4 = prior4yr_crop,
         crop_5 = prior5yr_crop,
         crop_6 = prior6yr_crop) 

## ============================================================================
## rotation_features_functional.R
## Multi-crop rotation features built from the RAW crop_1..crop_6 columns.
##
## Why raw, not prior*yr_functional: the existing functional recode collapses
## "Dbl Crop WinWht/Soybeans" (3.66% of slots) to "soybeans", erasing the wheat.
## Winter Wheat proper is only 0.85%, so the double crop IS the wheat signal --
## losing it makes wheat unanalyzable. We rebuild the recode from the raw names.
##
## Design: each slot is encoded as MULTI-HOT attribute flags (legume /
## small-grain / perennial) rather than one categorical, because a wheat-soy
## double crop is legitimately both a small grain and a legume in one year.
## For measures that need a single category (richness, diversity, transitions)
## the double crop gets its own class so it reads as a distinct rotational state.
##
## Requires data.table. crop_0 is always Corn in the estimation sample (you need
## a corn crop to observe a corn yield), so all of this describes the 6-year
## HISTORY; slots are ordered crop_1 = t-1 (most recent) ... crop_6 = t-6.
## ============================================================================

# ── Functional crosswalk from raw CDL names ───────────────────────────────────
# The primary "class" per raw name comes from cdl_functional_classify.R -- the
# project's established functional_group vocabulary (corn, soybeans, ley,
# annual_cereal, annual_legume, annual_broadleaf, fallow, or NA for
# non-agricultural / unclassified land) -- shared with cdl_functional_recode.R
# and rotation_setup_wa.R so every consumer agrees on what "ley" or
# "annual_cereal" means. Edit the name -> group rules there, not here.
#
# is_legume / is_small_grain / is_perennial are attribute flags LOCAL to this
# file's multi-hot rotation-feature design (not part of the shared
# classification): "ley" bundles alfalfa (legume, perennial) with grass hay /
# pasture (perennial, not legume), so is_legume treats the whole "ley" class
# as legume-bearing -- a coarser approximation than a raw-name-level split
# would give, but consistent with how rotation_setup_wa.R's legume_groups
# treats "ley". Grassland/Pasture always classifies to "ley" now (there is no
# separate noncrop option for it in the real vocabulary), so a class of NA
# means genuinely non-agricultural land; it's folded into a "noncrop" level
# below purely so every raw value lands in one of the categorical levels used
# for richness/diversity/transitions.
source("C:/Users/vf006/Box/crop_rotations_and_losses/code/cdl_functional_classify.R")

build_class_lut <- function(raw_vals) {
  lut <- classify_cdl_names(raw_vals)
  lut[, class := fifelse(is.na(class), "noncrop", class)]
  lut[, is_legume      := as.integer(class %in% c("soybeans", "annual_legume", "ley"))]
  lut[, is_small_grain := as.integer(class == "annual_cereal")]
  lut[, is_perennial   := as.integer(class == "ley")]
  lut[]
}

# ── Feature builder ───────────────────────────────────────────────────────────
# Vectorized over all rows; no per-row apply. Adds count/flag/gap/diversity
# columns describing the 6-year history.

add_rotation_features <- function(dt, crop_cols = paste0("crop_", 1:6)) {
  stopifnot(all(crop_cols %in% names(dt)))
  lut <- build_class_lut(unlist(lapply(crop_cols, function(c) dt[[c]])))

  # named lookups -> n x 6 matrices (columns = crop_1..crop_6 = lag1..lag6)
  flag_mat <- function(flag) {
    v <- lut[[flag]]; names(v) <- lut$raw
    m <- vapply(crop_cols, function(c) {
      x <- v[as.character(dt[[c]])]; fifelse(is.na(x), 0L, as.integer(x))
    }, integer(nrow(dt)))
    matrix(m, nrow = nrow(dt))
  }
  Lg <- flag_mat("is_legume")
  Sg <- flag_mat("is_small_grain")
  Pr <- flag_mat("is_perennial")

  # primary class code matrix
  cl <- lut$class; names(cl) <- lut$raw
  lev <- c("corn","soybeans","annual_cereal","annual_legume",
           "annual_broadleaf","ley","fallow","noncrop")
  code <- setNames(seq_along(lev), lev)
  nc   <- code[["noncrop"]]
  CM <- vapply(crop_cols, function(c) {
    k <- code[cl[as.character(dt[[c]])]]; fifelse(is.na(k), nc, as.integer(k))
  }, integer(nrow(dt)))
  CM <- matrix(CM, nrow = nrow(dt))

  crop_codes <- setdiff(unname(code), nc)              # classes that count as crops

  # counts / presence
  dt[, n_legume      := rowSums(Lg)]
  dt[, n_small_grain := rowSums(Sg)]
  dt[, n_perennial   := rowSums(Pr)]
  dt[, has_wheat     := as.integer(rowSums(Sg) > 0)]
  dt[, has_legume    := as.integer(rowSums(Lg) > 0)]
  dt[, noncrop_count := rowSums(CM == nc)]

  # years since most recent legume (1 = last year was a legume); NA if none
  ys <- max.col(Lg, ties.method = "first")
  ys[rowSums(Lg) == 0] <- NA_integer_
  dt[, yrs_since_legume := ys]

  # transitions (adjacent class changes across the 6 slots)
  dt[, n_transitions := rowSums(CM[, 1:5] != CM[, 2:6])]

  # crop richness (distinct crop classes, noncrop excluded)
  dt[, richness := Reduce(`+`, lapply(crop_codes,
        function(k) as.integer(rowSums(CM == k) > 0)))]

  # Shannon diversity over crop slots (noncrop excluded); NA if all noncrop
  counts <- lapply(crop_codes, function(k) rowSums(CM == k))
  tot <- Reduce(`+`, counts)
  sh  <- numeric(nrow(dt))
  for (cnt in counts) { p <- cnt / tot; sh <- sh - fifelse(p > 0, p * log(p), 0) }
  sh[tot == 0] <- NA_real_
  dt[, diversity := sh]

  dt[]
}

# ── Usage ─────────────────────────────────────────────────────────────────────
corn_data <- add_rotation_features(corn_data)          # adds the columns above
#
# Then build rot_index on OBSERVED histories only (fit the PCA on the unique
# feature rows that actually occur, not a factorial):
feat <- c("richness","diversity","n_transitions","n_legume",
          "n_small_grain","yrs_since_legume")

# confirm what's non-finite and why (should be all-noncrop / no-legume windows)
corn_data[, .N, by = .(all_noncrop = noncrop_count == 6L,
                       no_legume   = is.na(yrs_since_legume))]

X <- unique(
  corn_data[!is.na(corn_yield) & noncrop_count < 6L, ..feat]  # real rotations only
)
X[is.na(yrs_since_legume), yrs_since_legume := 7L]            # "no legume in 6 yrs"

stopifnot(all(is.finite(as.matrix(X))))                       # passes now
pca <- prcomp(X, scale. = TRUE)
print(summary(pca)$importance[, 1:3])
#   # CHECK PC1 SIGN before using it as a regressor (prcomp's sign is arbitrary).
#
# For the estimation sample the enrichment is entirely in these history features
# plus per-lag flags if you want them, e.g.:
#   for (k in 1:5) corn_data[, paste0("legume_lag", k) :=
#     as.integer(build_class_lut(get(paste0("crop_", k)))[
#       match(get(paste0("crop_", k)), raw), is_legume])]


feat <- c("richness","diversity","n_transitions","n_legume",
          "n_small_grain","yrs_since_legume")
 
# Unique feature rows on real rotations, carrying a grouping var for the biplot.
Xall <- unique(corn_data[!is.na(corn_yield) & noncrop_count < 6L,
                         c(feat, "has_wheat"), with = FALSE])
Xall[is.na(yrs_since_legume), yrs_since_legume := 7L]      # "no legume in window"
Xall[, wheat := factor(has_wheat, levels = c(0, 1),
                       labels = c("no wheat", "wheat in history"))]
 
pca_rot <- prcomp(Xall[, ..feat], scale. = TRUE)
cat("Variance explained by first three PCs:\n")
print(summary(pca_rot)$importance[, 1:3])
 
# ── Figure: rot_pca_plot — PCA biplot of multi-crop rotation features ──────────
# NB: check the PC1 sign against the loadings before reading rot_index as
# "more complex = higher"; prcomp's sign is arbitrary.
autoplot(pca_rot, data = Xall,
         colour = "wheat", alpha = 0.5, size = 1.6,
         loadings = TRUE,         loadings.colour       = "grey30",
         loadings.label = TRUE,   loadings.label.colour = "black",
         loadings.label.size = 4, loadings.label.repel  = TRUE) +
  scale_colour_manual(values = c("no wheat" = "#1f77b4",
                                 "wheat in history" = "#d62728"),
                      name = NULL) +
  labs(title   = "PCA of corn-soy-wheat rotation features",
       caption = "Each point is a distinct six-year rotation history observed in the data.") +
  theme_bw() +
  theme(legend.position = "bottom") ->
  rot_pca_plot
ggsave(paste0(fig_dir, "rot_pca_plot.png"), rot_pca_plot,
       width = 9, height = 7, dpi = 300)
 
# ── Figure: rot_scree_plot — variance explained per component ──────────────────
data.frame(PC  = seq_along(pca_rot$sdev),
           pve = pca_rot$sdev^2 / sum(pca_rot$sdev^2)) |>
  ggplot(aes(x = PC, y = pve)) +
  geom_col(fill = "grey70", width = 0.7) +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  geom_point(colour = "#1f77b4", size = 2.5) +
  scale_x_continuous(breaks = seq_along(pca_rot$sdev)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Principal component", y = "Share of variance explained",
       title   = "Variance explained — rotation-feature PCA",
       caption = "Six standardized features; corn-soy-wheat histories.") +
  theme_bw() ->
  rot_scree_plot
ggsave(paste0(fig_dir, "rot_scree_plot.png"), rot_scree_plot,
       width = 7, height = 5, dpi = 300)