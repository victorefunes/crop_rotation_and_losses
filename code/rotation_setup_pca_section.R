# ── PCA of rotation features ──────────────────────────────────────────────────
# Replaces the soy-centric feature block. Features are now defined over the CDL
# codes in the pattern strings, so wheat and alfalfa contribute distinctly:
#   legume years  = soybeans (5) + alfalfa (36) + WinWht/Soy double crop (26)
#   small-grain yrs = winter wheat (24) + WinWht/Soy double crop (26)
# The double crop (26) counts as BOTH, because it is literally both crops.
# (For alfalfa as its own axis, add n_perennial = count of 36.)
#
# NOTE: this renames the old feature helpers (count_soy_years, consecutive_soy,
# soy_gap, smallest_soy_gap) and the feature columns (n_soy -> n_legume, etc.).
# If the analysis files reference the old names, update them. The object that
# propagates downstream, rot_index, is preserved.

legume_codes      <- c(5L, 26L, 36L)   # soybeans, WinWht/Soy dbl crop, alfalfa
small_grain_codes <- c(24L, 26L)       # winter wheat, WinWht/Soy dbl crop

parse_rot   <- function(x) as.integer(strsplit(x, "-", fixed = TRUE)[[1]])
shannon_div <- function(v) { p <- table(v) / length(v); -sum(p * log(p)) }
longest_run <- function(v) max(rle(v)$lengths)
count_in    <- function(v, codes) sum(v %in% codes)
gap_mean    <- function(v, codes) { p <- which(v %in% codes); if (length(p) < 2) 0 else mean(diff(p)) }
gap_min     <- function(v, codes) { p <- which(v %in% codes); if (length(p) < 2) 0 else min(diff(p)) }

# One row of features per unique pattern.
build_rot_features <- function(patterns) {
  tibble(pattern = unique(patterns)) |>
    mutate(
      rot_vec    = lapply(pattern, parse_rot),
      n_legume   = sapply(rot_vec, count_in, codes = legume_codes),      # was n_soy
      n_grain    = sapply(rot_vec, count_in, codes = small_grain_codes), # new: wheat
      diversity  = sapply(rot_vec, shannon_div),                         # new: multi-crop
      max_run    = sapply(rot_vec, longest_run),
      no_mono    = 6 - max_run,                                          # was no_consec
      leg_gap    = sapply(rot_vec, gap_mean, codes = legume_codes),
      leg_mingap = sapply(rot_vec, gap_min,  codes = legume_codes),
      free_leg   = ifelse(leg_gap    == 0, 0, 1 / leg_gap),              # was free_soy
      tight_leg  = ifelse(leg_mingap == 0, 0, 1 / leg_mingap)            # was tight_soy
    ) |>
    select(pattern, n_legume, n_grain, diversity, no_mono, free_leg, tight_leg)
}

rot_pca_features <- c("n_legume", "n_grain", "diversity",
                      "no_mono", "free_leg", "tight_leg")

# Build features + fit PCA on a supplied set of patterns; returns the feature
# table with rot_index (PC1) and the fitted prcomp attached as an attribute.
build_rot_pca <- function(patterns) {
  feats <- build_rot_features(patterns)
  X <- as.data.frame(feats[, rot_pca_features])
  # Drop zero-variance columns: if a crop never appears in the supplied set
  # (e.g. observed data with no wheat), scale. = TRUE would divide by 0 and
  # prcomp's svd would error. Guarding here keeps the fit well-defined.
  keep <- vapply(X, function(z) stats::sd(z) > 0, logical(1))
  pca  <- prcomp(X[, keep, drop = FALSE], scale. = TRUE)
  feats$rot_index <- pca$x[, 1]
  attr(feats, "pca")      <- pca
  attr(feats, "pca_cols") <- names(keep)[keep]
  feats
}

# Default fit on the full candidate list, so this file still runs standalone and
# rot_features / pca_rot exist after sourcing.
rot_features <- build_rot_pca(corn_soy_patterns$pattern)
pca_rot      <- attr(rot_features, "pca")

cat("Rotation PCA — variance explained by first two PCs:\n")
print(summary(pca_rot)$importance[, 1:2])

# RECOMMENDED for the paper: refit on the rotations that actually occur, so the
# index reflects the empirical distribution rather than the 5^6 = 15,625 factorial.
# In each analysis file, AFTER loading the data:
#   rot_features <- build_rot_pca(unique(corn_df$rot_crop))
#   pca_rot      <- attr(rot_features, "pca")
# then join rot_features into the analysis data by pattern = rot_crop as before.
#
# CHECK THE PC1 SIGN against pca_rot$rotation before using rot_index as a
# regressor — prcomp's sign is arbitrary. Flip with rot_index <- -rot_index if
# "more complex" should be the high end.

# Figure: rot_pca_plot — PCA biplot (rotation_plots.R has the fuller version).
autoplot(pca_rot, data = rot_features,
         loadings = TRUE, loadings.label = TRUE, loadings.label.repel = TRUE) +
  labs(title   = "PCA of rotation features (corn, soy, wheat, alfalfa)",
       caption = "Each point is a distinct six-year rotation sequence.") +
  theme_bw() +
  theme(legend.position = "none") ->
  rot_pca_plot
ggsave(paste0(fig_dir, "rot_pca_plot.png"), rot_pca_plot,
       width = 9, height = 7, dpi = 300)
