## ==========================================================================
## CORN YIELD -- SLX ESTIMATION (final, stripped)
## ==========================================================================
## Reported model:
##   corn_yield ~ regime/RCI + WS + W*WS | COUNTY_FIPS + year
##     regime/RCI : nested -> per-regime RCI intercept & slope
##     WS         : weather/soil covariates, own-field
##     W*WS       : spatial lags of weather/soil only (spillover channel)
##   rho fixed at 0 (SLX). SEs: 25km spatial-grid clusters.
##
## Design choices (justification in the diagnostics companion, slx_corn_regime.R):
##   - rho = 0: SLM had strong first stage but rho_hat ~ 0; freeing rho in the
##     SDM needs W^2 X (weak, Wald ~10) and returns inadmissible rho_hat = 1.72.
##   - Rotation own-field only (not spatially lagged); weather/soil lagged.
##   - Nested regime/RCI so each slope is on its own RCI support.
##   - County x year FE (keeps regime + static spatial lags identified).
##   - No SARAR/Conley/field FE: within-year residual Moran ~ 0 (no spatial
##     error); ordinary clustered SEs suffice.
## Assumes the cleaned pipeline has produced: corn_df with corn_yield, covars,
## W_<ws> lags, X_c/Y_c; regime built (6 levels).
## ==========================================================================
library(arrow)
library(tidyverse)
library(statar)
library(fixest)
library(marginaleffects)
library(Matrix)
library(dbscan)
library(sf)
options(scipen = 999) 

setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

source("rotation_setup_wa.R")

# ── Load corn data ────────────────────────────────────────────────────────────
# Load, correct, and add degree-day variables in one pipeline.
# corn_df is kept as the raw data source; corn_jp_data is the analysis sample.

cat("Loading corn data...\n")
corn_df <- read_parquet(
  "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet") 

corn_df <- corn_df |>
  filter(STATE_ABBR == "IL") |> 
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         corn_yield = corn_yield / 62.77)  |> 
  arrange(tile_field_ID, year) 

cat("Corn raw rows:", nrow(corn_df), "\n")

## Add present year crop variable
setDT(corn_df)  

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

add_crop_year(corn_df)

corn_df <- corn_df |>
  rename(crop_0 = crop_year,
         crop_1 = prioryr_crop,
         crop_2 = prior2yr_crop,
         crop_3 = prior3yr_crop,
         crop_4 = prior4yr_crop,
         crop_5 = prior5yr_crop,
         crop_6 = prior6yr_crop) 

source("rci_vectorized.R")

corn_df <- corn_df |>
    mutate(RCI = rci(crop_0, crop_1, crop_2, crop_3, crop_4, crop_5))     

source("cdl_recode.R")

recode_cdl(corn_df, cols = paste0("crop_", 0:6))

corn_df <- corn_df |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-", 
          crop_2, "-", crop_1, "-", crop_0)) |>
  add_degree_days()

# ── Build analysis sample ─────────────────────────────────────────────────────
# Join rotation features, add lag dummies, recode to C/S labels,
# filter to complete cases on all controls.
 
corn_df <- corn_df |>
  left_join(
    rot_features |> select(pattern, rot_index),
    by = c("rot_crop" = "pattern")
  ) |>
  filter(!is.na(corn_yield)) |>
  filter(if_all(all_of(all_controls_cols), ~ !is.na(.))) |>
  mutate(vpd_name = case_when(
    vpdmax_7 >= 0   & vpdmax_7 < 1.9  ~ "normal",
    vpdmax_7 >= 1.9 & vpdmax_7 <= 2.1 ~ "somewhat dry",
    vpdmax_7 > 2.1                     ~ "dry",
    .default = NA_character_),
    vpd_name = factor(vpd_name, levels = c("normal", "somewhat dry", "dry")))
 
cat("Corn analysis sample:", nrow(corn_df), "rows\n")

xcols <- paste0("crop_", 0:5)

corn_code     <- 1L
soy_code      <- 5L
# Cereals broadly: wheat (all classes) + other small grains + sorghum,
# plus their double-crop combinations with corn/soy (MATCH has_cereal below)
cereal_codes  <- c(
  21L,  22L,  23L,  24L,  25L,  27L,  28L,  29L,  30L,  205L,  4L,   # standalone
  26L,  225L, 226L, 228L, 236L, 237L, 238L, 240L, 254L,               # dbl-crop w/ corn or soy
  234L, 235L                                                          # dbl-crop sorghum + cereal
)
alfalfa_codes <- 36L

# The residual "other" bucket turned out to be ~89% Grassland/Pasture (176) and
# Fallow/Idle Cropland (61) -- i.e. idled land, not genuine crop diversification.
# Split it: land taken out of active production vs. perennial forage vs. actual
# minor/specialty row crops. See corn_RCI.r conversation notes for the frequency
# breakdown that motivated this split.
fallow_pasture_codes <- c(176L, 61L)                # Grassland/Pasture, Fallow/Idle Cropland
forage_codes         <- c(36L, 37L, 58L)             # Alfalfa, Other Hay/Non-Alfalfa, Clover/Wildflowers

core_codes    <- c(corn_code, soy_code, cereal_codes, alfalfa_codes)

# non-crop land covers that exist in corn_df but NOT in corn_jp_data — must be excluded
nonag_codes <- c(0L, 63L, 64L, 65L, 81L, 82L, 83L, 87L, 88L, 111L, 112L,
                 121L, 122L, 123L, 124L, 131L, 141L, 142L, 143L, 152L, 190L, 195L)

all_codes        <- sort(unique(unlist(corn_df[, ..xcols], use.names = FALSE)))
other_crop_codes <- setdiff(all_codes, c(core_codes, nonag_codes, fallow_pasture_codes, forage_codes))

corn_df[, has_other := Reduce(`|`, lapply(.SD, \(v) v %in% other_crop_codes)),
        .SDcols = xcols]

corn_df[, has_cereal := Reduce(`|`, lapply(.SD, \(v) v %in% cereal_codes)),
        .SDcols = xcols]

corn_df[, has_alfalfa := Reduce(`|`, lapply(.SD, \(v) v %in% alfalfa_codes)),
        .SDcols = xcols]

corn_df[, has_fallow_pasture := Reduce(`|`, lapply(.SD, \(v) v %in% fallow_pasture_codes)),
        .SDcols = xcols]

corn_df[, has_forage := Reduce(`|`, lapply(.SD, \(v) v %in% forage_codes)),
        .SDcols = xcols]

corn_df[, has_soy := Reduce(`|`, lapply(.SD, \(v) v %in% soy_code)),
        .SDcols = xcols]

# Perfect corn/soy rotation = strict yearly alternation over the 6-yr window
# (crop_0 is always corn, so this means crop_1,3,5 = soy and crop_2,4 = corn)
corn_df[, perfect_cs := crop_0 == corn_code & crop_1 == soy_code &
                         crop_2 == corn_code & crop_3 == soy_code &
                         crop_4 == corn_code & crop_5 == soy_code]

# corn monoculture, corn/soy rotations (perfect vs other), corn/soy/cereal,
# other crops (forage + specialty merged), fallow/pasture.
# Precedence: forage/specialty checked before fallow/pasture, since perennial
# hay/alfalfa stands are frequently mis-coded as Grassland/Pasture (176) in
# establishment or dormant years -- when both appear, the crop is the more
# informative label than the incidental pasture code.
corn_df[, regime := fcase(
  has_forage | has_other,      "corn_other_crops",       # checked first — wins
  has_fallow_pasture,          "corn_fallow_pasture",
  has_cereal,                  "corn_soy_cereal",
  has_soy & perfect_cs,        "corn_soy_perfect",
  has_soy,                     "corn_soy_other",
  default =                    "corn_monoculture"
)]
corn_df[, regime := factor(regime, levels = c(
  "corn_monoculture","corn_soy_perfect","corn_soy_other","corn_soy_cereal",
  "corn_other_crops","corn_fallow_pasture"))]

corn_df[, .N, by = regime][order(-N)]
 
corn_df <- corn_df |>
  filter(!is.na(regime))

# Free raw data — no longer needed
gc() 

## --------------------------------------------------------------------------
## 0.  Settings
## --------------------------------------------------------------------------
corn_df <- corn_df |>
  mutate(pr_6_sq = pr_6^2,
         pr_7_sq = pr_7^2,
         pr_8_sq = pr_8^2)

corn_df <- corn_df[!is.na(RCI)] 

# Remove misclassified observations
corn_df <- corn_df |>
    filter(!(regime == "corn_monoculture" & RCI != 0))

## Build spatial lags
k <- 10L; CRS_M <- 5070
covars <- c("RCI",
            "pr_6","pr_7","pr_8","pr_6_sq","pr_7_sq","pr_8_sq",
            "cGDD_6m","cGDD_7m","cGDD_8m",
            "EDD_6","EDD_7","EDD_8",
            "vpd_6","vpd_7","vpd_8",
            "soil_6","soil_7","soil_8",
            "rootznaws_mean")
ws  <- setdiff(covars, "RCI")              # weather/soil = everything lagged
Wws <- paste0("W_", ws)

## --------------------------------------------------------------------------
## 1.  Reconcile keys and attach projected coordinates to corn_df
## --------------------------------------------------------------------------
## `fields`  : data.table with tile_field_ID, lat, lon
## `corn_df` : yield panel with field_id (bare num) + tile column, corn_yield, ...
 
## 1a. Rebuild corn's full key. ONLY prepend a constant prefix if `fields`
##     uses a single prefix; otherwise carry corn's real prefix through.

fields <- unique(corn_df[, .(tile_field_ID, lat, lon)])
setkey(fields, tile_field_ID)
stopifnot(!anyNA(fields), uniqueN(fields$tile_field_ID) == nrow(fields))

prefixes <- unique(sub("_h[0-9]+v[0-9]+_.*$", "", fields$tile_field_ID))
stopifnot(length(prefixes) == 1L)                       # guard: single prefix
if (!grepl("^T", corn_df$tile_field_ID[1]))             # not already prefixed
  corn_df[, tile_field_ID := paste0(prefixes, "_", tile_field_ID)]
 
## 1b. Project field coordinates lat/lon -> Albers meters
geo <- unique(fields[, .(tile_field_ID, lat, lon)])
pts <- st_transform(st_as_sf(geo, coords = c("lon", "lat"), crs = 4326), CRS_M)
xy  <- st_coordinates(pts)
geo[, `:=`(X_c = xy[, 1], Y_c = xy[, 2])]
 
## 1c. Join coords onto corn; every corn row must match
setkey(geo, tile_field_ID)
corn_df <- geo[, .(tile_field_ID, X_c, Y_c)][corn_df, on = "tile_field_ID"]
stopifnot(!anyNA(corn_df$X_c))                          # no unmatched fields
 

## --------------------------------------------------------------------------
## 2.  Row-standardized sparse W over the field set (one coord pair per field)
## --------------------------------------------------------------------------
fields2 <- unique(corn_df[, .(tile_field_ID, X_c, Y_c)]); setkey(fields2, tile_field_ID)
stopifnot(uniqueN(fields2$tile_field_ID) == nrow(fields2))
## sanity: coords in a sane Albers range, no (0,0) or wild outliers
stopifnot(all(is.finite(fields2$X_c)), all(is.finite(fields2$Y_c)))
 
coords <- as.matrix(fields2[, .(X_c, Y_c)])
n      <- nrow(fields2)
knn    <- dbscan::kNN(coords, k = k)                    # KD-tree; minutes at ~3e5
W      <- sparseMatrix(i = rep(seq_len(n), each = k),
                       j = as.vector(t(knn$id)),
                       x = 1 / k, dims = c(n, n))        # exactly k/row => row-stochastic
row_of <- setNames(seq_len(n), fields2$tile_field_ID)
## kNN W is asymmetric. For symmetry: W <- (W + t(W))/2; then re-standardize rows.


## 3. spatial lags of weather/soil (per year; W^2 never materialized) ----
spatial_lag <- function(data, var, W, row_of,
                        id_col = "tile_field_ID", time_col = "year", order = 1L) {
  ids <- as.character(data[[id_col]]); times <- data[[time_col]]
  out <- numeric(nrow(data))
  for (yr in unique(times)) {
    idx <- which(times == yr); r <- row_of[ids[idx]]
    if (anyNA(r)) stop("id/row_of mismatch")
    v  <- numeric(nrow(W)); v[r] <- data[[var]][idx]
    lv <- as.numeric(W %*% v)
    if (order == 2L) lv <- as.numeric(W %*% lv)
    out[idx] <- lv[r]
  }
  out
}
for (v in ws) corn_df[, (paste0("W_", v)) := spatial_lag(corn_df, v, W, row_of, order = 1L)]

stopifnot(all(Wws %in% names(corn_df)))    # ready for slx_corn_estimation.R    

## --- covariate blocks -----------------------------------------------------
## Weather/soil: enter own-field AND as spatial lags.
ws <- c("pr_6","pr_7","pr_8","pr_6_sq","pr_7_sq","pr_8_sq",
        "cGDD_6m","cGDD_7m","cGDD_8m",
        "EDD_6","EDD_7","EDD_8",
        "vpd_6","vpd_7","vpd_8",
        "soil_6","soil_7","soil_8",
        "rootznaws_mean")
Wws  <- paste0("W_", ws)
fe   <- c("COUNTY_FIPS", "year")
BASE <- "corn_monoculture"
 
corn_df[, regime := relevel(factor(regime), ref = BASE)]
corn_df[, grid_cell := paste(X_c %/% 25000, Y_c %/% 25000, sep = "_")]
stopifnot(all(Wws %in% names(corn_df)), uniqueN(corn_df$grid_cell) >= 50)

regime_dict <- setNames(levels(corn_df$regime), paste0("regime", levels(corn_df$regime)))
 
## --- fit (no collinearity drop now: constant-RCI regimes were only a problem
##         because of the interaction; as level effects they're identified) ----
f_reg <- as.formula(sprintf("corn_yield ~ regime + %s + %s | %s",
                            paste(ws, collapse = " + "),
                            paste(Wws, collapse = " + "),
                            paste(fe, collapse = " + ")))
m_reg <- feols(f_reg, data = corn_df, vcov = ~ grid_cell)
etable(m_reg, keep = "^corn", 
        dict = regime_dict, 
        tex = TRUE,
        file = paste0(tab_dir, "slx_corn_regime.tex"), 
        replace = TRUE)
 
## ==========================================================================
## RESULTS
## ==========================================================================
b <- coef(m_reg); V <- vcov(m_reg); se <- sqrt(diag(V))
 
## (1) Regime yield gaps = plain coefficient contrasts (no RCI to hold fixed) --
.d <- function(r) {
  L <- setNames(numeric(length(b)), names(b))
  if (r != BASE) L[paste0("regime", r)] <- 1
  L
}
regs <- levels(corn_df$regime)
cmb  <- combn(regs, 2)
gaps <- rbindlist(lapply(seq_len(ncol(cmb)), function(k) {
  a <- cmb[1, k]; z <- cmb[2, k]; L <- .d(a) - .d(z)
  est <- sum(L * b); s <- sqrt(as.numeric(t(L) %*% V %*% L))
  data.table(regime_a = a, regime_b = z, gap = est, se = s,
             z = est / s, p = 2 * pnorm(-abs(est / s)))
}))[order(-abs(gap))]
 
## regime levels vs baseline (just the dummy coefficients + SEs) --------------
reg_levels <- rbindlist(lapply(setdiff(regs, BASE), function(r) {
  k <- paste0("regime", r)
  data.table(regime = r, vs_base = b[k], se = se[k],
             z = b[k]/se[k], p = 2*pnorm(-abs(b[k]/se[k])))
}))[order(-vs_base)]
 
## (2) Weather/soil impacts -- exact SLX decomposition (rho = 0) --------------
impacts_ws <- rbindlist(lapply(seq_along(ws), function(i) {
  lk <- ws[i]; wk <- Wws[i]
  data.table(var = lk,
             direct = b[lk], direct_se = se[lk],
             indirect = b[wk], indirect_se = se[wk],
             total = b[lk] + b[wk],
             total_se = sqrt(V[lk,lk] + V[wk,wk] + 2*V[lk,wk]))
}))
impacts_ws[, total_p := 2 * pnorm(-abs(total/total_se))]
 
#options(scipen = 999)
cat("\n== Regime levels vs baseline ==\n");  print(reg_levels, digits = 4)
cat("\n== Pairwise regime gaps ==\n");        print(gaps,       digits = 4)
cat("\n== Weather/soil impacts ==\n");        print(impacts_ws, digits = 4)

## (3) Pairwise regime gaps -> LaTeX table (\input-able; matches house style) ---
##     writes tables/slx_corn_regime_gaps.tex : one row per regime pair, sorted
##     by |gap|, estimate + significance stars, delta-method SE beneath.
pretty_lab <- c(
  corn_monoculture    = "Corn monoculture",
  corn_soy_perfect    = "Perfect corn-soybean",
  corn_soy_other      = "Other corn-soybean",
  corn_soy_cereal     = "Corn-soybean-cereal",
  corn_other_crops    = "Corn + forage/specialty",
  corn_fallow_pasture = "Corn + fallow/pasture"
)
.star <- function(p) ifelse(p < .01, "$^{***}$",
                     ifelse(p < .05, "$^{**}$",
                     ifelse(p < .1,  "$^{*}$", "")))

.gaps_body <- unlist(lapply(seq_len(nrow(gaps)), function(i) {
  g <- gaps[i]
  c(sprintf("   %s $-$ %s & %.3f%s\\\\",
            pretty_lab[[g$regime_a]], pretty_lab[[g$regime_b]],
            g$gap, .star(g$p)),
    sprintf("   & (%.3f)\\\\", g$se))
}))

writeLines(c(
  "",
  "\\begingroup",
  "\\centering",
  "\\begin{tabular}{lc}",
  "   \\tabularnewline \\midrule \\midrule",
  "   Contrast & Difference (bu/acre)\\\\",
  "   \\midrule",
  .gaps_body,
  "   \\midrule \\midrule",
  "   \\multicolumn{2}{l}{\\emph{Delta-method standard errors, 25\\,km grid-cell clustered, in parentheses}}\\\\",
  "   \\multicolumn{2}{l}{\\emph{Signif. Codes: ***: 0.01, **: 0.05, *: 0.1}}\\\\",
  "\\end{tabular}",
  "\\par\\endgroup",
  ""
), paste0(tab_dir, "slx_corn_regime_gaps.tex"))
cat("wrote ", paste0(tab_dir, "slx_corn_regime_gaps.tex"), "\n", sep = "")

## ==========================================================================
## CHECK (RESOLVED): the hybrid is DEGENERATE -> regime-only is the answer.
## ==========================================================================
## We tested a hybrid (regime levels for all 6 + RCI slopes for the 4 regimes
## with RCI variation). It FAILED for two independent reasons, so its
## "significant" RCI slopes are artifacts, not evidence RCI matters:
##
##   (a) LEVEL/SLOPE COLLINEARITY. Within each varying regime RCI is nearly
##       constant (15-18 distinct values, sd 0.33-0.68), so the regime level
##       dummy and its regime_v:RCI slope are near-linearly dependent. Result:
##       regime level SEs exploded to ~51 (vs 0.24-0.61 in regime-only) while
##       the slopes looked tight -- the SUM is identified, the SPLIT is not.
##       The slopes are absorbing level differences, not a yield response to
##       rotation intensity.
##   (b) SAMPLE HALVING. regime_v is NA for the two constant-RCI regimes, and
##       regime_v:RCI drops those rows -> n fell 2,288,053 -> 1,136,956,
##       discarding all of corn_monoculture and corn_soy_perfect. The "hybrid"
##       silently became a 4-regime model on half the data.
##
## CONCLUSION: RCI is ~constant within regime, so it carries no identifiable
## within-regime response; its apparent signal in any interaction model is
## level differences leaking through the RCI-regime collinearity. The regime
## LEVEL effects (m_reg above) are where that variation belongs. Report m_reg.
## This also removes the interaction degeneracy that made fixest drop a slope
## arbitrarily (flipping corn_monoculture vs corn_soy_perfect across sessions).
##
## The block below is kept ONLY to reproduce the failure for the record; it is
## NOT a reported specification. Note the ~51 SEs and the halved Observations.
varying <- c("corn_soy_other","corn_soy_cereal","corn_fallow_pasture","corn_other_crops")

corn_df[, regime_v := factor(fifelse(as.character(regime) %in% varying,
                                     as.character(regime), NA_character_))]
regime_dict_v <- setNames(levels(corn_df$regime_v), paste0("regime_v", levels(corn_df$regime_v)))                                     
f_hyb <- as.formula(sprintf("corn_yield ~ regime + regime_v:RCI + %s + %s | %s",
                            paste(ws, collapse = " + "),
                            paste(Wws, collapse = " + "),
                            paste(fe, collapse = " + ")))
m_hyb <- feols(f_hyb, data = corn_df, vcov = ~ grid_cell)
etable(m_reg, m_hyb, keep = c("^corn", "RCI"),
       headers = c("regime-only", "hybrid"),
       dict = c(regime_dict, regime_dict_v), tex = TRUE,
       file = paste0(tab_dir, "slx_corn_regime_hybrid.tex"), replace = TRUE)
## Do NOT report m_hyb or a Wald test on its slopes -- the SEs beneath it are
## unstable (level/slope collinearity), so any joint test is uninterpretable.
 
## (A) Is there even a first stage? Regress RCI on the proposed instruments,
##     within the variation FE leave. If weak, no RCI IV can work.
first <- feols(RCI ~ nccpi3corn_mean + soc0_100_mean | COUNTY_FIPS + year,
               data = corn_df, vcov = ~ grid_cell)
fitstat(first, "f")            # partial F of the excluded instruments
## RCI is near-constant within regime, so expect this to be weak once regime-
## driven variation is accounted for.

## (B) The better use: land-quality controls to test the selection story.
f_ctrl <- as.formula(sprintf(
  "corn_yield ~ regime + nccpi3corn_mean + soc0_100_mean + %s + %s | %s",
  paste(ws, collapse=" + "), paste(Wws, collapse=" + "),
  paste(fe, collapse=" + ")))
m_ctrl <- feols(f_ctrl, data = corn_df, vcov = ~ grid_cell)
etable(m_reg, m_ctrl, keep = "^corn",
       headers = c("baseline", "+ land-quality controls"),
       dict = regime_dict, tex = TRUE,
       file = paste0(tab_dir, "slx_corn_regime_ctrl.tex"), replace = TRUE)
## Regime gaps STABLE -> ranking isn't land-quality sorting (strengthens paper).
## Regime gaps SHRINK -> land quality explained part of it (also a real finding).
