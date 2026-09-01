## ============================================================================
## corn_spatial.R
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

covars <- c("RCI", "pr_6", "pr_7", "pr_8", "pr_6_sq", "pr_7_sq", "pr_8_sq", "cGDD_6m",  # nolint: trailing_whitespace_linter.
            "cGDD_7m", "cGDD_8m", "EDD_6", "EDD_7", "EDD_8", "vpd_6", "vpd_7", "vpd_8", # nolint
            "soil_6", "soil_7", "soil_8", "rootznaws_mean")       # <- your exogenous regressors
fe     <- c("COUNTY_FIPS", "year")    # <- FE to absorb; use "0" for none
k      <- 6L                     # neighbours per field
CRS_M  <- 5070                   # CONUS Albers equal-area (meters)

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

## --------------------------------------------------------------------------
## 3.  Spatial lags (parameterized; fails loudly on any id mismatch)
## --------------------------------------------------------------------------
spatial_lag <- function(data, var, W, row_of,
                        id_col = "tile_field_ID", time_col = "year", order = 1L) {
  stopifnot(var %in% names(data), id_col %in% names(data))
  ids   <- as.character(data[[id_col]])
  times <- if (!is.null(time_col) && time_col %in% names(data)) data[[time_col]]
           else rep(1L, nrow(data))                      # cross-section fallback
  out <- numeric(nrow(data))
  for (yr in unique(times)) {
    idx <- which(times == yr)
    r   <- row_of[ids[idx]]
    if (anyNA(r))
      stop("Some ", id_col, " values aren't in row_of -- W/row_of and data ",
           "must share the same field ids (and type).")
    v  <- numeric(nrow(W)); v[r] <- data[[var]][idx]     # place into field order
    lv <- as.numeric(W %*% v)
    if (order == 2L) lv <- as.numeric(W %*% lv)          # W^2 v, matrix never formed
    out[idx] <- lv[r]
  }
  out
}
 
corn_df[, Wy := spatial_lag(corn_df, "corn_yield", W, row_of, order = 1L)]
for (v in covars) {
  corn_df[, (paste0("W_",  v)) := spatial_lag(corn_df, v, W, row_of, order = 1L)]
  corn_df[, (paste0("W2_", v)) := spatial_lag(corn_df, v, W, row_of, order = 2L)]
}

## --------------------------------------------------------------------------
## 4.  Estimate: FE absorbed + Wy instrumented, in one fixest IV call
## --------------------------------------------------------------------------
inst1 <- paste0("W_",  covars)      # WX   -- valid ONLY if the model is a pure SLM
inst2 <- paste0("W2_", covars)      # W^2 X
rhs   <- paste(covars, collapse = " + ")
fe_f  <- paste(fe, collapse = " + ")
 
## IDENTIFICATION:
##   Pure SLM    -> WX excludable; use inst1 (+ inst2 for efficiency).
##   Suspect SDM -> WX belongs in the equation, NOT a valid instrument:
##                  drop inst1, identify Wy off W^2 X (weaker -> check 1st-stage F).
f_slm <- as.formula(sprintf("corn_yield ~ %s | %s | Wy ~ %s",
                            rhs, fe_f, paste(c(inst1, inst2), collapse = " + ")))
 
m_slm <- feols(f_slm, data = corn_df, vcov = ~ COUNTY_FIPS)
## Spatial-HAC alt:  vcov = vcov_conley(...) with projected coords + cutoff.
 
etable(m_slm)                      # coefficient on 'fit_Wy' is rho
etable(m_slm, stage = 1)           # first stage -- confirm instruments are strong
fitstat(m_slm, "ivwald")            # weak-instrument check (critical on SDM path)

## --------------------------------------------------------------------------
## 5.  Average direct / indirect / total impacts  (large-n, trace approx.)
## --------------------------------------------------------------------------
## Pure SLM, row-stochastic W:
##   total_k    = beta_k / (1 - rho)
##   direct_k   = beta_k * (1/n) tr((I - rho W)^{-1}) = beta_k * sum_q rho^q (1/n) tr(W^q)
##   indirect_k = total_k - direct_k
## (1/n) tr(W^q) via Hutchinson probes -- no dense inverse, no W^q matrix.
mc_traces <- function(W, q_max = 120L, m = 100L, seed = 1L) {
  set.seed(seed); n <- nrow(W)
  tr <- numeric(q_max + 1L); tr[1] <- n                 # tr(W^0) = n
  for (l in seq_len(m)) {
    u <- sample(c(-1, 1), n, replace = TRUE); v <- u
    for (q in seq_len(q_max)) {
      v <- as.numeric(W %*% v)
      tr[q + 1L] <- tr[q + 1L] + sum(u * v)
    }
  }
  tr[-1] <- tr[-1] / m
  tr / n
}
 
rho   <- unname(coef(m_slm)["fit_Wy"])
q_max <- 120L                                            # raise if rho is large (~0.9)
tbar  <- mc_traces(W, q_max = q_max, m = 100L)
dir_mult <- sum(rho^(0:q_max) * tbar)
tot_mult <- 1 / (1 - rho)
 
b <- coef(m_slm)[covars]
impacts <- data.table(var = covars,
                      direct   = b * dir_mult,
                      indirect = b * (tot_mult - dir_mult),
                      total    = b * tot_mult)
print(impacts)

## ==========================================================================
## 6.  SLX  (rho = 0):  plain feols, no instrumenting, exact impacts
## ==========================================================================
lag_vars <- covars                     # WX terms to include as regressors

ws <- c("pr_6","pr_7","pr_8","pr_6_sq","pr_7_sq","pr_8_sq",
        "cGDD_6m","cGDD_7m","cGDD_8m",
        "EDD_6","EDD_7","EDD_8",
        "vpd_6","vpd_7","vpd_8",
        "soil_6","soil_7","soil_8",
        "rootznaws_mean")
Wws <- paste0("W_", ws)  
 
Wx  <- paste0("W_",  lag_vars)         # WX  -> SDM/SLX regressors
W2x <- paste0("W2_", covars)           # W^2 X -> excluded instruments for Wy
stopifnot(all(Wx  %in% names(corn_df)),
          all(W2x %in% names(corn_df)))   # lags must exist from the SLM script

rhs   <- paste(covars, collapse = " + ")
wxs   <- paste(Wx,     collapse = " + ")
fe_f  <- paste(fe,     collapse = " + ")

f_slx <- as.formula(sprintf("corn_yield ~ %s + %s | %s", rhs, wxs, fe_f))
m_slx <- feols(f_slx, data = corn_df, vcov = ~ COUNTY_FIPS)
## Conley (spatial-HAC) alternative -- carry lat/lon onto corn_df first, then:
##   vcov = vcov_conley(lat = "lat", lon = "lon", cutoff = 50)   # km
etable(m_slx, keep = c(covars, Wx))
 
## --- SLX impacts are EXACT (no (I-rhoW)^{-1}, because rho = 0) -------------
##   direct_k = beta_k ;  indirect_k = theta_k ;  total_k = beta_k + theta_k
b_slx  <- coef(m_slx)
impacts_slx <- data.table(
  var      = lag_vars,
  direct   = b_slx[lag_vars],
  indirect = b_slx[Wx],
  total    = b_slx[lag_vars] + b_slx[Wx]
)
print(impacts_slx, digits = 4)

## --------------------------------------------------------------------------
## 7.  Fit the SLX with nested regime/RCI
## --------------------------------------------------------------------------
f_slx <- as.formula(sprintf(
  "corn_yield ~ regime / RCI + %s + %s | %s",
  paste(ws,  collapse = " + "),
  paste(Wws, collapse = " + "),
  paste(fe,  collapse = " + ")
))
 
## Conley spatial-HAC SEs (needs lat/lon on corn_df; cutoff in km). If lat/lon
## aren't attached: corn_df <- geo[, .(tile_field_ID, lat, lon)][corn_df, on = "tile_field_ID"]
## Report sensitivity to a couple of cutoffs (25/50/100).
m_slx <- feols(f_slx, data = corn_df,
               vcov = vcov_conley(lat = "lat", lon = "lon", cutoff = 25))
 
## Proven fallback SE (within-county correlation across space and time):
m_slx_cl <- feols(f_slx, data = corn_df, vcov = ~ COUNTY_FIPS)
 
etable(m_slx, m_slx_cl,
       keep = c("regime", "RCI", ws, Wws),
       headers = c("SLX (Conley 50km)", "SLX (cluster county)"))
gc()

## --------------------------------------------------------------------------
## 8.  Joint tests
## --------------------------------------------------------------------------
wald(m_slx, keep = "^W_")            # are the weather/soil spillovers needed?
wald(m_slx, keep = "regime.*RCI|regime:RCI")   # does RCI slope differ by regime?
 
## --------------------------------------------------------------------------
## 8a. RCI impacts -- regime-specific, own-field only (no indirect: rotation
##     is not lagged, so direct = total for RCI within each regime).
## --------------------------------------------------------------------------
## NOTE: corn_soy_perfect has ~0 within-regime RCI variance, so fixest DROPS
## regimecorn_soy_perfect:RCI for collinearity -> that regime is intercept-only
## (no identifiable RCI response; report as such, NOT as a zero slope).
## Each remaining coefficient is that regime's ABSOLUTE slope on its own RCI
## support (regime/RCI drops the bare RCI term), not a contrast.
## RCI slope per regime = the regime:RCI coefficient. The model is LINEAR in
## RCI, so the "average slope" IS the coefficient -- do NOT run avg_slopes over
## 2.29M rows (it materializes per-row derivatives across the full frame and
## exhausts memory). Read coef + Conley SE directly; identical result, instant.
slope_tab <- function(model) {
  b <- coef(model); se <- sqrt(diag(vcov(model))); nm <- names(b)
  ks <- grep(":RCI$", nm, value = TRUE)
  rbindlist(lapply(ks, function(k) {
    data.table(regime = sub(":RCI$", "", sub("^regime", "", k)),
               slope = b[k], se = se[k],
               z = b[k] / se[k], p = 2 * pnorm(-abs(b[k] / se[k])))
  }))[order(slope)]
}
print(slope_tab(m_slx), digits = 4)
## corn_soy_perfect has no :RCI term (dropped) -> no row -> no identifiable slope.
 
## Report the RCI support each slope is estimated on -- the slopes are LOCAL
## derivatives on nearly non-overlapping RCI bands, not one global curve bending.
print(corn_df[, .(rci_min = min(RCI), rci_med = median(RCI),
                  rci_max = max(RCI), sd_RCI = sd(RCI), n = .N),
              by = regime][order(rci_med)])

## --------------------------------------------------------------------------
## 8b. Regime yield LEVELS at a common RCI
## --------------------------------------------------------------------------
## Do NOT read the regime dummies as yield gaps -- with the bare RCI term gone
## they are intercepts at RCI = 0, far outside support. Evaluate at a real RCI.
##
## (i) Regime yield levels at RCI = m, RELATIVE TO BASELINE, from coefficients.
##     avg_predictions with grid_type="counterfactual" replicates the full
##     2.29M-row frame per regime -> memory blowup; and FE uncertainty can't
##     propagate through predictions anyway. The baseline-relative level is a
##     coefficient contrast (absorbed FE/intercept cancel), so it's cheap AND
##     carries a valid Conley SE. index_r = D_r + S_r*m ; baseline = S_base*m.
m_rci <- median(corn_df$RCI)
regime_levels <- function(model, regimes, m, baseline = "corn_monoculture") {
  b <- coef(model); V <- vcov(model); nm <- names(b)
  load_reg <- function(r) {
    L <- setNames(numeric(length(nm)), nm)
    d <- paste0("regime", r); s <- paste0("regime", r, ":RCI")
    if (r != baseline && d %in% nm) L[d] <- 1
    if (s %in% nm)                  L[s] <- m
    L
  }
  Lb <- load_reg(baseline)
  rbindlist(lapply(regimes, function(r) {
    L <- load_reg(r) - Lb
    est <- sum(L * b); seL <- sqrt(as.numeric(t(L) %*% V %*% L))
    data.table(regime = r, level_vs_base = est, se = seL,
               z = est / seL, p = 2 * pnorm(-abs(est / seL)))
  }))[order(-level_vs_base)]
}
print(regime_levels(m_slx, levels(corn_df$regime), m_rci), digits = 4)
 
## (ii) Pairwise gaps WITH valid SEs. A regime-vs-regime gap is a coefficient
##      contrast; the absorbed intercept/FE cancel, so the FE-uncertainty issue
##      is irrelevant and the Conley vcov gives a proper SE.
##      Fitted yield (net of common FE) for regime r at RCI = m:
##         baseline: S_r * m           (dummy = 0)
##         other:    D_r + S_r * m
##      corn_soy_perfect has no slope term -> contributes only D_r (no RCI response).
regime_gaps <- function(model, regimes, m, baseline = "corn_monoculture") {
  b <- coef(model); V <- vcov(model); nm <- names(b)
  load_reg <- function(r) {
    L <- setNames(numeric(length(nm)), nm)
    d <- paste0("regime", r); s <- paste0("regime", r, ":RCI")
    if (r != baseline && d %in% nm) L[d] <- 1
    if (s %in% nm)                  L[s] <- m      # absent for perfect -> 0
    L
  }
  cmb <- combn(regimes, 2)
  rbindlist(lapply(seq_len(ncol(cmb)), function(k) {
    a <- cmb[1, k]; z <- cmb[2, k]
    L <- load_reg(a) - load_reg(z)
    est <- sum(L * b); seL <- sqrt(as.numeric(t(L) %*% V %*% L))
    data.table(regime_a = a, regime_b = z, gap = est, se = seL,
               z = est / seL, p = 2 * pnorm(-abs(est / seL)))
  }))[order(-abs(gap))]
}
## gap > 0 => regime_a yields more than regime_b at RCI = m, ceteris paribus.
gaps <- regime_gaps(m_slx, levels(corn_df$regime), m_rci)
print(gaps, digits = 4)
 
## --------------------------------------------------------------------------
## 8c. Weather/soil impacts -- EXACT SLX decomposition (rho = 0):
##     direct = beta, indirect = theta (= W_ coef), total = beta + theta.
##     total SE via Var(beta)+Var(theta)+2Cov(beta,theta).
## --------------------------------------------------------------------------
b <- coef(m_slx); V <- vcov(m_slx); se <- sqrt(diag(V))
impacts_ws <- rbindlist(lapply(seq_along(ws), function(i) {
  lk <- ws[i]; wk <- Wws[i]
  var_tot <- V[lk, lk] + V[wk, wk] + 2 * V[lk, wk]
  data.table(var = lk,
             direct = b[lk], direct_se = se[lk],
             indirect = b[wk], indirect_se = se[wk],
             total = b[lk] + b[wk], total_se = sqrt(var_tot))
}))
impacts_ws[, total_p := 2 * pnorm(-abs(total / total_se))]
options(scipen = 999); print(impacts_ws, digits = 4)
 
## --------------------------------------------------------------------------
## 9.  Residual spatial autocorrelation (any spatial error left => SARAR?)
## --------------------------------------------------------------------------
resid_moran <- function(model, data, W, row_of,
                        id_col = "tile_field_ID", time_col = "year") {
  e   <- resid(model, na.rm = FALSE)
  ids <- as.character(data[[id_col]]); tm <- data[[time_col]]
  ok  <- !is.na(e); den <- sum(e[ok]^2); num <- 0
  for (yr in unique(tm)) {
    idx <- which(tm == yr); r <- row_of[ids[idx]]
    ei <- e[idx]; ei[is.na(ei)] <- 0
    v  <- numeric(nrow(W)); v[r] <- ei
    We <- as.numeric(W %*% v)
    keep <- !is.na(e[idx])
    num  <- num + sum(e[idx][keep] * We[r][keep])
  }
  num / den
}
cat("SLX residual Moran-type coef:",
    round(resid_moran(m_slx, corn_df, W, row_of), 5), "\n")
 
## ==========================================================================
## Reporting: RCI effect is the four regime-specific slopes (Section 3a),
## caveated by the within-regime sd from Section 0. Weather/soil spillovers
## are the total-effect column of Section 3b. rho is fixed at 0 by the trail
## in the header. SARAR only if the residual Moran coef is non-trivial (and
## even then Conley SEs already protect inference).
## ==========================================================================
 

# 1. Is it soil-shaped? Check whether residuals correlate with a smooth
#    spatial trend or with finer soil data you may have held out.
corn_df[, e := resid(m_slx, na.rm = FALSE)]

build_W <- function(coords, k) {
  n   <- nrow(coords)
  knn <- dbscan::kNN(coords, k = k)
  sparseMatrix(i = rep(seq_len(n), each = k),
               j = as.vector(t(knn$id)),
               x = 1 / k, dims = c(n, n))       # row-stochastic
}

coords <- as.matrix(fields2[, .(X_c, Y_c)])            # same order as row_of
row_of <- setNames(seq_len(nrow(fields2)), fields2$tile_field_ID)

for (kk in c(3, 6, 12, 24)) {
  Wk <- build_W(coords, kk)
  cat("k =", kk, " Moran =",
      round(resid_moran(m_slx, corn_df, Wk, row_of), 4), "\n")
}

## --------------------------------------------------------------------------
## 8.  DIAGNOSIS + FIX: dynamic spatial-error process -> grid-cluster SEs
## --------------------------------------------------------------------------
## Not scale (Sec 5), not a smooth surface (Sec 6), not field heterogeneity
## (Sec 7) => a TIME-VARYING, spatially-correlated shock: sub-county/sub-monthly
## weather texture the monthly county-scale weather covariates miss. This is a
## genuine spatial-ERROR process, not a mean-model omission. beta/theta stay
## consistent; only inference needs fixing -> grid-cluster SEs (Conley-infeasible
## at this n; see the fit block).
##
## 8a. Confirm it's weather: residual autocorrelation should spike in extreme
##     weather years (drought/flood) and be mild in benign years.
W6 <- build_W(coords, 6)   

for (yr in sort(unique(corn_df$year))) {
  cat(yr, ":", round(resid_moran(m_ffe, corn_df[year == yr], W6, row_of), 4), "\n")
}
 
## 8b. Match the grid cell size to the (long) correlation range: sweep sizes and
##     report where SEs stabilize. This REPLACES the infeasible Conley-cutoff sweep.
for (km in c(25, 50, 100)) {
  corn_df[, gc_tmp := paste(X_c %/% (km * 1000), Y_c %/% (km * 1000), sep = "_")]
  m <- feols(f_slx, data = corn_df, vcov = ~ gc_tmp)             # SE-only refit: cheap
  cat("grid", km, "km (", uniqueN(corn_df$gc_tmp), "cells): se(mono RCI slope) =",
      round(se(m)["regimecorn_monoculture:RCI"], 4), "\n")
}
corn_df[, gc_tmp := NULL]