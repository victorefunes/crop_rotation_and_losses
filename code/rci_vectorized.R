# =============================================================================
# RCI (Rotational Complexity Index) — vectorized rewrite
# Drop-in replacement for rci() in RCI_func_char2.R
#
# Fixes vs. original:
#   1. Turnover is a CHANGE COUNT, not popcount of as.binary(diff(...)).
#      The original weighted each crop change by the Hamming weight of the
#      arithmetic difference of CDL codes, which inflated turnover for any
#      transition whose code gap is not a power of two (e.g. Alfalfa->Corn:
#      |36-1|=35 -> popcount 3 instead of 1). This version reproduces the
#      RCI_unique column exactly (0 mismatches over 14,865 sequences).
#   2. No apply(., 1, .): all row work is vectorized -> ~O(n) column ops.
#   3. No binaryLogic / car / sjlabelled dependencies.
#   4. Single source-of-truth code table; duplicate name->code entries
#      (Shrubland 64/152, Barren 65/131) collapsed to the first (as case_match did).
#   5. Input crop names are trimmed (fixes unmatched "Clouds/No Data " etc.).
# =============================================================================

# --- crop name -> CDL code (first-wins, matching the original case_match) -----
crop_codes <- c(
  "Background"=0L, "Corn"=1L, "Cotton"=2L, "Rice"=3L,
  "Sorghum"=4L, "Soybeans"=5L, "Sunflower"=6L, "Peanuts"=10L,
  "Tobacco"=11L, "Sweet Corn"=12L, "Pop or Orn Corn"=13L, "Mint"=14L,
  "Barley"=21L, "Durum Wheat"=22L, "Spring Wheat"=23L, "Winter Wheat"=24L,
  "Other Small Grains"=25L, "Dbl Crop WinWht/Soybeans"=26L, "Rye"=27L, "Oats"=28L,
  "Millet"=29L, "Speltz"=30L, "Canola"=31L, "Flaxseed"=32L,
  "Safflower"=33L, "Rape Seed"=34L, "Mustard"=35L, "Alfalfa"=36L,
  "Other Hay/Non Alfalfa"=37L, "Camelina"=38L, "Buckwheat"=39L, "Sugarbeets"=41L,
  "Dry Beans"=42L, "Potatoes"=43L, "Other Crops"=44L, "Sugarcane"=45L,
  "Sweet Potatoes"=46L, "Misc Vegs & Fruits"=47L, "Watermelons"=48L, "Onions"=49L,
  "Cucumbers"=50L, "Chick Peas"=51L, "Lentils"=52L, "Peas"=53L,
  "Tomatoes"=54L, "Caneberries"=55L, "Hops"=56L, "Herbs"=57L,
  "Clover/Wildflowers"=58L, "Sod/Grass Seed"=59L, "Switchgrass"=60L, "Fallow/Idle Cropland"=61L,
  "Forest"=63L, "Shrubland"=64L, "Barren"=65L, "Cherries"=66L,
  "Peaches"=67L, "Apples"=68L, "Grapes"=69L, "Christmas Trees"=70L,
  "Other Tree Crops"=71L, "Citrus"=72L, "Pecans"=74L, "Almonds"=75L,
  "Walnuts"=76L, "Pears"=77L, "Clouds/No Data"=81L, "Developed"=82L,
  "Water"=83L, "Wetlands"=87L, "Nonag/Undefined"=88L, "Aquaculture"=92L,
  "Open Water"=111L, "Perennial Ice/Snow"=112L, "Developed/Open Space"=121L, "Developed/Low Intensity"=122L,
  "Developed/Med Intensity"=123L, "Developed/High Intensity"=124L, "Deciduous Forest"=141L, "Evergreen Forest"=142L,
  "Mixed Forest"=143L, "Grassland/Pasture"=176L, "Woody Wetlands"=190L, "Herbaceous Wetlands"=195L,
  "Pistachios"=204L, "Triticale"=205L, "Carrots"=206L, "Asparagus"=207L,
  "Garlic"=208L, "Cantaloupes"=209L, "Prunes"=210L, "Olives"=211L,
  "Oranges"=212L, "Honeydew Melons"=213L, "Broccoli"=214L, "Avocados"=215L,
  "Peppers"=216L, "Pomegranates"=217L, "Nectarines"=218L, "Greens"=219L,
  "Plums"=220L, "Strawberries"=221L, "Squash"=222L, "Apricots"=223L,
  "Vetch"=224L, "Dbl Crop WinWht/Corn"=225L, "Dbl Crop Oats/Corn"=226L, "Lettuce"=227L,
  "Dbl Crop Triticale/Corn"=228L, "Pumpkins"=229L, "Dbl Crop Lettuce/Durum Wht"=230L, "Dbl crop Lettuce/Cantaloupe"=231L,
  "Dbl Crop Lettuce/Cotton"=232L, "Dbl Crop Lettuce/Barley"=233L, "Dbl Crop Durum Wht/Sorghum"=234L, "Dbl Crop Barley/Sorghum"=235L,
  "Dbl Crop WinWht/Sorghum"=236L, "Dbl Crop Barley/Corn"=237L, "Dbl Crop WinWht/Cotton"=238L, "Dbl Crop Soybeans/Cotton"=239L,
  "Dbl Crop Soybeans/Oats"=240L, "Dbl Crop Corn/Soybeans"=241L, "Blueberries"=242L, "Cabbage"=243L,
  "Cauliflower"=244L, "Celery"=245L, "Radishes"=246L, "Turnips"=247L,
  "Eggplants"=248L, "Gourds"=249L, "Cranberries"=250L, "Dbl Crop Barley/Soybeans"=254L
)

#' @param x0..x5 character vectors (year-0 .. year-5 crop names), equal length
#' @return numeric vector of RCI, NA where the sequence contains a non-ag year
#'         or an unmappable crop name.
rci <- function(x0, x1, x2, x3, x4, x5,
                nonag     = c(0L, 81L, 82L, 83L, 86L, 88L),
                perennial = c(36L, 37L, 61L, 176L, 224L)) {

  M <- cbind(x0, x1, x2, x3, x4, x5)                 # n x 6 character
  M[] <- trimws(M)                                   # fix stray whitespace
  code <- matrix(crop_codes[M], nrow = nrow(M),      # names -> codes (NA if unmapped)
                 dimnames = NULL)
  code[code %in% nonag] <- NA_integer_               # non-ag -> NA (was set_na)
  ok <- rowSums(is.na(code)) == 0L                   # keep complete cases only

  # distinct crops per row (first-occurrence trick; no rowwise apply)
  firstocc <- matrix(TRUE, nrow(code), 6L)
  for (j in 2:6) {
    eq_prev <- code[, j] == code[, 1L]
    if (j >= 3L) for (k in 2:(j - 1L)) eq_prev <- eq_prev | code[, j] == code[, k]
    firstocc[, j] <- !eq_prev
  }
  number <- rowSums(firstocc)

  # turnover as plain change counts
  t1 <- rowSums(code[, 2:6] != code[, 1:5])          # first-order (lag 1)
  t2 <- rowSums(code[, 3:6] != code[, 1:4])          # second-order (lag 2)

  # perennial correction = adjacent equal pairs that are both perennial
  # (equivalent to sum(run_length - 1) over perennial runs; non-perennials break runs)
  isper <- matrix(code %in% perennial, nrow(code))
  corr  <- rowSums((code[, 2:6] == code[, 1:5]) & isper[, 2:6])

  turnover <- (t2 + t1 + corr) / 2
  rci <- round(sqrt(turnover * number), 2)
  rci[!ok] <- NA_real_
  rci
}

# --- example: score rot_frequencies.csv on UNIQUE sequences, then join back ---
# library(data.table)
# d  <- fread("rot_frequencies.csv")
# u  <- unique(d[, .(rot_crop)])
# p  <- tstrsplit(u$rot_crop, "-", fixed = TRUE)          # 6 columns
# u[, RCI_new := rci(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]])]
# d  <- u[d, on = "rot_crop"]                              # 17,907 calcs, not 2.3M rows
