## ============================================================================
## cdl_recode.R
## Replace CDL land-cover NAMES with their numeric USDA NASS codes.
## Legend source: USDA NASS CDL data dictionary (official cross-reference list).
## Applies to the crop-name columns (crop_1..crop_6, crop_YYYY); the lowercase
## *_functional columns use a different vocabulary and are intentionally skipped.
## ============================================================================

library(data.table)

# ── Full CDL legend (code -> name) ────────────────────────────────────────────
cdl_legend <- data.table(
  code = c(1L,2L,3L,4L,5L,6L,10L,11L,12L,13L,14L,
           21L,22L,23L,24L,25L,26L,27L,28L,29L,30L,31L,32L,33L,34L,35L,36L,37L,
           38L,39L,41L,42L,43L,44L,45L,46L,47L,48L,49L,50L,51L,52L,53L,54L,55L,
           56L,57L,58L,59L,60L,61L,62L,63L,64L,65L,
           66L,67L,68L,69L,70L,71L,72L,74L,75L,76L,77L,
           81L,82L,83L,87L,88L,92L,
           111L,112L,121L,122L,123L,124L,131L,141L,142L,143L,152L,176L,190L,195L,
           204L,205L,206L,207L,208L,209L,210L,211L,212L,213L,214L,215L,216L,217L,
           218L,219L,220L,221L,222L,223L,224L,225L,226L,227L,228L,229L,230L,231L,
           232L,233L,234L,235L,236L,237L,238L,239L,240L,241L,242L,243L,244L,245L,
           246L,247L,248L,249L,250L,254L),
  name = c("Corn","Cotton","Rice","Sorghum","Soybeans","Sunflower","Peanuts",
           "Tobacco","Sweet Corn","Pop or Orn Corn","Mint",
           "Barley","Durum Wheat","Spring Wheat","Winter Wheat","Other Small Grains",
           "Dbl Crop WinWht/Soybeans","Rye","Oats","Millet","Speltz","Canola",
           "Flaxseed","Safflower","Rape Seed","Mustard","Alfalfa",
           "Other Hay/Non Alfalfa","Camelina","Buckwheat","Sugarbeets","Dry Beans",
           "Potatoes","Other Crops","Sugarcane","Sweet Potatoes","Misc Vegs & Fruits",
           "Watermelons","Onions","Cucumbers","Chick Peas","Lentils","Peas",
           "Tomatoes","Caneberries","Hops","Herbs","Clover/Wildflowers",
           "Sod/Grass Seed","Switchgrass","Fallow/Idle Cropland","Pasture/Grass",
           "Forest","Shrubland","Barren",
           "Cherries","Peaches","Apples","Grapes","Christmas Trees","Other Tree Crops",
           "Citrus","Pecans","Almonds","Walnuts","Pears",
           "Clouds/No Data","Developed","Water","Wetlands","Nonag/Undefined",
           "Aquaculture",
           "Open Water","Perennial Ice/Snow","Developed/Open Space",
           "Developed/Low Intensity","Developed/Med Intensity","Developed/High Intensity",
           "Barren","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrubland",
           "Grassland/Pasture","Woody Wetlands","Herbaceous Wetlands",
           "Pistachios","Triticale","Carrots","Asparagus","Garlic","Cantaloupes",
           "Prunes","Olives","Oranges","Honeydew Melons","Broccoli","Avocados",
           "Peppers","Pomegranates","Nectarines","Greens","Plums","Strawberries",
           "Squash","Apricots","Vetch","Dbl Crop WinWht/Corn","Dbl Crop Oats/Corn",
           "Lettuce","Dbl Crop Triticale/Corn","Pumpkins","Dbl Crop Lettuce/Durum Wht",
           "Dbl Crop Lettuce/Cantaloupe","Dbl Crop Lettuce/Cotton",
           "Dbl Crop Lettuce/Barley","Dbl Crop Durum Wht/Sorghum",
           "Dbl Crop Barley/Sorghum","Dbl Crop WinWht/Sorghum","Dbl Crop Barley/Corn",
           "Dbl Crop WinWht/Cotton","Dbl Crop Soybeans/Cotton","Dbl Crop Soybeans/Oats",
           "Dbl Crop Corn/Soybeans","Blueberries","Cabbage","Cauliflower","Celery",
           "Radishes","Turnips","Eggplants","Gourds","Cranberries",
           "Dbl Crop Barley/Soybeans")
)

# Two names map to two codes each: "Barren" (65 old / 131 NLCD) and
# "Shrubland" (64 old / 152 NLCD). Post-2012 CDL archives use the NLCD-derived
# codes, so we drop the retired 64/65 from the name->code inversion. If your
# extract predates that recode, change these two lines.
name2code <- cdl_legend[!code %in% c(64L, 65L)]
name2code <- setNames(name2code$code, name2code$name)

# ── In-place recode ───────────────────────────────────────────────────────────
# Converts each targeted column from character names to integer CDL codes.
# Unmatched names become NA and are reported so nothing is dropped silently.
recode_cdl <- function(dt, cols = grep("^crop_", names(dt), value = TRUE)) {
  cols <- cols[vapply(cols, function(c) is.character(dt[[c]]), logical(1))]
  vals <- unique(unlist(lapply(cols, function(c) dt[[c]]), use.names = FALSE))
  miss <- setdiff(vals[!is.na(vals)], names(name2code))
  if (length(miss))
    warning("Unmatched names left as NA: ", paste(miss, collapse = ", "))
  for (c in cols) set(dt, j = c, value = unname(name2code[dt[[c]]]))
  invisible(dt[])
}

# ── Usage ─────────────────────────────────────────────────────────────────────
# recode_cdl(corn_data)                       # all crop_* name columns, in place
# recode_cdl(corn_data, cols = paste0("crop_", 1:6))   # or a specific subset
#
# Spot-check (expect: Corn=1, Soybeans=5, Winter Wheat=24,
#                     Dbl Crop WinWht/Soybeans=26, Alfalfa=36, Grassland=176):
# corn_data[, .N, by = crop_1][order(-N)]
