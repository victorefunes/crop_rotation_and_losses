## ============================================================================
## cdl_functional_recode.R
## Replace CDL land-cover NAMES with functional-group categories (corn,
## soybeans, ley, annual_cereal, annual_legume, annual_broadleaf, fallow, or
## NA for non-agricultural / unclassified land) -- the categorical counterpart
## to cdl_recode.R's numeric-code recode.
##
## Classification rules live in cdl_functional_classify.R (shared with
## rotation_setup_wa.R and rotation_features_multicrop.R's build_class_lut())
## so all consumers agree on what "annual_cereal" or "ley" means. That
## function derives its lookup only from names it observes in the data; this
## one applies it to the full CDL legend (plus the functional-group labels
## that already appear as raw values for out-of-CDL-coverage lag years) so it
## can run standalone. Edit the classification rules in
## cdl_functional_classify.R, not here.
##
## Applies to the crop-name columns (crop_1..crop_6, crop_YYYY).
## ============================================================================

library(data.table)
source("C:/Users/vf006/Box/crop_rotations_and_losses/code/cdl_functional_classify.R")

# ── Full CDL legend (names only; see cdl_recode.R for the code cross-walk) ────
cdl_names <- unique(c(
  "Corn","Cotton","Rice","Sorghum","Soybeans","Sunflower","Peanuts",
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
  "Dbl Crop Barley/Soybeans",
  # functional-group labels used directly as fallback values for lag years
  # outside real CDL coverage (see cdl_functional_classify.R header)
  "corn","soybeans","ley","fallow","annual_cereal","annual_legume",
  "annual_broadleaf","noncrop"
))

# ── Functional classification (rules in cdl_functional_classify.R) ───────────
lut <- classify_cdl_names(cdl_names)
name2group <- setNames(lut$class, lut$raw)

# ── In-place recode ───────────────────────────────────────────────────────────
# Converts each targeted column from character CDL names (or already-fallback
# functional-group labels) to functional-group strings. Names absent from the
# legend entirely become NA and are reported so nothing is dropped silently;
# names present but classified NA on purpose (e.g. "Forest", "noncrop") do
# NOT warn -- that's an intentional non-agricultural classification, not a
# data gap.
recode_cdl_functional <- function(dt, cols = grep("^crop_", names(dt), value = TRUE)) {
  cols <- cols[vapply(cols, function(c) is.character(dt[[c]]), logical(1))]
  vals <- unique(unlist(lapply(cols, function(c) dt[[c]]), use.names = FALSE))
  miss <- setdiff(vals[!is.na(vals)], names(name2group))
  if (length(miss))
    warning("Unmatched names left as NA: ", paste(miss, collapse = ", "))
  for (c in cols) set(dt, j = c, value = unname(name2group[dt[[c]]]))
  invisible(dt[])
}

# ── Usage ─────────────────────────────────────────────────────────────────────
# recode_cdl_functional(corn_data)                     # all crop_* name columns, in place
# recode_cdl_functional(corn_data, cols = paste0("crop_", 1:6))   # or a specific subset
#
# Spot-check (expect: Corn=corn, Soybeans=soybeans, Winter Wheat=annual_cereal,
#                     Dbl Crop WinWht/Soybeans=soybeans, Alfalfa=ley,
#                     Dbl Crop Soybeans/Oats=annual_cereal, Forest=NA):
# corn_data[, .N, by = crop_1][order(-N)]
