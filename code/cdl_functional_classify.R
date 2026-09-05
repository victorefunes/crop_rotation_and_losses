## ============================================================================
## cdl_functional_classify.R
## Single source of truth for mapping raw CDL crop NAMES to the project's
## established functional-group vocabulary:
##   corn | soybeans | ley | annual_cereal | annual_legume | annual_broadleaf
##   | fallow | NA (non-agricultural land / perennial specialty crops / unclear)
##
## This vocabulary and the base-name assignments below were reverse-engineered
## from an observed CDL-name -> functional_group frequency table (Illinois
## field-year sample, ~19M crop-years, cross-tabulated by a human collaborator
## against a pre-existing `crop_previous_2_year` -> `functional_group` column).
## Rows marked VERIFIED were directly observed in that table. Rows marked
## EXTRAPOLATED do not occur in Illinois CDL data (e.g. citrus, almonds) and
## were assigned by analogy to the verified rows -- double check these if your
## sample ever includes them.
##
## Double-crop names ("Dbl Crop X/Y") are NOT listed individually. All 10
## double-crop rows observed in the frequency table follow one rule: the
## functional group is whatever the SECOND named crop (Y) would get on its
## own -- e.g. "Dbl Crop Soybeans/Oats" -> annual_cereal (Oats' group), NOT
## soybeans; "Dbl Crop WinWht/Corn" -> corn (Corn's group), NOT annual_cereal.
## dbl_crop_group() below implements that rule generally, so any "Dbl Crop
## X/Y" name -- not just the 10 observed ones -- resolves automatically as
## long as Y is a known base name.
##
## Used by: cdl_functional_recode.R, rotation_setup_wa.R, and
## rotation_features_multicrop.R's build_class_lut() (for the "class" field
## only -- that file layers its own is_legume/is_small_grain/is_perennial
## attribute flags on top locally, since those are a different, orthogonal
## concept from this table).
## ============================================================================

library(data.table)

# ── Base classification (single, non-double-crop CDL names) ──────────────────

base_group <- c(
  # corn -- VERIFIED
  "Corn" = "corn", "Sweet Corn" = "corn", "Pop or Orn Corn" = "corn",

  # soybeans -- VERIFIED
  "Soybeans" = "soybeans",

  # ley (perennial forage / pasture / hay) -- VERIFIED
  "Grassland/Pasture" = "ley", "Alfalfa" = "ley",
  "Other Hay/Non Alfalfa" = "ley", "Sod/Grass Seed" = "ley",
  "Pasture/Grass" = "ley",                              # EXTRAPOLATED (retired dup of Grassland/Pasture)

  # annual_cereal (small grains, cereals other than corn) -- VERIFIED
  "Winter Wheat" = "annual_cereal", "Spring Wheat" = "annual_cereal",
  "Durum Wheat" = "annual_cereal", "Other Small Grains" = "annual_cereal",
  "Oats" = "annual_cereal", "Barley" = "annual_cereal", "Rye" = "annual_cereal",
  "Sorghum" = "annual_cereal", "Rice" = "annual_cereal", "Millet" = "annual_cereal",
  "Triticale" = "annual_cereal", "Speltz" = "annual_cereal",
  "Switchgrass" = "annual_cereal",

  # annual_legume (grain/forage legumes other than soybeans) -- VERIFIED
  "Dry Beans" = "annual_legume", "Peas" = "annual_legume",
  "Lentils" = "annual_legume", "Clover/Wildflowers" = "annual_legume",
  "Peanuts" = "annual_legume", "Vetch" = "annual_legume",
  "Chick Peas" = "annual_legume",                       # EXTRAPOLATED (like Dry Beans/Peas/Lentils)

  # fallow -- VERIFIED
  "Fallow/Idle Cropland" = "fallow",

  # annual_broadleaf (other annual row / vegetable / industrial crops) -- VERIFIED
  "Cotton" = "annual_broadleaf", "Sunflower" = "annual_broadleaf",
  "Sugarbeets" = "annual_broadleaf", "Potatoes" = "annual_broadleaf",
  "Sweet Potatoes" = "annual_broadleaf", "Canola" = "annual_broadleaf",
  "Flaxseed" = "annual_broadleaf", "Safflower" = "annual_broadleaf",
  "Rape Seed" = "annual_broadleaf", "Mustard" = "annual_broadleaf",
  "Camelina" = "annual_broadleaf", "Buckwheat" = "annual_broadleaf",
  "Misc Vegs & Fruits" = "annual_broadleaf", "Watermelons" = "annual_broadleaf",
  "Onions" = "annual_broadleaf", "Cucumbers" = "annual_broadleaf",
  "Tomatoes" = "annual_broadleaf", "Herbs" = "annual_broadleaf",
  "Tobacco" = "annual_broadleaf", "Carrots" = "annual_broadleaf",
  "Cantaloupes" = "annual_broadleaf", "Broccoli" = "annual_broadleaf",
  "Peppers" = "annual_broadleaf", "Greens" = "annual_broadleaf",
  "Squash" = "annual_broadleaf", "Pumpkins" = "annual_broadleaf",
  "Cabbage" = "annual_broadleaf", "Cauliflower" = "annual_broadleaf",
  "Celery" = "annual_broadleaf", "Radishes" = "annual_broadleaf",
  "Turnips" = "annual_broadleaf", "Eggplants" = "annual_broadleaf",
  "Gourds" = "annual_broadleaf", "Mint" = "annual_broadleaf",
  "Strawberries" = "annual_broadleaf",
  "Honeydew Melons" = "annual_broadleaf",               # EXTRAPOLATED (like Cantaloupes)
  "Garlic" = "annual_broadleaf",                        # EXTRAPOLATED (annual vegetable)
  "Lettuce" = "annual_broadleaf",                       # EXTRAPOLATED (annual vegetable)

  # NA: non-agricultural land cover -- VERIFIED
  "Background" = NA_character_, "Clouds/No Data" = NA_character_,
  "Developed" = NA_character_, "Water" = NA_character_, "Wetlands" = NA_character_,
  "Nonag/Undefined" = NA_character_, "Aquaculture" = NA_character_,
  "Open Water" = NA_character_, "Developed/Open Space" = NA_character_,
  "Developed/Low Intensity" = NA_character_, "Developed/Med Intensity" = NA_character_,
  "Developed/High Intensity" = NA_character_, "Barren" = NA_character_,
  "Deciduous Forest" = NA_character_, "Evergreen Forest" = NA_character_,
  "Mixed Forest" = NA_character_, "Shrubland" = NA_character_,
  "Forest" = NA_character_, "Woody Wetlands" = NA_character_,
  "Herbaceous Wetlands" = NA_character_,
  "Perennial Ice/Snow" = NA_character_,                 # EXTRAPOLATED (non-ag)

  # NA: perennial specialty / tree / vine crops, and the catch-all
  # "Other Crops" -- VERIFIED ("Other Crops" and all perennial-tree entries
  # observed in the frequency table classify to NA, not annual_broadleaf)
  "Other Crops" = NA_character_,
  "Grapes" = NA_character_, "Christmas Trees" = NA_character_,
  "Other Tree Crops" = NA_character_, "Walnuts" = NA_character_,
  "Cranberries" = NA_character_, "Apples" = NA_character_, "Peaches" = NA_character_,
  "Asparagus" = NA_character_,
  "Sugarcane" = NA_character_, "Caneberries" = NA_character_, "Hops" = NA_character_,   # EXTRAPOLATED
  "Cherries" = NA_character_, "Citrus" = NA_character_, "Pecans" = NA_character_,        # EXTRAPOLATED
  "Almonds" = NA_character_, "Pears" = NA_character_, "Pistachios" = NA_character_,      # EXTRAPOLATED
  "Prunes" = NA_character_, "Olives" = NA_character_, "Oranges" = NA_character_,          # EXTRAPOLATED
  "Avocados" = NA_character_, "Pomegranates" = NA_character_, "Nectarines" = NA_character_, # EXTRAPOLATED
  "Plums" = NA_character_, "Apricots" = NA_character_, "Blueberries" = NA_character_       # EXTRAPOLATED
)

# ── Functional-group labels appearing directly as raw values ─────────────────
# For lag years outside real CDL coverage (deep positions in the 6-year crop
# history), the source data substitutes the functional-group label itself as
# a fallback proxy instead of a specific crop name -- e.g. crop_5 might hold
# "annual_cereal" directly rather than a named cereal. Recognize these as
# pass-through so recode_cdl_functional() doesn't flag them as unmatched.
# "noncrop" is this fallback vocabulary's own no-crop label, equivalent to our
# NA classification.
functional_passthrough <- c(
  "corn" = "corn", "soybeans" = "soybeans", "ley" = "ley",
  "fallow" = "fallow", "annual_cereal" = "annual_cereal",
  "annual_legume" = "annual_legume", "annual_broadleaf" = "annual_broadleaf",
  "noncrop" = NA_character_
)
base_group <- c(base_group, functional_passthrough)

# ── Double-crop names ─────────────────────────────────────────────────────────
# "Dbl Crop X/Y" -> functional group of Y (see file header). A few
# double-crop tokens abbreviate their base name; normalize those first.
dbl_token_alias <- c(
  "WinWht" = "Winter Wheat",
  "Durum Wht" = "Durum Wheat",
  "Cantaloupe" = "Cantaloupes"
)

dbl_crop_group <- function(names) {
  m <- regmatches(names, regexpr("(?<=/).+$", names, perl = TRUE))
  m[m == ""] <- NA_character_
  m <- ifelse(m %in% names(dbl_token_alias), dbl_token_alias[m], m)
  unname(base_group[m])
}

# ── Public API ─────────────────────────────────────────────────────────────
#' @param raw character vector of raw CDL crop names
#' @return data.table(raw, class) -- class is one of corn/soybeans/ley/
#'   annual_cereal/annual_legume/annual_broadleaf/fallow/NA. `raw` values not
#'   present in base_group and not a resolvable "Dbl Crop X/Y" name come back
#'   with class = NA but are NOT distinguished from a deliberate NA
#'   classification -- callers that need to warn on truly-unmatched names
#'   (see cdl_functional_recode.R) check membership in base_group /
#'   dbl-crop-resolvability themselves, not just is.na(class).
classify_cdl_names <- function(raw) {
  raw <- unique(as.character(raw))
  raw <- raw[!is.na(raw)]

  is_dbl <- grepl("^Dbl [Cc]rop ", raw)
  class <- character(length(raw))
  class[!is_dbl] <- unname(base_group[raw[!is_dbl]])
  class[is_dbl]  <- dbl_crop_group(raw[is_dbl])

  data.table(raw = raw, class = class)
}
