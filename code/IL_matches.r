library(arrow)
library(tidyverse)
library(statar)
library(fixest)
library(broom)
setwd("C:/Users/vf006/Box/crop_rotations_and_losses/code")

# ── Load corn data ────────────────────────────────────────────────────────────

cat("Loading corn data...\n")
corn_df <- read_parquet(
  "D:/Crop data/d_igis13_12_1_2025.with_rci.parquet") 

corn_df <- corn_df |>
  filter(STATE_ABBR == "IL") |> 
  mutate(tile_field_ID = paste0("T", STATE_FIPS, "_", tile, "_", field_id),
         corn_yield = corn_yield / 62.77)  |> 
  arrange(tile_field_ID, year) 

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
         crop_6 = prior6yr_crop) |>
  filter(!(is.na(crop_5) & is.na(crop_4) & is.na(crop_3)))    


corn_df <- corn_df |>
  rename(RCI = annual_RCI) |>
  mutate(rot_crop = paste0(crop_5, "-", crop_4, "-", crop_3, "-", 
          crop_2, "-", crop_1, "-", crop_0))    

rot_tab <- corn_df |>
    group_by(rot_crop, RCI) |>
    summarise(n = n()) |>
    arrange(rot_crop, RCI) |>
    ungroup()

rm(corn_df)

rot_tab <- rot_tab |>
    mutate(rot_crop = factor(rot_crop))

rot_unique <- rot_tab |>
    distinct(rot_crop)

library(tidyr)

rot_unique <- rot_unique |>
  separate(rot_crop, into = paste0("crop", 1:6), sep = "-", remove = FALSE)

source("rci_vectorized.R")

rot_unique <- rot_unique |>
    mutate(RCI = rci(crop1, crop2, crop3, crop4, crop5, crop6)) |>
    select(-starts_with("crop")) 

rot_tab <- rot_tab |>
    left_join(rot_unique |> 
                rename(RCI_unique = RCI), 
              by = "rot_crop") |>
    arrange(rot_crop, RCI)    

rot_tab <- rot_tab |>
    mutate(RCI_match = !is.na(RCI) & !is.na(RCI_unique) & RCI == RCI_unique)

write_csv(rot_tab, file = "../tables/rot_frequencies.csv")

rot_sum <- rot_tab |>
    group_by(rot_crop) |>
    summarise(n_match = sum(n[RCI_match], na.rm = TRUE),
              n_total = sum(n, na.rm = TRUE)) |>
    mutate(pct_match = n_match / n_total * 100) |>
    ungroup()

rot_sum |> 
    tab(pct_match)

rot_sum |>
    filter(rot_crop %in% c("Corn-Corn-Corn-Corn-Corn-Corn", 
    "Soybeans-Corn-Soybeans-Corn-Soybeans-Corn",
    "Corn-Soybeans-Corn-Corn-Soybeans-Corn",
    "Soybeans-Corn-Corn-Soybeans-Corn-Corn"))  