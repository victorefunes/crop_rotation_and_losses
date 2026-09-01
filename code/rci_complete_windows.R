# =============================================================================
# Identify fields usable for the era-to-era RCI decomposition.
# Input: long panel with columns id, year, crop, window_id, window_pos
#        (window_id / window_pos built as in the reshape step:
#           window_id  = (year - origin) %/% 6L
#           window_pos = (year - origin) %%  6L )
# =============================================================================
library(data.table)
setDT(corn_df)

stopifnot(anyDuplicated(corn_df, by = c("tile_field_ID", "year")) == 0L)   # dcast-safe check

# ---- 1. per (tile_field_ID, window) completeness summary ------------------------------
win <- corn_df[, .(
  n_yrs     = uniqueN(window_pos),        # distinct slots filled (max 6)
  has_start = 0L %in% window_pos,         # beginning of window observed
  has_end   = 5L %in% window_pos          # end of window observed
), by = .(tile_field_ID, window_id)]

win[, complete_full      := n_yrs == 6L]              # all six years -> RCI computable
win[, complete_endpoints := has_start & has_end]      # only beginning AND end observed

wins    <- sort(unique(win$window_id))
w_first <- wins[1L]; w_last <- wins[length(wins)]
n_win   <- length(wins)

# ---- 2. the set the DECOMPOSITION consumes ---------------------------------
# fields with a COMPLETE window at both the first and last era -> ΔRCI computable
ids_first <- win[complete_full & window_id == w_first, tile_field_ID]
ids_last  <- win[complete_full & window_id == w_last,  tile_field_ID]
usable_change <- intersect(ids_first, ids_last)       # <-- use THIS to filter

# ---- 3. stricter / looser alternatives -------------------------------------
# complete in EVERY window (fully balanced panel across all eras)
balanced_all <- win[complete_full == TRUE, .N, by = tile_field_ID][N == n_win, tile_field_ID]

# literal "observed at the beginning and end of each window" (endpoints only) --
# NB: weaker than RCI-computable; a field can pass this yet miss an interior
# year, so DON'T feed this to the decomposition. Kept for reference / coverage.
endpoints_every <- win[complete_endpoints == TRUE, .N, by = tile_field_ID][N == n_win, tile_field_ID]

# ---- 4. report + apply ------------------------------------------------------
cat(sprintf(
  "windows: %s | fields: %d total | change-usable (complete %d & %d): %d | balanced-all: %d\n",
  paste(range(wins), collapse = "-"), uniqueN(win$tile_field_ID),
  w_first, w_last, length(usable_change), length(balanced_all)))

# keep only field-windows for change-usable fields, both eras, RCI-ready:
panel_change <- corn_df[tile_field_ID %in% usable_change &
                          window_id %in% c(w_first, w_last)]
