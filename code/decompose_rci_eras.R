# =============================================================================
# decompose_rci_eras(): era-to-era Shapley decomposition of ΔRCI.
# Takes the LONG panel (id, year, integer CDL code) and does everything:
#   window mapping -> reshape to 6 columns per (id, window) -> completeness
#   filter -> keep fields complete in BOTH chosen eras -> component Shapley.
#
# Requires rci_components() and shapley_drci() from rci_shapley_decomp.R.
#   source("rci_shapley_decomp.R")
# =============================================================================
library(data.table)
if (!exists("rci_components") || !exists("shapley_drci"))
  stop("source('rci_shapley_decomp.R') first — needs rci_components() and shapley_drci().")

#' @param dt long panel, one row per field-year
#' @param id,year,crop column names; `crop` MUST be integer CDL codes
#' @param origin calendar year that starts window 0 (default = min year).
#'   Set explicitly (e.g. 2008L) to force specific eras like 2008-13 vs 2014-19.
#' @param anchor_A,anchor_B window_ids to compare (default = first & last complete)
#' @param merge_perennial TRUE -> 3-way split {number, t1c=t1+corr, t2}
#' @return one row per usable field: RCI_A, RCI_B, dRCI, shap_* , check(~0)
decompose_rci_eras <- function(dt,
                               id = "tile_field_ID", year = "year", crop = "crop",
                               origin = NULL, anchor_A = NULL, anchor_B = NULL,
                               merge_perennial = FALSE) {
  dt <- copy(as.data.table(dt))
  stopifnot(anyDuplicated(dt, by = c(id, year)) == 0L,
            is.numeric(dt[[crop]]))                       # codes, not name strings
  if (is.null(origin)) origin <- min(dt[[year]])

  dt[, `:=`(.wid = (get(year) - origin) %/% 6L,
            .pos = (get(year) - origin) %%  6L)]

  # reshape to one row per (id, window) with crop_0..crop_5
  f    <- as.formula(sprintf("%s + .wid ~ .pos", id))
  wide <- dcast(dt, f, value.var = crop)
  pos  <- as.character(0:5)
  miss <- setdiff(pos, names(wide)); if (length(miss)) wide[, (miss) := NA_integer_]
  setnames(wide, pos, paste0("crop_", 0:5))
  xcols <- paste0("crop_", 0:5)

  # keep only complete windows (all six years present)
  wide <- wide[rowSums(!is.na(as.matrix(wide[, ..xcols]))) == 6L]

  comp <- rci_components(as.matrix(wide[, ..xcols]))       # number,t1,t2,corr,RCI
  wide <- cbind(wide, comp)

  wins <- sort(unique(wide$.wid))
  if (is.null(anchor_A)) anchor_A <- wins[1L]
  if (is.null(anchor_B)) anchor_B <- wins[length(wins)]

  if (merge_perennial) { wide[, t1c := t1 + corr]; comps <- c("number","t1c","t2") }
  else                   comps <- c("number","t1","t2","corr")

  # fields complete in BOTH eras, aligned by id
  A <- wide[.wid == anchor_A]; B <- wide[.wid == anchor_B]
  keep <- intersect(A[[id]], B[[id]])
  A <- A[get(id) %in% keep]; B <- B[get(id) %in% keep]
  setorderv(A, id); setorderv(B, id)
  stopifnot(identical(A[[id]], B[[id]]))

  Ac <- A[, ..comps]; Bc <- B[, ..comps]
  contrib <- shapley_drci(A = Ac, comps = comps, B = Bc)

  rci_of <- function(S) sqrt(S$number *
                             rowSums(as.matrix(S[, setdiff(comps, "number"), with = FALSE])) / 2)
  res <- data.table(A[[id]], anchor_A, anchor_B, rci_of(Ac), rci_of(Bc))
  setnames(res, c(id, "window_A", "window_B", "RCI_A", "RCI_B"))
  res <- cbind(res, as.data.table(contrib))
  setnames(res, comps, paste0("shap_", comps))
  res[, dRCI  := RCI_B - RCI_A]
  res[, check := rowSums(.SD) - dRCI, .SDcols = patterns("^shap_")]   # ~1e-12
  res[]
}

# ---- usage ------------------------------------------------------------------
# source("rci_shapley_decomp.R")
# eras3 <- decompose_rci_eras(corn_df, id = "id", crop = "crop",
#                             origin = 2008L, merge_perennial = TRUE)   # 2008-13 vs 2014-19
# eras3[, .(max_abs_check = max(abs(check)))]                          # ~1e-12
# eras3[, lapply(.SD, function(s) sum(s)/sum(dRCI)), .SDcols = patterns("^shap_")]
