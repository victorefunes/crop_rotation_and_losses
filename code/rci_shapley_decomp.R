# =============================================================================
# Shapley decomposition of ΔRCI into its components (number, t1, t2, corr)
#
# RCI = sqrt( number * turnover ),  turnover = (t1 + t2 + corr) / 2
# Because RCI is nonlinear in a product under a sqrt, ΔRCI has no clean additive
# split. Shapley attributes ΔRCI between an "A" (before) and "B" (after) state to
# each component, averaging its marginal contribution over all orderings. The
# pieces sum to ΔRCI EXACTLY (telescoping), with no ordering dependence and no
# residual — verified.
#
# Fully vectorized: 24 permutations x 4 steps = ~96 vectorized RCI evaluations
# over the whole panel, regardless of N. Scales to millions of rows.
# =============================================================================

library(data.table)

PERENNIAL <- c(36L, 37L, 61L, 176L, 224L)
NONAG     <- c(0L, 63L, 64L, 65L, 81L, 82L, 83L, 87L, 88L, 111L, 112L,
               121L, 122L, 123L, 124L, 131L, 141L, 142L, 143L, 152L, 190L, 195L)

# ---- 1. components from the six windowed CDL-code columns -------------------
#' @param M integer matrix n x 6 (crop_0 .. crop_5)
#' @return data.table with number, t1, t2, corr, RCI (NA where window has non-ag)
rci_components <- function(M) {
  storage.mode(M) <- "double"
  M[M %in% NONAG] <- NA_real_
  ok <- rowSums(is.na(M)) == 0L

  # distinct crop count (first-occurrence trick)
  firstocc <- matrix(TRUE, nrow(M), 6L)
  for (j in 2:6) {
    eq <- M[, j] == M[, 1L]
    if (j >= 3L) for (k in 2:(j - 1L)) eq <- eq | M[, j] == M[, k]
    firstocc[, j] <- !eq
  }
  number <- rowSums(firstocc)
  t1   <- rowSums(M[, 2:6] != M[, 1:5])                       # lag-1 changes
  t2   <- rowSums(M[, 3:6] != M[, 1:4])                       # lag-2 changes
  isper<- matrix(M %in% PERENNIAL, nrow(M))
  corr <- rowSums((M[, 2:6] == M[, 1:5]) & isper[, 2:6])      # perennial run pairs

  RCI  <- sqrt(number * (t1 + t2 + corr) / 2)
  out  <- data.table(number, t1, t2, corr, RCI)
  out[!ok, c("number","t1","t2","corr","RCI") := NA]
  out[]
}

rci_from_comp <- function(number, t1, t2, corr)
  sqrt(number * (t1 + t2 + corr) / 2)

# ---- 2. vectorized Shapley over the whole panel ----------------------------
# A, B: data.tables/matrices with columns number, t1, t2, corr (same nrow)
# returns n x 4 matrix of contributions; rowSums == RCI(B) - RCI(A) exactly
shapley_drci <- function(A, B) {
  comps <- c("number","t1","t2","corr")
  A <- as.matrix(A[, ..comps]); B <- as.matrix(B[, ..comps])
  n <- nrow(A)

  permn <- function(x) if (length(x) == 1L) list(x) else
    unlist(lapply(seq_along(x),
      function(i) lapply(permn(x[-i]), function(p) c(x[i], p))), recursive = FALSE)
  perms <- permn(comps)

  contrib <- matrix(0, n, 4L, dimnames = list(NULL, comps))
  rc <- function(S) rci_from_comp(S[,"number"], S[,"t1"], S[,"t2"], S[,"corr"])

  for (p in perms) {
    cur  <- A                                   # start all-at-A
    prev <- rc(cur)
    for (k in p) {
      cur[, k] <- B[, k]                         # switch component k to B
      now <- rc(cur)
      contrib[, k] <- contrib[, k] + (now - prev)
      prev <- now
    }
  }
  contrib / length(perms)                        # average over the 24 orderings
}

# ---- 3. driver: per-field year-over-year decomposition ---------------------
# Default "change" = this year's window vs the same field's previous year.
# Swap the A/B construction for field-vs-baseline or observed-vs-counterfactual.
decompose_rci <- function(dt,
                          id = "tile_field_ID", time = "year",
                          xcols = paste0("crop_", 0:5)) {
  dt <- as.data.table(dt)
  setorderv(dt, c(id, time))

  comp <- rci_components(as.matrix(dt[, ..xcols]))
  dt   <- cbind(dt[, c(id, time), with = FALSE], comp)

  # B = current row; A = previous row within field
  lagcols <- c("number","t1","t2","corr")
  dt[, paste0("A_", lagcols) := shift(.SD), by = id, .SDcols = lagcols]
  A <- dt[, setNames(.SD, lagcols), .SDcols = paste0("A_", lagcols)]
  B <- dt[, ..lagcols]

  keep <- stats::complete.cases(A) & stats::complete.cases(B)
  contrib <- matrix(NA_real_, nrow(dt), 4L,
                    dimnames = list(NULL, paste0("shap_", lagcols)))
  contrib[keep, ] <- shapley_drci(A[keep], B[keep])

  res <- cbind(
    dt[, c(id, time), with = FALSE],
    RCI_prev = rci_from_comp(A$number, A$t1, A$t2, A$corr),
    RCI_now  = dt$RCI,
    as.data.table(contrib)
  )
  res[, dRCI  := RCI_now - RCI_prev]
  res[, check := shap_number + shap_t1 + shap_t2 + shap_corr - dRCI]  # ~0 by construction
  res[]
}

# ---- usage ------------------------------------------------------------------
# out <- decompose_rci(corn_df)                       # per field-year contributions
# out[!is.na(dRCI), .(max_abs_check = max(abs(check)))]        # should be ~1e-12
#
# # sample-level attribution: of the average RCI change, how much is each margin?
# out[dRCI > 0, lapply(.SD, mean), .SDcols = patterns("^shap_")]
#
# # share of total ΔRCI variation carried by each component:
# out[!is.na(dRCI), lapply(.SD, function(s) sum(s)/sum(dRCI)),
#     .SDcols = patterns("^shap_")]
