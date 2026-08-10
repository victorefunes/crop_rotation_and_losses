################################################################################
## Sequential updating method for maximum entropy densities
## Wu, X. (2003), "Calculation of maximum entropy densities with application
## to income distribution," Journal of Econometrics 115, 347-354.
##
## The maxent density subject to k arithmetic moment constraints is
##     p(x) = exp( -sum_{i=0}^{k} lambda_i * x^i )
## where lambda_0 is the normalizing (log-partition) constant.
##
## We solve for lambda = (lambda_1,...,lambda_k) by Newton's method on the
## dual objective, and use SEQUENTIAL UPDATING: the converged lambda for k
## moments, padded with a zero, seeds the Newton iteration for k+1 (here k+2).
################################################################################

## ----------------------------------------------------------------------------
## 0. Numerical integration grid
## ----------------------------------------------------------------------------
## All moments mu_i(lambda) = E[x^i] and the partition function Z are computed
## by quadrature on a fixed grid over the (rescaled) support. Wu rescales income
## by $100,000 so the support is roughly [0, ~22]; we standardize instead to
## keep x^k from overflowing for large k. Everything is done on a grid.

make_grid <- function(lo, hi, n = 4000L) {
  x  <- seq(lo, hi, length.out = n)
  dx <- x[2] - x[1]
  list(x = x, dx = dx)
}

## Trapezoidal integral of f over the grid
trap <- function(fx, dx) {
  dx * (sum(fx) - 0.5 * (fx[1] + fx[length(fx)]))
}

## ----------------------------------------------------------------------------
## 1. Core pieces of the dual problem
## ----------------------------------------------------------------------------
## Given lambda = (lambda_1,...,lambda_k) (NOT including lambda_0), compute on
## the grid: the unnormalized density w(x)=exp(-sum lambda_i x^i), the partition
## function Z, the normalized density p, and the predicted moments mu_i.

maxent_pieces <- function(lambda, grid, k) {
  x <- grid$x; dx <- grid$dx
  ## Powers matrix X[, i] = x^i for i = 1..k
  ## exponent = sum_i lambda_i x^i
  expo <- rep(0, length(x))
  for (i in seq_len(k)) expo <- expo + lambda[i] * x^i
  ## subtract max for numerical stability (cancels in normalization)
  m <- max(-expo)
  w <- exp(-expo - m)             # stabilized unnormalized weights
  Z <- trap(w, dx)                # stabilized partition function
  p <- w / Z                      # normalized density (independent of shift m)
  ## True Z = exp(m) * Z_stab, so log Z_true = m + log Z_stab.
  logZ <- log(Z) + m
  ## predicted moments mu_i = E[x^i], i = 1..k
  mu <- numeric(k)
  for (i in seq_len(k)) mu[i] <- trap((x^i) * p, dx)
  list(p = p, w = w, Z = Z, logZ = logZ, mu = mu, x = x, dx = dx)
}

## ----------------------------------------------------------------------------
## 2. Newton solver for a FIXED number of moments k
## ----------------------------------------------------------------------------
## Dual objective Gamma(lambda) = log Z + sum_i lambda_i * nu_i, whose
## stationarity conditions are mu_i(lambda) = nu_i (predicted = sample moments).
##   gradient_i = nu_i - mu_i(lambda)
##   Hessian_ij = mu_{i+j}(lambda) - mu_i(lambda) mu_j(lambda)     [Eq. (4)]
## The Hessian is positive definite, so the solution is unique.

solve_maxent_k <- function(nu, grid, k, lambda0 = NULL,
                           tol = 1e-7, maxit = 300L, verbose = FALSE) {
  if (is.null(lambda0)) lambda0 <- rep(0, k)
  lambda <- lambda0
  x <- grid$x; dx <- grid$dx

  for (it in seq_len(maxit)) {
    pc <- maxent_pieces(lambda, grid, k)
    ## gradient: nu_i - mu_i  (we minimize Gamma => step toward grad = 0)
    grad <- nu - pc$mu
    ## Build Hessian via moments up to order 2k
    mu_all <- numeric(2 * k)            # mu_1 .. mu_{2k}
    for (i in seq_len(2 * k)) mu_all[i] <- trap((x^i) * pc$p, dx)
    H <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(k)) {
      H[i, j] <- mu_all[i + j] - pc$mu[i] * pc$mu[j]
    }
    ## Newton step: lambda_new = lambda - H^{-1} grad   (Eq. 3)
    ## Solve H step = grad, with an adaptive ridge + pseudoinverse fallback,
    ## since the Hessian approaches singularity as k grows (Wu, p.350).
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      ridge <- 1e-12 * max(diag(H))
      ok <- FALSE
      for (r in ridge * 10^(0:8)) {
        step <- tryCatch(solve(H + diag(r, k), grad), error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) { ok <- TRUE; break }
      }
      if (!ok) {                                   # SVD pseudoinverse
        sv <- svd(H); dinv <- ifelse(sv$d > sv$d[1] * 1e-12, 1 / sv$d, 0)
        step <- sv$v %*% (dinv * (t(sv$u) %*% grad))
        step <- as.numeric(step)
      }
    }
    ## damped line search to keep Z finite and improve the gradient norm
    a <- 1
    gnorm0 <- sqrt(sum(grad^2))
    repeat {
      cand <- lambda - a * step
      pc2  <- maxent_pieces(cand, grid, k)
      if (is.finite(pc2$Z) && pc2$Z > 0) {
        gnew <- sqrt(sum((nu - pc2$mu)^2))
        if (gnew < gnorm0 || a < 1e-6) { lambda <- cand; break }
      }
      a <- a / 2
      if (a < 1e-10) { lambda <- cand; break }
    }
    gnorm <- sqrt(sum((nu - maxent_pieces(lambda, grid, k)$mu)^2))
    if (verbose) cat(sprintf("  k=%d it=%2d  |grad|=%.3e  step=%.2f\n",
                             k, it, gnorm, a))
    if (gnorm < tol) break
  }
  pc <- maxent_pieces(lambda, grid, k)
  lambda0_norm <- pc$logZ           # lambda_0 = log Z (normalizer)
  list(lambda = lambda, lambda0 = lambda0_norm, mu = pc$mu,
       p = pc$p, x = grid$x, converged = (gnorm < tol), iters = it,
       gradnorm = gnorm)
}

## ----------------------------------------------------------------------------
## 3. SEQUENTIAL UPDATING across k = k_seq (Wu's main contribution)
## ----------------------------------------------------------------------------
## Solve for the lowest k with zero initial values, then for each higher k use
## the previous solution padded with zeros for the new high-order coefficients.
## If a direct Newton solve from the seed fails, split the moment difference
## into segments and approach in steps (homotopy), exactly as Wu describes.

maxent_sequential <- function(x_std, grid, k_seq = c(4, 6, 8, 10, 12),
                              verbose = TRUE) {
  ## sample (raw) moments on the standardized scale: nu_i = mean(x^i)
  kmax <- max(k_seq)
  nu_all <- vapply(seq_len(kmax), function(i) mean(x_std^i), numeric(1))

  fits <- list()
  prev_lambda <- NULL
  for (idx in seq_along(k_seq)) {
    k  <- k_seq[idx]
    nu <- nu_all[seq_len(k)]
    if (is.null(prev_lambda)) {
      seed <- rep(0, k)                          # zeros for the first model
    } else {
      seed <- c(prev_lambda, rep(0, k - length(prev_lambda)))  # pad with 0s
    }
    if (verbose) cat(sprintf("Fitting k = %2d ... ", k))

    fit <- solve_maxent_k(nu, grid, k, lambda0 = seed, verbose = FALSE)

    ## Homotopy fallback: if not converged, walk from the seed's predicted
    ## moments to the true sample moments in M steps.
    if (!fit$converged) {
      if (verbose) cat("  direct solve failed; using staged homotopy\n")
      base_mu <- maxent_pieces(seed, grid, k)$mu
      M <- 10L
      lam <- seed
      for (m in seq_len(M)) {
        nu_m <- base_mu + (nu - base_mu) * (m / M)
        fit  <- solve_maxent_k(nu_m, grid, k, lambda0 = lam, verbose = FALSE)
        lam  <- fit$lambda
      }
      fit <- solve_maxent_k(nu, grid, k, lambda0 = lam, verbose = FALSE)
    }
    if (verbose) cat(sprintf("converged=%s  iters=%d  |grad|=%.2e\n",
                             fit$converged, fit$iters, fit$gradnorm))

    fits[[as.character(k)]] <- fit
    prev_lambda <- fit$lambda
  }
  list(fits = fits, nu_all = nu_all, k_seq = k_seq)
}

## ----------------------------------------------------------------------------
## 4. Goodness-of-fit / specification statistics (Table 1 of Wu 2003)
## ----------------------------------------------------------------------------
## For the exponential family, MoM = MLE, so the log-likelihood is
##   L = -N * sum_{i=0}^k lambda_i nu_i      (nu_0 = 1)
## LR test compares p_{k+2} vs p_k ~ chi^2(2).  ID index uses KL distance.
## KS compares the fitted CDF to the empirical CDF on the standardized data.

loglik_maxent <- function(fit, nu_all, N) {
  k <- length(fit$lambda)
  ## sum_{i=0}^k lambda_i * nu_i, with nu_0 = 1 and lambda_0 the normalizer
  s <- fit$lambda0 * 1 + sum(fit$lambda * nu_all[seq_len(k)])
  -N * s
}

maxent_cdf <- function(fit, grid) {
  cs <- cumsum(fit$p) * grid$dx
  cs / cs[length(cs)]
}

gof_table <- function(seq_out, x_std, grid) {
  N <- length(x_std)
  k_seq <- seq_out$k_seq
  nu_all <- seq_out$nu_all
  fits <- seq_out$fits

  L  <- sapply(k_seq, function(k) loglik_maxent(fits[[as.character(k)]], nu_all, N))
  names(L) <- k_seq

  ## LR: p_{k+2} vs p_k  (matched to the k that has a k+2 partner)
  LR <- rep(NA_real_, length(k_seq))
  for (i in seq_along(k_seq)) {
    kp <- k_seq[i] + 2
    if (as.character(kp) %in% names(fits)) {
      LR[i] <- 2 * (L[as.character(kp)] - L[as.character(k_seq[i])])
    }
  }

  ## KS statistic against empirical CDF
  Fx <- ecdf(x_std)
  KS <- sapply(k_seq, function(k) {
    fit <- fits[[as.character(k)]]
    cdf <- maxent_cdf(fit, grid)
    Femp <- Fx(grid$x)
    max(abs(cdf - Femp))
  })

  ## AIC / BIC per observation (penalty uses k params: lambda_1..lambda_k)
  AIC <- (-2 * L + 2 * k_seq) / N
  BIC <- (-2 * L + log(N) * k_seq) / N

  data.frame(k = k_seq, logLik = round(L, 1),
             LR_vs_kplus2 = round(LR, 1),
             KS = round(KS, 4),
             AIC = round(AIC, 4), BIC = round(BIC, 4),
             converged = sapply(k_seq, function(k) fits[[as.character(k)]]$converged))
}

################################################################################
## 5. APPLICATION TO CPS FAMILY INCOME
################################################################################
run_income_example <- function(csv_path,
                                income_col = "V7",
                                n_draw = 5000L,
                                k_seq = c(4, 6, 8, 10, 12),
                                seed = 1L) {
  d <- read.csv(csv_path, header = FALSE)
  inc <- d[[income_col]]
  inc <- inc[is.finite(inc) & inc > 0]           # positive family income only

  ## Draw 5,000 observations at random, as in Wu (2003)
  set.seed(seed)
  if (length(inc) > n_draw) inc <- sample(inc, n_draw)

  ## Wu rescales income by $100,000 and fits on roughly [0, 22]. With k up to
  ## 12 and a heavy right tail, the raw moments E[x^i] span >20 orders of
  ## magnitude (E[x^12] ~ 1e12 here), which makes the Hessian numerically
  ## singular. The standard conditioning fix is to map the support onto the
  ## UNIT INTERVAL u = x / x_max in [0,1], so that every moment E[u^i] lies in
  ## (0,1) and the powers u^i never overflow. The fitted maxent density g(u)
  ## is mapped back to the income scale by the change of variables
  ##     p_inc(x) = g(x / x_max) / x_max ,   and to $ by another / scale_unit.
  scale_unit <- 1e5
  x_inc <- inc / scale_unit                      # income in units of $100k
  xmax  <- max(x_inc) * 1.02                      # upper support, padded
  u     <- x_inc / xmax                           # fitting variable in [0,1)

  lo <- 0
  hi <- 1
  grid <- make_grid(lo, hi, n = 8000L)

  cat(sprintf("N = %d   $-scale = %.0f   x_max = %.3f ($100k)   u in [0,1]\n\n",
              length(u), scale_unit, xmax))

  seq_out <- maxent_sequential(u, grid, k_seq = k_seq, verbose = TRUE)

  cat("\n================ Specification / goodness-of-fit ================\n")
  tab <- gof_table(seq_out, u, grid)
  print(tab, row.names = FALSE)
  ks_crit <- 1.358 / sqrt(length(u))
  cat(sprintf("\nKS 5%% critical value = %.4f\n", ks_crit))

  invisible(list(seq_out = seq_out, grid = grid, u = u,
                 xmax = xmax, scale_unit = scale_unit))
}

## ---- Run it ----------------------------------------------------------------
res <- run_income_example("/mnt/user-data/uploads/cps_2014.csv",
                          k_seq = c(4, 6, 8, 10, 12))

## ---- Plot: histogram + fitted maxent density (k = 12), on income scale ------
png("/home/claude/maxent_income_fit.png", width = 900, height = 600, res = 110)
grid  <- res$grid
xmax <- res$xmax
fit12 <- res$seq_out$fits[["12"]]
## map u-density back to $100k income scale: x = xmax*u, p(x) = g(u)/xmax
x_axis <- xmax * grid$x
p_x    <- fit12$p / xmax
inc_100k <- res$u * xmax
hist(inc_100k, breaks = 80, freq = FALSE, col = "grey85", border = "white",
     xlim = c(0, 5), main = "CPS family income (Wu 2003 sequential maxent)",
     xlab = "Family income ($100,000)", ylab = "Density")
lines(x_axis, p_x, lwd = 2, col = "firebrick")
legend("topright", legend = "Maxent density (k = 12)",
       lwd = 2, col = "firebrick", bty = "n")
dev.off()
cat("\nSaved plot to maxent_income_fit.png\n")
