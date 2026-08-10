################################################################################
## Tack, J. (2013), "A Nested Test for Common Yield Distributions with
## Application to U.S. Corn," J. Agric. Resour. Econ. 38(1):64-77.
##
## This is a GENERALIZATION of the Wu (2003) sequential-updating implementation.
## Wu fits a "rank-k exponential" density exp(-sum lambda_j x^j) by sequential
## Newton updates. Tack notes that arbitrary characterizing functions phi_j(x)
## can play the same role -- e.g., phi_1 = ln(x), phi_2 = ln(1-x) gives the
## beta density. He then builds a HYBRID model that stacks both sets, so beta
## and the rank-k alternative are both nested:
##
##   f^{Hyb}(x) = exp[-lambda_0 - lambda_1 ln(x) - lambda_2 ln(1-x)
##                                    - sum_{j=1}^k lambda_{j+2} x^j],
##
## with restriction lambda_3=...=lambda_{k+2}=0 giving Beta and lambda_1=
## lambda_2=0 giving the rank-k exponential. A standard LR test compares
## restricted vs unrestricted, R = 2(loglik_J - loglik_m) ~ chi^2(J-m).
##
## This file implements: (1) the generalized ME solver, (2) Beta / rank-k /
## hybrid wrappers, (3) sequential warm-starting (Wu 2003), (4) the Table-2
## diagnostics (logL, AIC, BIC, ID index, LR test), (5) a Harri-style
## detrending pipeline, and (6) GRP premium-rate calculation.
################################################################################

## ---------------------------------------------------------------------------
## 0. Quadrature primitives
## ---------------------------------------------------------------------------
make_grid <- function(lo, hi, n = 8000L) {
  x  <- seq(lo, hi, length.out = n)
  list(x = x, dx = x[2] - x[1])
}
trap <- function(fx, dx) dx * (sum(fx) - 0.5 * (fx[1] + fx[length(fx)]))

## ---------------------------------------------------------------------------
## 1. Generic maxent density on the grid
## ---------------------------------------------------------------------------
## Given a list of characterizing functions phi_j (each a function of x), build
## the matrix Phi[i,j] = phi_j(x_i) once -- it's the same across iterations.

build_Phi <- function(phi_funs, x) {
  M <- do.call(cbind, lapply(phi_funs, function(f) f(x)))
  ## Guard against -Inf from log near edges (we use an interior grid, but
  ## just in case)
  M[!is.finite(M)] <- 0
  M
}

## Density pieces given lambda, with a max-shift for numerical stability.
maxent_pieces_g <- function(lambda, Phi, dx) {
  expo <- as.numeric(Phi %*% lambda)               # sum_j lambda_j phi_j(x)
  m    <- max(-expo)
  w    <- exp(-expo - m)
  Z    <- trap(w, dx)
  p    <- w / Z
  ## w_stab = exp(-expo - m)  =>  Z_true = exp(m) * Z_stab
  ## hence the true log-normalizer is  log Z_true = m + log Z_stab.
  logZ <- log(Z) + m
  J    <- ncol(Phi)
  mu   <- numeric(J)
  for (j in seq_len(J)) mu[j] <- trap(Phi[, j] * p, dx)
  list(p = p, Z = Z, logZ = logZ, mu = mu)
}

## ---------------------------------------------------------------------------
## 2. Newton solver for arbitrary characterizing functions
## ---------------------------------------------------------------------------
## Gradient_j  = mu_j_target - E_p[phi_j]
## Hessian_ij  = Cov_p[phi_i, phi_j] = E_p[phi_i phi_j] - E_p[phi_i] E_p[phi_j]
## Newton step lambda <- lambda - H^{-1} grad,  damped to keep Z finite.

solve_maxent_g <- function(phi_funs, mu_target, grid, lambda0 = NULL,
                            tol = 1e-6, maxit = 400L, verbose = FALSE,
                            Phi = NULL) {
  if (is.null(Phi)) Phi <- build_Phi(phi_funs, grid$x)
  J  <- ncol(Phi)
  dx <- grid$dx
  lambda <- if (is.null(lambda0)) rep(0, J) else lambda0
  gnorm <- Inf
  for (it in seq_len(maxit)) {
    pc <- maxent_pieces_g(lambda, Phi, dx)
    grad <- mu_target - pc$mu

    ## Build Hessian (symmetric)
    H <- matrix(0, J, J)
    for (i in seq_len(J)) for (j in i:J) {
      H[i, j] <- trap(Phi[, i] * Phi[, j] * pc$p, dx) - pc$mu[i] * pc$mu[j]
      H[j, i] <- H[i, j]
    }

    ## Solve H step = grad with ridge / SVD-pseudoinverse fallback.
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      base_ridge <- 1e-12 * max(diag(H), 1)
      ok <- FALSE
      for (r in base_ridge * 10^(0:9)) {
        step <- tryCatch(solve(H + diag(r, J), grad), error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) { ok <- TRUE; break }
      }
      if (!ok) {
        sv   <- svd(H)
        dinv <- ifelse(sv$d > sv$d[1] * 1e-12, 1 / sv$d, 0)
        step <- as.numeric(sv$v %*% (dinv * (t(sv$u) %*% grad)))
      }
    }

    ## Damped line search: keep Z finite and reduce |grad|.
    a <- 1; gnorm0 <- sqrt(sum(grad^2))
    repeat {
      cand <- lambda - a * step
      pc2  <- maxent_pieces_g(cand, Phi, dx)
      if (is.finite(pc2$Z) && pc2$Z > 0) {
        gnew <- sqrt(sum((mu_target - pc2$mu)^2))
        if (gnew < gnorm0 || a < 1e-6) { lambda <- cand; break }
      }
      a <- a / 2
      if (a < 1e-10) { lambda <- cand; break }
    }
    gnorm <- sqrt(sum((mu_target - maxent_pieces_g(lambda, Phi, dx)$mu)^2))
    if (verbose) cat(sprintf("  it=%2d |grad|=%.3e step=%.3f\n", it, gnorm, a))
    if (gnorm < tol) break
  }
  pc <- maxent_pieces_g(lambda, Phi, dx)
  list(lambda = lambda, lambda0 = pc$logZ, mu = pc$mu, p = pc$p,
       converged = gnorm < tol, iters = it, gradnorm = gnorm,
       phi_funs = phi_funs, Phi = Phi)
}

## ---------------------------------------------------------------------------
## 3. Characterizing-function constructors (Tack 2013)
## ---------------------------------------------------------------------------
## Beta on (0,1): phi = (ln x, ln(1-x))                 [Tack p.67]
## Rank-k exponential: phi = (x, x^2, ..., x^k)         [Wu 2003]
## Hybrid: stack the two -- phi = (ln x, ln(1-x), x, x^2, ..., x^k)  [Eq. 11]

make_beta_funs <- function() list(
  function(x) log(x),
  function(x) log1p(-x)            # log(1 - x) with better precision near x=1
)

make_rank_k_funs <- function(k) {
  lapply(seq_len(k), function(i) { force(i); function(x) x^i })
}

make_hybrid_funs <- function(k) c(make_beta_funs(), make_rank_k_funs(k))

## Sample moments for an arbitrary phi list.
sample_moments_g <- function(x, phi_funs) {
  vapply(phi_funs, function(f) mean(f(x)), numeric(1))
}

## ---------------------------------------------------------------------------
## 4. Wu-style sequential warm-starting for the rank-k and hybrid families
## ---------------------------------------------------------------------------
## Tack uses k = 3, 5, 7 with the first THREE raw moments always included
## (yield skewness; see his footnote 7). We sequence the fits: each higher k
## is seeded with the previous lambda padded by zeros, exactly as in Wu.

fit_beta_me <- function(x, grid) {
  phi_funs <- make_beta_funs()
  mu  <- sample_moments_g(x, phi_funs)
  fit <- solve_maxent_g(phi_funs, mu, grid)
  fit$mu_target <- mu; fit$family <- "Beta"; fit
}

fit_rank_k <- function(x, k, grid, k_seq = NULL) {
  if (is.null(k_seq)) k_seq <- seq(3, k, by = 2)   # 3, 5, 7, ...
  prev_lambda <- NULL; fit <- NULL
  for (kk in k_seq) {
    phi_funs <- make_rank_k_funs(kk)
    mu       <- sample_moments_g(x, phi_funs)
    seed <- if (is.null(prev_lambda)) rep(0, kk)
            else c(prev_lambda, rep(0, kk - length(prev_lambda)))
    fit <- solve_maxent_g(phi_funs, mu, grid, lambda0 = seed)
    prev_lambda <- fit$lambda
    fit$mu_target <- mu
  }
  fit$family <- sprintf("Rank-%d exp", k); fit
}

fit_hybrid <- function(x, k, grid, beta_fit = NULL, rank_fit = NULL,
                        k_seq = NULL) {
  ## The hybrid has 2 + k Lagrange multipliers, which is hard to fit cold.
  ## We use a *two-source* warm start: the beta fit gives the two log-moment
  ## multipliers, and the rank-k fit gives the k power-moment multipliers.
  ## Stacking them is a much better seed than padding with zeros, because each
  ## marginal piece is already near its conditional optimum.
  ##
  ## Because rank-k is nested in the hybrid (set lambda_1 = lambda_2 = 0), the
  ## hybrid's log-likelihood must be >= rank-k's at the true optimum. We pick
  ## across seeds by ACHIEVED log-likelihood (the convex dual's actual value)
  ## rather than gradient norm, which is a safer convergence proxy when the
  ## 9-dim Hessian is ill-conditioned.
  if (is.null(beta_fit)) beta_fit <- fit_beta_me(x, grid)
  if (is.null(rank_fit)) rank_fit <- fit_rank_k(x, k, grid)

  phi_funs_final <- make_hybrid_funs(k)
  mu_final       <- sample_moments_g(x, phi_funs_final)
  N              <- length(x)
  needed         <- 2L + k

  ## Seeds tried (in order of expected quality):
  seeds <- list(
    stacked  = c(beta_fit$lambda, rank_fit$lambda),
    rank_pad = c(0, 0, rank_fit$lambda),               # nesting-feasible start
    beta_pad = c(beta_fit$lambda, rep(0, k)),
    zeros    = rep(0, needed)
  )

  ## Solve once from each seed, optionally polished, and keep the one with the
  ## highest log-likelihood -- which equals -N * (lambda_0 + sum lambda * mu).
  best <- NULL; best_ll <- -Inf
  for (s_name in names(seeds)) {
    fit_try <- try(solve_maxent_g(phi_funs_final, mu_final, grid,
                                  lambda0 = seeds[[s_name]],
                                  tol = 1e-9, maxit = 600L),
                   silent = TRUE)
    if (inherits(fit_try, "try-error")) next
    fit_try$mu_target <- mu_final
    ll <- loglik_me(fit_try, N)
    if (ll > best_ll) { best <- fit_try; best_ll <- ll }
  }

  ## Polish: warm-restart once more from the best lambda with a fresh Newton.
  pol <- try(solve_maxent_g(phi_funs_final, mu_final, grid,
                            lambda0 = best$lambda,
                            tol = 1e-10, maxit = 600L), silent = TRUE)
  if (!inherits(pol, "try-error")) {
    pol$mu_target <- mu_final
    if (loglik_me(pol, N) > best_ll) best <- pol
  }
  best$family <- sprintf("Hybrid (k=%d)", k)
  best
}

## ---------------------------------------------------------------------------
## 5. Likelihood, information criteria, ID index, LR test
## ---------------------------------------------------------------------------
## For exponential-family ME, MLE = MoM, so the log-likelihood at the fitted
## lambdas is
##     loglik = sum_t log f(x_t; lambda)
##            = sum_t [ -lambda_0 - sum_j lambda_j phi_j(x_t) ]
##            = -N [ lambda_0 + sum_j lambda_j * mu_j_hat ],
## where mu_j_hat is the sample moment and lambda_0 is the normalizer.

loglik_me <- function(fit, N) {
  -N * (fit$lambda0 + sum(fit$lambda * fit$mu_target))
}
aic_me <- function(fit, N) -2 * loglik_me(fit, N) + 2 * length(fit$lambda)
bic_me <- function(fit, N) -2 * loglik_me(fit, N) + log(N) * length(fit$lambda)

## ID index of Soofi-Ebrahimi-Habibullah (1995): 1 - exp(-KL(p1 : p2))
id_index <- function(p1, p2, dx) {
  keep <- p1 > 1e-300
  K <- trap(p1[keep] * log(p1[keep] / pmax(p2[keep], 1e-300)), dx)
  1 - exp(-K)
}

## LR test of restricted (m parameters) vs unrestricted (J parameters)
lr_test <- function(loglik_unr, loglik_res, df) {
  R <- 2 * (loglik_unr - loglik_res)
  list(R = R, df = df,
       crit_5 = qchisq(0.95, df),
       pval   = pchisq(R, df, lower.tail = FALSE))
}

## ---------------------------------------------------------------------------
## 6. Convenience: full Tack Table 2 for one sample
## ---------------------------------------------------------------------------
tack_table2 <- function(x, grid = NULL, k_max = 7L) {
  if (is.null(grid)) {
    eps <- max(1e-4, 0.5 * min(min(x), 1 - max(x)))
    grid <- make_grid(eps, 1 - eps, n = 12000L)
  }
  N <- length(x)
  fit_beta  <- fit_beta_me(x, grid)
  fit_r3    <- fit_rank_k(x, 3, grid, k_seq = 3)
  fit_r5    <- fit_rank_k(x, 5, grid, k_seq = c(3, 5))
  fit_r7    <- fit_rank_k(x, 7, grid, k_seq = c(3, 5, 7))
  fit_hyb   <- fit_hybrid(x, k_max, grid, beta_fit = fit_beta, rank_fit = fit_r7)
  fits <- list(Beta = fit_beta, rank3 = fit_r3, rank5 = fit_r5,
               rank7 = fit_r7, Hybrid = fit_hyb)

  logL <- sapply(fits, loglik_me, N = N)
  AIC  <- sapply(fits, aic_me,   N = N) / N           # per obs, as Tack reports
  BIC  <- sapply(fits, bic_me,   N = N) / N

  ## ID indices: rank3 vs Beta, rank5 vs rank3, rank7 vs rank5, then
  ##             Hybrid vs rank7 and Hybrid vs Beta.
  IDx <- c(
    Beta   = NA_real_,
    rank3  = id_index(fits$rank3$p,  fits$Beta$p,  grid$dx),
    rank5  = id_index(fits$rank5$p,  fits$rank3$p, grid$dx),
    rank7  = id_index(fits$rank7$p,  fits$rank5$p, grid$dx),
    Hybrid = id_index(fits$Hybrid$p, fits$rank7$p, grid$dx)
  )
  IDx_HybBeta <- id_index(fits$Hybrid$p, fits$Beta$p, grid$dx)

  ## LR tests (Tack column 5)
  ##  - rank5 vs rank3  : df = 2
  ##  - rank7 vs rank5  : df = 2
  ##  - Hybrid vs rank7 : df = 2  (adds ln x, ln(1-x))
  ##  - Hybrid vs Beta  : df = k_max = 7
  LR <- c(
    Beta   = NA_real_,
    rank3  = NA_real_,
    rank5  = lr_test(logL["rank5"],  logL["rank3"], 2)$R,
    rank7  = lr_test(logL["rank7"],  logL["rank5"], 2)$R,
    Hybrid = lr_test(logL["Hybrid"], logL["rank7"], 2)$R
  )
  LR_HybBeta <- lr_test(logL["Hybrid"], logL["Beta"], k_max)$R

  tab <- data.frame(
    Model  = names(fits),
    logL   = round(logL, 1),
    AIC    = round(AIC,  4),
    BIC    = round(BIC,  4),
    ID     = signif(IDx, 3),
    LR     = round(LR,   2),
    conv   = sapply(fits, function(f) f$converged),
    row.names = NULL
  )
  list(table = tab, fits = fits, grid = grid,
       extras = list(ID_HybVsBeta = IDx_HybBeta,
                     LR_HybVsBeta = LR_HybBeta,
                     crit_chi2_2  = qchisq(0.95, 2),
                     crit_chi2_k  = qchisq(0.95, k_max)))
}

## ---------------------------------------------------------------------------
## 7. Yield detrending pipeline (Harri et al. 2011, used by Tack)
## ---------------------------------------------------------------------------
## y_ist = alpha_i + beta_{1s} t + beta_{2s} t^2 + e_it          (state-specific
##                                                                quadratic
##                                                                trend, county FE)
## v_it    = e_it * (yhat_isT / yhat_ist)    (proportional-variance rescaling)
## ytilde  = yhat_isT + v_it
## y_norm  = ytilde / (1.5 * max(ytilde))    (push onto (0, 2/3) so we are
##                                            safely inside (0,1))

detrend_yields_proportional <- function(df) {
  ## df columns: county, state, year, yield (numeric).  Returns df augmented
  ## with `yield_norm` in (0,1), state-by-state.
  stopifnot(all(c("county","state","year","yield") %in% names(df)))
  out_list <- list()
  for (st in unique(df$state)) {
    sub <- df[df$state == st, , drop = FALSE]
    sub$county <- factor(sub$county)
    fit <- lm(yield ~ county + I(year) + I(year^2), data = sub)
    yhat <- fitted(fit); resid <- residuals(fit)
    T <- max(sub$year)
    ## Predict at year=T holding county fixed
    pred_T <- predict(fit, newdata = transform(sub, year = T))
    v       <- resid * (pred_T / yhat)
    ytilde  <- pred_T + v
    sub$yield_detrended <- ytilde
    sub$yield_norm      <- ytilde / (1.5 * max(ytilde))
    out_list[[st]] <- sub
  }
  do.call(rbind, out_list)
}

## ---------------------------------------------------------------------------
## 8. GRP premium-rate calculation (Tack Eqs 14-16)
## ---------------------------------------------------------------------------
grp_rate <- function(fit, grid, coverage) {
  p <- fit$p; x <- grid$x; dx <- grid$dx
  Ey <- trap(x * p, dx)
  vapply(coverage, function(cov) {
    y_trig <- cov * Ey
    ind    <- pmax(y_trig - x, 0) / cov         # disappearing deductible
    Eind   <- trap(ind * p, dx)
    Eind / y_trig
  }, numeric(1))
}

################################################################################
## 9. Sanity check: maxent Beta vs closed-form MLE on simulated Beta data
################################################################################
sanity_beta <- function(alpha = 17, beta = 20, N = 5000L, seed = 1L) {
  set.seed(seed)
  x   <- rbeta(N, alpha, beta)
  eps <- max(1e-4, 0.5 * min(min(x), 1 - max(x)))
  g   <- make_grid(eps, 1 - eps, n = 10000L)
  me  <- fit_beta_me(x, g)
  ## Recover (alpha, beta) from lambda: f = exp(-lambda_0 -lambda_1 ln x - lambda_2 ln(1-x))
  ##                                      = x^{-lambda_1} (1-x)^{-lambda_2} / B
  ##   => alpha_hat = 1 - lambda_1, beta_hat = 1 - lambda_2
  alpha_hat_me <- 1 - me$lambda[1]
  beta_hat_me  <- 1 - me$lambda[2]
  ## Closed-form MLE for beta via optim
  nll <- function(par) -sum(dbeta(x, par[1], par[2], log = TRUE))
  o   <- optim(c(2, 2), nll, method = "L-BFGS-B", lower = c(0.01, 0.01))
  alpha_hat_ml <- o$par[1]; beta_hat_ml <- o$par[2]
  cat(sprintf("Sanity check (true alpha=%.2f, beta=%.2f, N=%d):\n", alpha, beta, N))
  cat(sprintf("  Maxent beta:   alpha_hat=%.4f   beta_hat=%.4f\n",
              alpha_hat_me, beta_hat_me))
  cat(sprintf("  MLE beta:      alpha_hat=%.4f   beta_hat=%.4f\n",
              alpha_hat_ml, beta_hat_ml))
  invisible(list(me = me, alpha_me = alpha_hat_me, beta_me = beta_hat_me,
                 alpha_ml = alpha_hat_ml, beta_ml = beta_hat_ml))
}

################################################################################
## 10. Demonstration: nested test correctly identifies beta misspecification
################################################################################
demo_tack_table2 <- function() {

  cat("\n=== Sanity: ME-Beta vs MLE-Beta on simulated Beta(17,20) ===\n")
  sanity_beta()

  ## ---- Case A: data really ARE Beta(17, 20). LR should NOT reject. ----
  cat("\n=== Case A: true DGP is Beta(17, 20), N=5000 ===\n")
  set.seed(7); xA <- rbeta(5000, 17, 20)
  resA <- tack_table2(xA, k_max = 7L)
  cat("Per-model diagnostics:\n"); print(resA$table, row.names = FALSE)
  cat(sprintf("\nID(Hybrid : Beta)   = %.4g\n", resA$extras$ID_HybVsBeta))
  cat(sprintf("LR(Hybrid vs Beta)  = %.2f   (chi^2_7 5%% crit = %.2f)\n",
              resA$extras$LR_HybVsBeta, resA$extras$crit_chi2_k))
  cat(sprintf("LR(Hybrid vs rank7) = %.2f   (chi^2_2 5%% crit = %.2f)\n",
              resA$table$LR[resA$table$Model == "Hybrid"],
              resA$extras$crit_chi2_2))
  cat("Expected: Beta is the true DGP, so neither LR test should reject.\n")

  ## ---- Case B: mixture - bulk like beta, but extra mass in the left tail
  ##                       that beta cannot match. LR SHOULD reject beta.   ----
  cat("\n=== Case B: 0.9*Beta(20,18) + 0.1*Beta(2,12), N=5000 ===\n")
  set.seed(7)
  N  <- 5000; w <- rbinom(N, 1, 0.9)
  xB <- ifelse(w == 1, rbeta(N, 20, 18), rbeta(N, 2, 12))
  resB <- tack_table2(xB, k_max = 7L)
  cat("Per-model diagnostics:\n"); print(resB$table, row.names = FALSE)
  cat(sprintf("\nID(Hybrid : Beta)   = %.4g\n", resB$extras$ID_HybVsBeta))
  cat(sprintf("LR(Hybrid vs Beta)  = %.2f   (chi^2_7 5%% crit = %.2f)\n",
              resB$extras$LR_HybVsBeta, resB$extras$crit_chi2_k))
  cat(sprintf("LR(Hybrid vs rank7) = %.2f   (chi^2_2 5%% crit = %.2f)\n",
              resB$table$LR[resB$table$Model == "Hybrid"],
              resB$extras$crit_chi2_2))
  cat("Expected: beta is misspecified, so LR(Hybrid vs Beta) >> 14.07.\n")

  ## ---- Plot Case B (mimicking Tack's Figure 2) ----
  png("/home/claude/tack_demo_case_B.png", width = 950, height = 620, res = 110)
  par(mar = c(4.1, 4.3, 3, 1))
  hist(xB, breaks = 60, freq = FALSE, col = "grey88", border = "white",
       xlim = c(0, 1), ylim = c(0, 6),
       main = "Demonstration: 0.9*Beta(20,18) + 0.1*Beta(2,12), N = 5000",
       xlab = "x", ylab = "Density")
  ## kernel density for the underlying data (Tack's dashed kernel line)
  kd <- density(xB, bw = "SJ", from = 0, to = 1)
  lines(kd, lwd = 1.5, lty = 2, col = "grey30")
  lines(resB$grid$x, resB$fits$Beta$p,  lwd = 2, col = "steelblue")
  lines(resB$grid$x, resB$fits$rank7$p, lwd = 2, col = "firebrick")
  legend("topright",
         legend = c("Kernel density", "Beta fit", "Rank-7 exp fit"),
         lwd = c(1.5, 2, 2), lty = c(2, 1, 1),
         col = c("grey30", "steelblue", "firebrick"), bty = "n")
  dev.off()

  ## ---- GRP rate ratios (Beta over Rank-7), Tack's Figure 3 panel ----
  cov_grid <- seq(0.5, 0.9, by = 0.01)
  rate_B   <- grp_rate(resB$fits$Beta,  resB$grid, cov_grid)
  rate_R7  <- grp_rate(resB$fits$rank7, resB$grid, cov_grid)
  png("/home/claude/tack_demo_grp_ratio.png", width = 750, height = 520, res = 110)
  par(mar = c(4.1, 4.3, 3, 1))
  plot(cov_grid, rate_B / rate_R7, type = "o", pch = 1, lwd = 2,
       xlab = "Coverage level", ylab = "Rate ratio (Beta / Rank-7)",
       ylim = c(0, 1.2),
       main = "GRP premium-rate ratio under beta vs rank-7 (Case B)")
  abline(h = 1, lty = 3, col = "grey50")
  dev.off()

  invisible(list(case_A = resA, case_B = resB))
}

## ============================== Run it ==============================
out <- demo_tack_table2()
cat("\nSaved demo plots to tack_demo_case_B.png and tack_demo_grp_ratio.png\n")
