################################################################################
# Bandwidth selection
################################################################################
if (!requireNamespace("DirStatsOld", quietly = TRUE)) {
  install.packages("DirStatsOld_0.1.5.tar.gz", repos = NULL, type = "source")
}
library(DirStatsOld)
source("init.R")

#===============================================================================
# CROSS-VALIDATION APPROACH (for independent data)
#===============================================================================

##------Cross-Validation----------------------------------------
# Computes leave-one-out cross-validation score for bandwidth h.
# CV(h) = sum_i [ Y_i - m_hat_{h,p,-i}(X_i) ]^2 
cv_loo <- function(X, Y, h, p) {
  n <- nrow(X)
  res <- numeric(n)
  for (i in seq_len(n)) {
    idx_train <- setdiff(seq_len(n), i)
    
    # Fit model on data excluding point i and evaluate at point i
    yhat_i <- loc.directional.linear(x = X[i, , drop = FALSE],
                                     data.dir = X[idx_train, , drop = FALSE],
                                     data.lin = Y[idx_train],
                                     h = h, p = p)$Yhat
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res) 
}

# Computes simplified cross-validation score for bandwidth h.
# CV(h) = (1/n) * sum_i [ (Y_i - m_hat_{h,p}(X_i)) / (1 - S_{ii}) ]^2 
cv <- function(X, Y, h, p) {
  # Fit full model to extract fitted values and smoother matrix S
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                   h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  Yhat <- as.numeric(fit$Yhat)
  
  denom <- 1 - diag(S)
  if (any(denom <= 1e-5)) return(Inf)
  
  mean(((Y - Yhat) / denom)^2, na.rm = TRUE)
}

##------Generalized Cross-Validation ----------------------------------------
# Computes generalized cross-validation score for bandwidth h.
# --- GCV(h) = (1/n) * sum_i [ (Y_i - m_hat_{h,p}(X_i)) / (1 - tr(S)/n) ]^2 ---

gcv <- function(X, Y, h, p) {
  n <- nrow(X)
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  Yhat <- as.numeric(fit$Yhat)
  
  denom <- 1 - sum(diag(S)) / n
  if (!is.finite(denom) || denom <= 0) return(Inf)
  
  mean(((Y - Yhat) / denom)^2, na.rm = TRUE)
}

#===============================================================================
# MODIFIED CROSS-VALIDATION APPROACH (for dependent data)
#===============================================================================

##------ Modified  cross-validation -----------------------------------------
# MCV(h) = sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
# N(i) = { j : theta(X_j, X_i) <= ell }
# For S^2, N(i) is a spherical cap of geodesic radius ell around X_i.

mcv_loo <- function(X, Y, h, p, ell, D = NULL) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi  # (n x n) geodesic distances
  res <- numeric(n)
  
  for (i in seq_len(n)) {
    in_nbhd <- which(D[i, ] <= ell)
    idx_train <- setdiff(seq_len(n), in_nbhd)
    
    # Skip if too few points remain after removing N(i) to fit a 
    # degree-p polynomial
    if (length(idx_train) < (p + 2)) {
      res[i] <- NA
      next
    }
    yhat_i <- loc.directional.linear(x = X[i, , drop = FALSE],
                                     data.dir = X[idx_train, , drop = FALSE],
                                     data.lin = Y[idx_train],
                                     h = h, p = p)$Yhat
    
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res, na.rm = TRUE)
}

##------ Modified generalized cross-validation --------------------------------
#  MGCV(h) = (1/n) sum_i [ (Y_i - m_hat_{h,p}(X_i)) / (1 - tr(S_p R)/n) ]^2

mgcv <- function(X, Y, h, p, R = NULL, D = NULL, h_pilot = NULL,
                 J = 20, tol = 0.005, min_pairs = 30) {
  # Arguments:
  # X: Matrix of directional covariates on the unit sphere
  # Y: Vector of scalar responses
  # h: Bandwidth parameter to evaluate
  # p: Polynomial degree for local directional fitting
  # R: n x n error correlation matrix (computed via correlation_matrix if NULL)
  # D: n x n normalized geodesic distance matrix
  # h_pilot: Fixed pilot bandwidth used to estimate R (required if R is NULL)
  # J: Number of distance bins used for empirical variogram estimation
  # tol: Half-width tolerance for pairing points into distance bins
  # min_pairs: Minimum required pair count per distance bin
  
  n <- nrow(X)
  
  # Estimate error correlation matrix R if not provided
  if (is.null(R)) {
    if (is.null(D)) D <- geodesic_dist(X) / pi
    if (is.null(h_pilot))
      stop("Supply h_pilot: using the h being scored as the pilot makes R change with h.")
    R <- correlation_matrix(X, Y, h = h_pilot, p = p, J = J, tol = tol,
                            D = D, min_pairs = min_pairs)$R
  }
  
  # Fit full model to extract fitted values and smoother matrix S
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]      # Fortran branch: drop h index
  Yhat <- as.numeric(fit$Yhat)
  
  trSR  <- sum(S * t(R))                        # tr(S_p R)
  denom <- 1 - (1/n) * trSR
  if (!is.finite(denom) || denom <= 0) return(Inf)
  
  mean(((Y - Yhat) / denom)^2, na.rm = TRUE)
}

# Empirical semivariogram of the pilot residuals on the grid built by
# gamma_hat(d_j) = 1/(2 n(d_j,t)) * sum_{(i,j) in S(d_j,t)} (eps_i - eps_j)^2.

empirical_variogram <- function(eps, dist) {
  # eps: Vector of residuals from the pilot fit, length n
  # dist: Output of dist_est(), supplying the lags and the pairs in each bin
  
  v <- (eps[dist$S_set$i] - eps[dist$S_set$j])^2
  rs <- rowsum(v, dist$S_set$bin)              # one row per non-empty bin
  
  num <- numeric(dist$J)
  num[as.integer(rownames(rs))] <- rs[, 1]   # place each sum at its bin index
  
  # a data frame with one row per lag: the distance d, the semivariogram
  # estimate gamma, and the number of pairs n that produced it.
  data.frame(d = dist$d, gamma = num / (2 * dist$counts), n = dist$counts)
}


# Construct truncated quantile distance grid and the index sets S(d_j, t)
# Builds the lag distance grid \eqn{d_1, \dots, d_J} and assigns observation pairs 
# to distance bins for empirical semivariogram estimation on the unit sphere

dist_est <- function(X, D = NULL, J = 20, tol = 0.005, max_lag_prop = 0.3) {
  
  # X: Matrix of directional covariates on the unit sphere
  # D: n x n normalized geodesic distance matrix 
  # J: Number of distance bins (lags) in the grid
  # tol: Half-width of the distance bins; capped at half the smallest lag
  # spacing so the bins do not overlap
  # max_lag_prop : Proportion of shortest pairwaise distances to retain. 
  # Truncates the distance domain to focus on short lags and prevent binning on
  # the variogram sill.
  if (is.null(D)) D <- geodesic_dist(X) / pi

  dvec <- D[upper.tri(D)]
  pairs <- which(upper.tri(D), arr.ind = TRUE)
  
  # Filter out non-finite or zero (self-distance) entries
  keep <- is.finite(dvec) & dvec > 0
  dvec <- dvec[keep]
  pairs <- pairs[keep, , drop = FALSE]
  if (!length(dvec)) stop("No positive pairwise distances.")

  # Truncated quantile grid xonstruction ------------------------------------
  # Determine max distance threshold capturing the bottom `max_lag_prop` of pairs
  max_d <- quantile(dvec, probs = max_lag_prop, names = FALSE)
  short_dvec <- dvec[dvec <= max_d]
  
  # Generate J interior quantiles uniformly across the truncated short-lag distribution
  probs <- seq_len(J) / (J + 1)
  d <- as.numeric(quantile(short_dvec, probs = probs, names = FALSE))
  
  ## --- tolerance: cannot exceed half the smallest spacing, or the bins overlap
  tol <- min(tol, min(diff(d)) / 2)
  
  ## --- bin label: which d_j each pair belongs to, NA if none --------------
  bin <- rep(NA_integer_, length(dvec))
  for (j in seq_len(J)){
    bin[dvec >= d[j] - tol & dvec < d[j] + tol] <- j
  }

  S_set <- data.frame(i = pairs[, 1], j = pairs[, 2], dist = dvec, bin = bin)
  S_set <- S_set[!is.na(S_set$bin), ]
  
  # the lags d, the tolerance, the data frame S_set with
  # one row per retained pair (its two indices, its distance and its bin), the
  # pair count per bin, and the full vector of pairwise distances.
  list(d = d, probs = probs, tol = tol, J = J,
       S_set = S_set, counts = tabulate(S_set$bin, nbins = J), dvec = dvec)
}



# Estimates the n x n error correlation matrix R required by the MGCV
# criterion. Residuals from a pilot fit give sigma2 and an empirical
# semivariogram; the exponential correlation model is then fitted to that
# semivariogram by nonlinear least squares (NLS), and R is rebuilt from the fitted
# range parameter.

correlation_matrix <- function(X, Y, h, p, J = 20, tol = 0.005,
                               D = NULL, min_pairs = 30, max_lag_prop = 0.3) {
  
  # X: Matrix of directional covariates on the unit sphere
  # Y: Vector of scalar responses
  # h: Pilot bandwidth used for the fit whose residuals feed the semivariogram
  # p: Polynomial degree for the pilot fit (0 = NW, 1 = local linear)
  # J: Number of distance bins used for the empirical semivariogram
  # tol: Half-width tolerance for pairing points into distance bins
  # D: n x n normalized geodesic distance matrix (computed from X if NULL)
  # min_pairs: Minimum pair count for a bin to enter the fit
  # max_lag_prop : Proportion of shortest pairwaise distances to retain. 
  
  n <- nrow(X)
  q <- ncol(X) - 1
  if (is.null(D)) D <- geodesic_dist(X) / pi

  ## --- pilot fit and residuals --------------------------------------------
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  eps <- Y - as.numeric(fit$Yhat)

  # DF correction for variance
  S <- fit$weights
  trS <- if(!is.null(dim(S))) sum(diag(S[,,1])) else sum(diag(S))
  sigma2 <- sum(eps^2) / max(1, (n - trS))

  # Truncated quantile grid and empirical variogram
  g <- dist_est(X, D = D, J = J, tol = tol, max_lag_prop = max_lag_prop)
  emp_var <- empirical_variogram(eps, g)

  # Keep finite bins with enough pairs
  ok <- is.finite(emp_var$gamma) & emp_var$n >= min_pairs & emp_var$d > 0
  if (sum(ok) < 3) stop("Cannot estimate alpha: too few usable bins.")

  d_j <- emp_var$d[ok]
  emp_j <- emp_var$gamma[ok]
  counts_j <- emp_var$n[ok]

  # Weighted NLS
  # gamma(d) = sigma2 * (1 - exp(-n^(1/q) d / alpha))
  # alpha-hat = argmin_a sum_j [ gamma-hat(d_j) - sigma2 (1 - exp(-n^(1/q) d_j / a)) ]^2
  # optimized over log(a) so that a > 0 holds automatically
  obj <- function(log_alpha) {
    a <- exp(log_alpha)
    sum(counts_j * (emp_j - sigma2 * (1 - exp(-n^(1/q) * d_j / a)))^2)
  }

  alpha <- exp(optimize(obj, interval = log(c(1e-4, 100)))$minimum)

  R <- exp(- D * n^(1/q) / alpha)
  diag(R) <- 1

  list(R = R, alpha = alpha, sigma2 = sigma2, variogram = emp_var)
}



#===============================================================================
# CROSS-VALIDATION APPROACH (for independent data)
#===============================================================================
library(MASS)
library(DirStats) 

test <- FALSE

if (test) {
  set.seed(1)
  n <- 400 
  alpha_true <- 1.2 
  q <- 2
  
  # 1. Data Generation
  X <- unif_sphere(n, q) 
  D <- geodesic_dist(X) / pi
  R_true <- build_Sigma(X, alpha_true, sigma2 = 1)
  
  # m(X) = X_1 with spatially correlated errors
  Y <- X[, 1] + as.numeric(MASS::mvrnorm(1, mu = rep(0, n), Sigma = R_true))
  
  ## -------------------------------------------------------------------
  ## TEST 1: CV and CV simplified yield same result
  ## -------------------------------------------------------------------
  hs <- c(0.2, 0.4, 0.8)
  res_cv <- cbind(
    slow = sapply(hs, function(h) cv_loo(X, Y, h, p = 1)),
    fast = sapply(hs, function(h) cv(X, Y, h, p = 1))
  )
  print(res_cv)
  
  ## -------------------------------------------------------------------
  ## TEST 2: MGCV with R = I reproduces standard GCV
  ## -------------------------------------------------------------------
  res_gcv <- cbind(
    gcv  = sapply(hs, function(h) gcv(X, Y, h, p = 1)),
    mgcv = sapply(hs, function(h) mgcv(X, Y, h, p = 1, R = diag(n)))
  )
  print(res_gcv)
  
  ## -------------------------------------------------------------------
  ## TEST 3: Spatial Parameter Recovery (alpha and sigma2)
  ## -------------------------------------------------------------------
  h_pilot <- bw_dir_rot(X)
  cm <- correlation_matrix(X, Y, h = h_pilot, p = 1, D = D, max_lag_prop = 0.3)
  
  res_params <- c(
    alpha_hat  = cm$alpha, 
    alpha_true = alpha_true, 
    sigma2_hat = cm$sigma2,
    sigma2_true = 1.0
  )
  print(round(res_params, 4))
  
  ## -------------------------------------------------------------------
  ## TEST 4: Fitted Variogram vs Empirical Points
  ## -------------------------------------------------------------------
  ev <- cm$variogram
  
  # Filter for valid/used bins to prevent plotting NA or empty bins
  valid_bins <- is.finite(ev$gamma) & ev$n > 0
  
  if (sum(valid_bins) > 0) {
    plot(ev$d[valid_bins], ev$gamma[valid_bins], pch = 19, 
         xlab = "Normalized Geodesic Distance (d)", 
         ylab = expression(hat(gamma)(d)),
         xlim = c(0, max(ev$d[valid_bins], na.rm = TRUE) * 1.05),
         ylim = c(0, max(cm$sigma2 * 1.2, max(ev$gamma[valid_bins], na.rm = TRUE))))
    
    # Superimpose fitted theoretical exponential model
    curve(cm$sigma2 * (1 - exp(- (n^(1/q)) * x / cm$alpha)), 
          add = TRUE, col = "firebrick", lwd = 2)
    
    # Superimpose estimated sill level (sigma^2)
    abline(h = cm$sigma2, lty = 2, col = "blue")
    legend("bottomright", legend = c("Empirical bins", "Fitted model", "Estimated sill"),
           col = c("black", "firebrick", "blue"), pch = c(19, NA, NA), 
           lty = c(NA, 1, 2), lwd = c(NA, 2, 1))
  } else {
    warning("No valid variogram bins available to plot.")
  }
  
  ## -------------------------------------------------------------------
  ## TEST 5: MGCV Bandwidth Selection Comparison
  ## -------------------------------------------------------------------
  hg <- seq(0.15, 1.5, length.out = 30) 
  
  # 1. Optimal h using estimated correlation matrix R_hat
  mgcv_Rhat_scores <- sapply(hg, function(h) mgcv(X, Y, h, p = 1, R = cm$R))
  h_Rhat <- hg[which.min(mgcv_Rhat_scores)]
  
  # 2. Optimal h using true correlation matrix R_true
  mgcv_Rtrue_scores <- sapply(hg, function(h) mgcv(X, Y, h, p = 1, R = R_true))
  h_Rtrue <- hg[which.min(mgcv_Rtrue_scores)]
  
  # 3. Standard CV (Ignorant of spatial correlation)
  cv_scores <- sapply(hg, function(h) cv(X, Y, h, p = 1))
  h_CV <- hg[which.min(cv_scores)]
  
  # ASE (Average Squared Error against true signal X[, 1])
  ase_scores <- sapply(hg, function(h) {
    fit_h <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = 1)
    mean((as.numeric(fit_h$Yhat) - X[, 1])^2)
  })
  h_ASE <- hg[which.min(ase_scores)]
  
  res_h <- c(
    h_Rhat  = h_Rhat,  # MGCV with estimated R
    h_Rtrue = h_Rtrue, # MGCV with true R
    h_CV    = h_CV,    # Standard CV (typically undersmoothes)
    h_ASE   = h_ASE    # bandwidth minimizing true risk
  )
  print(round(res_h, 4))
  
}