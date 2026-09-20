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
# dist_est():  gamma-hat(d_j) = 1/(2 n(d_j,t)) * sum_{(i,j) in S(d_j,t)}
# (eps_i - eps_j)^2.

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





# Builds the distance grid and the index sets S(d_j, t) for the empirical
# semivariogram. The J lags are interior quantiles of the pairwise geodesic
# distance sample, and each lag collects the pairs lying within +/- tol of it.

dist_est <- function(X, D = NULL, J = 20, tol = 0.005) {
  # X: Matrix of directional covariates on the unit sphere
  # D: n x n normalized geodesic distance matrix 
  # J: Number of distance bins (lags) in the grid
  # tol: Half-width of the distance bins; capped at half the smallest lag
  # spacing so the bins do not overlap
  
  if (is.null(D)) D <- geodesic_dist(X) / pi 
  
  ## --- pairwise geodesic distances, i < j ---------------------------------
  dvec <- D[upper.tri(D)]
  pairs <- which(upper.tri(D), arr.ind = TRUE)
  
  keep <- is.finite(dvec) & dvec > 0
  dvec <- dvec[keep]
  pairs <- pairs[keep, , drop = FALSE]
  if (!length(dvec)) stop("No positive pairwise distances.")
  
  ## --- grid: J interior quantiles of the distance sample ------------------
  probs <- seq_len(J) / (J + 1)
  d <- as.numeric(quantile(dvec, probs = probs, names = FALSE))
  
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
# semivariogram by nonlinear least squares, and R is rebuilt from the fitted
# range parameter.

correlation_matrix2 <- function(X, Y, h, p, J = 20, tol = 0.005,
                               D = NULL, min_pairs = 30) {
  # X: Matrix of directional covariates on the unit sphere
  # Y: Vector of scalar responses
  # h: Pilot bandwidth used for the fit whose residuals feed the semivariogram
  # p: Polynomial degree for the pilot fit (0 = NW, 1 = local linear)
  # J: Number of distance bins used for the empirical semivariogram
  # tol: Half-width tolerance for pairing points into distance bins
  # D: n x n normalized geodesic distance matrix (computed from X if NULL)
  # min_pairs: Minimum pair count for a bin to enter the fit
  n <- nrow(X)
  q <- ncol(X) - 1 
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  ## --- pilot fit and residuals --------------------------------------------
  yhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                 h = h, p = p)$Yhat
  eps <- Y - as.numeric(yhat)
  sigma2 <- mean(eps^2)                               
  
  ## --- empirical semivariogram on the quantile grid -----------------------
  g <- dist_est(X, D = D, J = J, tol = tol)
  emp_var <- empirical_variogram(eps, g)    
  
  # Keep every bin with enough pairs. The least-squares method does not need
  # gamma < sigma2
  ok <- is.finite(emp_var$gamma) & emp_var$n >= min_pairs & emp_var$d > 0
  if (sum(ok) < 3) stop("Too few usable bins to fit the variogram.")
  
  d_j <- emp_var$d[ok]
  emp_j <- emp_var$gamma[ok]
  
  ## --- alpha calculated by least squares ------------------------------------
  ## gamma(d) = sigma2 * (1 - exp(-n^(1/q) d / alpha))
  ## alpha-hat = argmin_a sum_j [ gamma-hat(d_j) - sigma2 (1 - exp(-n^(1/q) d_j / a)) ]^2
  ## optimised over log(a) so that a > 0 holds automatically
  obj <- function(log_alpha) {
    a <- exp(log_alpha)
    sum((emp_j - sigma2 * (1 - exp(-n^(1/q) * d_j / a)))^2)
  }
  alpha <- exp(optimize(obj, interval = log(c(1e-3, 1e3)))$minimum)
  
  ## --- Estimation of R = R(theta-hat) --------------------------------------
  R <- exp(- D * n^(1/q) / alpha)
  diag(R) <- 1     
  
  
  # the estimated correlation matrix R, the fitted range alpha, the
  # residual variance sigma2, the empirical semivariogram, the logical vector of
  # bins used, and the residual sum of squares of the fit.
  list(R = R, alpha = alpha, sigma2 = sigma2,
       variogram = emp_var, used = ok, fit_ss = obj(log(alpha)))
}


correlation_matrix <- function(X, Y, h, p, J = 20, tol = 0.005,
                               D = NULL, min_pairs = 30) {
  
  n <- nrow(X)
  q <- ncol(X) - 1 
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  ## --- pilot fit and residuals --------------------------------------------
  yhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                 h = h, p = p)$Yhat
  eps <- Y - as.numeric(yhat)
  sigma2 <- mean(eps^2)                                          
  
  ## --- empirical semivariogram on the quantile grid -----------------------
  g <- dist_est(X, D = D, J = J, tol = tol)
  emp_var <- empirical_variogram(eps, g)                         # (10)
  
  # gamma <- pmin(emp_var$gamma, 0.95 * sigma2)
  ## --- alpha_j, method of moments -----------------------------------------
  ## rho(d) = exp(-d/alpha)  =>  alpha_j = d_j / [ln(sigma2) - ln(sigma2 - gamma(d_j))]
  check <- is.finite(emp_var$gamma) & emp_var$gamma < sigma2 &
    emp_var$d > 0 & emp_var$n >= min_pairs
  if (!any(check))
    stop("Cannot estimate alpha: no bin has gamma-hat below sigma2-hat.")
  
  # Invert rho(d) = exp(-d / alpha) -> alpha = d / (ln(sigma2) - ln(sigma2 - gamma))
  alpha_j <- rep(NA_real_, J)
  alpha_j[check] <- (emp_var$d[check] * n^(1/q)) / 
    (log(sigma2) - log(sigma2 - emp_var$gamma[check]))
  
  alpha <- mean(alpha_j[check], na.rm = TRUE)
  
  ## --- Estimation of R --------
  R <- exp(- D * n^(1/q) / alpha)
  diag(R) <- 1     
  
  list(R = R, alpha = alpha, alpha_j = alpha_j, sigma2 = sigma2,
       variogram = emp_var, used = check)
}



# 
# ## UPDATED: Truncated quantile distance grid
# dist_est <- function(X, D = NULL, J = 20, tol = 0.005, max_lag_prop = 0.3) {
#   if (is.null(D)) D <- geodesic_dist(X) / pi 
#   
#   dvec <- D[upper.tri(D)]
#   pairs <- which(upper.tri(D), arr.ind = TRUE)
#   
#   keep <- is.finite(dvec) & dvec > 0
#   dvec <- dvec[keep]
#   pairs <- pairs[keep, , drop = FALSE]
#   if (!length(dvec)) stop("No positive pairwise distances.")
#   
#   # Filter for short lags, then calculate J quantiles
#   max_d <- quantile(dvec, probs = max_lag_prop, names = FALSE)
#   short_dvec <- dvec[dvec <= max_d]
#   
#   probs <- seq_len(J) / (J + 1)
#   d <- as.numeric(quantile(short_dvec, probs = probs, names = FALSE))
#   tol <- min(tol, min(diff(d)) / 2)
#   
#   bin <- rep(NA_integer_, length(dvec))
#   for (j in seq_len(J)){
#     bin[dvec >= d[j] - tol & dvec < d[j] + tol] <- j
#   }
#   
#   S_set <- data.frame(i = pairs[, 1], j = pairs[, 2], dist = dvec, bin = bin)
#   S_set <- S_set[!is.na(S_set$bin), ]
#   
#   list(d = d, probs = probs, tol = tol, J = J,
#        S_set = S_set, counts = tabulate(S_set$bin, nbins = J), dvec = dvec)
# }
# 
# # UPDATED: Correlation matrix with weighted NLS
# correlation_matrix <- function(X, Y, h, p, J = 20, tol = 0.005,
#                                D = NULL, min_pairs = 30, max_lag_prop = 0.3) {
#   n <- nrow(X)
#   q <- ncol(X) - 1 
#   if (is.null(D)) D <- geodesic_dist(X) / pi
#   
#   # Pilot fit
#   fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
#   eps <- Y - as.numeric(fit$Yhat)
#   
#   # DF correction for variance
#   S <- fit$weights
#   trS <- if(!is.null(dim(S))) sum(diag(S[,,1])) else sum(diag(S))
#   sigma2 <- sum(eps^2) / max(1, (n - trS)) 
#   
#   # Truncated quantile grid and empirical variogram
#   g <- dist_est(X, D = D, J = J, tol = tol, max_lag_prop = max_lag_prop)
#   emp_var <- empirical_variogram(eps, g)    
#   
#   # Keep finite bins with enough pairs
#   ok <- is.finite(emp_var$gamma) & emp_var$n >= min_pairs & emp_var$d > 0
#   if (sum(ok) < 3) stop("Cannot estimate alpha: too few usable bins.")
#   
#   d_j <- emp_var$d[ok]
#   emp_j <- emp_var$gamma[ok]
#   counts_j <- emp_var$n[ok]
#   
#   # Weighted NLS
#   obj <- function(log_alpha) {
#     a <- exp(log_alpha)
#     sum(counts_j * (emp_j - sigma2 * (1 - exp(-n^(1/q) * d_j / a)))^2)
#   }
#   
#   alpha <- exp(optimize(obj, interval = log(c(1e-4, 100)))$minimum)
#   
#   R <- exp(- D * n^(1/q) / alpha)
#   diag(R) <- 1    
#   
#   list(R = R, alpha = alpha, sigma2 = sigma2, variogram = emp_var)
# }

#===============================================================================
# CROSS-VALIDATION APPROACH (for independent data)
#===============================================================================


test <- FALSE
if (test){
  set.seed(1)
  n <- 400; 
  alpha_true <- 1.2; 
  q <- 2
  X <- unif_sphere(n, q); 
  D <- geodesic_dist(X)/pi
  Y <- X[,1] + as.numeric(MASS::mvrnorm(1, rep(0,n), build_Sigma(X, alpha_true, 1)))
  
  ## TEST 1: CV and CV simplified yield same result
  hs <- c(0.2, 0.4, 0.8)
  cbind(slow = sapply(hs, function(h) cv_loo(X, Y, h, 1)),
        fast = sapply(hs, function(h) cv(X, Y, h, 1)))
  
  ## TEST 2: MGCV with R = I reproduces gcv
  cbind(gcv  = sapply(hs, function(h) gcv(X, Y, h, 1)),
        mgcv = sapply(hs, function(h) mgcv(X, Y, h, 1, R = diag(n))))
  
  ## 3. alpha-hat recupera alpha  (lo importante)
  cm <- correlation_matrix(X, Y, h =  bw_dir_rot(X), p = 1, D = D)
  c(alpha_hat = cm$alpha, alpha_true = alpha_true, sigma2 = cm$sigma2, used = sum(cm$used))
  
  ## 4. el variograma ajustado sigue a los puntos empíricos
  ev <- cm$variogram
  plot(ev$d, ev$gamma, pch = 19, xlab = "d", ylab = expression(hat(gamma)(d)))
  curve(cm$sigma2 * (1 - exp(-n^(1/q)*x/cm$alpha)), add = TRUE, col = "firebrick", lwd = 2)
  abline(h = cm$sigma2, lty = 2)
  
  ## 5. MGCV con R estimada frente a R verdadera y al oráculo
  R0 <- build_Sigma(X, alpha_true, 1)      # sigma2 = 1, luego esto ES R
  hg <- seq(0.05, 2, length.out = 40)
  c(h_Rhat  = hg[which.min(sapply(hg, function(h) mgcv(X, Y, h, 1, R = cm$R)))],
    h_Rtrue = hg[which.min(sapply(hg, function(h) mgcv(X, Y, h, 1, R = R0)))],
    h_CV    = hg[which.min(sapply(hg, function(h) cv(X, Y, h, 1)))],
    h_ASE   = hg[which.min(sapply(hg, function(h)
      mean((as.numeric(loc.directional.linear(x = X, data.dir = X,
                                              data.lin = Y, h = h, p = 1)$Yhat) - X[,1])^2)))])
}
