################################################################################
# Simulation study: Linear-spherical regression under spatial dependence
# Bandwidth selection: standard CV, GCV, vs. Non-Parametric MGCV (npsp)
################################################################################
rm(list = ls())

if (!requireNamespace("DirStatsOld", quietly = TRUE)) {
  install.packages("DirStatsOld_0.1.5.tar.gz", repos = NULL, type = "source")
}

library(DirStatsOld)
library(DirStats)
library(MASS)
library(npsp)
library(foreach)
library(progressr)
library(future)
library(doRNG)
library(doFuture)

################################################################################
# Initial functions
################################################################################

##------ Spherical distances -------------------------------------------------
geodesic_dist <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  ip <- X_norm %*% t(X_norm)
  ip <- pmax(pmin(ip, 1), -1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
}
##------ Regression functions on S^2 -----------------------------------------
m_funs <- list(
  m1 = function(X) X[, 1],
  m2 = function(X) sin(pi * X[, 1]) * X[, 2],
  m3 = function(X, a = 1, b = 1.5) a * sin(2 * pi * X[, 2]) + b * cos(2 * pi * X[, 1])
)

##------ Uniform sample on S^d -----------------------------------------------
unif_sphere <- function(n, d) {
  X <- matrix(rnorm(n * (d + 1)), nrow = n, ncol = d + 1)
  return(X / sqrt(rowSums(X^2)))
}


################################################################################
##------ Cross-Validation & MGCV Criteria ------------------------------------
################################################################################
cv_loo <- function(X, Y, h, p) {
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  Yhat <- as.numeric(fit$Yhat)
  
  denom <- 1 - diag(S)
  if (any(denom <= 1e-5)) return(Inf)
  if (any(!is.finite(denom))) return(Inf)
  mean(((Y - Yhat) / denom)^2)
}

gcv <- function(X, Y, h, p) {
  n <- nrow(X)
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  Yhat <- as.numeric(fit$Yhat)
  
  denom <- 1 - sum(diag(S)) / n
  if (!is.finite(denom) || denom <= 0) return(Inf)
  
  mean(((Y - Yhat) / denom)^2)
}

################################################################################ 
# Isotropic variogram estimation from pairwise distances 
################################################################################
svariso.from.pairs <- function(dist, gamma, maxlag = NULL, nlags = 101,
                               h = NULL, degree = 1,
                               hat.bin = TRUE) {
  
  dist <- as.numeric(dist)
  gamma <- as.numeric(gamma)
  
  correct <- is.finite(dist) & is.finite(gamma)
  dist <- dist[correct]
  gamma <- gamma[correct]
  
  if (length(dist) == 0L){
    stop("No valid pairs available for variogram estimation")
  }
  
  # If maxlag is not supplied, use 55% of the maximum pairwise distance
  if (is.null(maxlag))  maxlag <- 0.55 * max(dist, na.rm = TRUE)
  if (!is.numeric(maxlag) ||  length(maxlag) != 1L || !is.finite(maxlag) ||
      maxlag <= 0) {
    stop("'maxlag' must be a positive finite number")
  } 
  
  # Keep only pairs below maximum lag 
  ok <- dist <= maxlag
  dist <- dist[ok]
  gamma <- gamma[ok]
  
  if (length(dist) < 10L) warning("Too few pairs available below 'maxlag'")
  
  # Linear binning with npsp
  bin <- npsp::binning(x = dist, y = gamma, nbin = nlags,type = "linear",
                       set.NA = TRUE)
  
  # Make it a svar.bin object
  bin$svar <- list( type = "isotropic", estimator = "classical")
  class(bin) <- c("svar.bin", class(bin))
  
  h.cv.result <- NULL
  h.cv.value <- NA_real_
  
  #Bandwidth selection
  if (is.null(h)) {
    # Automatic bandwidth selection by MRSE
    h.cv.result <- npsp::h.cv(bin, loss = "MRSE", degree = degree, ncv = 1)
    # Selected bandwidth 
    h <- as.numeric(h.cv.result$h)
    # MRSE criterion at the selected bandwidth
    h.cv.value <- h.cv.result$value
  }
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) { 
    stop("'h' must be a positive finite number") }
  
  # Local linear smoothing with npsp
  svar <- npsp::locpol(bin, h = h, degree = degree,  drv = FALSE,
                       hat.bin = hat.bin, ncv = 0)
  
  if (is.null(svar$est) || !any(is.finite(svar$est))) { 
    stop("The estimated variogram contains no finite values") }
  
  svar$directional <- list(maxlag = maxlag,nlags = nlags,
                           npairs = length(dist), 
                           h = h, h.cv.value = h.cv.value, h.cv = h.cv.result)
  
  return(svar)
}

################################################################################
# Covariance matrix from the directional variogram with the geodesic distance
################################################################################
varcov.directional <- function(svar, D, sill = NULL) {
  
  if (is.null(sill)) sill <- max(svar$est, na.rm = TRUE)
  if (!is.finite(sill) || sill <= 0) stop("'sill' must be positive and finite")
  n <- nrow(D)
  dists <- D[lower.tri(D)]
  covs <- numeric(length(dists))              # 0 más allá de maxlag (taper)
  idx <- dists <= svar$grid$max
  covs[idx] <- npsp::covar(svar, dists[idx], sill = sill)
  C <- matrix(0, n, n)
  C[lower.tri(C)] <- covs
  C <- C + t(C)
  diag(C) <- sill
  C
}

################################################################################
# Complete directional variogram estimation
################################################################################
variogram.est <- function(X, Y, D = NULL, h_reg = NULL, h_var = NULL, maxlag = NULL, 
                          max_iter = 15, nlags = 101, tol = 0.05){
  
  stopifnot(!missing(X), !missing(Y))
  n <- nrow(X)
  if (length(Y) != n) stop("'X' and 'Y' have incompatible dimensions")
  
  if (is.null(D)) D <- geodesic_dist(X) / pi
  if (is.null(maxlag)) maxlag <- 0.55 * max(D)
  
  # Initial bandwidth assuming independence
  if(is.null(h_reg)){
    h_reg <- DirStatsOld::bw.pi.loc(data.dir = X, data.lin = Y, p = 1)$h.opt
  }
  
  # Initial estimate regression function
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                h = h_reg, p = 1)
  
  # Extract smoother matrix S
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L] # Ensure S is a 2D matrix
  
  # Obtain residuals
  Yhat <- as.numeric(fit$Yhat)
  residuals <- Y - Yhat
  
  rm(fit,Yhat)
  
  # Paiswise distances i < j
  ind <- upper.tri(D)
  dist <- D[ind]
  
  # Classical semivariogram:
  # gamma_ij = 1/2 (residuals_i - residuals_j)^2
  ii <- row(D)[ind]
  jj <- col(D)[ind]
  gamma <- 0.5 * (residuals[ii] - residuals[jj])^2
  
  correct <- is.finite(dist) & is.finite(gamma)
  
  dist <- dist[correct]
  gamma <- gamma[correct]
  ii <- ii[correct]
  jj <- jj[correct]
  
  if (length(dist) == 0L){
    stop("No valid pairs available for variogram estimation")
  }
  
  
  # Initial nonparametric variogram; h_var = NULL, selected  minimizing MRSE.
  svar <- svariso.from.pairs(dist = dist, gamma = gamma, maxlag = maxlag,
                             nlags = nlags, h = h_var, degree = 1, 
                             hat.bin = TRUE)
  h_var <- svar$directional$h
  maxlag <- svar$directional$maxlag
  
  # --------------------------------
  # Iterative bias correction
  # --------------------------------
  error <- Inf
  
  for (iter in 1:max_iter) {
    
    # Covariance matrix estimated from the current nonparametric variogram
    C_hat <- varcov.directional(svar = svar, D = D, sill = NULL)
    
    SC <- S %*% C_hat
    
    # Bias matrix
    B_hat <- SC %*% t(S) - SC - t(SC)
    rm(SC, C_hat)
    
    # Pairwise bias correction
    b_diag <- diag(B_hat)
    bias_pair <- 0.5 * (b_diag[ii] + b_diag[jj]) - B_hat[cbind(ii, jj)]
    
    # Empirical residual semivariogram
    gamma.resid <- 0.5 * ( residuals[ii] - residuals[jj])^2
    
    # Bias-corrected pairwise semivariogram
    gamma.corrected <- gamma.resid - bias_pair
  
    rm(B_hat, b_diag, gamma.resid, bias_pair)
    
    
    ok <- is.finite(dist) & is.finite(gamma.corrected) 
    dist.corrected <- dist[ok]
    gamma.corrected <- gamma.corrected[ok]
    if (length(dist.corrected) < 10L) { stop("Too few valid pairs after bias correction") }
    
    # Re-estimate the variogram from corrected pairs
    # If h_var was given, keep it fixed.
    # If h_var was NULL, select a new bandwidth by MRSE.
    
    svar.new <- svariso.from.pairs(dist = dist.corrected,
                                   gamma = gamma.corrected, maxlag = maxlag,
                                   nlags = nlags, h = h_var, degree = 1, 
                                   hat.bin = TRUE)
    
    # Relative squared error between consecutive variograms
    #error <- mean((svar$est / svar.new$est - 1)^2, na.rm = TRUE)
    denom <- max(abs(svar$est), na.rm = TRUE)
    
    error <- sqrt(mean((svar.new$est - svar$est)^2, na.rm = TRUE)) / max(denom, 1e-8)
    
    # Updates
    svar <- svar.new
    h_var <- svar$directional$h
    rm(svar.new, dist.corrected, gamma.corrected)
    
    # Convergence
    if (is.finite(error) && error < tol) break
    
  }
  
  # Final Shapiro–Botha variogram model
  svm <- npsp::fitsvar.sb.iso(svar, dk = 0)
  if (!is.finite(svm$sill) || svm$sill <= 0) { 
    stop("The fitted Shapiro-Botha sill is not positive and finite") }
  
  # return(
  #   list(fit = fit, Yhat = Yhat, residuals = residuals, S = S, D = D,
  #        svar = svar,svm = svm, h_reg = h_reg, h_var = h_var, maxlag = maxlag,
  #        nlags = nlags, iter = iter,error = error)
  # )
  return(list(svm = svm, h_reg = h_reg, h_var = h_var, maxlag = maxlag, iter = iter, 
              error = error, sill = svm$sill, svar = svar))
  
}

################################################################################
# Covariance matrix from a fitted Shapiro-Botha model
################################################################################
varcov.svm.directional <- function(svm, D) {
  
  n <- nrow(D)
  
  if (ncol(D) != n) stop("'D' must be a square matrix")
  if (!isTRUE(all.equal(D, t(D), tolerance = 1e-10))){
    stop("'D' must be a symmetric distance matrix")}
  
  # Distances corresponding to the lower triangle
  dists <- D[lower.tri(D)]
  
  # Covariance from the fitted variogram model
  covs <- npsp::covar(svm, dists, sill = svm$sill)
  if (any(!is.finite(covs))) { stop("The fitted covariance contains non-finite values") }
  
  # Build covariance matrix
  Sigma <- matrix(0,nrow = n, ncol = n)
  Sigma[lower.tri(Sigma)] <- covs
  Sigma <- Sigma + t(Sigma)
  
  # Variance at zero
  diag(Sigma) <- svm$sill
  
  return(Sigma)
}

################################################################################
# Modified Generalized Cross-Validation for directional regression
################################################################################

mgcv <- function(X, Y, h, p, R) {
  
  n <- nrow(X)
  
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  
  Yhat <- as.numeric(fit$Yhat)
  trSR <- sum(S * R)
  denom <- 1 - (1/n) * trSR
  if (!is.finite(denom) || denom <= 0) return(Inf)
  
  res <- ((Y - Yhat) / denom)^2
  if (any(!is.finite(res))) return(Inf)
  mean(res)
}

################################################################################
# Simulation Functions
################################################################################
one_rep <- function(n, alpha, sigma2, m_fun, h_grid, d = 2,
                    h_pilot = NULL, do_gcv = FALSE) {
  
  
  X <- unif_sphere(n, d)
  m_vals <- m_fun(X)
  
  # Normalized geodesic distance
  D <- geodesic_dist(X) / pi
  
  Sigma <- sigma2 * exp(-D / alpha)
  eps <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
  Y <- m_vals + eps
  
  # TRUE correlation matrix
  R_true <- Sigma / sigma2
  
  rm(Sigma, eps)
  
  
  if (is.null(h_pilot)) {
    h_pilot <- DirStatsOld::bw.pi.loc(data.dir = X, data.lin = Y, p = 1)$h.opt }
  
  # ASE function
  ase <- function(h, p) {
    yhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                   h = h, p = p)$Yhat
    mean((yhat - m_vals)^2, na.rm = TRUE)
  }

  out <- c()  
  # Estimate spatial variogram and Shapiro-Botha model
  vario <- tryCatch({
    v <- variogram.est(X = X, Y = Y, D = D, h_reg = h_pilot, h_var = NULL,
                       maxlag = 0.55 * max(D), max_iter = 15, nlags = 101, tol = 0.05)
    v$R <- varcov.svm.directional(v$svm, D) / v$svm$sill
    if (any(!is.finite(v$R))) stop("R no finita")
    v
  }, error = function(e) NULL)

  svm <- NULL
  R <- NULL
  
  if (!is.null(vario)) {
    
    svm <- vario$svm
    Sigma <- varcov.svm.directional(svm = svm, D = D)
    R <- Sigma / svm$sill
    
    # Check estimated correlation matrix
    if (any(!is.finite(R))) {
      stop("Estimated correlation matrix contains non-finite values")
    }
    
    if (max(abs(diag(R) - 1)) > 1e-8) {
      warning("Estimated correlation matrix does not have unit diagonal")
    }
    
    if (!isTRUE(all.equal(R, t(R), tolerance = 1e-10))) {
      warning("Estimated correlation matrix is not exactly symmetric")
    }
    
    R_error <- sqrt(mean(( R[lower.tri(R)] - R_true[lower.tri(R)])^2))
    
    } else {
      
    R <- NULL
    R_error <- NA_real_
  }
  
  for (p in c(0, 1)) {
    tag <- if (p == 0) "nw" else "ll"
    
    # ASE Benchmark
    ase_vals <- sapply(h_grid, function(h) ase(h, p))
    idx_case <- which.min(ase_vals)
    out[paste0(tag, "_h_case")] <- h_grid[idx_case]
    out[paste0(tag, "_ase_case")] <- ase_vals[idx_case]
    
    # LOO CV
    cv_vals <- sapply(h_grid, function(h) cv_loo(X, Y, h, p))
    idx_cv<- which.min(cv_vals)
    out[paste0(tag, "_h_cv")] <- h_grid[idx_cv]
    out[paste0(tag, "_ase_cv")] <- ase_vals[idx_cv]
    
    # GCV
    if (do_gcv) {
      gcv_vals <- sapply(h_grid, function(h) gcv(X, Y, h, p))
      idx_gcv  <- which.min(gcv_vals)
      out[paste0(tag, "_h_gcv")]  <- h_grid[idx_gcv]
      out[paste0(tag, "_ase_gcv")] <- ase_vals[idx_gcv]
    }
    
    
    # Non-parametric MGCV
    
    # If variogram estimation fails
    if (is.null(vario)) {
      out[paste0(tag, "_h_mgcv")] <- NA_real_
      out[paste0(tag, "_ase_mgcv")] <- NA_real_
      out[paste0(tag, "_h_mgcv_true")] <- NA_real_
      out[paste0(tag, "_ase_mgcv_true")] <- NA_real_
      out[paste0(tag, "_sigma2_hat")] <- NA_real_
      out[paste0(tag, "_R_error")] <- NA_real_
      next
    }
    # Estimated still
    out[paste0(tag, "_sigma2_hat")] <- svm$sill

    # Error of estimated correlation matrix
    out[paste0(tag, "_R_error")] <- R_error
    
    # Evaluate MGCV over h_grid
    mgcv_vals <- sapply(h_grid, function(h) {
      tryCatch( mgcv( X = X, Y = Y, h = h, p = p, R = R),
                error = function(e) Inf)}
    )
    
    # Select MGCV bandwidth
    if (all(!is.finite(mgcv_vals))) {
      out[paste0(tag, "_h_mgcv")] <- NA_real_
      out[paste0(tag, "_ase_mgcv")] <- NA_real_
    } else {
      idx_mgcv <- which.min(mgcv_vals)
      out[paste0(tag, "_h_mgcv")] <-h_grid[idx_mgcv]
      out[paste0(tag, "_ase_mgcv")] <- ase_vals[idx_mgcv]
    }
    
    # MGCV using TRUE R
    mgcv_true_vals <- sapply( h_grid,function(h) 
      { tryCatch( mgcv( X = X, Y = Y, h = h, p = p, R = R_true),
          error = function(e) Inf)
      }
    )
    
    if (all(!is.finite(mgcv_true_vals))) {
      out[paste0(tag, "_h_mgcv_true")] <- NA_real_
      out[paste0(tag, "_ase_mgcv_true")] <- NA_real_
    } else {
      idx_mgcv_true <- which.min(mgcv_true_vals)
      out[paste0(tag, "_h_mgcv_true")] <- h_grid[idx_mgcv_true]
      out[paste0(tag, "_ase_mgcv_true")] <- ase_vals[idx_mgcv_true]
    }
  }
  
  return(out)
}

run_simulation <- function(MC, n_values, alpha_vals, sigma2, m_idx,
                           h_grid, d = 2, cores = 1,
                           h_pilot = NULL, do_gcv = FALSE) {
  
  doFuture::registerDoFuture()
  future::plan(future::multisession(), workers = cores)
  
  handlers(handler_progress(
    format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
    clear  = FALSE
  ))
  
  results <- list()
  
  for (n in n_values) {
    for (alpha_idx in seq_along(alpha_vals)) {
      alpha <- alpha_vals[alpha_idx]
      
      cat(sprintf("\n--- m%d | n = %d | alpha = %.2f ---\n", m_idx, n, alpha))
      
      progressr::with_progress({
        prog <- progressr::progressor(along = seq_len(MC))
        
        reps <- foreach(k = seq_len(MC), .inorder = TRUE,
                        .packages = c("DirStatsOld", "DirStats", "MASS", "npsp"),
                        .export = c("one_rep", "cv_loo", "gcv", "mgcv",
                                    "geodesic_dist", "unif_sphere",
                                    "m_funs","varcov.svm.directional",
                                    "svariso.from.pairs","varcov.directional", 
                                    "variogram.est","h_grid", "do_gcv", "sigma2",
                                    "d", "h_pilot", "m_idx")) %dorng% {
                                      prog()
                                      one_rep(n = n, alpha = alpha, sigma2 = sigma2,
                                              m_fun   = m_funs[[m_idx]],
                                              h_grid  = h_grid,
                                              d = d, do_gcv = do_gcv,
                                              h_pilot = h_pilot)
                                    }
      })
      
      mat <- do.call(rbind, reps)
      
      key <- paste0("n", n, "_a", alpha_idx)
      results[[key]] <- list(
        n = n, 
        alpha = alpha,
        mat = mat,
        
        mean_ase_nw_case = mean(mat[, "nw_ase_case"], na.rm = TRUE),
        mean_ase_ll_case = mean(mat[, "ll_ase_case"], na.rm = TRUE),
        sd_ase_nw_case   = sd(mat[, "nw_ase_case"], na.rm = TRUE),
        sd_ase_ll_case   = sd(mat[, "ll_ase_case"], na.rm = TRUE),
        mean_h_nw_case   = mean(mat[, "nw_h_case"], na.rm = TRUE),
        mean_h_ll_case   = mean(mat[, "ll_h_case"], na.rm = TRUE),
        
        mean_ase_nw_cv   = mean(mat[, "nw_ase_cv"], na.rm = TRUE),
        mean_ase_ll_cv   = mean(mat[, "ll_ase_cv"], na.rm = TRUE),
        sd_ase_nw_cv     = sd(mat[, "nw_ase_cv"], na.rm = TRUE),
        sd_ase_ll_cv     = sd(mat[, "ll_ase_cv"], na.rm = TRUE),
        mean_h_nw_cv     = mean(mat[, "nw_h_cv"], na.rm = TRUE),
        mean_h_ll_cv     = mean(mat[, "ll_h_cv"], na.rm = TRUE),
        
        
        # MGCV using the TRUE correlation matrix
        mean_ase_nw_mgcv = mean(mat[, "nw_ase_mgcv"], na.rm = TRUE),
        mean_ase_ll_mgcv = mean(mat[, "ll_ase_mgcv"], na.rm = TRUE),
        sd_ase_nw_mgcv   = sd(mat[, "nw_ase_mgcv"],   na.rm = TRUE),
        sd_ase_ll_mgcv   = sd(mat[, "ll_ase_mgcv"],   na.rm = TRUE),
        mean_h_nw_mgcv = mean(mat[, "nw_h_mgcv"], na.rm = TRUE),
        mean_h_ll_mgcv = mean(mat[, "ll_h_mgcv"], na.rm = TRUE),
        mean_h_nw_mgcv_true = mean(mat[, "nw_h_mgcv_true"], na.rm = TRUE),
        mean_h_ll_mgcv_true = mean(mat[, "ll_h_mgcv_true"], na.rm = TRUE),
        
        # Error in estimated correlation matrix
        mean_R_error_nw =  mean(mat[, "nw_R_error"], na.rm = TRUE),
        mean_R_error_ll = mean(mat[, "ll_R_error"], na.rm = TRUE),
        sd_R_error_nw = sd(mat[, "nw_R_error"], na.rm = TRUE),
        sd_R_error_ll = sd(mat[, "ll_R_error"], na.rm = TRUE),
        
        med_sigma2_hat_nw = median(mat[, "nw_sigma2_hat"], na.rm = TRUE),
        med_sigma2_hat_ll = median(mat[, "ll_sigma2_hat"], na.rm = TRUE)
        
      )
      
      if (do_gcv) {
        results[[key]]$mean_ase_nw_gcv <- mean(mat[, "nw_ase_gcv"])
        results[[key]]$mean_ase_ll_gcv <- mean(mat[, "ll_ase_gcv"])
        results[[key]]$sd_ase_nw_gcv   <- sd(mat[, "nw_ase_gcv"])
        results[[key]]$sd_ase_ll_gcv   <- sd(mat[, "ll_ase_gcv"])
        results[[key]]$mean_h_nw_gcv   <- mean(mat[, "nw_h_gcv"])
        results[[key]]$mean_h_ll_gcv   <- mean(mat[, "ll_h_gcv"])
      }
    }
  }
  
  future::plan(future::sequential())
  return(results)
}


################################################################################
# Print Tables
################################################################################

make_table <- function(results, print_h = TRUE) {
  all_n     <- sort(unique(sapply(results, `[[`, "n")))
  all_alpha <- sort(unique(sapply(results, `[[`, "alpha")))
  has_gcv   <- !is.null(results[[1]]$mean_ase_nw_gcv)
  
  fmt   <- function(m, s) sprintf("%.5f (%.5f)", m, s)
  fmt_h <- function(h)    sprintf("%.5f", h)
  
  rows <- list()
  
  for (ai in seq_along(all_alpha)) {
    for (n in all_n) {
      key <- paste0("n", n, "_a", ai)
      r <- results[[key]]
      if (is.null(r)) next
      
      row <- data.frame(
        alpha   = all_alpha[ai],
        n       = n,
        NW_CV   = fmt(r$mean_ase_nw_cv,   r$sd_ase_nw_cv),
        NW_MGCV = fmt(r$mean_ase_nw_mgcv, r$sd_ase_nw_mgcv),
        NW_CASE = fmt(r$mean_ase_nw_case, r$sd_ase_nw_case),
        LL_CV   = fmt(r$mean_ase_ll_cv,   r$sd_ase_ll_cv),
        LL_MGCV = fmt(r$mean_ase_ll_mgcv, r$sd_ase_ll_mgcv),
        LL_CASE = fmt(r$mean_ase_ll_case, r$sd_ase_ll_case),
        stringsAsFactors = FALSE
      )
      
      if (has_gcv) {
        row$NW_GCV <- fmt(r$mean_ase_nw_gcv, r$sd_ase_nw_gcv)
        row$LL_GCV <- fmt(r$mean_ase_ll_gcv, r$sd_ase_ll_gcv)
      }
      
      if (print_h) {
        row$NW_h_CV   <- fmt_h(r$mean_h_nw_cv)
        row$NW_h_MGCV <- fmt_h(r$mean_h_nw_mgcv)
        row$NW_h_CASE <- fmt_h(r$mean_h_nw_case)
        row$LL_h_CV   <- fmt_h(r$mean_h_ll_cv)
        row$LL_h_MGCV <- fmt_h(r$mean_h_ll_mgcv)
        row$LL_h_CASE <- fmt_h(r$mean_h_ll_case)
        
        if (has_gcv) {
          row$NW_h_GCV <- fmt_h(r$mean_h_nw_gcv)
          row$LL_h_GCV <- fmt_h(r$mean_h_ll_gcv)
        }
      }
      
      rows <- c(rows, list(row))
    }
  }
  
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  
  nw_ase <- c("NW_CV", if (has_gcv) "NW_GCV", "NW_MGCV", "NW_CASE")
  ll_ase <- c("LL_CV", if (has_gcv) "LL_GCV", "LL_MGCV", "LL_CASE")
  
  if (print_h) {
    nw_h <- c("NW_h_CV", if (has_gcv) "NW_h_GCV", "NW_h_MGCV", "NW_h_CASE")
    ll_h <- c("LL_h_CV", if (has_gcv) "LL_h_GCV", "LL_h_MGCV", "LL_h_CASE")
    col_order <- c("alpha", "n", nw_ase, nw_h, ll_ase, ll_h)
  } else {
    col_order <- c("alpha", "n", nw_ase, ll_ase)
  }
  
  tab[, col_order]
}

