################################################################################
# Simulation study: Linear-spherical regression under spatial dependence
################################################################################
rm(list=ls())
if (!requireNamespace("DirStatsOld", quietly = TRUE)) {
  install.packages("DirStatsOld_0.1.5.tar.gz", repos = NULL, type = "source")
}
library(DirStatsOld)
library(MASS)
library(foreach)
library(progressr)
library(future)
library(doRNG)
library(doFuture)

################################################################################
##------ Spherical distances & Covariance ------------------------------------
geodesic_dist <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  ip <- X_norm %*% t(X_norm)
  ip <- pmax(pmin(ip, 1), -1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
}

m_funs <- list(
  m1 = function(X) X[, 1],
  m2 = function(X) sin(pi * X[, 1]) * X[, 2],
  m3 = function(X, a = 1, b = 1.5) a * sin(2 * pi * X[, 2]) + b * cos(2 * pi * X[, 1])
)

unif_sphere <- function(n, d) {
  X <- matrix(rnorm(n * (d + 1)), nrow = n, ncol = d + 1)
  return(X / sqrt(rowSums(X^2)))
}

################################################################################
# Estimators & Cross-Validation
################################################################################
cv <- function(X, Y, h_grid, p, plot = FALSE) {
  
  cv_error <- sapply(h_grid, function(h){ 
    fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                  h = h, p = p)
    S <- fit$weights
    if (length(dim(S)) == 3L) S <- S[, , 1L]
    
    denom <- 1 - diag(S)
    if (any(denom <= 1e-5)) return(Inf)
    
    res <- ((Y - as.numeric(fit$Yhat)) / denom)^2
    if (any(!is.finite(res))) return(Inf)
    mean(res)
  })
  
  # Bandwidth minimizing CV (NA if CV is not finite for any h)
  if (all(!is.finite(cv_error))) return(list(h = NA_real_, cv_min = NA_real_))
  idx <- which.min(cv_error)
  h <- h_grid[idx]
  cv_min <- cv_error[idx]
  
  if (plot) {
    plot(h_grid, cv_error, type = "b", pch = 19, xlab = "Bandwidth h",
         ylab = "CV(h)")
    points(h,cv_min,pch = 19,cex = 1.5)
    abline(v = h,lty = 2)
  }
  return(list(h = h, cv_min = cv_min))
}

##------ Modified cross-validation (used only to choose the pilot) -----------
# MCV(h) = (1/n) sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
#  N(i) = { j : theta(X_j, X_i) / pi <= ell }
mcv <- function(X, Y, h_grid, p, ell, D = NULL) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi
  min_train <- p * (ncol(X) - 1) + 1
  res <- matrix(NA_real_, nrow = n, ncol = length(h_grid))

  for (i in seq_len(n)) {
    # N(i) always contains i itself (D[i, i] = 0)
    idx_train <- which(D[i, ] > ell)
    if (length(idx_train) < min_train) next

    yhat_i <- loc.directional.linear(x = X[i, , drop = FALSE],
                                     data.dir = X[idx_train, , drop = FALSE],
                                     data.lin = Y[idx_train], h = h_grid,
                                     p = p)$Yhat
    res[i, ] <- (Y[i] - drop(yhat_i))^2
  }
  mcv_error <- colMeans(res)
  mcv_error[!is.finite(mcv_error)] <- Inf

  if (all(!is.finite(mcv_error))) return(list(h = NA_real_, mcv_min = NA_real_))
  idx <- which.min(mcv_error)
  list(h = h_grid[idx], mcv_min = mcv_error[idx])
}

##------ Data-driven pilot bandwidth ------------------------------------------
# h_pilot = c_pilot * h_MCV. MCV removes the correlated neighbours, so h_MCV
# adapts to the curvature of m without collapsing like CV; the factor c_pilot
# oversmooths so the pilot residuals keep the spatial correlation.
pilot_bandwidth <- function(X, Y, h_grid, p, D = NULL, ell_pilot = 0.1,
                            c_pilot = 2) {
  h_mcv <- mcv(X, Y, h_grid, p, ell = ell_pilot, D = D)$h
  if (is.na(h_mcv)) return(NA_real_)
  min(c_pilot * h_mcv, max(h_grid))
}

mgcv <- function(X, Y, h_grid, p, R = NULL, D = NULL, h_pilot = NULL,
                 J = 20, min_pairs = 5, plot = FALSE) {
  n <- nrow(X)
  if (is.null(R)) {
    if (is.null(D)) D <- geodesic_dist(X) / pi
    if (is.null(h_pilot)) stop("Supply h_pilot for MGCV.")
    R <- correlation_matrix(X, Y, h = h_pilot, p = p, J = J,
                            D = D, min_pairs = min_pairs)$R
  }
  mgcv_error <- sapply(h_grid, function(h){ 
    fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                  h = h, p = p)
    S <- fit$weights
    if (length(dim(S)) == 3L) S <- S[, , 1L]
    
    trSR  <- sum(S * t(R))
    denom <- 1 - (1/n) * trSR
    if (!is.finite(denom) || denom <= 0) return(Inf)
    
    res <- ((Y - as.numeric(fit$Yhat)) / denom)^2
    if (any(!is.finite(res))) return(Inf)
    mean(res)
  })
  # Bandwidth minimizing MGCV (NA if MGCV is not finite for any h)
  if (all(!is.finite(mgcv_error))) return(list(h = NA_real_, mgcv_min = NA_real_))
  idx <- which.min(mgcv_error)
  h <- h_grid[idx]
  mgcv_min <- mgcv_error[idx]
  
  if (plot) {
    plot(h_grid, mgcv_error, type = "b", pch = 19, xlab = "Bandwidth h",
         ylab = "MGCV(h)")
    points(h, mgcv_min, pch = 19,cex = 1.5)
    abline(v = h,lty = 2)
  }
  return(list(h = h, mgcv_min = mgcv_min))
}

################################################################################
# Distance Binning & Variogram Estimation
################################################################################
# Classifies pairwise geodesic distances into J logarithmically spaced bins
dist_est <- function(X, D = NULL, J = 20, plot = FALSE) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  dvec <- D[upper.tri(D)]
  pairs <- which(upper.tri(D), arr.ind = TRUE)
  
  keep <- is.finite(dvec) & dvec > 0
  dvec <- dvec[keep]
  pairs <- pairs[keep, , drop = FALSE]
  if (!length(dvec)) stop("No positive pairwise distances.")
  
  # Define the range of distances used for binning
  d_min <- max(quantile(dvec, 0.005, names = FALSE), 0.005)
  d_max <- quantile(dvec, 0.95, names = FALSE)
  
  # Define logarithmically spaced reference distances
  d <- exp(seq(log(d_min), log(d_max), length.out = J))
  
  # Differences between consecutive reference distances.
  diffs <- diff(d)
  # half-width of each bin
  tol <- c(diffs / 2, diffs[length(diffs)] / 2)
  
  # Assign each pairwise distance to a bin
  bin <- rep(NA_integer_, length(dvec))
  for (j in seq_len(J)) {
    # A distance belongs to bin j if it falls within the interval
    # centred at d[j] with half-width tol[j].
    bin[dvec >= (d[j] - tol[j]) & dvec < (d[j] + tol[j])] <- j
  }
  
  S_set <- data.frame(i = pairs[, 1], j = pairs[, 2], dist = dvec, bin = bin)
  S_set <- S_set[!is.na(S_set$bin), ]
  counts <- as.vector(tabulate(S_set$bin, nbins = J))
  
  if (plot) {
    par(mfrow = c(1, 2))
    # Histogram of pairwise distances
    hist(dvec, breaks = 30, xlab = "Normalized distance",
         ylab = "Number of pairs")
    # Bin centres
    abline(v = d, lty = 2)
    
    # Number of pairs per bin
    plot(d, counts, type = "b", pch = 19, log = "x", xlab = "Bin centre",
      ylab = "Number of pairs", main = "Pairs per bin")
    par(mfrow = c(1, 1))
  }
  
  list(d = d, tol = tol, J = J, S_set = S_set, counts = counts, dvec = dvec)
}


# Computes classical empirical semivariogram from binned regression residuals
empirical_variogram <- function(eps, g) {
  # Pairwise distances and bin assignments
  S_set <- g$S_set
  J <- g$J
  
  gamma_val <- rep(NA_real_, J)
  n_pairs   <- rep(0L, J)
  
  # Compute semivariances for all pairs
  if (!is.null(S_set) && nrow(S_set) > 0) {
    # Pairwise semivariance: 0.5 * (e_i - e_j)^2
    sq_diff <- 0.5 * (eps[S_set$i] - eps[S_set$j])^2
    
    # Sum and count semivariances within each bin
    ag_sum <- tapply(sq_diff, S_set$bin, sum)
    ag_count <- tapply(sq_diff, S_set$bin, length)
    
    
    bins_present <- as.integer(names(ag_sum))
    # Average semivariance within each bin
    gamma_val[bins_present] <- as.numeric(ag_sum) / as.numeric(ag_count)
    # Number of pairs per bin
    n_pairs[bins_present] <- as.integer(ag_count)
  }
  
  list(d = g$d, gamma = gamma_val, n = n_pairs)
}

# Estimates spatial correlation matrix R by fitting an exponential variogram via NLS
correlation_matrix <- function(X, Y, h, p, J = 20, D = NULL, min_pairs = 5) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  fit_pilot <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  yhat <- as.numeric(fit_pilot$Yhat)
  eps <- Y - yhat
  
  #Distance bins and empirical variogram
  g <- dist_est(X, D = D, J = J)
  emp_var <- empirical_variogram(eps, g)
  
  check <- is.finite(emp_var$gamma) & emp_var$d > 0 & emp_var$n >= min_pairs
  if (sum(check) < 3) stop("Too few usable bins for NLS estimation.")
  
  d_j <- emp_var$d[check]
  emp_j <- emp_var$gamma[check]
  
  # Define the exponential variogram objective
  obj_nls <- function(par) {
    # Define the exponential variogram objective
    s2 <- exp(par[1]) #sigma2
    a <- exp(par[2])  #alpga
    # Log-scale parametrization ensures positive parameters
    gamma_teorico <- s2 * (1 - exp(- d_j / a))
    # Theoretical exponential variogram
    val <- sum((emp_j - gamma_teorico)^2)
    if (!is.finite(val)) return(1e10)
    val
  }
  
  # Initial values optimization
  init_s2 <- log(pmax(quantile(emp_j, 0.9, na.rm = TRUE, names = FALSE), 1e-3))
  init_a <- log(pmax(median(d_j) / 2, 0.1))
  
  # Fit the exponential variogram
  opt <- tryCatch(
    optim(par = c(init_s2, init_a), fn = obj_nls, method = "L-BFGS-B",
          lower = c(log(1e-4), log(1e-4)), upper = c(log(100), log(10))),
    error = function(e) {
      optim(par = c(init_s2, init_a), fn = obj_nls, method = "Nelder-Mead")
    }
  )
  
  sigma2_hat <- exp(opt$par[1])
  alpha_hat  <- exp(opt$par[2])
  
  R <- exp(- D / alpha_hat)
  diag(R) <- 1
  
  list(R = R, alpha = alpha_hat, sigma2 = sigma2_hat, variogram = emp_var)
}

################################################################################
# Simulation Core Loop
################################################################################

# h_pilot = NULL: data-driven pilot (c_pilot * h_MCV with ell_pilot), computed
# separately for p = 0 and p = 1. A positive number fixes the pilot instead
# (useful for sensitivity analyses).
one_rep <- function(n, alpha, sigma2, m_fun, h_grid, h_pilot = NULL,
                    ell_pilot = 0.1, c_pilot = 2, d = 2) {

  if (!is.null(h_pilot) && (!is.numeric(h_pilot) || length(h_pilot) != 1L ||
                            !is.finite(h_pilot) || h_pilot <= 0)) {
    stop("h_pilot must be NULL or a single positive number.")
  }
 
  X <- unif_sphere(n, d)               # (n x 3) points on S^2
  m_vals <- m_fun(X)                   # true regression values at X
  D <- geodesic_dist(X) / pi
  
  Sigma <- sigma2 * exp(-D/ alpha) 
  
  if (inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    stop(sprintf(" Sigma is not positive definite for n = %d and alpha = %.3f", 
                 n, alpha))
  }
  
  eps <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
  Y <- m_vals + eps
  
  out <- c()
  
  for (p in c(0, 1)) {
    tag <- if (p == 0) "nw" else "ll"
    
    ## --- ASE : true ase over the whole grid (one fit for all h) ---
    mhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                   h = h_grid, p = p)$Yhat
    ase_vals <- colMeans((mhat - m_vals)^2)
    ase_vals[!is.finite(ase_vals)] <- Inf
    
    idx_ase <- if (all(!is.finite(ase_vals))) {NA} else {which.min(ase_vals)}
    out[paste0(tag, "_h_ase")] <- h_grid[idx_ase]
    out[paste0(tag, "_mse_ase")] <- ase_vals[idx_ase]
    
    ## --- CV --- (ASE of the selected h read from ase_vals; NA if CV failed)
    h_cv <- cv(X, Y, h_grid, p)$h
    out[paste0(tag, "_h_cv")] <- h_cv
    out[paste0(tag, "_mse_cv")] <- if (is.na(h_cv)) NA_real_ else
      ase_vals[match(h_cv, h_grid)]

    # MGCV
    h_pil <- if (is.null(h_pilot)) {
      pilot_bandwidth(X, Y, h_grid, p, D = D, ell_pilot = ell_pilot,
                      c_pilot = c_pilot)
    } else h_pilot
    out[paste0(tag, "_h_pilot")] <- h_pil

    cm <- if (is.na(h_pil)) NULL else tryCatch(
      correlation_matrix(X, Y, h = h_pil, p = p, J = 20, D = D, min_pairs = 5),
      error = function(e) NULL
    )
    
    Rhat <- if (is.null(cm)) NULL else cm$R
    out[paste0(tag, "_alpha_hat")] <- if (is.null(cm)) NA_real_ else cm$alpha
    
    if (is.null(Rhat)) {
      out[paste0(tag, "_h_mgcv")]   <- NA_real_
      out[paste0(tag, "_mse_mgcv")] <- NA_real_
    } else {
      h_mgcv <- mgcv(X, Y, h_grid, p, R = Rhat)$h
      
      out[paste0(tag, "_h_mgcv")] <- h_mgcv
      
      if (is.na(h_mgcv)) {
        out[paste0(tag, "_mse_mgcv")] <- NA_real_
      } else {
        out[paste0(tag, "_mse_mgcv")] <-
          ase_vals[match(h_mgcv, h_grid)]
      }
    }
  }
  
  return(list(results = out, X = X, Y = Y))
}

run_simulation <- function(MC, n_values, alpha_vals, sigma2, m_idx,
                           h_grid, h_pilot = NULL, ell_pilot = 0.1,
                           c_pilot = 2, d = 2, cores = 1) {
  
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
                        .packages = c("DirStatsOld", "MASS"),
                        .export = c("one_rep", "cv", "mgcv", "mcv",
                                    "pilot_bandwidth",
                                    "correlation_matrix", "dist_est",
                                    "empirical_variogram", "geodesic_dist",
                                    "unif_sphere",
                                    "m_funs", "h_grid", "sigma2",
                                    "d", "h_pilot", "ell_pilot", "c_pilot",
                                    "m_idx")) %dorng% {
                                      prog()
                                      one_rep(n = n, alpha = alpha,
                                              sigma2 = sigma2,
                                              m_fun   = m_funs[[m_idx]],
                                              h_grid  = h_grid,
                                              h_pilot = h_pilot,
                                              ell_pilot = ell_pilot,
                                              c_pilot = c_pilot, d = d)
                                    }
      })
      
      mat <- do.call(rbind, lapply(reps, `[[`, "results"))
      samples <- lapply(reps, function(r) list(X = r$X, Y = r$Y))
      
      key <- paste0("n", n, "_a", alpha_idx)
      results[[key]] <- list(
        n = n,
        alpha = alpha,
        h_pilot = h_pilot,
        ell_pilot = ell_pilot,
        c_pilot = c_pilot,
        mat = mat,
        samples = samples,
         
        mean_mse_nw_cv = mean(mat[, "nw_mse_cv"], na.rm = TRUE),
        mean_mse_nw_mgcv = mean(mat[, "nw_mse_mgcv"], na.rm = TRUE),
        mean_mse_nw_ase = mean(mat[, "nw_mse_ase"], na.rm = TRUE),
        mean_mse_ll_cv  = mean(mat[, "ll_mse_cv"], na.rm = TRUE),
        mean_mse_ll_mgcv = mean(mat[, "ll_mse_mgcv"], na.rm = TRUE),
        mean_mse_ll_ase = mean(mat[, "ll_mse_ase"], na.rm = TRUE),
         
        median_mse_nw_cv  = median(mat[, "nw_mse_cv"], na.rm = TRUE),
        median_mse_nw_mgcv  = median(mat[, "nw_mse_mgcv"], na.rm = TRUE),
        median_mse_nw_ase = median(mat[, "nw_mse_ase"], na.rm = TRUE),
        median_mse_ll_cv  = median(mat[, "ll_mse_cv"], na.rm = TRUE),
        median_mse_ll_mgcv  = median(mat[, "ll_mse_mgcv"], na.rm = TRUE),
        median_mse_ll_ase = median(mat[, "ll_mse_ase"], na.rm = TRUE),
         
        sd_mse_nw_cv = sd(mat[, "nw_mse_cv"], na.rm = TRUE),
        sd_mse_nw_mgcv = sd(mat[, "nw_mse_mgcv"], na.rm = TRUE),
        sd_mse_nw_ase = sd(mat[, "nw_mse_ase"], na.rm = TRUE),
        sd_mse_ll_cv = sd(mat[, "ll_mse_cv"], na.rm = TRUE),
        sd_mse_ll_mgcv = sd(mat[, "ll_mse_mgcv"], na.rm = TRUE),
        sd_mse_ll_ase = sd(mat[, "ll_mse_ase"], na.rm = TRUE),
                               
        mean_h_nw_cv = mean(mat[, "nw_h_cv"], na.rm = TRUE),
        mean_h_nw_mgcv = mean(mat[, "nw_h_mgcv"], na.rm = TRUE),
        mean_h_nw_ase = mean(mat[, "nw_h_ase"], na.rm = TRUE),
        mean_h_ll_cv = mean(mat[, "ll_h_cv"], na.rm = TRUE),
        mean_h_ll_mgcv = mean(mat[, "ll_h_mgcv"], na.rm = TRUE),
        mean_h_ll_ase = mean(mat[, "ll_h_ase"], na.rm = TRUE),
         
        median_h_nw_cv = median(mat[, "nw_h_cv"], na.rm = TRUE),
        median_h_nw_mgcv = median(mat[, "nw_h_mgcv"], na.rm = TRUE),
        median_h_nw_ase = median(mat[, "nw_h_ase"], na.rm = TRUE),
        median_h_ll_cv = median(mat[, "ll_h_cv"], na.rm = TRUE),
        median_h_ll_mgcv = median(mat[, "ll_h_mgcv"], na.rm = TRUE),
        median_h_ll_ase = median(mat[, "ll_h_ase"], na.rm = TRUE),
                               
    
        n_fail_nw_mgcv = sum(is.na(mat[, "nw_mse_mgcv"])),
        n_fail_ll_mgcv = sum(is.na(mat[, "ll_mse_mgcv"])),
        med_alpha_hat_nw = median(mat[, "nw_alpha_hat"], na.rm = TRUE),
        med_alpha_hat_ll = median(mat[, "ll_alpha_hat"], na.rm = TRUE),
        med_h_pilot_nw = median(mat[, "nw_h_pilot"], na.rm = TRUE),
        med_h_pilot_ll = median(mat[, "ll_h_pilot"], na.rm = TRUE)
      )
    }
  }
  
  future::plan(future::sequential())
  return(results)
}

make_table <- function(results, med = FALSE, print_h = TRUE) {
  
  stat <- if (med) "median" else "mean"
  fmt <- function(x, s) {sprintf("%.5f (%.5f)", x, s)}
  fmt_h <- function(h) {sprintf("%.5f", h)}
  
  rows <- list()
  
  for (r in results) {
    
    row <- data.frame(alpha = r$alpha, n = r$n,
      NW_CV = fmt(r[[paste0(stat, "_mse_nw_cv")]], r$sd_mse_nw_cv),
      NW_MGCV = fmt(r[[paste0(stat, "_mse_nw_mgcv")]], r$sd_mse_nw_mgcv),
      NW_ASE = fmt(r[[paste0(stat, "_mse_nw_ase")]], r$sd_mse_nw_ase),
      
      LL_CV = fmt(r[[paste0(stat, "_mse_ll_cv")]], r$sd_mse_ll_cv),
      LL_MGCV = fmt(r[[paste0(stat, "_mse_ll_mgcv")]], r$sd_mse_ll_mgcv),
      LL_ASE = fmt(r[[paste0(stat, "_mse_ll_ase")]], r$sd_mse_ll_ase),
      stringsAsFactors = FALSE)
    
    if (print_h) {
      
      row$NW_h_CV <- fmt_h(r[[paste0(stat, "_h_nw_cv")]])
      row$NW_h_MGCV <- fmt_h(r[[paste0(stat, "_h_nw_mgcv")]])
      row$NW_h_ASE <- fmt_h(r[[paste0(stat, "_h_nw_ase")]])
      row$LL_h_CV <- fmt_h(r[[paste0(stat, "_h_ll_cv")]])
      row$LL_h_MGCV <- fmt_h(r[[paste0(stat, "_h_ll_mgcv")]])
      row$LL_h_ASE <- fmt_h(r[[paste0(stat, "_h_ll_ase")]])
    }
    rows <- c(rows, list(row))
  }
  
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  
  # Error / ASE columns
  nw_mse <- c("NW_CV", "NW_MGCV", "NW_ASE")
  ll_mse <- c( "LL_CV","LL_MGCV","LL_ASE")

  if (print_h) {
    nw_h <- c("NW_h_CV", "NW_h_MGCV","NW_h_ASE")
    ll_h <- c("LL_h_CV","LL_h_MGCV","LL_h_ASE" )
    col_order <- c("alpha","n",nw_mse,nw_h, ll_mse,ll_h)
  } else {
    col_order <- c( "alpha","n",nw_mse,ll_mse)
  }
  
  tab <- tab[, col_order, drop = FALSE]
  # Order by alpha and then sample size
  tab <- tab[order(tab$alpha, tab$n), ]
  
  rownames(tab) <- NULL
  
  tab
}

alpha_table <- function(results) {
  tab <- do.call(rbind, lapply(results, function(r)
    data.frame(alpha = r$alpha, n = r$n,
               alpha_hat_NW = round(r$med_alpha_hat_nw, 3),
               alpha_hat_LL = round(r$med_alpha_hat_ll, 3),
               h_pilot_NW = round(r$med_h_pilot_nw, 3),
               h_pilot_LL = round(r$med_h_pilot_ll, 3))))
  rownames(tab) <- NULL
  tab[order(tab$alpha, tab$n), ]
}