################################################################################
# Simulation study: Linear-spherical regression under spatial dependence
# Bandwidth selection: standard CV vs. modified GCV (MGCV)
################################################################################

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
##------ Spherical distances -------------------------------------------------
geodesic_dist <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  ip <- X_norm %*% t(X_norm)
  ip <- pmax(pmin(ip, 1), -1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
}

##------ Covariance structure for the error term -----------------------------
# Exponential: Sigma_ij = sigma2 * exp(-d_geo(X_i, X_j) * n^(1/d) / alpha)
build_Sigma_noN <- function(X, alpha, sigma2) {
  D <- geodesic_dist(X) / pi     # normalized to [0, 1]
  sigma2 * exp(-D / alpha)
}


build_Sigma<- function(X, alpha, sigma2) {
  D <- geodesic_dist(X) / pi     # normalized to [0, 1]
  n <- nrow(X)
  q <- ncol(X) - 1               # q = 2 for S^2 embedded in R^3
  sigma2 * exp(- (n^(1/q) * D) / alpha)
}

##------ Regression functions on S^2 -----------------------------------------
m_funs <- list(
  m1 = function(X) X[, 1],
  m2 = function(X) sin(pi * X[, 1]) * X[, 2],
  m3 = function(X, a = 1, b = 1.5) a * sin(2 * pi * X[, 2]) + b * 
    cos(2 * pi * X[, 1])
)

##------ Uniform sample on S^d -----------------------------------------------
unif_sphere <- function(n, d) {
  X <- matrix(rnorm(n * (d + 1)), nrow = n, ncol = d + 1)
  return(X / sqrt(rowSums(X^2)))
}

################################################################################
# Bandwidth selection
################################################################################

##------Projected local estimator fit -----------------------------------------
lpe <- function(eval_pts, dir_data, lin_data, h, p) {
  # x: evaluation points
  # data.dir: directional training data
  # data.lin: linear training response
  # h: bandwith parameter
  # p: polynomial order
  fit <- loc.directional.linear(x = eval_pts, data.dir = dir_data,
                                data.lin = lin_data, h = h, p = p)
  return(fit$Yhat)
}


##------Cross-Validation (leave-one-out)----------------------------------------
# 1) CV(h) = sum_i [ Y_i - m_hat_{h,p,-i}(X_i) ]^2

cv_loo_old <- function(X, Y, h, p) {
  n <- nrow(X)
  res <- numeric(n)
  for (i in seq_len(n)) {
    idx_train <- setdiff(seq_len(n), i)
    yhat_i <- lpe( eval_pts = X[i, , drop = FALSE],
                   dir_data = X[idx_train, , drop = FALSE],
                   lin_data = Y[idx_train],
                   h = h, p = p
    )
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res) 
}

cv_loo <- function(X, Y, h, p) {
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  S <- fit$weights
  if (length(dim(S)) == 3L) S <- S[, , 1L]
  Yhat <- as.numeric(fit$Yhat)
  
  denom <- 1 - diag(S)
  if (any(denom <= 1e-5)) return(Inf)
  
  mean(((Y - Yhat) / denom)^2, na.rm = TRUE)
}

##------Generalized Cross-Validation ----------------------------------------
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

##------ Modified  cross-validation -----------------------------------------
# 2) MCV(h) = sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
#  N(i) = { j : theta(X_j, X_i) <= ell }
#  For S^2 (sphere), N(i) is a spherical cap of geodesic radius ell around X_i.
#  NOTE: kept for reference, not used in the current simulation.

mcv_loo <- function(X, Y, h, p, ell, D = NULL) {
  n <- nrow(X)
  if (is.null(D)) {
    D <- geodesic_dist(X) / pi
  } # (n x n) geodesic distances
  res <- numeric(n)
  
  for (i in seq_len(n)) {
    in_nbhd <- which(D[i, ] <= ell)
    idx_train <- setdiff(seq_len(n), in_nbhd)
    
    # Skip if too few points remain after removing N(i) to fit a 
    #degree-p polynomial
    if (length(idx_train) < (p + 2)) {
      res[i] <- NA
      next
    }
    
    yhat_i <- lpe( eval_pts = X[i, , drop = FALSE],
                   dir_data = X[idx_train, , drop = FALSE],
                   lin_data = Y[idx_train],
                   h = h, p= p
    )
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res, na.rm = TRUE)
}

##------ Modified generalized cross-validation --------------------------------
# 3) MGCV(h) = (1/n) sum_i [ (Y_i - m_hat_{h,p}(X_i)) / (1 - tr(S_p R)/n) ]^2
#  S_p = smoother matrix (row i = S_{X_i;p}), R = error correlation matrix.
#  No observation is left out: the factor 1 - tr(S_p R)/n plays that role.

mgcv <- function(X, Y, h, p, R = NULL, D = NULL, h_pilot = NULL,
                 J = 20, tol = 0.05, min_pairs = 30) {
  n <- nrow(X)
  
  if (is.null(R)) {
    if (is.null(D)) D <- geodesic_dist(X) / pi
    if (is.null(h_pilot))
      stop("Supply h_pilot: using the h being scored as the pilot makes R change with h.")
    R <- correlation_matrix(X, Y, h = h_pilot, p = p, J = J, tol = tol,
                            D = D, min_pairs = min_pairs)$R
  }
  
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

# dist_est <- function(X, D = NULL, J = 20, tol = 0.005) {
# 
#   if (is.null(D)) D <- geodesic_dist(X) / pi
# 
#   ## --- pairwise geodesic distances, i < j ---------------------------------
#   dvec <- D[upper.tri(D)]
#   pairs <- which(upper.tri(D), arr.ind = TRUE)
# 
#   keep <- is.finite(dvec) & dvec > 0
#   dvec <- dvec[keep]
#   pairs <- pairs[keep, , drop = FALSE]
#   if (!length(dvec)) stop("No positive pairwise distances.")
# 
#   ## --- grid: J interior quantiles of the distance sample ------------------
#   probs <- seq_len(J) / (J + 1)
#   d <- as.numeric(quantile(dvec, probs = probs, names = FALSE))
# 
#   ## --- tolerancia: no puede exceder medio espaciado, o los bins se solapan
#   tol <- min(tol, min(diff(d)) / 2)
# 
#   ## --- bin label: which d_j each pair belongs to, NA if none --------------
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

empirical_variogram <- function(eps, dist) {
  ## eps : residuals from the pilot fit, length n
  ## dist : output of dist_est()
  
  v <- (eps[dist$S_set$i] - eps[dist$S_set$j])^2
  rs <- rowsum(v, dist$S_set$bin)              # one row per non-empty bin
  
  num <- numeric(dist$J)
  num[as.integer(rownames(rs))] <- rs[, 1]   # place each sum at its bin index
  
  data.frame(d = dist$d, gamma = num / (2 * dist$counts), n = dist$counts)
}

# 
# correlation_matrix <- function(X, Y, h, p, J = 20, tol = 0.05,
#                                D = NULL, min_pairs = 30) {
# 
#   n <- nrow(X)
#   q <- ncol(X) - 1
#   if (is.null(D)) D <- geodesic_dist(X) / pi
# 
#   ## --- pilot fit and residuals --------------------------------------------
#   yhat <- lpe(eval_pts = X, dir_data = X, lin_data = Y, h = h, p = p)
#   eps <- Y - as.numeric(yhat)
#   sigma2 <- mean(eps^2)
# 
#   ## --- empirical semivariogram on the quantile grid -----------------------
#   g <- dist_est(X, D = D, J = J, tol = tol)
#   emp_var <- empirical_variogram(eps, g)                         # (10)
# 
#   # gamma <- pmin(emp_var$gamma, 0.95 * sigma2)
#   ## --- alpha_j, method of moments -----------------------------------------
#   ## rho(d) = exp(-d/alpha)  =>  alpha_j = d_j / [ln(sigma2) - ln(sigma2 - gamma(d_j))]
#   check <- is.finite(emp_var$gamma) & emp_var$gamma < sigma2 &
#     emp_var$d > 0 & emp_var$n >= min_pairs
#   if (!any(check))
#     stop("Cannot estimate alpha: no bin has gamma-hat below sigma2-hat.")
# 
#   # Invert rho(d) = exp(-d / alpha) -> alpha = d / (ln(sigma2) - ln(sigma2 - gamma))
#   alpha_j <- rep(NA_real_, J)
#   alpha_j[check] <- (emp_var$d[check] * n^(1/q)) /
#     (log(sigma2) - log(sigma2 - emp_var$gamma[check]))
# 
#   alpha <- mean(alpha_j[check], na.rm = TRUE)
# 
#   ## --- Estimation of R --------
#   R <- exp(- D * n^(1/q) / alpha)
#   diag(R) <- 1
# 
#   list(R = R, alpha = alpha, alpha_j = alpha_j, sigma2 = sigma2,
#        variogram = emp_var, used = check)
# }
# 



# correlation_matrixNLS <- function(X, Y, h, p, J = 20, tol = 0.05,
#                                    D = NULL, min_pairs = 30) {
#   
#   n <- nrow(X)
#   q <- ncol(X) - 1 
#   if (is.null(D)) D <- geodesic_dist(X) / pi
#   
#   ## --- 1. Ajuste piloto y residuos -----------------------------------
#   yhat <- lpe(eval_pts = X, dir_data = X, lin_data = Y, h = h, p = p)
#   eps <- Y - as.numeric(yhat)
#   sigma2 <- mean(eps^2)                                          
#   
#   ## --- 2. Variograma empírico con dist_est original ------------------
#   g <- dist_est(X, D = D, J = J, tol = tol)
#   emp_var <- empirical_variogram(eps, g)                         
#   
#   # Bins válidos (ya NO requerimos que gamma < sigma2)
#   check <- is.finite(emp_var$gamma) & emp_var$d > 0 & emp_var$n >= min_pairs
#   if (sum(check) < 3)
#     stop("Cannot estimate alpha: too few usable bins.")
#   
#   d_j <- emp_var$d[check]
#   emp_j <- emp_var$gamma[check]
#   counts_j <- emp_var$n[check] # Número de pares por bin como peso
#   
#   ## --- 3. Estimación por Mínimos Cuadrados No Lineales (NLS) ---------
#   # Modelo teórico: gamma(d) = sigma2 * (1 - exp(-n^(1/q) * d / alpha))
#   # Se optimiza sobre log(alpha) para garantizar alpha > 0
#   obj_nls <- function(log_alpha) {
#     a <- exp(log_alpha)
#     gamma_teorico <- sigma2 * (1 - exp(- (n^(1/q) * d_j) / a))
#     
#     # Suma de residuos cuadráticos ponderados por pares por bin
#     sum(counts_j * (emp_j - gamma_teorico)^2)
#   }
#   
#   opt <- optimize(obj_nls, interval = log(c(1e-4, 100)))
#   alpha <- exp(opt$minimum)
#   
#   ## --- 4. Estimación de R --------------------------------------------
#   R <- exp(- D * n^(1/q) / alpha)
#   diag(R) <- 1    
#   
#   list(R = R, alpha = alpha, sigma2 = sigma2,
#        variogram = emp_var, used = check)
# }




dist_est <- function(X, D = NULL, J = 20, d_max) {
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  dvec  <- D[upper.tri(D)]
  pairs <- which(upper.tri(D), arr.ind = TRUE)
  keep  <- is.finite(dvec) & dvec > 0 & dvec <= d_max
  dvec  <- dvec[keep]
  pairs <- pairs[keep, , drop = FALSE]
  if (!length(dvec)) stop("No pairs below d_max.")
  
  ## J lags equiespaciados en (0, d_max]; bins contiguos sin solape
  d   <- seq_len(J) * d_max / J
  tol <- d_max / (2 * J)
  
  bin <- rep(NA_integer_, length(dvec))
  for (j in seq_len(J)) bin[dvec >= d[j] - tol & dvec < d[j] + tol] <- j
  
  S_set <- data.frame(i = pairs[, 1], j = pairs[, 2], dist = dvec, bin = bin)
  S_set <- S_set[!is.na(S_set$bin), ]
  
  list(d = d, tol = tol, J = J, S_set = S_set,
       counts = tabulate(S_set$bin, nbins = J))
}

correlation_matrix <- function(X, Y, h, p, J = 20, D = NULL,
                               min_pairs = 30, smax = 6, ...) {
  n <- nrow(X); q <- ncol(X) - 1
  if (is.null(D)) D <- geodesic_dist(X) / pi
  
  ## Ajuste piloto, residuos y sigma2 con corrección de grados de libertad
  fit <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, h = h, p = p)
  eps <- Y - as.numeric(fit$Yhat)
  S <- fit$weights; if (length(dim(S)) == 3L) S <- S[, , 1L]
  sigma2 <- sum(eps^2) / max(1, n - sum(diag(S)))
  
  ## Variograma solo en la zona informativa: sqrt(n) * d <= smax
  g <- dist_est(X, D = D, J = J, d_max = smax / n^(1/q))
  emp_var <- empirical_variogram(eps, g)
  
  ok <- is.finite(emp_var$gamma) & emp_var$n >= min_pairs
  if (sum(ok) < 3) stop("Cannot estimate alpha: too few usable bins.")
  d_j <- emp_var$d[ok]; emp_j <- emp_var$gamma[ok]; w_j <- emp_var$n[ok]
  
  ## NLS con el sill fijado a sigma2
  obj <- function(la) {
    a <- exp(la)
    sum(w_j * (emp_j - sigma2 * (1 - exp(-n^(1/q) * d_j / a)))^2)
  }
  alpha <- exp(optimize(obj, interval = log(c(1e-3, 50)))$minimum)
  
  R <- exp(-D * n^(1/q) / alpha); diag(R) <- 1
  list(R = R, alpha = alpha, sigma2 = sigma2, variogram = emp_var)
}
################################################################################
# Simulation functions
################################################################################

## ------ Single Monte Carlo iteration ----------------------------------------
## Returns a named vector with ase (and selected h) for:
##   NW  (p=0): CV, MGCV, CASE
##   LL  (p=1): CV, MGCV, CASE
one_rep <- function(n, alpha, sigma2, m_fun, h_grid, d = 2,
                    h_pilot = NULL, do_gcv = FALSE) {
  
  ## Data generation
  X <- unif_sphere(n, d)               # (n x 3) points on S^2
  m_vals <- m_fun(X)                   # true regression values at X
  
  Sigma  <- build_Sigma(X, alpha, sigma2)
  eps <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
  Y <- m_vals + eps
  
  ## Compute standardized distances
  D <- geodesic_dist(X) / pi
  
  if (is.null(h_pilot)) h_pilot <- bw_dir_rot(X)
  
  ##  error for a given  bandwidth and p at the data points
  ase <- function(h, p) {
    yhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                   h = h, p = p)$Yhat
    
    mean((yhat - m_vals)^2)
  }
  
  ## Bandwidth selection for each estimator (p = 0 NW, p = 1 LL)
  out <- c()
  
  for (p in c(0, 1)) {
    tag <- if (p == 0) "nw" else "ll"
    
    ## --- CASE : minimize true ase over the grid ---
    ase_vals <- sapply(h_grid, function(h) ase(h, p))
    idx_case <- which.min(ase_vals)
    out[paste0(tag, "_h_case")] <- h_grid[idx_case]
    out[paste0(tag, "_ase_case")] <- ase_vals[idx_case]
    
    ## --- CV ---
    cv_vals <- sapply(h_grid, function(h) cv_loo(X, Y, h, p))
    idx_cv <- which.min(cv_vals)
    out[paste0(tag, "_h_cv")] <- h_grid[idx_cv]
    out[paste0(tag, "_ase_cv")] <- ase_vals[idx_cv]
    
    ## --- GCV (uncorrelated) ---
    if (do_gcv) {
      gcv_vals <- sapply(h_grid, function(h) gcv(X, Y, h, p))
      idx_gcv <- which.min(gcv_vals)
      out[paste0(tag, "_h_gcv")] <- h_grid[idx_gcv]
      out[paste0(tag, "_ase_gcv")] <- ase_vals[idx_gcv]
    }
    
    ## --- MGCV ---
    ## R is estimated ONCE per replicate with the fixed pilot bandwidth:
    ## re-estimating it inside the h_grid sweep would score each h under a
    ## different correlation model.
    
    Rhat <- if (!is.finite(h_pilot)) NULL else tryCatch(
      correlation_matrix(X, Y, h = h_pilot, p = p, J = 20, tol = 0.05,
                         D = D, min_pairs = 30)$R,
      error = function(e) NULL)
    
    if (is.null(Rhat)) {
      out[paste0(tag, "_h_mgcv")]   <- NA_real_
      out[paste0(tag, "_ase_mgcv")] <- NA_real_
    } else {
      mgcv_vals <- sapply(h_grid, function(h)
        tryCatch(mgcv(X, Y, h, p, R = Rhat), error = function(e) Inf))
      if (all(!is.finite(mgcv_vals))) {
        out[paste0(tag, "_h_mgcv")]   <- NA_real_
        out[paste0(tag, "_ase_mgcv")] <- NA_real_
      } else {
        idx_mgcv <- which.min(mgcv_vals)
        out[paste0(tag, "_h_mgcv")]   <- h_grid[idx_mgcv]
        out[paste0(tag, "_ase_mgcv")] <- ase_vals[idx_mgcv]
      }
    }
    
  }
  
  return(list(results = out, X = X, Y = Y))
}


## ------ Simulation over all the MC iterations -------------------------------
run_simulation <- function(MC, n_values, alpha_vals, sigma2, m_idx,
                           h_grid, d = 2, cores = 1,
                           h_pilot = NULL, do_gcv = FALSE) {
  
  ## Register parallel backend
  doFuture::registerDoFuture()
  future::plan(future::multisession(), workers = cores)
  
  ## Progress bar style
  handlers(handler_progress(
    format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
    clear  = FALSE
  ))
  
  results <- list()
  
  for (n in n_values) {
    for (alpha_idx in seq_along(alpha_vals)) {
      alpha <- alpha_vals[alpha_idx]
      
      cat(sprintf("\n--- m%d | n = %d | alpha = %.2f ---\n",
                  m_idx, n, alpha))
      
      progressr::with_progress({
        prog <- progressr::progressor(along = seq_len(MC))
        
        # Run MC replications in parallel
        # Each replication returns list(results = ..., X = ..., Y = ...)
        reps <- foreach(k = seq_len(MC), .inorder = TRUE,
                        .packages = c("DirStatsOld", "MASS"),
                        .export = c("one_rep", "cv_loo", "gcv", "mgcv", 
                                    "correlation_matrix", "dist_est", 
                                    "empirical_variogram", "geodesic_dist",
                                    "build_Sigma", "unif_sphere", "lpe", 
                                    "m_funs", "h_grid", "do_gcv", "sigma2", 
                                    "d", "h_pilot", "m_idx")) %dorng% {
                                      prog()
                                      one_rep(n = n, alpha = alpha, 
                                              sigma2 = sigma2,
                                              m_fun   = m_funs[[m_idx]],
                                              h_grid  = h_grid,
                                              d = d, do_gcv = do_gcv,
                                              h_pilot = h_pilot)
                                    }
      })
      
      # Save results of h and ASE per method into an MC x 12 matrix
      mat <- do.call(rbind, lapply(reps, `[[`, "results"))
      
      # Store the data (X, Y) from each replication as a list of length MC
      samples <- lapply(reps, function(r) list(X = r$X, Y = r$Y))
      
      
      # Summary for this (n, alpha) scenario
      key <- paste0("n", n, "_a", alpha_idx)
      results[[key]] <- list(
        n = n, 
        alpha = alpha,
        mat = mat,    # MC x 12, saves h and ASE for each replicate
        samples = samples, # saves the MC samples
        
        # Mean bandwidth for each method 
        mean_ase_nw_case = mean(mat[, "nw_ase_case"]),
        mean_ase_ll_case = mean(mat[, "ll_ase_case"]),
        sd_ase_nw_case = sd(mat[, "nw_ase_case"]),
        sd_ase_ll_case = sd(mat[, "ll_ase_case"]),
        mean_h_nw_case = mean(mat[, "nw_h_case"]),
        mean_h_ll_case = mean(mat[, "ll_h_case"]),
        
        mean_ase_nw_cv  = mean(mat[, "nw_ase_cv"]),
        mean_ase_ll_cv  = mean(mat[, "ll_ase_cv"]),
        sd_ase_nw_cv = sd(mat[, "nw_ase_cv"]),
        sd_ase_ll_cv = sd(mat[, "ll_ase_cv"]),
        mean_h_nw_cv = mean(mat[, "nw_h_cv"]),
        mean_h_ll_cv = mean(mat[, "ll_h_cv"]),
        
        
        # MGCV summaries (na.rm: a replicate may fail to estimate R)
        mean_ase_nw_mgcv = mean(mat[, "nw_ase_mgcv"], na.rm = TRUE),
        mean_ase_ll_mgcv = mean(mat[, "ll_ase_mgcv"], na.rm = TRUE),
        sd_ase_nw_mgcv   = sd(mat[, "nw_ase_mgcv"],   na.rm = TRUE),
        sd_ase_ll_mgcv   = sd(mat[, "ll_ase_mgcv"],   na.rm = TRUE),
        mean_h_nw_mgcv   = mean(mat[, "nw_h_mgcv"],   na.rm = TRUE),
        mean_h_ll_mgcv   = mean(mat[, "ll_h_mgcv"],   na.rm = TRUE),
        n_fail_nw_mgcv   = sum(is.na(mat[, "nw_ase_mgcv"])),
        n_fail_ll_mgcv   = sum(is.na(mat[, "ll_ase_mgcv"]))
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
# Print tables
################################################################################

make_table <- function(results, print_h = TRUE) {
  
  all_n <- sort(unique(sapply(results, `[[`, "n")))
  all_alpha <- sort(unique(sapply(results, `[[`, "alpha")))
  
  ## GCV columns only if the simulation was run with do_gcv = TRUE
  has_gcv <- !is.null(results[[1]]$mean_ase_nw_gcv)
  
  fmt   <- function(m, s) sprintf("%.5f (%.5f)", m, s)
  fmt_h <- function(h)    sprintf("%.5f", h)
  
  rows <- list()
  
  for (ai in seq_along(all_alpha)) {
    for (n in all_n) {
      
      key <- paste0("n", n, "_a", ai)
      r <- results[[key]]
      
      if (is.null(r)) {
        warning(sprintf("Key '%s' not found — skipping.", key))
        next
      }
      
      row <- data.frame(
        alpha   = all_alpha[ai],
        n = n,
        NW_CV   = fmt(r$mean_ase_nw_cv,   r$sd_ase_nw_cv),
        NW_MGCV = fmt(r$mean_ase_nw_mgcv, r$sd_ase_nw_mgcv),
        NW_CASE = fmt(r$mean_ase_nw_case, r$sd_ase_nw_case),
        LL_CV = fmt(r$mean_ase_ll_cv,   r$sd_ase_ll_cv),
        LL_MGCV = fmt(r$mean_ase_ll_mgcv, r$sd_ase_ll_mgcv),
        LL_CASE = fmt(r$mean_ase_ll_case, r$sd_ase_ll_case),
        stringsAsFactors = FALSE
      )
      
      if (has_gcv) {
        row$NW_GCV <- fmt(r$mean_ase_nw_gcv, r$sd_ase_nw_gcv)
        row$LL_GCV <- fmt(r$mean_ase_ll_gcv, r$sd_ase_ll_gcv)
      }
      
      if (print_h) {
        row$NW_h_CV <- fmt_h(r$mean_h_nw_cv)
        row$NW_h_MGCV <- fmt_h(r$mean_h_nw_mgcv)
        row$NW_h_CASE <- fmt_h(r$mean_h_nw_case)
        row$LL_h_CV <- fmt_h(r$mean_h_ll_cv)
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
