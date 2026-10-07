################################################################################
# Simulation study: Linear-spherical regression under spatial dependence
# Bandwidth selection: standard CV vs. modified CV (MCV)
################################################################################
if (!requireNamespace("DirStatsOld", quietly = TRUE)) {
  install.packages("DirStatsOld_0.1.5.tar.gz", repos = NULL, type = "source")
}
library(DirStats)
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
  ip <- pmin(pmax(ip, -1), 1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
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

##------ Modified  cross-validation -----------------------------------------
# MCV(h) = (1/n) sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
#  N(i) = { j : theta(X_j, X_i) / pi <= ell }
#  For S^2 (sphere), N(i) is a spherical cap of geodesic radius pi * ell around
#  X_i (ell is on the standardized scale [0, 1]; ell = 0.5 is a hemisphere).
#  loc.directional.linear accepts a vector of bandwidths, so each leave-out fit
#  is computed once for the whole h_grid.

mcv <- function(X, Y, h_grid, p, ell, D = NULL, plot = FALSE) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi # (n x n) geodesic distances
  min_train <- p * (ncol(X) - 1) + 1
  res <- matrix(NA_real_, nrow = n, ncol = length(h_grid))

  for (i in seq_len(n)) {
    # N(i) always contains i itself (D[i, i] = 0)
    idx_train <- which(D[i, ] > ell)

    # Skip if too few points remain after removing N(i) to fit a
    # degree-p polynomial
    if (length(idx_train) < min_train) next

    yhat_i <- loc.directional.linear(x = X[i, , drop = FALSE],
                                     data.dir = X[idx_train, , drop = FALSE],
                                     data.lin = Y[idx_train], h = h_grid,
                                     p = p)$Yhat
    res[i, ] <- (Y[i] - drop(yhat_i))^2
  }
  mcv_error <- colMeans(res)
  mcv_error[!is.finite(mcv_error)] <- Inf

  # No valid bandwidth in the grid: return NA instead of h_grid[1]
  if (all(!is.finite(mcv_error))) return(list(h = NA_real_, mcv_min = NA_real_))

  # Bandwidth minimizing MCV
  idx <- which.min(mcv_error)
  h <- h_grid[idx]
  mcv_min <- mcv_error[idx]

  if (plot) {
    plot(h_grid, mcv_error, type = "b", pch = 19, xlab = "Bandwidth h",
         ylab = "MCV(h)")
    points(h, mcv_min,pch = 19,cex = 1.5)
    abline(v = h,lty = 2)
  }
  return(list(h = h, mcv_min = mcv_min))
}

##------Cross-Validation (leave-one-out)----------------------------------------
# CV(h) = (1/n) sum_i [ Y_i - m_hat_{h,p,-i}(X_i) ]^2
# Computed as MCV with ell = 0 (exact leave-one-out). The hat-matrix shortcut
# (Y_i - Yhat_i) / (1 - S_ii) breaks down for small h, where S_ii rounds to 1.

cv <- function(X, Y, h_grid, p, D = NULL, plot = FALSE) {
  out <- mcv(X, Y, h_grid, p, ell = 0, D = D, plot = plot)
  return(list(h = out$h, cv_min = out$mcv_min))
}


################################################################################
# Simulation functions
################################################################################

## ------ Single Monte Carlo iteration ----------------------------------------
## Returns a named vector with ase (and selected h) for:
##   NW  (p=0): CV, MCV, ASE
##   LL  (p=1): CV, MCV, ASE

one_rep <- function(n, alpha, sigma2, m_fun, h_grid, ell_vals, d = 2) {
  
  ## Data generation
  X <- unif_sphere(n, d)               # (n x 3) points on S^2
  m_vals <- m_fun(X)                   # true regression values at X
  
  ## Compute standardized distances
  D <- geodesic_dist(X) / pi
  
  Sigma <- sigma2 * exp(-D / alpha) 
  
  if (inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    stop(sprintf(" Sigma is not positive definite for n = %d and alpha = %.3f", 
                 n, alpha))
  }
  
  eps <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
  #eps <- rnorm(n, 0, sqrt(sigma2)) #Caso independiente
  Y <- m_vals + eps
  
  
  ## Bandwidth selection for each estimator (p = 0 NW, p = 1 LL)
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
    h_cv <- cv(X, Y, h_grid, p, D = D)$h
    out[paste0(tag, "_h_cv")] <- h_cv
    out[paste0(tag, "_mse_cv")] <- ase_vals[match(h_cv, h_grid)]

    ## --- MCV ---
    for (b in seq_along(ell_vals)) {
      h_mcv <- mcv(X, Y, h_grid, p, ell_vals[b], D = D)$h
      out[paste0(tag, "_h_mcv", b)] <- h_mcv
      out[paste0(tag, "_mse_mcv", b)] <- ase_vals[match(h_mcv, h_grid)]
    }
    
  }
  
  return(list(results = out, X = X, Y = Y))
}


## ------ Simulation over all the MC iterations -------------------------------
run_simulation <- function(MC, n_values, alpha_vals, sigma2, m_idx,
                           h_grid, ell_vals, d = 2, cores) {
  
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
                        .export   = c("one_rep", "cv", "mcv", 
                                      "geodesic_dist", "unif_sphere", 
                                      "m_funs")) %dorng% { prog()
                                        one_rep(n = n, alpha = alpha, 
                                                sigma2 = sigma2,
                                                m_fun   = m_funs[[m_idx]],
                                                h_grid  = h_grid,
                                                ell_vals = ell_vals, d = d)
                                      }
      })
      
      # Save results of h and ASE per method into a matrix
      mat <- do.call(rbind, lapply(reps, `[[`, "results"))
      
      # Store the data (X, Y) from each replication as a list of length MC
      samples <- lapply(reps, function(r) list(X = r$X, Y = r$Y))
      
      # Summary for this (n, alpha) scenario
      key <- paste0("n", n, "_a", alpha_idx)
      results[[key]] <- list(n = n, alpha = alpha, mat = mat, samples = samples,
                              # Replications where a selector failed (NA h)
                              n_na = colSums(is.na(mat[, grep("_h_", colnames(mat)),
                                                       drop = FALSE])),
                              # Mean bandwidth for each method
                              mean_mse_nw_cv = mean(mat[, "nw_mse_cv"], na.rm = TRUE),
                              mean_mse_nw_ase = mean(mat[, "nw_mse_ase"], na.rm = TRUE),
                              mean_mse_ll_cv  = mean(mat[, "ll_mse_cv"], na.rm = TRUE),
                              mean_mse_ll_ase = mean(mat[, "ll_mse_ase"], na.rm = TRUE),
                              
                              median_mse_nw_cv  = median(mat[, "nw_mse_cv"], na.rm = TRUE),
                              median_mse_nw_ase = median(mat[, "nw_mse_ase"], na.rm = TRUE),
                              median_mse_ll_cv  = median(mat[, "ll_mse_cv"], na.rm = TRUE),
                              median_mse_ll_ase = median(mat[, "ll_mse_ase"], na.rm = TRUE),
                              
                              sd_mse_nw_cv = sd(mat[, "nw_mse_cv"], na.rm = TRUE),
                              sd_mse_nw_ase = sd(mat[, "nw_mse_ase"], na.rm = TRUE),
                              sd_mse_ll_cv = sd(mat[, "ll_mse_cv"], na.rm = TRUE),
                              sd_mse_ll_ase = sd(mat[, "ll_mse_ase"], na.rm = TRUE),
                              
                              mean_h_nw_cv = mean(mat[, "nw_h_cv"], na.rm = TRUE),
                              mean_h_nw_ase = mean(mat[, "nw_h_ase"], na.rm = TRUE),
                              mean_h_ll_cv = mean(mat[, "ll_h_cv"], na.rm = TRUE),
                              mean_h_ll_ase = mean(mat[, "ll_h_ase"], na.rm = TRUE),
                              
                              median_h_nw_cv = median(mat[, "nw_h_cv"], na.rm = TRUE),
                              median_h_nw_ase = median(mat[, "nw_h_ase"], na.rm = TRUE),
                              median_h_ll_cv = median(mat[, "ll_h_cv"], na.rm = TRUE),
                              median_h_ll_ase = median(mat[, "ll_h_ase"], na.rm = TRUE) 
                              )
      
      # MCV summaries for each neighborhood size
      for (b in seq_along(ell_vals)) {
        results[[key]][[paste0("mean_mse_nw_mcv", b)]] <-
          mean(mat[, paste0("nw_mse_mcv", b)], na.rm = TRUE)
        results[[key]][[paste0("mean_mse_ll_mcv", b)]] <-
          mean(mat[, paste0("ll_mse_mcv", b)], na.rm = TRUE)
        
        results[[key]][[paste0("median_mse_nw_mcv", b)]] <-
          median(mat[, paste0("nw_mse_mcv", b)], na.rm = TRUE)
        results[[key]][[paste0("median_mse_ll_mcv", b)]] <-
          median(mat[, paste0("ll_mse_mcv", b)], na.rm = TRUE)
        
        results[[key]][[paste0("sd_mse_nw_mcv", b)]] <-
          sd(mat[, paste0("nw_mse_mcv", b)], na.rm = TRUE)
        results[[key]][[paste0("sd_mse_ll_mcv", b)]] <-
          sd(mat[, paste0("ll_mse_mcv", b)], na.rm = TRUE)
        
        results[[key]][[paste0("mean_h_nw_mcv", b)]] <-
          mean(mat[, paste0("nw_h_mcv", b)], na.rm = TRUE)
        results[[key]][[paste0("mean_h_ll_mcv", b)]] <-
          mean(mat[, paste0("ll_h_mcv", b)], na.rm = TRUE)
        
        results[[key]][[paste0("median_h_nw_mcv", b)]] <-
          median(mat[, paste0("nw_h_mcv", b)], na.rm = TRUE)
        results[[key]][[paste0("median_h_ll_mcv", b)]] <-
          median(mat[, paste0("ll_h_mcv", b)], na.rm = TRUE)
      }
    }
  }
  
  future::plan(future::sequential())
  return(results)
}


################################################################################
# Print tables
################################################################################
make_table <- function(results, ell_vals, med = FALSE, print_h = TRUE) {
  
  stat <- if (med) "median" else "mean"
  
  fmt <- function(x, s) sprintf("%.5f (%.5f)", x, s)
  fmt_h <- function(h) sprintf("%.5f", h)
  rows <- list()

  for (r in results) {
    row <- data.frame(alpha = r$alpha,n = r$n,
    NW_CV = fmt( r[[paste0(stat, "_mse_nw_cv")]], r$sd_mse_nw_cv),
    NW_ASE = fmt(r[[paste0(stat, "_mse_nw_ase")]], r$sd_mse_nw_ase),
    LL_CV = fmt(r[[paste0(stat, "_mse_ll_cv")]], r$sd_mse_ll_cv),
    LL_ASE = fmt(r[[paste0(stat, "_mse_ll_ase")]], r$sd_mse_ll_ase),
    stringsAsFactors = FALSE)
    
    for (b in seq_along(ell_vals)) {
      row[[paste0("NW_MCV_b", b)]] <- fmt(r[[paste0(stat, "_mse_nw_mcv", b)]],
                                          r[[paste0("sd_mse_nw_mcv", b)]])
      row[[paste0("LL_MCV_b", b)]] <- fmt(r[[paste0(stat, "_mse_ll_mcv", b)]],
                                          r[[paste0("sd_mse_ll_mcv", b)]])
    }
    
    if (print_h) {
      
      row$NW_h_CV <- fmt_h(r[[paste0(stat, "_h_nw_cv")]])
      row$NW_h_ASE <- fmt_h(r[[paste0(stat, "_h_nw_ase")]])
      row$LL_h_CV <- fmt_h(r[[paste0(stat, "_h_ll_cv")]])
      row$LL_h_ASE <- fmt_h(r[[paste0(stat, "_h_ll_ase")]])
      
      for (b in seq_along(ell_vals)) {
        row[[paste0("NW_h_MCV_b", b)]] <- fmt_h(
          r[[paste0(stat, "_h_nw_mcv", b)]])
        row[[paste0("LL_h_MCV_b", b)]] <- fmt_h(
          r[[paste0(stat, "_h_ll_mcv", b)]]
        )
      }
    }
    rows <- c(rows, list(row))
  }
  
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  
  nw_mse <- c("NW_CV", paste0("NW_MCV_b", seq_along(ell_vals)),"NW_ASE")
  ll_mse <- c("LL_CV", paste0("LL_MCV_b", seq_along(ell_vals)), "LL_ASE")
  
  
  if (print_h) {
    nw_h <- c("NW_h_CV",paste0("NW_h_MCV_b", seq_along(ell_vals)),"NW_h_ASE")
    ll_h <- c("LL_h_CV", paste0("LL_h_MCV_b", seq_along(ell_vals)),"LL_h_ASE")
    col_order <- c("alpha", "n", nw_mse, nw_h, ll_mse, ll_h)
  } else {
    col_order <- c("alpha","n",nw_mse,ll_mse)
  }
  
  tab <- tab[, col_order, drop = FALSE]
  
  ## Order by alpha first, then by n
  tab <- tab[order(tab$alpha, tab$n), ]
  
  rownames(tab) <- NULL
  
  tab
}
