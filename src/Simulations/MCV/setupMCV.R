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

##------Cross-Validation (leave-one-out)----------------------------------------
# 1) CV(h) = sum_i [ Y_i - m_hat_{h,p,-i}(X_i) ]^2

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
  
  # Bandwidth minimizing CV
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

##------ Modified  cross-validation -----------------------------------------
# 2) MCV(h) = sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
#  N(i) = { j : theta(X_j, X_i) <= ell }
#  For S^2 (sphere), N(i) is a spherical cap of geodesic radius ell around X_i.

mcv <- function(X, Y, h_grid, p, ell, D = NULL, plot = FALSE) {
  n <- nrow(X)
  if (is.null(D)) D <- geodesic_dist(X) / pi # (n x n) geodesic distances
  res <- numeric(n)
  
  mcv_error <- sapply(h_grid, function(h){ 
    for (i in seq_len(n)) {
      in_nbhd <- which(D[i, ] <= ell)
      idx_train <- setdiff(seq_len(n), in_nbhd)
      
      # Skip if too few points remain after removing N(i) to fit a 
      # degree-p polynomial
      if (length(idx_train) < p * (ncol(X) - 1) + 1) {
        res[i] <- NA_real_
        next
      }
      
      yhat_i <- loc.directional.linear(x = X[i, , drop = FALSE], 
                                       data.dir = X[idx_train, , drop = FALSE],
                                       data.lin = Y[idx_train], h = h, p = p)$Yhat
      res[i] <- (Y[i] - yhat_i)^2
    }
    if (any(!is.finite(res))) return(Inf)
    mean(res)
  })
  # Bandwidth minimizing CV
  idx <- which.min(mcv_error)
  h <- h_grid[idx]
  mcv_min <- mcv_error[idx]
  
  if (plot) {
    plot(h_grid, mcv_error, type = "b", pch = 19, xlab = "Bandwidth h",
         ylab = "CV(h)")
    points(h, mcv_min,pch = 19,cex = 1.5)
    abline(v = h,lty = 2)
  }
  return(list(h = h, mcv_min = mcv_min))
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
  
  ##  Error for a given  bandwidth and p at the data points
  ase <- function(h, p) {
   
  }
  
  ## Bandwidth selection for each estimator (p = 0 NW, p = 1 LL)
  out <- c()
  
  for (p in c(0, 1)) {
    tag <- if (p == 0) "nw" else "ll"
    
    
    ## --- ASE : minimize true ase over the grid ---
    ase_vals <- sapply(h_grid, function(h) {
      mhat <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                     h = h, p = p)$Yhat
      mean((mhat - m_vals)^2)
    })
    
    idx_ase <- if (all(!is.finite(ase_vals))) {NA} else {which.min(ase_vals)}
    
    out[paste0(tag, "_h_ase")] <- h_grid[idx_ase]
    out[paste0(tag, "_mse_ase")] <- ase_vals[idx_ase]
    
    ## --- CV ---
    h_cv <- cv(X, Y, h_grid, p)$h
    mhat_cv <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                   h = h_cv, p = p)$Yhat
    mse_cv <- mean((mhat_cv - m_vals)^2)
    out[paste0(tag, "_h_cv")] <- h_cv
    out[paste0(tag, "_mse_cv")] <- mse_cv
    
    ## --- MCV ---
    for (b in seq_along(ell_vals)) {
      ell <- ell_vals[b]
      h_mcv <- mcv(X, Y, h_grid, p, ell, D = D)$h
      
      mhat_mcv <- loc.directional.linear(x = X, data.dir = X, data.lin = Y, 
                                        h = h_mcv, p = p)$Yhat
      mse_mcv <- mean((mhat_mcv - m_vals)^2)
      out[paste0(tag, "_h_mcv", b)] <- h_mcv
      out[paste0(tag, "_mse_mcv", b)] <- mse_mcv
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
                              # Mean bandwidth for each method 
                              mean_mse_nw_cv = mean(mat[, "nw_mse_cv"]),
                              mean_mse_nw_ase = mean(mat[, "nw_mse_ase"]),
                              mean_mse_ll_cv  = mean(mat[, "ll_mse_cv"]),
                              mean_mse_ll_ase = mean(mat[, "ll_mse_ase"]),
                              
                              median_mse_nw_cv  = median(mat[, "nw_mse_cv"]),
                              median_mse_nw_ase = median(mat[, "nw_mse_ase"]),
                              median_mse_ll_cv  = median(mat[, "ll_mse_cv"]),
                              median_mse_ll_ase = median(mat[, "ll_mse_ase"]),
                              
                              sd_mse_nw_cv = sd(mat[, "nw_mse_cv"]),
                              sd_mse_nw_ase = sd(mat[, "nw_mse_ase"]),
                              sd_mse_ll_cv = sd(mat[, "ll_mse_cv"]),
                              sd_mse_ll_ase = sd(mat[, "ll_mse_ase"]),
                              
                              mean_h_nw_cv = mean(mat[, "nw_h_cv"]),
                              mean_h_nw_ase = mean(mat[, "nw_h_ase"]),
                              mean_h_ll_cv = mean(mat[, "ll_h_cv"]),
                              mean_h_ll_ase = mean(mat[, "ll_h_ase"]),
                              
                              median_h_nw_cv = median(mat[, "nw_h_cv"]),
                              median_h_nw_ase = median(mat[, "nw_h_ase"]),
                              median_h_ll_cv = median(mat[, "ll_h_cv"]),
                              median_h_ll_ase = median(mat[, "ll_h_ase"]) 
                              )
      
      # MCV summaries for each neighborhood size
      for (b in seq_along(ell_vals)) {
        results[[key]][[paste0("mean_mse_nw_mcv", b)]] <-
          mean(mat[, paste0("nw_mse_mcv", b)])
        results[[key]][[paste0("mean_mse_ll_mcv", b)]] <-
          mean(mat[, paste0("ll_mse_mcv", b)])
        
        results[[key]][[paste0("median_mse_nw_mcv", b)]] <-
          median(mat[, paste0("nw_mse_mcv", b)])
        results[[key]][[paste0("median_mse_ll_mcv", b)]] <-
          median(mat[, paste0("ll_mse_mcv", b)])
        
        results[[key]][[paste0("sd_mse_nw_mcv", b)]] <-
          sd(mat[, paste0("nw_mse_mcv", b)])
        results[[key]][[paste0("sd_mse_ll_mcv", b)]] <-
          sd(mat[, paste0("ll_mse_mcv", b)])
        
        results[[key]][[paste0("mean_h_nw_mcv", b)]] <-
          mean(mat[, paste0("nw_h_mcv", b)])
        results[[key]][[paste0("mean_h_ll_mcv", b)]] <-
          mean(mat[, paste0("ll_h_mcv", b)])
        
        results[[key]][[paste0("median_h_nw_mcv", b)]] <-
          median(mat[, paste0("nw_h_mcv", b)])
        results[[key]][[paste0("median_h_ll_mcv", b)]] <-
          median(mat[, paste0("ll_h_mcv", b)])
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
