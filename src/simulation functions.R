################################################################################
# Simulation study: Linear-spherical regression under spatial dependence
# Bandwidth selection: standard CV vs. modified CV (MCV)
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

# Exponential: Sigma_ij = sigma2 * exp(-alpha * n^(1 / d) * d_geo(X_i, X_j))
build_Sigma <- function(X, alpha, sigma2, d = 2) {
  n <- nrow(X)
  D <- geodesic_dist(X) / pi
  
  sigma2 * exp(-alpha * n^(1 / d) * D)
}

##------ Regression functions on S^2 -----------------------------------------

m_funs <- list(
  m1 = function(X) X[, 1],
  m2 = function(X) X[, 1]^2 - X[, 3],
  m3 = function(X) sin(pi * X[, 1]) * X[, 2],
  m4 = function(X, a = 1, b = 1.5) a * sin(2 * pi * X[, 2]) + b * 
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

cv_loo <- function(X, Y, h, p) {
  n   <- nrow(X)
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

##------ Modified  cross-validation -----------------------------------------
# 2) MCV(H) = sum_i [ Y_i - m_hat_{h,p,-N(i)}(X_i) ]^2
#  N(i) = { j : theta(X_j, X_i) <= ell }
#  For S^2 (sphere), N(i) is a spherical cap of geodesic radius ell around X_i.

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

################################################################################
# Simulation functions
################################################################################

## ------ Single Monte Carlo iteration ----------------------------------------
## Returns a named vector with ase (and selected h) for:
##   NW  (p=0): CV, MCV, CASE
##   LL  (p=1): CV, MCV, CASE

one_rep <- function(n, alpha, sigma2, m_fun, h_grid, ell_vals, d = 2) {
  
  ## Data generation
  X <- unif_sphere(n, d)               # (n x 3) points on S^2
  m_vals <- m_fun(X)                   # true regression values at X
  Sigma  <- build_Sigma(X, alpha, sigma2)
  eps <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
  Y <- m_vals + eps
  
  ## Compute standardized distances
  D <- geodesic_dist(X) / pi
  
  ##  error for a given  bandwidth and p at the data points
  ase <- function(h, p) {
    yhat <- lpe(eval_pts = X, dir_data = X, lin_data = Y, h = h, p = p)
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
    
    ## --- MCV ---
    for (b in seq_along(ell_vals)) {
      ell <- ell_vals[b]
      
      mcv_vals <- sapply(h_grid, function(h) {
        mcv_loo(X, Y, h, p, ell, D = D)
      })
      
      idx_mcv <- which.min(mcv_vals)
      
      out[paste0(tag, "_h_mcv", b)] <- h_grid[idx_mcv]
      out[paste0(tag, "_ase_mcv", b)] <- ase_vals[idx_mcv]
    }
    
  }
  
  return(out)
}


## ------ Simulation over all the MC iterations -------------------------------

run_simulation <- function(MC, n_values, alpha_vals, sigma2, m_idx,
                           h_grid, ell_vals, d = 2, cores = 1) {
  
  ## Register parallel backed
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
      
      cat(sprintf("\n--- m%d | n = %d | alpha = %.1f ---\n",
                  m_idx, n, alpha))
      
      progressr::with_progress({
        prog <- progressr::progressor(along = seq_len(MC))
        mat <- foreach(k = seq_len(MC), .inorder   = TRUE, 
                       .combine   = "rbind", 
                       .packages  = c("DirStatsOld", "MASS"),
                       .export = c("one_rep", "lpe", "cv_loo",
                                   "mcv_loo", "geodesic_dist",
                                   "build_Sigma", "unif_sphere", "m_funs",
                                   "h_grid", "ell_vals", 
                                   "sigma2", "d")) %dorng% {
                                     
                                     # Progress bar
                                     prog()
                                     # select the regression function
                                     m_function <- m_funs[[m_idx]] 
                                     # Computations
                                     one_rep(n = n, alpha = alpha, 
                                             sigma2 = sigma2,
                                             m_fun = m_function, 
                                             h_grid = h_grid,
                                             ell_vals = ell_vals, d = d)
                                   }
      })
      
      key <- paste0("n", n, "_a", alpha_idx)
      results[[key]] <- list(
        n = n,
        alpha = alpha,
        mat = mat,
        
        mean_ase_nw_cv = mean(mat[, "nw_ase_cv"]),
        mean_ase_nw_case = mean(mat[, "nw_ase_case"]),
        mean_ase_ll_cv = mean(mat[, "ll_ase_cv"]),
        mean_ase_ll_case = mean(mat[, "ll_ase_case"]),
        
        sd_ase_nw_cv = sd(mat[, "nw_ase_cv"]),
        sd_ase_nw_case = sd(mat[, "nw_ase_case"]),
        sd_ase_ll_cv = sd(mat[, "ll_ase_cv"]),
        sd_ase_ll_case = sd(mat[, "ll_ase_case"]),
        
        mean_h_nw_cv = mean(mat[, "nw_h_cv"]),
        mean_h_nw_case = mean(mat[, "nw_h_case"]),
        mean_h_ll_cv = mean(mat[, "ll_h_cv"]),
        mean_h_ll_case = mean(mat[, "ll_h_case"])
      )
      
      for (b in seq_along(ell_vals)) {
        results[[key]][[paste0("mean_ase_nw_mcv", b)]] <- 
          mean(mat[, paste0("nw_ase_mcv", b)])
        results[[key]][[paste0("sd_ase_nw_mcv", b)]] <- 
          sd(mat[, paste0("nw_ase_mcv", b)])
        results[[key]][[paste0("mean_ase_ll_mcv", b)]] <- 
          mean(mat[, paste0("ll_ase_mcv", b)])
        results[[key]][[paste0("sd_ase_ll_mcv", b)]] <- 
          sd(mat[, paste0("ll_ase_mcv", b)])
        
        results[[key]][[paste0("mean_h_nw_mcv", b)]] <- 
          mean(mat[, paste0("nw_h_mcv", b)])
        results[[key]][[paste0("mean_h_ll_mcv", b)]] <- 
          mean(mat[, paste0("ll_h_mcv", b)])
      }
    }
  }
  
  future::plan(future::sequential())
  return(results)
}


################################################################################
# Print tables
################################################################################

make_table <- function(results, ell_vals, print_h = FALSE) {
  
  all_n     <- sort(unique(sapply(results, `[[`, "n")))
  all_alpha <- sort(unique(sapply(results, `[[`, "alpha")))
  
  fmt   <- function(m, s) sprintf("%.4f (%.4f)", m, s)
  fmt_h <- function(h)    sprintf("%.4f", h)
  
  rows <- list()
  
  for (ai in seq_along(all_alpha)) {
    for (n in all_n) {
      
      key <- paste0("n", n, "_a", ai)
      r   <- results[[key]]
      
      if (is.null(r)) {
        warning(sprintf("Key '%s' not found — skipping.", key))
        next
      }
      
      row <- data.frame(
        alpha   = all_alpha[ai],
        n       = n,
        NW_CV   = fmt(r$mean_ase_nw_cv,   r$sd_ase_nw_cv),
        NW_CASE = fmt(r$mean_ase_nw_case, r$sd_ase_nw_case),
        LL_CV   = fmt(r$mean_ase_ll_cv,   r$sd_ase_ll_cv),
        LL_CASE = fmt(r$mean_ase_ll_case, r$sd_ase_ll_case),
        stringsAsFactors = FALSE
      )
      
      for (b in seq_along(ell_vals)) {
        row[[paste0("NW_MCV_b", b)]] <- fmt(r[[paste0("mean_ase_nw_mcv", b)]],
                                            r[[paste0("sd_ase_nw_mcv",  b)]])
        row[[paste0("LL_MCV_b", b)]] <- fmt(r[[paste0("mean_ase_ll_mcv", b)]],
                                            r[[paste0("sd_ase_ll_mcv",  b)]])
      }
      
      if (print_h) {
        row$NW_h_CV   <- fmt_h(r$mean_h_nw_cv)
        row$NW_h_CASE <- fmt_h(r$mean_h_nw_case)
        row$LL_h_CV   <- fmt_h(r$mean_h_ll_cv)
        row$LL_h_CASE <- fmt_h(r$mean_h_ll_case)
        for (b in seq_along(ell_vals)) {
          row[[paste0("NW_h_MCV_b", b)]] <- fmt_h(r[[paste0("mean_h_nw_mcv", b)]])
          row[[paste0("LL_h_MCV_b", b)]] <- fmt_h(r[[paste0("mean_h_ll_mcv", b)]])
        }
      }
      
      rows <- c(rows, list(row))
    }
  }
  
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  
  nw_ase <- c("NW_CV", paste0("NW_MCV_b", seq_along(ell_vals)), "NW_CASE")
  ll_ase <- c("LL_CV", paste0("LL_MCV_b", seq_along(ell_vals)), "LL_CASE")
  
  if (print_h) {
    nw_h <- c("NW_h_CV", paste0("NW_h_MCV_b", seq_along(ell_vals)), "NW_h_CASE")
    ll_h <- c("LL_h_CV", paste0("LL_h_MCV_b", seq_along(ell_vals)), "LL_h_CASE")
    col_order <- c("alpha", "n", nw_ase, nw_h, ll_ase, ll_h)
  } else {
    col_order <- c("alpha", "n", nw_ase, ll_ase)
  }
  
  tab[, col_order]
}
