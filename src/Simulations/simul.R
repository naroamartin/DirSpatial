library(foreach)
library(progressr)
library(future)
library(doRNG)
library(doFuture)

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
                        .export = c("one_rep", "cv_loo", "gcv",
                                    "mgcv", "correlation_matrix",
                                    "dist_est", "empirical_variogram",
                                    "geodesic_dist",
                                    "build_Sigma", "unif_sphere",
                                    "m_funs", "h_grid","do_gcv",
                                    "sigma2", "d", "h_pilot", "m_idx")) %dorng% {
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
