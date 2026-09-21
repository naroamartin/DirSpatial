################################################################################
# Simulation
################################################################################
rm(list=ls())
source("bandwidth.R")
source("simul.R")

## ------ Simulation parameters  ---------------------------------------------
MC <- 500
d <- 2
n_values <- c(100, 200, 400)
alpha_vals <- c(0.6)   #weak / medium / strong
sigma2 <- 1
m_idx_vals  <- c(1)
h_grid  <- seq(0.03, 1.0, length.out = 40)
h_pilot <- 0.33
do_gcv <- FALSE
cores <- parallel::detectCores() - 1

output_dir <- "sim_workspaces"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- FALSE

## ------ Simulation loop over m  --------------------------------------------
set.seed(123, kind = "Mersenne-Twister")

for (m_idx in m_idx_vals) {
  
  res <- run_simulation(MC = MC, n_values = n_values, alpha_vals = alpha_vals,
                        sigma2 = sigma2, m_idx = m_idx, h_grid = h_grid,
                        d = d, cores = cores, h_pilot = h_pilot,
                        do_gcv = do_gcv)
  
  for (res_key in names(res)) {
    alpha_idx <- as.integer(sub(".*_a", "", res_key))   # extract 1, 2, 3
    res[[res_key]]$alpha_scaled <- res[[res_key]]$alpha
    res[[res_key]]$alpha<- alpha_vals[alpha_idx]
  }
  
  obj_name <- paste0("all_results_m", m_idx)
  assign(obj_name, res)
  
  if (save_results) {
    fname <- file.path(output_dir, sprintf("all_results_m%d.RData", m_idx))
    save(list = obj_name, file = fname)
    cat(sprintf("Saved %s\n", fname))
  }
}


m_idx <- 1
results <- get(paste0("all_results_m", m_idx))

tab <- make_table(results)
print(tab, row.names = FALSE)


## ------ Obtaining median ASE per selector -----------------
get_ase_summary <- function(results, key, est) {
  
  if (is.null(results[[key]])) {
    stop("Key '", key, "' not found. Available: ",
         paste(names(results), collapse = ", "))
  }
  mat <- results[[key]]$mat
  
  ase_cols     <- c(paste0(est, "_ase_cv"), paste0(est, "_ase_mgcv"))
  method_names <- c("CV", "MGCV")
  
  missing <- setdiff(c(ase_cols, paste0(est, "_ase_case")), colnames(mat))
  if (length(missing)) {
    stop("Columns not in mat: ", paste(missing, collapse = ", "))
  }
  
  # Median ASE for CV and MGCV (raw scale)
  med_ase <- vapply(ase_cols,
                    function(col) median(mat[, col], na.rm = TRUE),
                    numeric(1))
  names(med_ase) <- method_names
  
  best_idx <- which.min(med_ase)
  
  list(median_ase_cv   = unname(med_ase["CV"]),
       median_ase_mgcv = unname(med_ase["MGCV"]),
       median_ase_case = median(mat[, paste0(est, "_ase_case")], na.rm = TRUE),
       n_fail_mgcv     = sum(is.na(mat[, paste0(est, "_ase_mgcv")])),
       best_method     = method_names[best_idx],
       median_ase_best = unname(med_ase[best_idx]))
}

m_idx <- m_idx_vals[1]
results <- get(paste0("all_results_m", m_idx))

n_sel <- 100   
alpha_sel <- 3     
key <- paste0("n", n_sel, "_a", alpha_sel)

cat(sprintf("\nScenario: m%d, n = %d, alpha = %.2f\n",
            m_idx, n_sel, alpha_vals[alpha_sel]))

nw_summary <- get_ase_summary(results, key, "nw")
ll_summary <- get_ase_summary(results, key, "ll")

cat(sprintf(paste0("NW: median ASE \n CV = %.5f | MGCV = %.5f | CASE = %.5f",
                   " | best (%s) = %.5f  [MGCV failures: %d]\n"),
            nw_summary$median_ase_cv, nw_summary$median_ase_mgcv,
            nw_summary$median_ase_case,
            nw_summary$best_method, nw_summary$median_ase_best,
            nw_summary$n_fail_mgcv))
cat(sprintf(paste0("LL: median ASE \n CV = %.5f | MGCV = %.5f | CASE = %.5f",
                   " | best (%s) = %.5f  [MGCV failures: %d]\n"),
            ll_summary$median_ase_cv, ll_summary$median_ase_mgcv,
            ll_summary$median_ase_case,
            ll_summary$best_method, ll_summary$median_ase_best,
            ll_summary$n_fail_mgcv))

