################################################################################
# Simulation
################################################################################
rm(list=ls())
source("simulation functions.R")

## ------ Simulation parameters  ---------------------------------------------
MC <- 500
d <- 2
n_values <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1.0, 1.5)
alpha_vals  <- pi * alpha_base_vals   # internal scaled values
sigma2 <- 1
m_idx_vals  <- c(1, 2, 3)
h_grid  <- seq(0.03, 1.0, length.out = 40)
ell_vals <- c(0.01, 0.05, 0.1, 0.2, 0.25)

cores <- parallel::detectCores() - 1

output_dir <- "sim_workspaces"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- FALSE


#----- Expected proportion of observations removed by the MCV neighborhood.---
q <- function(ell) {
  (1 - cos(pi * ell)) / 2
}

tab_ell <- expand.grid( ell = ell_vals, n = n_values)

tab_ell$q_ell <- q(tab_ell$ell)
tab_ell$deleted_obs <- 1 + (tab_ell$n - 1) * tab_ell$q_ell
tab_ell$deleted_prop <- tab_ell$deleted_obs / tab_ell$n

tab_ell$deleted_perce <- (tab_ell$deleted_obs / tab_ell$n)*100

tab_ell <- tab_ell[order(tab_ell$ell, tab_ell$n), ]

tab_ell$q_ell <- sprintf("%.4f", tab_ell$q_ell)
tab_ell$deleted_obs <- sprintf("%.2f", tab_ell$deleted_obs)
tab_ell$deleted_perce<- sprintf("%.4f", tab_ell$deleted_perce)

print(tab_ell)

## ------ Simulation loop over m  --------------------------------------------
set.seed(123, kind = "Mersenne-Twister")

for (m_idx in m_idx_vals) {
  
  res <- run_simulation(MC = MC, n_values = n_values, alpha_vals = alpha_vals,
                        sigma2 = sigma2, m_idx = m_idx, h_grid = h_grid,
                        ell_vals = ell_vals, d = d, cores = cores)
  
  # Fix $alpha: convert from scaled (pi * base) back to base (0.5, 1.0, 1.5)
  for (res_key in names(res)) {
    alpha_idx <- as.integer(sub(".*_a", "", res_key))   # extract 1, 2, 3
    res[[res_key]]$alpha_scaled <- res[[res_key]]$alpha
    res[[res_key]]$alpha<- alpha_base_vals[alpha_idx]
  }
  
  obj_name <- paste0("all_results_m", m_idx)
  assign(obj_name, res)
  
  if (save_results) {
    fname <- file.path(output_dir, sprintf("all_results_m%d.RData", m_idx))
    save(list = obj_name, file = fname)
    cat(sprintf("Saved %s\n", fname))
  }
}

## ------ Print results  -----------------------------------------------------
m_idx_vals  <- c(3) # Put the number of the one you want
print_h <- FALSE 
m_obj_names <- paste0("all_results_m", m_idx_vals)
tables <- lapply(m_obj_names, function(obj_name) {
  res_obj <- get(obj_name)
  m_label <- sub("all_results_", "", obj_name)
  tab <- make_table(res_obj, ell_vals = ell_vals, print_h = print_h)
  cbind(m = m_label, tab)
})
final_table <- do.call(rbind, tables)
rownames(final_table) <- NULL
print(final_table)

## ------ Obtaining median h_ -----------------------------------------------------

get_ase_summary <- function(results, key, est, ell_vals) {
  mat <- results[[key]]$mat
  
  mcv_ase_cols <- paste0(est, "_ase_mcv", seq_along(ell_vals))
  ase_cols  <- c(paste0(est, "_ase_cv"), mcv_ase_cols)
  method_names <- c("CV", paste0("MCV(", ell_vals, ")"))
  
  # Median ASE for CV and each MCV (raw scale)
  med_ase <- sapply(ase_cols, function(col) median(mat[, col], na.rm = TRUE))
  names(med_ase) <- method_names
  
  best_idx <- which.min(med_ase)
  best_method <- method_names[best_idx]
  
  list(
    median_ase_cv = med_ase["CV"],
    median_ase_case = median(mat[, paste0(est, "_ase_case")], na.rm = TRUE),
    best_method = best_method,
    median_ase_best = med_ase[best_idx]
  )
}

m_idx_vals  <- c(3)
results <- get(paste0("all_results_m", m_idx_vals))
key <- "n400_a3"

nw_summary <- get_ase_summary(results, key, "nw", ell_vals)
ll_summary <- get_ase_summary(results, key, "ll", ell_vals)

cat(sprintf("NW: median ASE \n CV = %.5f | CASE = %.5f | best (%s) = %.5f\n",
            nw_summary$median_ase_cv, nw_summary$median_ase_case,
            nw_summary$best_method, nw_summary$median_ase_best))
cat(sprintf("LL: median ASE \n CV = %.5f | CASE (oracle) = %.5f | best (%s) = %.5f\n",
            ll_summary$median_ase_cv, ll_summary$median_ase_case,
            ll_summary$best_method, ll_summary$median_ase_best))

