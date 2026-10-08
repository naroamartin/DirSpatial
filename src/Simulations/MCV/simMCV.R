################################################################################
# Simulation
################################################################################
rm(list = ls())
source("setupMCV_n.R")

## ------ Simulation parameters -----------------------------------------------
MC <- 500
d <- 2
n_values <- c(100, 200, 400)
alpha_vals <- c(0.1, 0.3, 0.6)
sigma2 <- 1
m_idx_vals <- c(1)
h_grid <- seq(0.03, 1.0, length.out = 40)
ell_vals <- c(0.05, 0.1, 0.15, 0.2, 0.25)

cores <- parallel::detectCores() - 3

output_dir <- "sim_workspaces"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- FALSE

## ------ Simulation loop over m  --------------------------------------------
set.seed(123, kind = "Mersenne-Twister")

for (m_idx in m_idx_vals) {
  
  results <- run_simulation(MC = MC, n_values = n_values, alpha_vals = alpha_vals,
                            sigma2 = sigma2, m_idx = m_idx, h_grid = h_grid,
                            ell_vals = ell_vals, d = d, cores = cores)
  obj_name <- paste0("all_results_m", m_idx)
  assign(obj_name, results)
  
  if (save_results) {
    fname <- file.path(output_dir, sprintf("all_results_m%d.RData", m_idx))
    save(list = obj_name, file = fname)
    cat(sprintf("Saved %s\n", fname))
  }
}

tab <- make_table(results, ell_vals, print_h = TRUE)
print(tab, row.names = FALSE)
