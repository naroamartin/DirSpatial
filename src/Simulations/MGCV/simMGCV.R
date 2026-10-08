################################################################################
# Simulation
################################################################################
rm(list = ls())
source("src/Simulations/MGCV/setupMGCV.R")

MC <- 100
d <- 2
n_values <- c(100, 200, 400)
# Effective correlation range is alpha / n^{1/d} (same design as MCV)
alpha_vals <- c(0.1, 0.3, 0.6)
sigma2 <- 1
m_idx_vals <- c(1)

h_grid <- seq(0.03, 2.5, length.out = 60)
# Pilot bandwidth for the residuals used in the variogram fit (MGCV only).
# NULL: data-driven, h_pilot = c_pilot * h_MCV(ell_pilot), chosen per sample
# and per estimator. Set a number to fix it instead (sensitivity analysis).
h_pilot <- NULL
ell_pilot <- 0.1
c_pilot <- 2

cores <- parallel::detectCores() - 3

output_dir <- "sim_workspaces"

if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- FALSE

set.seed(123, kind = "Mersenne-Twister")


for (m_idx in m_idx_vals) {

  res <- run_simulation(MC = MC, n_values = n_values, alpha_vals = alpha_vals,
    sigma2 = sigma2, m_idx = m_idx, h_grid = h_grid, h_pilot = h_pilot,
    ell_pilot = ell_pilot, c_pilot = c_pilot, d = d, cores = cores)

  obj_name <- paste0("all_results_m", m_idx)
  assign(obj_name, res)


  ## Optional saving
  if (save_results) {
    fname <- file.path(output_dir, sprintf("all_results_m%d.RData", m_idx))
    save(list = obj_name,file = fname )
    cat(sprintf("Saved %s\n", fname))
  }
}

m_idx <- 1
results <- get(paste0("all_results_m", m_idx))

tab <- make_table(results,med = FALSE,print_h = TRUE)
print(tab, row.names = FALSE)

print(alpha_table(results),row.names = FALSE)
