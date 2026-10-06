################################################################################
# Simulation
################################################################################
rm(list = ls())
source("setupMGCV.R")

## ------ Simulation parameters  ---------------------------------------------
MC <- 300
d <- 2
n_values <- c(100, 200, 400)
alpha_vals <- c(0.1, 0.3, 0.6)     #weak / medium / strong
sigma2 <- 1
m_idx_vals  <- c(1)
h_grid  <- seq(0.03, 1.0, length.out = 40)
h_pilot <- NULL
do_gcv <- FALSE
cores <- parallel::detectCores() - 3


output_dir <- "sim_workspaces"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- FALSE

## ------ Simulation loop over m ----------------------------------------------
set.seed(123)

for (m_idx in m_idx_vals) {
  
  res <- run_simulation(
    MC = MC, 
    n_values = n_values, 
    alpha_vals = alpha_vals,
    sigma2 = sigma2, 
    m_idx = m_idx, 
    h_grid = h_grid,
    d = d, 
    cores = cores, 
    h_pilot = h_pilot,
    do_gcv = do_gcv
  )
  
  obj_name <- paste0("all_results_m", m_idx)
  assign(obj_name, res)
  
  if (save_results) {
    fname <- file.path(output_dir, sprintf("all_results_m%d.RData", m_idx))
    save(list = obj_name, file = fname)
    cat(sprintf("Saved %s\n", fname))
  }
}

## ------ Output Summaries -----------------------------------------------------
m_idx <- 1
results <- get(paste0("all_results_m", m_idx))

# Print bandwidth and ASE comparison table
tab <- make_table(results)
print(tab, row.names = FALSE)

