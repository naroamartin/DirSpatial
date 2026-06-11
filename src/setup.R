################################################################################
# Simulation
################################################################################
source("simulation functions.R")

## ------ Simulation parameters  ---------------------------------------------
MC <- 500
d <- 2
n_values <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1.0, 1.5)
alpha_vals  <- pi * alpha_base_vals   # internal scaled values
sigma2 <- 1
m_idx_vals <- c(1)
# m_idx_vals   <- c(1, 2, 3, 4)
h_grid   <- seq(0.03, 1.0, length.out = 40)
ell_vals <- c(0.1, 0.2, 0.3, 0.4)
cores <- parallel::detectCores() - 1

output_dir <- "sim_workspaces_mod"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_results <- TRUE

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



