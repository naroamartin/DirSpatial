###############################################################################
# Simulation 
###############################################################################
source("src/simulation functions.R")

## ------ Simulation parameters  ----------------------------------------------
MC <- 500
d <- 2
n_values <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1.0, 1.5)
alpha_vals <- pi * alpha_base_vals # Scale alpha 
sigma2 <- 1
m_idx_vals <- c(1)
#m_idx_vals <- c(1, 2, 3, 4)   
h_grid <- seq(0.03, 1.0, length.out = 40)
ell_vals <- c(0.1, 0.2, 0.3, 0.4) 
cores <- parallel::detectCores() - 1


output_dir <- "sim_workspaces_mod"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

all_results <- list()
save_results <- TRUE

## ------ Simulation loop over m ----------------------------------------------

set.seed(123, kind = "Mersenne-Twister")

for (m_idx in m_idx_vals) {
  
  res <- run_simulation(MC = MC,  n_values = n_values, alpha_vals = alpha_vals,
                        sigma2 = sigma2, m_idx = m_idx, h_grid = h_grid,
                        ell_vals = ell_vals, d = d, cores = cores)
  
  # Relabel alpha in the result objects back to the base/original scale
  # so tables show alpha = 0.5 instead of alpha = pi * 0.5.
  for (res_key in names(res)) {
    res[[res_key]]$alpha_scaled <- res[[res_key]]$alpha
    res[[res_key]]$alpha <- res[[res_key]]$alpha / pi
  }
  
  names(res) <- sub( pattern = "_a.*$", 
                     replacement = paste0("_a", alpha_base_vals[1]),
                     x = names(res))
  
  key <- paste0("m", m_idx)
  all_results[[key]] <- res
  
  if (save_results) {
    for (res_key in names(res)) {
      
      r <- res[[res_key]]
      n_val <- r$n
      alpha_val <- r$alpha
      mat <- r$mat
      m_label <- key
      
      fname <- sprintf("%s_n%d_a%s.RData", m_label, n_val,
        gsub("\\.", "", format(alpha_val, nsmall = 1)))
      
      save(m_label,n_val, alpha_val, mat, r, 
           file = file.path(output_dir, fname))
    }
    
    save(all_results,file = "sim_workspace.RData")
    cat(sprintf("Workspace saved after m%d\n", m_idx))
  }
}

## ------ Print results  -----------------------------------------------------
print_h <- TRUE
tables <- lapply(names(all_results), function(key) {
  tab <- make_table(all_results[[key]], ell_vals = ell_vals, print_h = print_h)
  cbind(m = key, tab)
})

final_table <- do.call(rbind, tables)
rownames(final_table) <- NULL

print(final_table)






