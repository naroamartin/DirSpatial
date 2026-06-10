source("simulation functions.R")
MC <- 5
n_values <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1, 1.5)
alpha_vals <- pi * alpha_base_vals
sigma2 <- 1
m_idx_vals <- c(1)
h_grid <- seq(0.03, 1.0, length.out = 60)

# Normalized geodesic neighbourhood radii.
# These correspond to 10%, 20%, and 30% of the maximum geodesic distance.
ell_vals <- c(0.1, 0.2, 0.3, 0.4)

d <- 2
cores <- parallel::detectCores() - 1

output_dir <- "sim_workspaces_mod"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

all_results <- list()
save_results <- FALSE

set.seed(123, kind = "Mersenne-Twister")

for (m_idx in m_idx_vals) {
  res <- run_simulation(
    MC = MC,
    n_values = n_values,
    alpha_vals = alpha_vals,
    sigma2 = sigma2,
    m_idx = m_idx,
    h_grid = h_grid,
    ell_vals = ell_vals,
    d = d,
    cores = cores
  )
  
  # Relabel alpha in the result objects back to the base/original scale
  # so tables show alpha = 0.5 instead of alpha = pi * 0.5.
  for (res_key in names(res)) {
    res[[res_key]]$alpha_scaled <- res[[res_key]]$alpha
    res[[res_key]]$alpha <- res[[res_key]]$alpha / pi
  }
  
  names(res) <- sub(
    pattern = "_a.*$",
    replacement = paste0("_a", alpha_base_vals[1]),
    x = names(res)
  )
  
  key <- paste0("m", m_idx)
  all_results[[key]] <- res
  
  if (save_results) {
    for (res_key in names(res)) {
      r <- res[[res_key]]
      
      n_val <- r$n
      alpha_val <- r$alpha
      mat <- r$mat
      m_label <- key
      
      fname <- sprintf(
        "%s_n%d_a%s.RData",
        m_label,
        n_val,
        gsub("\\.", "", format(alpha_val, nsmall = 1))
      )
      
      save(
        m_label,
        n_val,
        alpha_val,
        mat,
        r,
        file = file.path(output_dir, fname)
      )
    }
    
    saveRDS(all_results, file = "sim_all_results_checkpoint.rds")
    cat(sprintf("Checkpoint saved after m%d\n", m_idx))
  }
}

print_h <- FALSE
tables <- lapply(names(all_results), function(key) {
  tab   <- make_table(all_results[[key]], ell_vals = ell_vals, print_h = print_h)
  tab$m <- key
  tab
})

final_table <- do.call(rbind, tables)
b_labels    <- paste0("b", seq_along(ell_vals))

if (print_h) {
  nw_ase_cols <- c("NW_CV",   paste0("NW_MCV_b",   seq_along(ell_vals)), "NW_CASE")
  ll_ase_cols <- c("LL_CV",   paste0("LL_MCV_b",   seq_along(ell_vals)), "LL_CASE")
  nw_h_cols   <- c("NW_h_CV", paste0("NW_h_MCV_b", seq_along(ell_vals)), "NW_h_CASE")
  ll_h_cols   <- c("LL_h_CV", paste0("LL_h_MCV_b", seq_along(ell_vals)), "LL_h_CASE")
  col_order   <- c("m", "alpha", "n", nw_ase_cols, nw_h_cols, ll_ase_cols, ll_h_cols)
} else {
  nw_cols   <- c("NW_CV", paste0("NW_MCV_", b_labels), "NW_CASE")
  ll_cols   <- c("LL_CV", paste0("LL_MCV_", b_labels), "LL_CASE")
  col_order <- c("m", "alpha", "n", nw_cols, ll_cols)
}

final_table <- final_table[, col_order]
rownames(final_table) <- NULL
print(final_table)

