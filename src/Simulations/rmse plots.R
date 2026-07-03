################################################################################
# RMSE fields
################################################################################
# RMSE field for a particular sceanrio: using  NW-MCV2 and alpha=0.5
rm(list = ls())
library(polykde)
library(viridis)
library(DirStats)
source("hammer functions.R")
source("simulation functions.R")

## Grid over S^2
m_grid_size <- 501
mu <- c(0, -1, 0)
gr <- expand.grid(
  theta = seq(0, 2 * pi, length.out = m_grid_size + 1)[-(m_grid_size + 1)],
  phi   = seq(0, pi,     length.out = m_grid_size / 2)
)
s <- DirStats::to_sph(th = gr$theta, ph = gr$phi)
G <- nrow(s) # rows of the grid 

if (!dir.exists("rmse_fields")) dir.create("rmse_fields")

## Fixed parameters
m_idx <- 3  #change this to obtain the plots for each regression function
method <- "nw_h_mcv3" 
p  <- 0
alpha_idx <- 1  # fix alpha = 0.5

# Load results
load(paste0("sim_workspaces/all_results_m", m_idx, ".RData"))

# True function on grid (fixed, same for all n)
m_true <- m_funs[[m_idx]](s)

# ---------- Compute all RMSE fields and find global range ----------
rmse_fields <- list()  # strore minum and maximum values 
rmse_min <- Inf
rmse_max <- -Inf

set.seed(123, kind = "Mersenne-Twister")

for (n in c(100, 200, 400)) { 
  
  key <- paste0("n", n, "_a", alpha_idx)
  res <- get(paste0("all_results_m", m_idx))[[key]]
  mat <- res$mat
  samples <- res$samples
  MC <- nrow(mat)
  
  cat(sprintf("Computing RMSE: %s | MC = %d\n", key, MC))
  
  sum_sq <- numeric(G)
  for (j in seq_len(MC)) {
    X_j <- samples[[j]]$X
    Y_j <- samples[[j]]$Y
    h_j <- mat[j, method]
    yhat_j <- as.numeric(lpe(eval_pts = s, dir_data = X_j, lin_data = Y_j, 
                             h = h_j, p = p))
    sum_sq <- sum_sq + (yhat_j - m_true)^2
  }
  
  rmse_field <- as.numeric(sqrt(sum_sq / MC))
  rmse_fields[[key]] <- rmse_field
  
  rmse_min <- min(rmse_min, min(rmse_field)) 
  rmse_max <- max(rmse_max, max(rmse_field))
}

cat(sprintf("Global RMSE range: [%.4f, %.4f]\n", rmse_min, rmse_max))

# ---------- Plot all the RMSE fields using the global range ----------
rmse_brks  <- seq(rmse_min, rmse_max, length.out = 100)
rmse_ticks <- pretty(c(rmse_min, rmse_max), n = 4)

for (n in c(100, 200, 400)) {
  key <- paste0("n", n, "_a", alpha_idx)
  rmse_field <- rmse_fields[[key]]
  
  rmse_cols <- col_cuts(rmse_field,
                        pal = colorRampPalette(c("white", "deepskyblue1")),
                        breaks = rmse_brks)
  
  fname <- sprintf("rmse_fields/m%d_n%d_a%d_rmse_sphere.jpeg",
                   m_idx, n, alpha_idx)
  
  jpeg(fname, width = 12, height = 6, units = "in", res = 200)
  par(mar = c(0, 0, 2, 0))
  hammer_plot(x = s, cols = rmse_cols)
  points(sph_to_hammer(mu), pch = 16, cex = 2.25)
  dev.off()
  cat(sprintf("  Saved: %s\n", fname))
}

# ---------- General colorbar ----------
jpeg(sprintf("rmse_fields/m%d_a%d_rmse_colorbar.jpeg", m_idx, alpha_idx),
     width = 10, height = 1.25, units = "in", res = 200)
color_bar(pal = colorRampPalette(c("white", "deepskyblue1")), 
          breaks = rmse_brks, ticks = rmse_ticks)
dev.off()



