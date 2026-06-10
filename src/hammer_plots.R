################################################################################
#---------------- Hammer projection plot ----------------
################################################################################
library(polykde)
library(viridis)
library(DirStats)
hammer_plot <- function(x, cols, pch = 16, cex = 1.5, fct = 1.012,
                        lwd_border = 5, lwd_grid = 1, ...) {
  
  # Check input
  stopifnot(nrow(x) == length(cols))
  
  # Heatmap
  plot(sph_to_hammer(x), col = cols, pch = pch, cex = cex, axes = FALSE,
       xlab = "", ylab = "", ...)
  
  # Spherical grids
  theta <- seq(0, 2 * pi, l = 201)
  phi <- seq(0, pi, l = 200)
  for (phi_i in seq(0, pi, l = 10)[-c(1, 10)]) {
    
    ang <- cbind(phi_i, theta)
    lines(sph_to_hammer(angles_to_sph(ang)[, 3:1]), lwd = lwd_grid)
    
  }
  for (theta_i in seq(0, 2 * pi, l = 10)[-c(1, 10)]) {
    
    ang <- cbind(phi, theta_i)
    lines(sph_to_hammer(angles_to_sph(ang)[, 3:1]), lwd = lwd_grid)
    
  }
  for (theta_i in c(0, 2 * pi)) {
    
    ang <- cbind(phi, theta_i)
    lines(fct * sph_to_hammer(angles_to_sph(ang)[, 3:1]), lwd = lwd_border)
    
  }
  
}

#---------------- Color cuts ----------------
col_cuts <- function(x, pal = colorRampPalette(c("blue", "white", "red")),
                     breaks = NULL) {
  #intervals have equal numbers of observations in each bin
  if (is.null(breaks)) { 
    
    breaks <- quantile(x, probs = seq(0, 1, length.out = 20))
    
  }
  # give a number and returns that many colour
  cols <- pal(length(breaks))
  # assigns each value of x to a interval defined by the brakpoints
  cuts <- cut(x, breaks = breaks)
  #wiith the bin index pick correpoding color
  return(cols[cuts])
  
}

#---------------- Simple color bar ----------------
color_bar <- function(breaks, ticks,
                      pal = colorRampPalette(c("blue", "white", "red")),
                      ...) {
  
  # Create horizontal colorbar
  par(mar = c(2, 3, 2, 3))
  n_cols <- length(breaks) - 1
  cols <- pal(n_cols)
  image(x = seq_along(cols), y = 1, z = matrix(seq_along(cols), ncol = 1),
        col = cols, axes = FALSE, xlab = "", ylab = "", ...)
  
  # Map tick values to positions along the colorbar
  tick_positions <- (ticks - min(breaks)) / (max(breaks) - min(breaks)) *
    (n_cols - 1) + 1
  axis(1, at = tick_positions, labels = ticks, lwd = 4, cex.axis = 2)
  box(lwd = 4)
  
}
