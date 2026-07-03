library(DirStats)
library(doFuture)
library(foreach)
library(doRNG)
library(progressr)
library(polykde)
library(viridis)
################################################################################
# Functions
################################################################################

# ------ Hammer plots -------------------------------------------------
hammer_plot <- function(x, cols, pch = 16, cex = 1.5, fct = 1.012,
                        lwd_border = 5, lwd_grid = 1, ...) {
  
  stopifnot(nrow(x) == length(cols))
  
  plot(sph_to_hammer(x), col = cols, pch = pch, cex = cex, axes = FALSE,
       xlab = "", ylab = "", ...)
  
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

col_cuts <- function(x, pal = colorRampPalette(c("blue", "white", "red")),
                     breaks = NULL) {
  if (is.null(breaks)) {
    breaks <- quantile(x, probs = seq(0, 1, length.out = 20))
  }
  cols <- pal(length(breaks))
  cuts <- cut(x, breaks = breaks)
  return(cols[cuts])
}

color_bar <- function(breaks, ticks,
                      pal = colorRampPalette(c("blue", "white", "red")),
                      ...) {
  par(mar = c(2, 3, 2, 3))
  n_cols <- length(breaks) - 1
  cols <- pal(n_cols)
  image(x = seq_along(cols), y = 1, z = matrix(seq_along(cols), ncol = 1),
        col = cols, axes = FALSE, xlab = "", ylab = "", ...)
  tick_positions <- (ticks - min(breaks)) / (max(breaks) - min(breaks)) *
    (n_cols - 1) + 1
  axis(1, at = tick_positions, labels = ticks, lwd = 4, cex.axis = 2)
  box(lwd = 4)
}

# ------ Bandwidth selection -------------------------------------------------

## Spherical distances
geodesic_dist <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  ip <- X_norm %*% t(X_norm)
  ip <- pmax(pmin(ip, 1), -1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
}

## Projected local estimator fit 
lpe <- function(eval_pts, dir_data, lin_data, h, p) {
  # x: evaluation points
  # data.dir: directional training data
  # data.lin: linear training response
  # h: bandwith parameter
  # p: polynomial order
  fit <- loc.directional.linear(x = eval_pts, data.dir = dir_data,
                                data.lin = lin_data, h = h, p = p)
  return(fit$Yhat)
}

## Cross-Validation (leave-one-out)
cv_loo <- function(X, Y, h, p) {
  n   <- nrow(X)
  res <- numeric(n)
  for (i in seq_len(n)) {
    idx_train <- setdiff(seq_len(n), i)
    yhat_i <- lpe( eval_pts = X[i, , drop = FALSE],
                   dir_data = X[idx_train, , drop = FALSE],
                   lin_data = Y[idx_train],
                   h = h, p = p
    )
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res) 
}

## Modified  cross-validation 
mcv_loo <- function(X, Y, h, p, ell, D = NULL) {
  n <- nrow(X)
  if (is.null(D)) {
    D <- geodesic_dist(X) / pi
  } # (n x n) geodesic distances
  res <- numeric(n)
  
  for (i in seq_len(n)) {
    in_nbhd <- which(D[i, ] <= ell)
    idx_train <- setdiff(seq_len(n), in_nbhd)
    
    # Skip if too few points remain after removing N(i) to fit a 
    #degree-p polynomial
    if (length(idx_train) < (p + 2)) {
      res[i] <- NA
      next
    }
    
    yhat_i <- lpe( eval_pts = X[i, , drop = FALSE],
                   dir_data = X[idx_train, , drop = FALSE],
                   lin_data = Y[idx_train],
                   h = h, p= p
    )
    res[i] <- (Y[i] - yhat_i)^2
  }
  mean(res, na.rm = TRUE)
}

