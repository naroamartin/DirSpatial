rm(list = ls())
source("hammer functions.R")
library(polykde)
library(viridis)
library(DirStats)

################################################################################
# Regression functions m(x)
################################################################################

## Grid on S^2
m_grid <- 501
mu <- c(0, -1, 0)
gr <- expand.grid(
  theta = seq(0, 2 * pi, length.out = m_grid + 1)[-(m_grid + 1)],
  phi   = seq(0, pi,     length.out = m_grid / 2)
)
s <- DirStats::to_sph(th = gr$theta, ph = gr$phi)

## Regression functions
m1 <- function(X) X[, 1]
m2 <- function(X) sin(pi * X[, 1]) * X[, 2]
m3 <- function(X, a = 1, b = 3/2) a * sin(2 * pi * X[, 2]) + b * cos(2 * pi * X[, 1])

reg_list <- list(
  list(fun = m1, title = expression(m[1](x) == x[1])),
  list(fun = m2, title = expression(m[2](x) == sin(pi*x[1])*x[2])),
  list(fun = m3, title = expression(m[3](x) == sin(2*pi*x[2]) + (3/2)*cos(2*pi*x[1])))
)

if (!dir.exists("fields")) dir.create("fields")

pal <- viridis::viridis

## Compute global v_lim across all functions BEFORE the loop
v_lim_global <- max(sapply(reg_list, function(r) {
  vals <- r$fun(s)
  vals <- (vals - mean(vals)) / sd(vals)
  max(abs(vals))
}))

brks <- seq(-v_lim_global, v_lim_global, length.out = 101)

for (k in seq_along(reg_list)) {
  
  vals <- reg_list[[k]]$fun(s)
  vals <- (vals - mean(vals)) / sd(vals)
  
  cols <- col_cuts(vals, pal = pal, breaks = brks)
  
  jpeg(paste0("fields/m", k, "_sphere.jpeg"),
       width = 12, height = 6, units = "in", res = 200)
  par(mar = c(0, 0, 2, 0))
  hammer_plot(x = s, cols = cols)
  points(sph_to_hammer(mu), pch = 16, cex = 2.25)
  dev.off()
}

## Shared colorbar
jpeg("fields/colorbar_shared.jpeg",
     width = 10, height = 1.25, units = "in", res = 200)
color_bar(breaks = brks,
          ticks  = pretty(c(-v_lim_global, v_lim_global), n = 5),
          pal    = pal)
dev.off()

################################################################################
# Correlation fields
################################################################################

geo_dist_to_t <- function(X, t) {
  X_norm <- X / sqrt(rowSums(X^2))
  t_norm <- t / sqrt(sum(t^2))
  ip <- pmax(pmin(X_norm %*% t_norm, 1), -1)
  return(as.numeric(acos(ip)))
}

t <- c(0, 0, 1)
alpha_vals <- c(0.5, 1, 1.5)
geo <- geo_dist_to_t(s, t)

if (!dir.exists("corr")) dir.create("corr")

for (alpha in alpha_vals) {
  
  cor_field <- exp(-alpha * geo)
  
  cols <- col_cuts(cor_field,
                   pal    = colorRampPalette(c("white", "red")),
                   breaks = seq(0, 1, length.out = 100))
  
  jpeg(paste0("corr/cor_alpha", alpha, ".jpeg"),
       width = 12, height = 6, units = "in", res = 200)
  par(mar = c(0, 0, 0, 0))
  hammer_plot(x = s, cols = cols)
  points(sph_to_hammer(t), pch = 18, cex = 4)
  dev.off()
}

jpeg("corr/cor_colorbar.jpeg",
     width = 10, height = 1.25, units = "in", res = 200)
color_bar(pal    = colorRampPalette(c("white", "red")),
          breaks = seq(0, 1, length.out = 100),
          ticks  = seq(0, 1, by = 0.25))
dev.off()
