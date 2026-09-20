################################################################################
# Initial functions
################################################################################
rm(list = ls())
library(MASS)

##------ Spherical distances -------------------------------------------------
geodesic_dist <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  ip <- X_norm %*% t(X_norm)
  ip <- pmax(pmin(ip, 1), -1)
  D <- acos(ip) 
  diag(D) <- 0 
  return(D)
}

##------ Covariance structure for the error term -----------------------------
# Exponential: Sigma_ij = sigma2 * exp(-d_geo(X_i, X_j) * n^(1/d) / alpha)
build_Sigma_noN <- function(X, alpha, sigma2) {
  D <- geodesic_dist(X) / pi     # normalized to [0, 1]
  sigma2 * exp(-D / alpha)
}


build_Sigma<- function(X, alpha, sigma2) {
  D <- geodesic_dist(X) / pi     # normalized to [0, 1]
  n <- nrow(X)
  q <- ncol(X) - 1               # q = 2 for S^2 embedded in R^3
  sigma2 * exp(- (n^(1/q) * D) / alpha)
}

##------ Regression functions on S^2 -----------------------------------------
m_funs <- list(
  m1 = function(X) X[, 1],
  m2 = function(X) sin(pi * X[, 1]) * X[, 2],
  m3 = function(X, a = 1, b = 1.5) a * sin(2 * pi * X[, 2]) + b * 
    cos(2 * pi * X[, 1])
)

##------ Uniform sample on S^d -----------------------------------------------
unif_sphere <- function(n, d) {
  X <- matrix(rnorm(n * (d + 1)), nrow = n, ncol = d + 1)
  return(X / sqrt(rowSums(X^2)))
}
