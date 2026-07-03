################################################################################
##               RHEA CRATERS REAL DATA APPLICATION                          
################################################################################
rm(list=ls())
if (!requireNamespace("DirStatsOld", quietly = TRUE)) {
  install.packages("DirStatsOld_0.1.5.tar", repos = NULL, type = "source")
}
library(DirStatsOld)
library(sphunif)
library(rsample)
source("functions.R")


#-------------------------------------------------------------------------------
# Load dataset
#-------------------------------------------------------------------------------
# rhea: locations and diameters for all known craters in Rhea (Saturn's moon)

data(rhea)
cat("Columns:", names(rhea), "\n")
cat("Number of craters: n =", nrow(rhea), "\n")
cat("Diameter range:", round(range(rhea$diameter), 1))

#-------------------------------------------------------------------------------
# Spherical covariate X and response Y
#-------------------------------------------------------------------------------

# Spherical covariate
# x=cos(phi)*cos(theta)
# y=cos(phi)*sin(theta)
# z=sin(phi)
# theta = longitude in [0, 2*pi) and phi  = latitude  in [-pi/2, pi/2]  

X <- to.sphere(theta = rhea$theta, phi = rhea$phi, type = "sin")
colnames(X)<-c("x1","x2","x3")

# Log-diameter : Response
Y <- rhea$diameter

# sum(Y > 15) / length(Y)
# keep <- Y > 15
# X <- X[keep, ]      # matrix: keep rows; use X[keep] only if X is a vector
# Y <- log(Y[keep])
# length(Y)
#Transform Y since its right skewed

Y <- log10(rhea$diameter)
sum(Y > 2.5 & Y < 3.5) / length(Y)

cat("\nSpherical covariate , dimension =", nrow(X), "x", ncol(X), "\n") #3596 x 3 
cat("Range of log(diameter):", round(range(Y), 3), "\n") # 2.303 6.109
cat("Mean  of log(diameter):", round(mean(Y), 3), "\n") # 2.953

#-------------------------------------------------------------------------------
# Train / test split (80% / 20%), stratified on the response
#-------------------------------------------------------------------------------
set.seed(4567)

split_df <- data.frame(.row = seq_len(nrow(X)), Y = Y)
rhea_split <- initial_split(split_df, prop = 0.8, strata = Y, breaks = 5)

# Recover original row indices from each side of the split
train_idx <- training(rhea_split)$.row
test_idx <- testing(rhea_split)$.row

# Subset into matrix (X) and vector (Y) form
X_train <- X[train_idx, ]
Y_train <- Y[train_idx]
X_test <- X[test_idx, ]
Y_test <- Y[test_idx]

cat(sprintf("\nTrain/test split: n_train = %d  n_test = %d\n", 
            length(train_idx), length(test_idx)))
#-------------------------------------------------------------------------------
# Bandwidth selection by CV and MCV (on TRAINING data only)
#-------------------------------------------------------------------------------
h_grid <- 10^seq(log10(0.05), log10(5), l = 40)
p <- 0   # Nadaraya-Watson 

# -------- 1) CV --------
cv_vals <- sapply(h_grid, function(h) cv_loo(X_train, Y_train, h, p))
h_cv <- h_grid[which.min(cv_vals)]
cat(sprintf("CV-optimal bandwidth: h_cv = %.3f\n", h_cv)) 
# CV-optimal bandwidth: h_cv = 0.129

if (!dir.exists("Rhea/Plots/Cross-Validation")) {
  dir.create("Rhea/Plots/Cross-Validation", recursive = TRUE)
}
# CV curve plot
jpeg("Rhea/Plots/Cross-Validation/rhea_cv_curve.jpeg", width = 9, height = 6,
     units = "in", res = 200)
par(mar = c(4, 4, 2, 1))
plot(h_grid, cv_vals, type = "b", pch = 19, col = "steelblue",
     xlab = "h", ylab = "CV(h)")
abline(v = h_cv, col = "firebrick", lwd = 2, lty = 2)
legend("topright", legend = paste0("h_CV = ", round(h_cv, 2)),
       col = "firebrick", lwd = 2, lty = 2, bty = "n")
dev.off()

# -------- 2) MCV --------
D_train <- geodesic_dist(X_train) / pi
ell_vals <- seq(from = 0.01, to = 0.2, by = 0.02)

cat("n_train =", nrow(X), ":\n")
for (ell in ell_vals) {
  # Count how many points are within ell for each i, on average
  in_nbhd <- rowSums(D_train <= ell)   # how many neighbours excluded
  pct_removed <- mean(in_nbhd / nrow(X)) * 100
  n_remaining <- round(nrow(X) * (1 - mean(in_nbhd / nrow(X))))
  cat(sprintf("  ell = %.2f -> avg removed = %.1f%%  avg remaining = %d points\n",
              ell, pct_removed, n_remaining))
}

mcv_vals_list <- lapply(ell_vals, function(ell) {
  cat(sprintf("  ell = %.2f...\n", ell))
  sapply(h_grid, function(h) mcv_loo(X_train, Y_train, h, p, ell, D = D_train))
})
h_mcv <- sapply(mcv_vals_list, function(v) h_grid[which.min(v)])

for (j in seq_along(ell_vals)) {
  cat(sprintf("  ell = %.2f  ->  h_mcv = %.3f\n", ell_vals[j], h_mcv[j]))
}

# ell = 0.01  ->  h_mcv = 0.145
# ell = 0.03  ->  h_mcv = 0.294
# ell = 0.05  ->  h_mcv = 0.331
# ell = 0.07  ->  h_mcv = 0.372
# ell = 0.09  ->  h_mcv = 0.372
# ell = 0.11  ->  h_mcv = 0.372
# ell = 0.13  ->  h_mcv = 0.294
# ell = 0.15  ->  h_mcv = 0.372
# ell = 0.17  ->  h_mcv = 0.372
# ell = 0.19  ->  h_mcv = 0.331


# MCV curves
for (j in seq_along(ell_vals)) {
  fname <- sprintf("Rhea/Plots/Cross-Validation/rhea_mcv_curve_ell%03.0f.jpeg", 
                   ell_vals[j] * 100)
  jpeg(fname, width = 9, height = 6, units = "in", res = 200)
  par(mar = c(4, 4, 2, 1))
  plot(h_grid, mcv_vals_list[[j]], type = "b", pch = 19, col = "steelblue",
       xlab = "h",
       ylab = sprintf("MCV(h),  ell = %.2f", ell_vals[j]))
  abline(v = h_mcv[j], col = "firebrick", lwd = 2, lty = 2)
  legend("topright", legend = paste0("h_MCV = ", round(h_mcv[j], 2)),
         col = "firebrick", lwd = 2, lty = 2, bty = "n")
  dev.off()
}

#-------------------------------------------------------------------------------
# Evaluate prediction error on TEST data (MSE)
#-------------------------------------------------------------------------------
# -------- CV --------
pred_cv <- lpe(eval_pts = X_test, dir_data = X_train,
               lin_data = Y_train, h = h_cv, p = p)
ase_cv <- mean((Y_test - as.vector(pred_cv))^2)
cat(sprintf("Test ASE  (h_cv  = %.3f): %.4f\n", h_cv, ase_cv)) 
# Test ASE  (h_cv  = 0.129): 0.2682

# -------- MCV --------
ase_mcv <- numeric(length(ell_vals))
for (j in seq_along(ell_vals)) {
  pred_j <- lpe(eval_pts = X_test, dir_data = X_train,
                lin_data = Y_train, h = h_mcv[j], p = p)
  ase_mcv[j] <- mean((Y_test - as.vector(pred_j))^2)
  cat(sprintf("Test ASE  (h_mcv = %.3f, ell = %.2f): %.4f\n",
              h_mcv[j], ell_vals[j], ase_mcv[j]))
}

# Test ASE  (h_mcv = 0.14, ell = 0.01): 0.2678
# Test ASE  (h_mcv = 0.29, ell = 0.03): 0.2684
# Test ASE  (h_mcv = 0.33, ell = 0.05): 0.2688
# Test ASE  (h_mcv = 0.37, ell = 0.07): 0.2692
# Test ASE  (h_mcv = 0.37, ell = 0.09): 0.2692
# Test ASE  (h_mcv = 0.37, ell = 0.11): 0.2692
# Test ASE  (h_mcv = 0.29, ell = 0.13): 0.2684
# Test ASE  (h_mcv = 0.37, ell = 0.15): 0.2692
# Test ASE  (h_mcv = 0.37, ell = 0.17): 0.2692
# Test ASE  (h_mcv = 0.33, ell = 0.19): 0.2688



#-------------------------------------------------------------------------------
# Plot Error difference
#-------------------------------------------------------------------------------
all_ase <- c(ase_cv, ase_mcv)
all_h <- c(h_cv, h_mcv)
all_labels <- c(0,ell_vals)

if (!dir.exists("Rhea/Plots/Cross-Validation")) {
  dir.create("Rhea/Plots/Cross-Validation", recursive = TRUE)
}

jpeg("Rhea/Plots/Cross-Validation/rhea_test_ase_comparison.jpeg", width = 7, 
     height = 5, units = "in", res = 200)
par(mar = c(5, 5, 2, 1))
plot(seq_along(all_ase), all_ase, pch = 8, cex = 1.0, xaxt = "n", xlab = "",
     ylab = "Prediction errors", xlim = c(0.5, length(all_ase) + 0.5),
     panel.first = grid(nx = NA, ny = NULL, lty = 1, col = "grey90"))
axis(1, at = seq_along(all_ase), labels = all_labels, las = 2, cex.axis = 0.75)    
dev.off()

#-------------------------------------------------------------------------------
# Best prediction error
#-------------------------------------------------------------------------------
best_idx <- which.min(all_ase)
h_best <- all_h[best_idx]

# Select the bandwidth with minimum test ase
cat(sprintf("\nBest selector: %s\n",all_labels[best_idx])) # 0.21
cat(sprintf("Best bandwidth: h_best = %.2f\n", h_best)) # h_best = 0.23
cat(sprintf("Best test MSE:  %.4f\n",  all_ase[best_idx])) # 0.2677

#-------------------------------------------------------------------------------
# NW estimator with h_cv
#-------------------------------------------------------------------------------
fit_cv <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                 h = h_cv, p = 0)
nw_cv <- as.vector(fit_cv$Yhat)
cat(" Estimation range:", round(range(nw_cv, na.rm = TRUE), 3), "\n")
#Estimation range: 2.7 3.327 

#-------------------------------------------------------------------------------
# Compute residuals
#-------------------------------------------------------------------------------
residuals_cv <- Y - nw_cv
cat(sprintf("Residuals: mean = %.4f   sd = %.4f\n",mean(residuals_cv), 
            sd(residuals_cv)))
#Residuals: mean = 0.0056   sd = 0.4645

#-------------------------------------------------------------------------------
# Test for spatial dependence
#-------------------------------------------------------------------------------

moran_I <- function(residuals, W) {
  # W: symmetric n×n weight matrix with zero diagonal
  n <- length(residuals)
  e <- residuals - mean(residuals)
  W_sum <- sum(W)
  I <- (n / W_sum) * (sum(W * outer(e, e)) / sum(e^2))
  return(I)
}

# ---- Weight matrix 1: k-nearest neighbors ----
# w_ij = 1 if j is among the k nearest neighbors of i, 0 otherwise

knn_weights <- function(D, k) {
  n <- nrow(D)
  W <- matrix(0, n, n)
  for (i in 1:n) {
    nn_idx <- order(D[i, ])[2:(k + 1)]  # exclude self (rank 1)
    W[i, nn_idx] <- 1
  }
  W <- (W + t(W)) / 2   # w_ij = w_ji
  diag(W) <- 0          # ensure w_ii = 0
  stopifnot(all(diag(D) == 0), !any(duplicated(X)))
  return(W)
}

# Geodesic distances on the full dataset (normalised by pi so range is [0,1])
D <- geodesic_dist(X) / pi
k <- 10
W_knn <- knn_weights(D, k)
B <- 999   # number of permutations

#--------- 1) Moran's I -------------------------------------------------------
# I > 0: positive autocorrelation 
# I ~ 0: no autocorrelation 
# I < 0: negative autocorrelation 

# -- kNN weights --
set.seed(123)
I_knn <- moran_I(residuals_cv, W_knn) #-0.01436599
perm_knn <- replicate(B, moran_I(sample(residuals_cv), W_knn))

p_left_knnI <- (sum(perm_knn <= I_knn) + 1) / (B + 1)  # negative autocorrelation
p_right_knnI  <- (sum(perm_knn >= I_knn) + 1) / (B + 1)  # positive autocorration
p_two_sided_knnI <- 2 * min(p_left_knnI, p_right_knnI)


cat(sprintf("Moran's I (kNN): I = %.5f\n", I_knn))
cat(sprintf("  Left p:  %.3f  |  Right p:  %.3f  |  Two-sided p:  %.3f\n",
            p_left_knnI, p_right_knnI, p_two_sided_knnI))
#Moran's I (kNN):    I = -0.01437
#Left p:  0.021  |  Right p:  0.980  |  Two-sided p:  0.042

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
cat(sprintf("Weight matrices: kNN (k=%d)\n", k))
cat(sprintf("Moran's I: kNN  I=%.5f  p(two)=%.3f\n", I_knn,  p_two_sided_knnI))

cat(sprintf("Moran's I (kNN): I = %.5f\n", I_knn))
cat(sprintf("  Left p:  %.3f  |  Right p:  %.3f  |  Two-sided p:  %.3f\n",
            p_left_knnI, p_right_knnI, p_two_sided_knnI))


# Residuals from h_best
fit_best <- loc.directional.linear(x = X, data.dir = X, data.lin = Y,
                                     h = h_best, p = 0)
nw_best <- as.numeric(fit_best$Yhat)
residuals_best <- Y - nw_best

set.seed(123)
I_knn_best  <- moran_I(residuals_best, W_knn)
perm_knn_best <- replicate(B, moran_I(sample(residuals_best), W_knn))
p_left_best  <- (sum(perm_knn_best <= I_knn_best) + 1) / (B + 1)
p_right_best <- (sum(perm_knn_best >= I_knn_best) + 1) / (B + 1)
p_two_best   <- 2 * min(p_left_best, p_right_best)

#Left p:  0.031  |  Right p:  0.970  |  Two-sided p:  0.062
cat(sprintf("Moran's I (h_best): I = %.5f  p(two) = %.3f\n",
            I_knn_best, p_two_best))

#-------------------------------------------------------------------------------
# Save workspace
#-------------------------------------------------------------------------------
save.image("Rhea/rhea_results_log10.RData")

