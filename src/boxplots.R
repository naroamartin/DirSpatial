
################################################################################
# Boxplots: ASE distribution by method, estimator, n, alpha
################################################################################
rm(list = ls())

m_idx_vals <- 1
n_vals <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1.0, 1.5)
a_vals <- 1:3
ell_vals <- c(0.1, 0.2, 0.3, 0.4)

load(paste0("sim_workspaces/all_results_m", m_idx_vals, ".RData"))
results <- get(paste0("all_results_m", m_idx_vals))

col_cv <- "#4292C6"
col_mcv <- "#4292C6"
col_case <- "#CD6090"

cols <- c(col_cv, rep(col_mcv, length(ell_vals)), col_case)
method_names <- c("CV", paste0("MCV(", ell_vals, ")"), "CASE")

if (!dir.exists("boxplots")) {
  dir.create("boxplots")
}

at_pos <- c(0.5, 0.8, 1.1, 1.4, 1.7, 2)

ase_cols <- c("nw_ase_cv","nw_ase_mcv1", "nw_ase_mcv2","nw_ase_mcv3",
              "nw_ase_mcv4", "nw_ase_case","ll_ase_cv","ll_ase_mcv1",
              "ll_ase_mcv2", "ll_ase_mcv3","ll_ase_mcv4","ll_ase_case"
)

all_ase <- numeric(0)

for (ai in a_vals) {
  for (n in n_vals) {
    key <- paste0("n", n, "_a", ai)
    r <- results[[key]]
    
    if (is.null(r)) {
      next
    }
    
    all_ase <- c(all_ase, as.numeric(r$mat[, ase_cols]))
  }
}

all_ase <- all_ase[is.finite(all_ase)]

y_max_orig <- max(all_ase, na.rm = TRUE)
y_max_orig <- ceiling(y_max_orig * 10) / 10

ylim_sqrt <- sqrt(c(0, y_max_orig))

y_breaks_orig <- pretty(c(0, y_max_orig), n = 5)
y_breaks_orig <- y_breaks_orig[y_breaks_orig >= 0 & y_breaks_orig <= y_max_orig]
y_breaks_sqrt <- sqrt(y_breaks_orig)

for (ai in a_vals) {
  for (n in n_vals) {
    key <- paste0("n", n, "_a", ai)
    r <- results[[key]]
    
    if (is.null(r)) {
      next
    }
    
    mat <- r$mat
    
    box_nw <- data.frame(CV = sqrt(mat[, "nw_ase_cv"]),
                         MCV1 = sqrt(mat[, "nw_ase_mcv1"]),
                         MCV2 = sqrt(mat[, "nw_ase_mcv2"]),
                         MCV3 = sqrt(mat[, "nw_ase_mcv3"]),
                         MCV4 = sqrt(mat[, "nw_ase_mcv4"]),
                         CASE = sqrt(mat[, "nw_ase_case"])
    )
    
    box_ll <- data.frame(CV = sqrt(mat[, "ll_ase_cv"]),
                         MCV1 = sqrt(mat[, "ll_ase_mcv1"]),
                         MCV2 = sqrt(mat[, "ll_ase_mcv2"]),
                         MCV3 = sqrt(mat[, "ll_ase_mcv3"]),
                         MCV4 = sqrt(mat[, "ll_ase_mcv4"]),
                         CASE = sqrt(mat[, "ll_ase_case"])
    )
    
    for (est in c("nw", "ll")) {
      box_data <- if (est == "nw") box_nw else box_ll
      
      jpeg(paste0("boxplots/m", m_idx_vals, "_", est, "_a", ai, "_n", n, ".jpeg"),
           width = 4, height = 3, units = "in", res = 200)
      
      par(mar = c(4, 2.5, 0.5, 0.5))
      
      boxplot(box_data, at = at_pos, names = method_names,col = cols,
              outline = FALSE, ylim = ylim_sqrt, xlim = c(0.3, 2.2),
              boxwex = 0.2, cex.axis = 0.7, yaxt = "n", ylab = "", las = 2)
      
      axis(2, at = y_breaks_sqrt,labels = format(y_breaks_orig, trim = TRUE),
        cex.axis = 0.7, las = 1)
      
      dev.off()
    }
  }
}

