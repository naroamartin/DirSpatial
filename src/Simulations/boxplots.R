################################################################################
# Boxplots: ASE distribution by method, estimator, n, alpha — ggplot2 version
################################################################################
rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)

# Parameters
m_idx_vals <- 3
n_vals <- c(100, 200, 400)
alpha_base_vals <- c(0.5, 1.0, 1.5)
a_vals <- 1:3
ell_vals <- c(0.01, 0.05, 0.10, 0.20, 0.30)

# Load saved workspace
setwd("/Users/naroamartin/Desktop/TFM/Simulations")
load(paste0("sim_workspaces/all_results_m", m_idx_vals, "_final.RData"))
results <- get(paste0("all_results_m", m_idx_vals))

# Derived counts
n_mcv <- length(ell_vals)
n_methods <- n_mcv + 2

# Method names and ordering
method_names  <- c("CV", paste0("MCV(", ell_vals, ")"), "CASE")
method_levels <- method_names

eligible_methods <- c("CV", paste0("MCV(", ell_vals, ")"))

# Fill colours by role:
#   other = non-best CV/MCV (medium blue), best = best CV/MCV (dark blue),
#   CASE  = pink
fill_cols <- c(other = "lightskyblue", best = "blue2", CASE = "hotpink2")

# Column names in mat
nw_mcv_cols <- paste0("nw_ase_mcv", seq_len(n_mcv))
ll_mcv_cols <- paste0("ll_ase_mcv", seq_len(n_mcv))
ase_cols    <- c("nw_ase_cv", nw_mcv_cols, "nw_ase_case",
                 "ll_ase_cv", ll_mcv_cols, "ll_ase_case")

# ------ Global y upper limit -----------------------------------------------
global_whisker_max <- 0
for (ai in a_vals) {
  for (n in n_vals) {
    key <- paste0("n", n, "_a", ai)
    mat <- results[[key]]$mat
    for (est in c("nw", "ll")) {
      mcv_cols <- paste0(est, "_ase_mcv", seq_len(n_mcv))
      all_cols <- c(paste0(est, "_ase_cv"), mcv_cols, paste0(est, "_ase_case"))
      for (col in all_cols) {
        w <- boxplot.stats(log1p(mat[, col]))$stats[5]
        global_whisker_max <- max(global_whisker_max, w, na.rm = TRUE)
      }
    }
  }
}
y_upper <- global_whisker_max * 1.05

# y-axis breaks on the original ASE scale, restricted to [0, y_upper]
all_ase  <- numeric(0)
for (ai in a_vals)
  for (n in n_vals)
    all_ase <- c(all_ase, as.numeric(results[[paste0("n", n, "_a", ai)]]$mat[, ase_cols]))
y_max_raw          <- max(all_ase, na.rm = TRUE)
y_breaks           <- pretty(c(0, y_max_raw), n = 5)
y_breaks_panel     <- y_breaks[log1p(y_breaks) <= y_upper]
y_breaks_log_panel <- log1p(y_breaks_panel)

if (!dir.exists("boxplots/boxplots final")) {
  dir.create("boxplots/boxplots final")
}

# --- Plot loop -------------------------------------------------------------
for (est in c("nw", "ll")) {
  for (ai in a_vals) {
    for (n in n_vals) {
      key      <- paste0("n", n, "_a", ai)
      mat      <- results[[key]]$mat
      mcv_cols <- paste0(est, "_ase_mcv", seq_len(n_mcv))
      all_cols <- c(paste0(est, "_ase_cv"), mcv_cols, paste0(est, "_ase_case"))
      
      df <- as.data.frame(mat[, all_cols])
      colnames(df) <- method_names
      df_long <- pivot_longer(df, cols = everything(),
                              names_to = "method", values_to = "ase")
      df_long$method <- factor(df_long$method, levels = method_levels)
      
      # Best median among CV and MCVs only
      stats <- df_long %>%
        group_by(method) %>%
        summarise(med = median(log1p(ase), na.rm = TRUE), .groups = "drop")
      stats_eligible <- filter(stats, method %in% eligible_methods)
      best_method    <- stats_eligible$method[which.min(stats_eligible$med)]
      
      # Fill role per box
      df_long$fillcat <- ifelse(df_long$method == "CASE", "CASE",
                                ifelse(df_long$method == best_method, "best", "other"))
      df_long$fillcat <- factor(df_long$fillcat, levels = c("other", "best", "CASE"))
      
      p <- ggplot(df_long, aes(x = method, y = log1p(ase), fill = fillcat)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.35,
                     colour = "grey30") +
        scale_fill_manual(values = fill_cols) +
        scale_y_continuous(
          breaks = y_breaks_log_panel,
          labels = format(round(y_breaks_panel, 4), trim = TRUE)
        ) +
        coord_cartesian(ylim = c(0, y_upper)) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 9) +
        theme(
          legend.position    = "none",
          axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y        = element_text(size = 7),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin        = margin(3, 5, 3, 3)
        )
      
      ggsave(
        paste0("boxplots/boxplots final/m", m_idx_vals, "_", est, "_a", ai, "_n", n, ".jpeg"),
        plot = p, width = 4, height = 3, units = "in", dpi = 200
      )
    }
  }
}