rm(list = ls())

library(sphunif)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

d <- 4
n <- 5

# Variograms 
gamma_sph <- function(dist, c0, c1, a) {
  (c0 + c1 * (1.5 * dist / a - 0.5 * (dist / a)^3)) * (dist <= a) +
    (c0 + c1) * (dist > a)
}

gamma_rq <- function(dist, c0, c1, a) {
  (c0 + c1 * (dist^2 / ((1 / 19) * a^2 + dist^2))) * (dist != 0)
}

gamma_exp <- function(dist, c0, c1, a) {
  (c0 + c1 * (1 - exp(-3 * dist / a))) * (dist != 0)
}

gamma_funs <- list( 
  spherical = gamma_sph,
  rational_quadratic = gamma_rq,
  exponential = gamma_exp
)

# Correlations
cov_sph <- function(dist, c0, c1, a, gamma_fun) {
  1 - gamma_fun(dist, c0, c1, a)
}

# Using connection (7) in Gneitting
cov_th <- function(th, gamma_fun) cov_sph(dist = 2 * sin(th / 2), 
                               c0 = 0, c1 = 1, a = 0.8, gamma_fun)


df <- imap_dfr(gamma_funs, function(gamma_fun, model_name) {
  map_dfr(1:6, function(dim) {
    
    psi_fun <- function(th) {cov_th(th,gamma_fun)}
    
    coefs <- Gegen_coefs(k = 1:20, p = dim + 1, psi = psi_fun, 
                         N = 1280, Gauss = FALSE)
    
    tibble(model = model_name, d = dim, k = 1:20, coef = coefs
    )
  })
})

# Check if there is any negative term
df %>%
  filter(model == "rational_quadratic") %>%
  select(d, k, coef) %>%
  print(n = Inf)

df %>%
  filter(model == "exponential") %>%
  select(d, k, coef) %>%
  print(n = Inf)

##----------- Rename the model labels before plotting
df$model_label <- factor(df$model,
                         levels = c("exponential", "spherical", "rational_quadratic"),
                         labels = c("Exponential", "Spherical", "Rational quadratic")
)

ggplot(df, aes(x = k, y = coef, color = factor(d))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~model_label, scales = "free_y") +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 1e-3),
    limits = range(df$coef)
  ) +
  scale_x_continuous(breaks = 1:20) +
  labs(
    x     = "k",
    y     = "Gegenbauer coefficient",
    color = "d"
  ) +
  theme_bw() +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 13),  # larger title
    panel.grid.minor  = element_line(linewidth = 0.2),
    panel.grid.major  = element_line(linewidth = 0.3),
    legend.position   = "right",
    axis.text.x       = element_text(angle = 60, hjust = 1),
    axis.ticks.y      = element_line(),
    panel.spacing.x   = unit(1, "lines"),
    panel.border      = element_rect(colour = "black", fill = NA),
    axis.title.x      = element_text(size = 12),   # x label size
    axis.title.y      = element_text(size = 12)    # y label size
  )
