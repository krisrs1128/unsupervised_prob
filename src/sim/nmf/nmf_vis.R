#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Visualize a combined experiment across various NMF / fitting
## parameterizations.
##
## author: kriss1@stanford.edu

## ---- libraries ----
## assumed running from NMF directory
library("data.table")
base_dir <- file.path("..", "src", "sim", "nmf")
source(file.path(base_dir, "nmf_utils.R"))

## ---- theta-reshape ----
## extract theta information from the fits
fits_dir <- file.path("..", "data", "fits")
fits <- list.files(fits_dir, "fit-*", full.names = TRUE)
expers <- fromJSON(
  file.path(base_dir, "config.json"),
  simplifyVector = FALSE
)

theta_fits <- reshape_all_samples(
  fits,
  file.path(base_dir, "config.json"),
  "theta",
  c("i", "k")
)
theta_fits$method <- basename(as.character(theta_fits$method))

## ---- visualizethetas-prep -----
## Visualize the fitted thetas, according to a few different simulation properties
plot_opts <- list(
  "x" = "value_1",
  "y" = "value_2",
  "fill" = "log(..level..)",
  "fill_type" = "gradient",
  "facet_terms" = c("N", "inference", "P"),
  "group" = "i",
  "alpha" = 0.05,
  "h" = 0.1,
  "mean_col" = "#e34a33",
  "x_lim" = c(0, 4),
  "y_lim" = c(0, 5),
  "text_size" = 2,
  "panel_border" = 0.1
)

## first, visualization in the non-zero-inflated case
gamma_pois_data <- theta_fits %>%
  filter(zero_inf_prob == 0, method == "nmf_gamma_poisson.stan")

theta_plots <- scores_contours(gamma_pois_data, plot_opts)

## ---- visualizethetas ----
#theta_plots$grouped +
  labs(
    "x" = expression(theta[1]),
    "y" = expression(theta[2])
  )

## ---- visualizethetashist ----
theme_set(min_theme(list(border_size = 0.5)))
mgamma_pois_data <- gamma_pois_data %>%
  melt(
    variable.name = "dimension",
    value.name = "estimate",
    measure.vars = c("value_1", "value_2")
  ) %>%
  melt(
    variable.name = "truth_dimension",
    measure.vars = c("truth_1", "truth_2"),
    value.name = "truth"
  )

ggplot(mgamma_pois_data) +
  geom_histogram(
    aes(x = sqrt(truth) - sqrt(estimate), fill = dimension),
    position = "identity", alpha = 0.7, bins = 75
  ) +
  facet_grid(formula(paste(plot_opts$facet_terms, collapse = "~"))) +
  scale_fill_manual(values = c("#d95f02", "#7570b3"))

## ---- visualize-zinf-thetas ----
zinf_data <- theta_fits %>%
  filter(zero_inf_prob != 0, P == 75, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
theta_plots <- scores_contours(zinf_data, plot_opts)
##theta_plots$grouped

## ---- visualize-betas ----
beta_fits <- reshape_all_samples(
  fits,
  file.path("batch", "config.json"),
  "beta",
  c("v", "k")
)
beta_fits$method <- basename(as.character(beta_fits$method))

plot_opts$facet_terms <- c("N", "inference", "P")
plot_opts$group <- "v"

gamma_pois_data <- beta_fits %>%
  filter(zero_inf_prob == 0, method == "nmf_gamma_poisson.stan")

beta_plots <- scores_contours(gamma_pois_data, plot_opts)
ggsave("~/test.png", beta_plots$grouped)

zinf_data <- beta_fits %>%
  filter(zero_inf_prob != 0, P == 75, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
theta_plots <- scores_contours(zinf_data, plot_opts)
##theta_plots$grouped
