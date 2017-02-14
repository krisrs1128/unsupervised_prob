#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("feather")
library("phyloseq")
library("treelapse")
library("ggscaffold")
set.seed(11242016)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

## ---- vis-times ----
raw_times <- sample_data(abt)$time
X <- t(asinh(get_taxa(abt)))
X[] <- as.integer(round(X, 0) * 1)

times <- 4 * round(raw_times / 4)
times_mapping <- match(times, unique(times))
times <- unique(times)

## ----  run-model ----
N <- nrow(X)
V <- ncol(X)
T <- length(times)

m <- stan_model("../src/stan/dtm.stan")
stan_data <- list(
  N = N,
  V = V,
  T = T,
  K = 2,
  sigma_hyper = c(0.5, 0.5),
  delta_hyper = c(0.5, 0.5),
  times = times,
  times_mapping = times_mapping,
  X = X
)
stan_fit <- vb(m, data = stan_data)
samples <- rstan::extract(stan_fit)

## ---- prepare-theta ----
softmax <- function(mu) {
  exp(mu) / sum(exp(mu))
}

theta_hat <- apply(samples$alpha, c(1, 2), softmax) %>%
  melt(
    varnames = c("cluster", "iteration", "time"),
    value.name = "theta"
  )
theta_hat$time <- times[theta_hat$time]

cur_samples <- data.frame(sample_data(abt))
cur_samples$time <- 4 * round(cur_samples$time / 4)
cur_samples <- cur_samples[c(1, which(diff(cur_samples$time) != 0)), ]

theta_hat <- cur_samples %>%
  right_join(theta_hat)

plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "theta",
  "fill" = "as.factor(cluster)",
  "col" = "as.factor(cluster)",
  "facet_terms" = c("cluster", "condition"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)

## ---- visualize-theta ----
ggboxplot(theta_hat, plot_opts) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    fill = "Cluster",
    x = "time"
  ) +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.2)
  )

## ---- prepare-beta ----
beta_hat <- apply(samples$beta, c(1, 2, 3), softmax) %>%
  melt(
    varnames = c("rsv_ix", "iteration", "time", "cluster")
  )
beta_hat$time <- times[beta_hat$time]

## merge in taxonomic information (for labeling evolutionary families)
taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- factor(rownames(tax_table(abt)))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
sorted_taxa <- names(sort(table(taxa$Taxon_5), decreasing = TRUE))
taxa$Taxon_5 <- factor(taxa$Taxon_5, levels = sorted_taxa)

beta_hat$rsv <- taxa$rsv[beta_hat$rsv_ix]
beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  filter(
    Taxon_5 %in% sorted_taxa[1:8],
    time <= 16,
    time >= 8
  )

plot_opts <- list(
  "x" = "rsv",
  "y" = "sqrt(value)",
  "col" = "Taxon_5",
  "fill" = "Taxon_5",
  "outlier.shape" = NA,
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)

## ---- visualize-beta ----
ggboxplot(beta_hat, plot_opts) +
  scale_y_continuous(limits = c(0.053, 0.0541), expand = c(0, 0)) +
  facet_grid(
    time ~ cluster + Taxon_5,
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_rect(fill = "transparent", size = 0.2),
    legend.position = "bottom"
  )
