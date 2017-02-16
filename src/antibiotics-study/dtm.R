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

times <- 4 * floor(raw_times / 4)
times <- raw_times
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
rm(stan_fit)

## ---- prepare-theta ----
alpha <- samples$alpha
for (i in seq_len(dim(alpha)[2])) {
  alpha[, i, ] <- alpha[, i, ] - mean(alpha[, i, ])
}

alpha_hat <- alpha %>%
  melt(
    varnames = c("iteration", "time", "cluster"),
    value.name = "alpha"
  )
alpha_hat$time <- times[alpha_hat$time]

cur_samples <- data.frame(sample_data(abt))

#cur_samples$time <- 4 * floor(cur_samples$time / 4)

alpha_hat <- cur_samples %>%
  unique() %>%
  right_join(alpha_hat)

plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "alpha",
  "fill" = "as.factor(cluster)",
  "col" = "as.factor(cluster)",
  "fill_colors" = brewer.pal(3, "Set2"),
  "col_colors" = brewer.pal(3, "Set2"),
  "facet_terms" = c("cluster", "condition"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)

## ---- visualize-theta ----
ggboxplot(alpha_hat, plot_opts) +
  labs(
    fill = "Topic",
    x = "time"
  ) +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.2),
    legend.position = "bottom"
  )

## ---- prepare-beta ----
mu <- samples$mu
for (k in seq_len(dim(mu)[3])) {
  mu[,, k, ] <- mu[,, k, ] - mean(mu[,, k, ])
}

cur_times <- c(11, 12, 13)
mu_hat <- mu[, cur_times,, ] %>%
  melt(
    varnames = c("iteration", "time", "cluster", "rsv_ix")
  )
mu_hat$time <- cur_times[mu_hat$time]

## merge in taxonomic information (for labeling evolutionary families)
taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- factor(rownames(tax_table(abt)))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
sorted_taxa <- names(sort(table(taxa$Taxon_5), decreasing = TRUE))
taxa$Taxon_5 <- factor(taxa$Taxon_5, levels = sorted_taxa)

mu_hat$rsv <- taxa$rsv[mu_hat$rsv_ix]
mu_hat <- mu_hat %>%
  left_join(taxa) %>%
  filter(
    Taxon_5 %in% sorted_taxa[1:8]
  )

plot_opts <- list(
  "x" = "rsv",
  "y" = "value",
  "col" = "Taxon_5",
  "fill" = "Taxon_5",
  "outlier.shape" = NA,
  "col_cols" = brewer.pal(3, "Set2"),
  "fill_cols" = brewer.pal(3, "Set2")
)

## ---- visualize-mu ----
ggboxplot(mu_hat %>% filter(rsv_ix < 100), plot_opts) +
  facet_grid(
    time ~ cluster + Taxon_5,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(limits = c(-3.5, 8), oob = scales::rescale_none) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.3) +
  labs(
    "col" = "Taxonomic Family",
    "fill" = "Taxonomic Family"
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_rect(fill = "transparent", size = 0.2),
    legend.position = "bottom"
  )
