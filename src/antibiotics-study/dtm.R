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

## ---- run-model ----
times <- sample_data(abt)$time
x <- t(get_taxa(abt))
dimnames(x) <- NULL

m <- stan_model("../stan/dtm.stan")
stan_data <- list(
  N = nrow(x),
  V = ncol(x),
  T = length(times),
  K = 2,
  sigma_hyper = c(0.5, 0.5),
  delta_hyper = c(0.5, 0.5),
  times = times,
  times_mapping = times,
  x = x
)

## Fit and save variational bayes model
stan_fit <- vb(m, data = stan_data)
save(
  stan_fit,
  file = sprintf("../../data/fits/dtm-%s.rda", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- prepare-alpha ----
## center the alpha
alpha <- samples$alpha
for (d in seq_len(stan_data$N)) {
  for (i in seq_len(dim(alpha)[1])) {
    alpha[i, d, ] <- alpha[i, d, ] - mean(alpha[i, d, ])
  }
}

alpha_hat <- alpha %>%
  melt(
    varnames = c("iteration", "time", "topic"),
    value.name = "alpha"
  )

alpha_hat$time <- times[alpha_hat$time]
alpha_hat <- sample_data(abt) %>%
  data.frame() %>%
  unique() %>%
  right_join(alpha_hat)

alpha_plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "alpha",
  "fill" = "as.factor(topic)",
  "col" = "as.factor(topic)",
  "fill_colors" = brewer.pal(3, "Set1"),
  "col_colors" = brewer.pal(3, "Set1"),
  "facet_terms" = c("topic", "condition"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x",
  "theme_opts" = list(border_size = 0.7)
)

## ---- prepare-mu ----
mu <- samples$mu
for (k in seq_len(dim(mu)[3])) {
  for (i in seq_len(stan_data$K)) {
    mu[i,, k, ] <- mu[i,, k, ] - mean(mu[i,, k, ])
  }
}

cur_times <- times[times %in% seq(10, 20, by = 3)]
mu_hat <- mu[, cur_times,, ] %>%
  melt(
    varnames = c("iteration", "time", "topic", "rsv_ix")
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
    Taxon_5 %in% sorted_taxa[1:4]
  )

mu_plot_opts <- list(
  "x" = "rsv",
  "y" = "value",
  "col" = "Taxon_5",
  "fill" = "Taxon_5",
  "outlier.shape" = NA,
  "fill_colors" = brewer.pal(4, "Set2"),
  "col_colors" = brewer.pal(4, "Set2"),
  "theme_opts" = list(border_size = 0.7)
)

## ---- visualizealpha ----
p <- ggboxplot(alpha_hat, alpha_plot_opts) +
  labs(
    fill = "Topic",
    x = "time"
  ) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  theme(
    legend.position = "bottom"
  )
ggsave("../../doc/figure/visualize-alpha-1.pdf", p)

## ---- visualizemu ----
p <- ggboxplot(mu_hat, mu_plot_opts) +
  facet_grid(
    time ~ topic + Taxon_5,
    scales = "free_x",
    space = "free_x"
  ) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  labs(
    "col" = "Taxonomic Family",
    "fill" = "Taxonomic Family"
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave("../../doc/figure/visualize-mu-1.pdf", p)

## ---- posterior-checks ----
counts_data_checker(x, samples$n_sim, "../../doc/figure/dtm_post_checks")
