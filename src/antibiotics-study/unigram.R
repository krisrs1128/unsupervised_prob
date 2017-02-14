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
library("phyloseq")
library("treelapse")
library("feather")
library("ggscaffold")
set.seed(11242016)

theme_set(
  min_theme(list(text_size = 8, subtitle_size = 12))
)

## ---- utils ----
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

# Code Block -------------------------------------------------------------------
## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

## ---- vis-times ----
raw_times <- sample_data(abt)$time
X <- asinh(t(otu_table(abt)@.Data))
X[] <- as.integer(round(X, 0) * 1)

times <- 4 * floor(raw_times / 4)
times_mapping <- match(times, unique(times))
times <- unique(times)

## ----  run-model ----
N <- nrow(X)
V <- ncol(X)
T <- length(times)

stan_data <- list(
  N = N,
  V = V,
  T = T,
  times = times,
  times_mapping = times_mapping,
  X = X,
  a0 = 0.5,
  b0 = 0.5
)

m <- stan_model("../src/stan/unigram.stan")
stan_fit <- vb(m, data = stan_data)
samples <- rstan::extract(stan_fit)

## ---- prepare-beta ----
taxa <- data.table(
  rsv = rownames(tax_table(abt)),
  tax_table(abt)@.Data
)
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- samples$beta %>%
  melt(
    varnames = c("iteration", "time", "rsv_ix"),
    value.name = "beta"
  )
beta_hat$rsv <- rownames(otu_table(abt))[beta_hat$rsv_ix]
beta_hat$time <- times[beta_hat$time]
beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  group_by(time) %>%
  mutate(prob = softmax(beta))

group_order <- sort(table(taxa$Taxon_5), decreasing = TRUE)
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = names(group_order))
beta_hat$rsv <- factor(taxa[beta_hat$rsv_ix]$rsv, levels = rownames(tax_table(abt)))

## ---- unigram-series ----
plot_opts <- list(
  "x" = "time",
  "y" = "sqrt(mean_prob)",
  "col" = "Taxon_5",
  "alpha" = 0.4,
  "group" = "rsv",
  "facet_terms" = c("Taxon_5", ".")
)
gglines(
  beta_hat %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8]) %>%
  group_by(rsv, time) %>%
  summarise(mean_prob = mean(prob), Taxon_5 = Taxon_5[1]) %>%
  as.data.frame(),
  plot_opts
) +
  scale_y_continuous(breaks = scales::pretty_breaks(2)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme(
    strip.text.y = element_blank(),
    legend.position = "bottom"
  )

## ---- unigram-beta-boxplots ----
plot_opts <- list(
  "x" = "rsv",
  "y" = "sqrt(prob)",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "outlier.shape" = NA,
  "alpha" = 0.4,
  "facet_terms" = c("time", "Taxon_5"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)
ggboxplot(
  beta_hat %>%
  filter(
    Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8]
  ) %>%
  as.data.frame(),
  plot_opts
) +
  scale_y_continuous(breaks = scales::pretty_breaks(2)) +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
