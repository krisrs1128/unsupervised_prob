#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

# Code block ------------------------------------------------------------------

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("phyloseq")
library("treelapse")
library("feather")
library("ggscaffold")
set.seed(11242016)

 theme_set(
  min_theme(list(text_size = 8, subtitle_size = 12))
)

softmax <- function(x) {
  exp(x) / sum(exp(x))
}

## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

## ---- vis-times ----
times <- sample_data(abt)$time
X <- t(get_taxa(abt))

## ---- run-model ----
stan_data <- list(
  N = nrow(X),
  V = ncol(X),
  T = length(times),
  times = times,
  times_mapping = times,
  X = X,
  a0 = 0.5,
  b0 = 0.5
)

m <- stan_model("../stan/unigram.stan")
stan_fit <- vb(m, data = stan_data)
save(
  stan_fit,
  file = sprintf("../../data/fits/unigram-%s.rda", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- prepare-beta ----
taxa <- data.table(
  rsv = rownames(tax_table(abt)),
  tax_table(abt)@.Data
)
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

## center the betas
beta <- samples$beta
for (i in seq_len(stan_data$T)) {
  beta[, i,] <- beta[, i, ] - mean(beta[, i,])
}

beta_hat <- samples$beta %>%
  melt(
    varnames = c("iteration", "time", "rsv_ix"),
    value.name = "beta"
  )

beta_hat$rsv <- rownames(otu_table(abt))[beta_hat$rsv_ix]
beta_hat$time <- times[beta_hat$time]
beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  left_join(sample_data(abt)[, c("time", "condition")]) %>%
  group_by(time) %>%
  mutate(prob = softmax(beta))

group_order <- sort(table(taxa$Taxon_5), decreasing = TRUE)
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = names(group_order))
beta_hat$rsv <- factor(
  taxa[beta_hat$rsv_ix]$rsv,
  levels = rownames(tax_table(abt))
)

## ---- unigramseries ----
plot_opts <- list(
  "x" = "time",
  "y" = "mean_beta",
  "col" = "Taxon_5",
  "facet_terms" = c("Taxon_5", "."),
  "alpha" = 0.4,
  "group" = "rsv"
)
p <- gglines(
  beta_hat %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:4]) %>%
  group_by(rsv, time) %>%
  summarise(mean_beta = mean(beta), Taxon_5 = Taxon_5[1]) %>%
  as.data.frame(),
  plot_opts
) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme(
    strip.text.y = element_blank(),
    legend.position = "bottom"
  )
ggsave("../../doc/figure/unigramseries.png", p)

## ---- unigramboxplots ----
plot_opts <- list(
  "x" = "rsv",
  "y" = "beta",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "outlier.shape" = NA,
  "alpha" = 0.4,
  "col_colors" = brewer.pal(8, "Set2"),
  "fill_colors" = brewer.pal(8, "Set2"),
  "theme_opts" = list(border_size = 0.7)
)
p <- ggboxplot(
  beta_hat %>%
  filter(
    Taxon_5 %in% levels(beta_hat$Taxon_5)[1:4],
    time %in% seq(10, 20, by = 3)
  ) %>%
  as.data.frame(),
  plot_opts
) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  facet_grid(condition + time ~ Taxon_5, scales = "free_x", space = "free_x") +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave("../../doc/figure/unigramboxplots.png", p, width = 10, height = 8)
