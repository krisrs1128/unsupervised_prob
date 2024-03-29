#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Generate data for multinomial mixture model, and run rstan code for it
# Based on
# https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf
# page 157

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phyloseq")
library("treelapse")
library("RColorBrewer")
library("ggscaffold")
library("feather")
source("./posterior_checks.R")
set.seed(11242016)
theme_set(min_theme())

# Code Block -------------------------------------------------------------------
## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

releveled_sample_data <- abt %>%
  sample_data %>%
  mutate(
    condition = revalue(
      condition,
      c("Pre Cp" = "Pre",
        "1st Cp" = "1st Course",
        "1st WPC" = "1st Course",
        "2nd Cp" = "2nd Course",
        "2nd WPC" = "2nd Course",
        "Post Cp" = "Post")
    )
  )
rownames(releveled_sample_data) <- abt %>%
  sample_data %>%
  rownames
sample_data(abt) <- releveled_sample_data

## ---- histograms ---
transformed_counts <- data.frame(
  count = c(get_taxa(abt), asinh(get_taxa(abt))),
  transformation = c(
    rep("original", ntaxa(abt) * nsamples(abt)),
    rep("asinh", ntaxa(abt) * nsamples(abt))
  )
)

p <- ggplot(transformed_counts) +
  geom_histogram(aes(x = count)) +
  facet_grid(. ~ transformation, scale = "free_x") +
  min_theme(list(text_size = 8, subtitle_size = 12))
ggsave("../../doc/figure/histograms-1.pdf", p)

## ---- heatmaps ----
x_order <- names(sort(taxa_sums(abt)))
y_order <- names(sort(sample_sums(abt)))
ordered_map <- function(x) {
  ggheatmap(
    x %>%
    melt(value.name = "fill", varnames = c("x", "y")),
    list("x_order" = x_order, "y_order" = y_order)
  ) +
    min_theme(list(text_size = 0)) +
    labs(x = "Sample", y = "Microbe")
}

p <- ordered_map(get_taxa(abt)) + ggtitle("Raw")
ggsave("../../doc/figure/heatmaps-1.pdf", p)

p <- ordered_map(asinh(get_taxa(abt))) + ggtitle("asinh")
ggsave("../../doc/figure/heatmaps-2.pdf", p)

## ---- lda ----
x <- t(get_taxa(abt))
dimnames(x) <- NULL
stan_data <- list(
  K = 4,
  V = ncol(x),
  D = nrow(x),
  n = x,
  alpha = rep(1, 4),
  gamma = rep(0.5, ncol(x))
)

m <- stan_model(file = "../stan/lda_counts.stan")
n_iter <- 1000
stan_fit <- vb(m, stan_data, iter = 2 * n_iter)
save(
  stan_fit,
  file = sprintf("../../data/fits/lda-%s.rda", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- extract_beta ----
# underlying RSV distributions
beta_logit <- samples$beta

for (i in seq_len(n_iter)) {
  for (k in seq_len(stan_data$K)) {
    beta_logit[i, k, ] <- log(beta_logit[i, k, ])
    beta_logit[i, k, ] <- beta_logit[i, k, ] - mean(beta_logit[i, k, ])
  }
}

beta_hat <- beta_logit %>%
  melt() %>%
  setnames(c("iterations", "topic", "rsv_ix", "beta_logit"))

beta_hat$rsv <- rownames(tax_table(abt))[beta_hat$rsv_ix]
taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  mutate(
    topic = paste("Topic", topic),
    Taxon_5 = stringr::str_extract(Taxon_5, "[^_]+")
  )

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)

## ---- extract_theta ----
theta_logit <- samples$theta
for (i in seq_len(n_iter)) {
  for (d in seq_len(stan_data$D)) {
    theta_logit[i, d, ] <- log(theta_logit[i, d, ])
    theta_logit[i, d, ] <- theta_logit[i, d, ] - mean(theta_logit[i, d, ])
  }
}

theta_hat <- theta_logit %>%
  melt(
    varnames = c("iteration", "sample", "topic"),
    value.name = "theta_logit"
  )

theta_hat$sample <- sample_names(abt)[theta_hat$sample]
sample_info <- sample_data(abt)
sample_info$sample <- rownames(sample_info)
theta_hat$topic <- paste("Topic", theta_hat$topic)

theta_hat <- theta_hat %>%
  left_join(sample_info, by = "sample")

## ---- visualize_lda_theta_heatmap ----
plot_opts <- list(
  "x" = "time",
  "y" = "topic",
  "fill" = "mean_theta",
  "y_order" = paste("Topic", stan_data$K:1)
)

p <- ggheatmap(
  theta_hat %>%
  group_by(topic, time) %>%
  summarise(mean_theta = mean(theta_logit, na.rm = TRUE)) %>%
  as.data.frame(),
  plot_opts
) +
  labs(fill = "g(theta)")
ggsave("../../doc/figure/visualize_lda_theta_heatmap-1.pdf", p, width = 7, height = 0.9)

## ---- visualize_lda_theta_boxplot ----
p <- ggplot(theta_hat) +
  geom_boxplot(
    aes(x = as.factor(time), y = theta_logit),
    fill = "#C9C9C9",
    outlier.size = 0.05,
    size = 0.1,
    notchwidth = 0.1
  ) +
scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  min_theme(list(border_size = 0.7, text_size = 10, subtitle_size = 11)) +
  facet_grid(topic ~ condition, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  labs(x = "Time", y = expression(paste("g(", theta[k], ")"))) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks = seq(1, 60, by = 10) - 1)
ggsave("../../doc/figure/visualize_lda_theta_boxplot-1.pdf", p, width = 6, height = 2.9)

## ---- visualize_lda_beta ----
plot_opts <- list("x" = "rsv",
  "y" =  "beta_logit",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "facet_terms" = c("topic", "Taxon_5"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x",
  "fill_colors" = brewer.pal(stan_data$K, "Set2"),
  "col_colors" = brewer.pal(stan_data$K, "Set2"),
  "theme_opts" = list(border_size = 0.7, text_size = 10, subtitle_size = 11),
  "outlier.shape" = NA,
  "alpha" = 1,
  "size" = 0.4
)
p <- ggboxplot(
  beta_hat %>%
  data.frame() %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:4]),
  plot_opts
) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-9, 7)) +
  labs(x = "Species", y = expression(paste("g(", beta[k], ")")), fill = "Family") +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave("../../doc/figure/visualize_lda_beta-1.pdf", p, width = 6, height = 3.5)

## ---- posterior-checks ----
counts_data_checker(x, samples$n_sim, "lda_post_checks")
