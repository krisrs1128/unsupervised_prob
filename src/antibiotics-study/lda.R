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
library("ggscaffold")
library("feather")
set.seed(11242016)

# Code Block -------------------------------------------------------------------
## ---- get-data ----
data(abt)
n_samples <- ncol(otu_table(abt))
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .4 * n_samples, prune = TRUE) %>%
  subset_samples(ind == "D")

## ---- histograms ---
transformed_counts <- data.frame(
  count = c(get_taxa(abt), asinh(get_taxa(abt))),
  transformation = c(
    rep("original", ntaxa(abt) * nsamples(abt)),
    rep("asinh", ntaxa(abt) * nsamples(abt))
  )
)

ggplot(transformed_counts) +
  geom_histogram(aes(x = count)) +
  facet_grid(. ~ transformation, scale = "free_x") +
  min_theme(list(text_size = 8, subtitle_size = 12))

## ---- heatmaps ----
y_order <- names(sort(taxa_sums(abt)))
x_order <- names(sort(sample_sums(abt)))
ordered_map <- function(x) {
  ggheatmap(
    x %>%
    melt(value.name = "fill", varnames = c("y", "x")),
    list("x_order" = x_order, "y_order" = y_order)
  ) +
    min_theme(list(text_size = 0))
}

ggheatmap(
  asinh(get_taxa(abt)) %>%
  melt(value.name = "fill", varnames = c("y", "x")),
  list("x_order" = x_order, "y_order" = y_order)
) +
  min_theme(list(text_size = 0))

## ----  lda ----
m <- stan_model(file = "lda_counts.stan")
X <- asinh(t(otu_table(abt)@.Data))
X[] <- as.integer(round(X, 2) * 100)

stan_data <- list(
  K = 4,
  V = ncol(X),
  D = nrow(X),
  n = X,
  alpha = rep(1e-10, 4)
)

stan_fit <- vb(m, stan_data, eta = .1, adapt_engaged = FALSE, grad_samples = 1, iter = 2000)
samples <- rstan::extract(stan_fit)

## ---- extract_beta ----
# underlying RSV distributions
beta_hat <- samples$beta %>%
  melt() %>%
  setnames(c("iterations", "cluster", "rsv_ix", "beta"))
beta_hat$rsv <- rownames(tax_table(abt))[beta_hat$rsv_ix]

taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- beta_hat %>%
  left_join(taxa)

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)

## ---- visualize_beta ----
# might want to set prior for more extreme decay
plot_opts <- list(
  "x" = "rsv",
  "y" =  "beta",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "facet_terms" = c("cluster", "Taxon_5"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x",
  "theme_opts" = list(border_size = 0)
)
p2 <- ggboxplot(
  beta_hat %>%
  data.frame() %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8]),
  plot_opts
) +
  labs(y = "beta", fill = "Family") +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_blank()
  )

## ---- extract_theta ----
theta_hat <- samples$theta %>%
  melt(
    varnames = c("iteration", "sample", "cluster"),
    value.name = "theta"
  )

theta_hat$sample <- rownames(X)[theta_hat$sample]
sample_info <- sample_data(abt)
sample_info$sample <- rownames(sample_info)
theta_hat$cluster <- paste("Cluster", theta_hat$cluster)

theta_hat <- theta_hat %>%
  left_join(sample_info, by = "sample")

## ---- visualize_theta ----
plot_opts <- list(
  "x" = "time",
  "y" = "cluster",
  "fill" = "mean(theta)"
)
p1 <- ggheatmap(theta_hat, plot_opts)

plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "theta",
  "fill" = "cluster",
  "col" = "cluster",
  "fill_colors" = RColorBrewer::brewer.pal(stan_data$K, "Set2"),
  "col_colors" = RColorBrewer::brewer.pal(stan_data$K, "Set2"),
  "facet_terms" = c("cluster", "."),
  "theme_opts" = list("panel_border" = 0.7)
)
p1 <- ggboxplot(data.frame(theta_hat), plot_opts) +
  labs(x = "Time") +
  theme(legend.position = "none")

## ---- save_results ----
dir.create("results")
write_feather(theta_hat, path = "theta.feather")
write_feather(beta_hat, path = "beta.feather")
