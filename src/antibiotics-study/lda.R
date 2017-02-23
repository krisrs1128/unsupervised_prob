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
set.seed(11242016)
theme_set(min_theme())

# Code Block -------------------------------------------------------------------
## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

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
orde#79B5B7_map <- function(x) {
  ggheatmap(
    x %>%
    melt(value.name = "fill", varnames = c("x", "y")),
    list("x_order" = x_order, "y_order" = y_order)
  ) +
    min_theme(list(text_size = 0)) +
    labs(x = "Sample", y = "Microbe")
}

p <- orde#79B5B7_map(get_taxa(abt)) + ggtitle("Raw")
ggsave("../../doc/figure/heatmaps-1.pdf", p)

p <- orde#79B5B7_map(asinh(get_taxa(abt))) + ggtitle("asinh")
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

## ---- sample-test ----
i <- 100
samples$beta[i,, ]
samples$theta[i,, ]

x_p#79B5B7 <- array(dim = c(length(test_ix), stan_data$V, n_iter))
for (j in seq_along(test_ix)) {
  for (i in seq_len(n_iter)) {
    theta <- rgamma(4, 1)
    theta <- theta / sum(theta)
    n_j <- sum(X[test_ix[j], ])
    x_p#79B5B7[j, , i] <- rmultinom(1, size = n_j, prob = t(samples$beta[i,, ]) %*% theta)
  }
}

dim(x_p#79B5B7)
rowSums(x_p#79B5B7[,,1])
rowSums(X[test_ix, ])
hist(x_p#79B5B7[1,1, ])

x_err <- x_p#79B5B7
for (j in 1:1000) {
  x_err[,, j] <- asinh(x_p#79B5B7[,, j]) - asinh(X[test_ix, ])
}

hist(x_err)
mean((x_err) ^ 2)
mean((asinh(X[test_ix, ]) - mean(asinh(X[test_ix, ]))) ^ 2)

dmultinom(c(1, 2, 3), 6, prob = rep(1, 3) / 3, log = TRUE)

## ---- extract_beta ----
# underlying RSV distributions
beta_logit <- samples$beta

for (i in seq_len(n_iter / 2)) {
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
  left_join(taxa)

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)

## ---- extract_theta ----
theta_logit <- samples$theta
for (i in seq_len(n_iter / 2)) {
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
plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "theta_logit",
  "fill" = "topic",
  "col" = "topic",
  "fill_colors" = brewer.pal(stan_data$K, "Set1"),
  "col_colors" = brewer.pal(stan_data$K, "Set1"),
  "facet_terms" = c("topic", "."),
  "theme_opts" = list(border_size = 0.7)
)
p <- ggboxplot(data.frame(theta_hat), plot_opts) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  labs(x = "Time", y = "g(theta)") +
  theme(legend.position = "none")
ggsave("../../doc/figure/visualize_lda_theta_boxplot-1.pdf", p)

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
  "theme_opts" = list(border_size = 0.7),
  "outlier.shape" = NA
)
p <- ggboxplot(
  beta_hat %>%
  data.frame() %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:4]),
  plot_opts
) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5, col = "#999999") +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  labs(y = "g(beta)", fill = "Family") +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )
ggsave("../../doc/figure/visualize_lda_beta-1.pdf", p)

## ---- posterior-data ----
x_sim <- samples$n_sim %>%
  melt(
    varnames = c("iteration", "sample", "rsv"),
    value.name = "sim_value"
  ) %>%
  as.data.table()

mx <- x %>%
  melt(
    varnames = c("sample", "rsv"),
    value.name = "truth"
  ) %>%
  as.data.table()

## ---- posterior-hists ----
overall_hists <- rbind(
  cbind(
    mx %>%
    rename(value = truth),
    iteration = NA,
    type = "truth"
  ),
  cbind(
    x_sim %>%
    filter(iteration %in% round(seq(1, 1000, length.out = 4))) %>%
    rename(value = sim_value),
    type = "simulated"
  )
)

ggplot(overall_hists) +
  geom_histogram(aes(x = asinh(value), fill = type), bins = 100) +
  facet_grid(iteration ~ .) +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
theme(
  panel.border = element_rect(fill = "transparent", size = 0.5)
)

## ---- posterior-quantiles ----
q_probs <- seq(0, 1, 0.01)
quantiles_comp <- x_sim %>%
  group_by(iteration) %>%
  do(
    data.frame(
      type = "sim", 
      q_ix = q_probs,
      q = quantile(asinh(.$sim_value), q_probs)
    )
  )

ggplot(quantiles_comp) +
  geom_step(
    aes(x = q, y = q_ix, group = iteration),
    alpha = 0.1, position = position_jitter(h = 0.005),
  ) +
  geom_step(
    data = data.frame(
      q_ix = q_probs,
      q = quantile(asinh(mx$truth), q_probs)
    ),
    aes(x = q, y = q_ix),
    col = "#79B5B7",
    size = 0.5
  ) +
  labs(
    "x" = "x",
    "y" = "Pr(asinh(count) < x)"
  )

## ---- col-margins ----
rsv_totals <- mx %>%
  group_by(rsv) %>%
  summarise(rsv_total = sum(asinh(truth)))
rsv_totals$rank <- rank(rsv_totals$rsv_total)

sim_rsv_totals <- x_sim %>%
  group_by(iteration, rsv) %>%
  summarise(sim_total = sum(asinh(sim_value))) %>%
  left_join(rsv_totals)

p <- ggplot() +
  geom_point(
    data = sim_rsv_totals,
    aes(y = rank, x = sim_total),
    alpha = 0.1, size = 0.5
  ) +
  geom_step(
    data = sim_rsv_totals %>% filter(iteration == 1),
    aes(y = rank, x = rsv_total),
    col = "#79B5B7"
  ) +
  labs(
    "x" = "x",
    "y" = "Prob(microbe sum < x)"
  )

p <- ggplot() +
  geom_boxplot(
    data = sim_rsv_totals,
    aes(y = as.factor(rank), x = sim_total),
    alpha = 0.1, size = 0.1
  ) +
  geom_step(
    data = sim_rsv_totals %>% filter(iteration == 1),
    aes(y = rank, x = rsv_total),
    col = "#79B5B7"
  ) +
  labs(
    "x" = "x",
    "y" = "Prob(microbe sum < x)"
  ) +
  theme(
    axis.text.y = element_blank()
  )

## ---- time-series ----
mx_samples <- mx
mx_samples$sample_id  <- sample_names(abt)[mx_samples$sample]
mx_samples <- mx_samples %>%
  left_join(
    data.frame(
      sample_id = sample_names(abt),
      sample_data(abt)
    )
  ) %>%
  filter(rsv %in% sample(seq_len(ntaxa(abt)), 12)) %>%
  left_join(x_sim)

ggplot() +
  geom_point(
    data = mx_samples,
    aes(x = time, y = asinh(sim_value), group = interaction(iteration, rsv)),
    alpha = 0.01, size = 0.1
  ) +
  geom_line(
    data = mx_samples %>% filter(iteration == 1),
    aes(x = time, y = asinh(truth), group = rsv),
    size = 0.5, col = "#79B5B7"
  ) +
  facet_wrap(~rsv, scales = "free", ncol = 4)
