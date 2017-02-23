#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Evaluate LDA performance on held out data

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
source("./eval.R")
set.seed(11242016)
theme_set(min_theme())

## ---- get-data ----
data(abt)
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .45 * nsamples(abt), prune = TRUE) %>%
  subset_samples(ind == "F")

## ---- lda-eval ----
test_ix <- seq(1, nrow(X), by = 5)
stan_data <- list(
  K = 4,
  V = ncol(X),
  D = nrow(X[-test_ix, ]),
  n = X[-test_ix, ],
  alpha = rep(1, 4),
  gamma = rep(0.5, ncol(X))
)

stan_fit <- vb(m, stan_data, iter = 2 * n_iter)
save(
  stan_fit,
  file = sprintf("../../data/fits/lda-%s.rda", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- sample-test ----
X_test_hat <- lda_pred(
  samples$beta,
  rowSums(X[test_ix, ]),
  stan_data$alpha
) %>%
  melt(
    varnames = c("iter", "ix", "v"),
    value.name = "predicted"
  )

X_test_hat$ix <- test_ix[X_test_hat$ix]
X_test <- X[test_ix, ] %>%
  melt(
    varnames = c("ix", "v"),
    value.name = "truth"
  )
X_test$ix <- test_ix[X_test$ix]

X_test_hat <- X_test_hat %>%
  left_join(X_test)

## ---- visualize-preds ----
ggplot(
  X_test_hat %>%
  filter(v %in% sample(1:ncol(X), 8)),
) +
  geom_point(
    aes(x = asinh(truth), y = asinh(predicted), col = as.factor(v)),
    alpha = 0.05, size = 0.4,
    position = position_jitter(w = 0.2, h = 0.2)
  ) +
  scale_color_brewer(palette = "Set3") +
  guides(
    col = guide_legend(override.aes = list(size = 2, alpha = 1))
  ) +
  coord_fixed()

ggplot(X_test_hat) +
  geom_histogram(
    aes(asinh(predicted - truth)),
    bins = 100
  )

ggplot(X_test_hat) +
  geom_histogram(
    aes(asinh(mean(truth) - truth)),
    bins = 100
  )

ggplot(X_test_hat) +
  geom_histogram(
    aes(asinh(median(truth) - truth)),
    bins = 100
  )

## ---- posterior-likelihoods ----
ll_test <- posterior_lda_loglik(
  X[test_ix, ],
  aperm(samples$beta, c(1, 3, 2)),
  stan_data$alpha
)

ll_train <- posterior_lda_loglik(
  X[-test_ix, ],
  aperm(samples$beta, c(1, 3, 2)),
  stan_data$alpha
)
