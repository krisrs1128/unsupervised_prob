
rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), alpha)
  x / sum(x)
}

#' Get LDA predictions, when theta[i] is unknown
#'
#' This generates predictions for a new sample using posterior beta values and
#' also sampling new dirichlet theta[i]'s. I wouldn't expect this approach to do
#' particularly well, since it doesn't use any topic information.
#'
#' @param beta_samples [n_samples x K x V array] An array containing posterior
#'   samples of the beta topic distributions.
#' @param n [int] The total count in the sample that we are trying to make
#'   predictions for.
#' @param alpha [numeric vector] The alpha parameter used to simulate dirichlet
#'   theta[i]'s.
#' @return preds [n_samples x V matrix] A matrix whose columns are associated
#'   with microbes are rows are predictions given different sampled betas.
#' @examples
#' lda_pred(samples$beta, sum(X[test_ix[1], ]), stan_data$alpha)
lda_pred_sample <- function(beta_samples, n, alpha) {
  n_samples <- nrow(beta_samples)
  V <- dim(beta_samples)[3]
  preds <- matrix(0, n_samples, V)

  for (i in seq_len(n_samples)) {
    theta <- rdirichlet(alpha)
    preds[i, ] <- rmultinom(1, n, prob = t(beta_samples[i,, ]) %*% theta)
  }
  preds
}

#' Get LDA predictions across all samples
#'
#' This wraps lda_pred_sample, giving predictions over multiple samples.
#'
#' @param beta_samples [n_samples x K x V array] An array containing posterior
#'   samples of the beta topic distributions.
#' @param ns [length n vector] The total counts across samples that we are
#'   trying to make predictions for.
#' @param alpha [numeric vector] The alpha parameter used to simulate dirichlet
#'   theta[i]'s.
#' @return preds [n_samples x n x V matrix] A matrix whose columns are
#'   associated with microbes are rows are predictions given different sampled
#'   betas.
#' @examples
#' preds <- lda_pred(samples$beta, rowSums(X[test_ix, ]), stan_data$alpha)
#' preds[1,,] <- preds[1,, ] + runif(prod(dim(preds[1,,])), max = 0.5)
#' X[test_ix,] <- X[test_ix, ] + runif(prod(dim(X[test_ix, ])), max = 0.5)
#' plot(asinh(preds[1,,]), asinh(X[test_ix, ]))
lda_pred <- function(beta_samples, ns, alpha) {
  preds <- array(0, c(nrow(beta_samples), length(ns), dim(beta_samples)[3]))
  for (j in seq_along(ns)) {
    preds[, j, ] <- lda_pred_sample(beta_samples, ns[j], alpha)
  }
  preds
}

#' LDA loglikelihood with fixed beta and theta
#'
#' This computes the loglikelihood for the LDA model, given theta and beta.
#' Specifically, it takes
#' log((sum(n) choose n[1], ..., n[V])) + sum_{i, v} n[iv] * log(theta[i]^T beta[v, ])
#'
#' @param n [integer vector] A vector of counts n[i] for the sample on which to
#'   evaluate the loglikelihood.
#' @param beta [V x K numeric matrix] The topic probabilities distribution.
#' @param theta [length K vector] The topic mixture probabilities for the current document
#' @return loglik_sample The loglikelihood for the current sample, with the
#'   specified parameters.
lda_loglik_sample <- function(n, beta, theta) {
  multinom_coef <- dmultinom(n, sum(n), p = rep(1 / length(n), length(n)), log = TRUE)
  multinom_coef + sum(n * (beta %*% theta))
}

#' LDA loglikelihood with fixed beta and theta, across a dataset
#'
#' This wraps lda_loglik_sample to evaluate it across the entire dataset.
#' @param n [integer vector] A vector of counts n[i] for the sample on which to
#'   evaluate the loglikelihood.
#' @param beta [V x K numeric matrix] The topic probabilities distribution.
#' @param theta [length K vector] The topic mixture probabilities for the current document
#' @return loglik_sample The loglikelihood for the current sample, with the
#'   specified parameters.
lda_loglik <- function(N, beta, theta) {
  logliks <- vector(length = nrow(N))
  for (i in seq_len(nrow(N))) {
    logliks[i] <- lda_loglik_sample(N[i, ], beta, theta)
  }
  logliks
}

#' Compute the loglikelihood across posterior samples, for a single data point
#'
#' @param n [integer vector] A vector of counts n[i] for the sample on which to
#'   evaluate the loglikelihood.
#' @param beta_posterior [n_samples x V x K array] The topic probabilities
#'   distribution.
#' @param alpha [length K vector] The hyperparameter in the underlying topics.
#' @return logliks [length n_samples vector] The loglikelihood for the current
#'   sample, across samples from beta.
posterior_lda_loglik_sample <- function(n, beta_posterior, alpha) {
  n_posterior <- nrow(beta_posterior)
  logliks <- vector(length = n_posterior)
  for (i in seq_len(n_posterior)) {
    theta <- rdirichlet(alpha)
    logliks[i] <- lda_loglik_sample(n, beta_posterior[i,, ], rdirichlet(alpha))
  }
  logliks
}

#' Compute LDA likelihoods across samples
#'
#' This wraps posterior lda_loglik to evaluate likelihoods across posterior
#' samples of beta, but across multiple data points.
#'
#' @param n [length n data points integer vector] A vector of counts n[i] for
#'   the sample on which to evaluate the loglikelihood.
#' @param beta_posterior [n_samples x V x K array] The topic probabilities
#'   distribution.
#' @param alpha [length K vector] The hyperparameter in the underlying topics.
#' @return logliks [n_samples x n data points matrix] The loglikelihood for the
#'   current sample, across samples from beta.
#' @examples
#' logliks <- posterior_lda_loglik(
#'   X[-test_ix, ],
#'   aperm(samples$beta, c(1, 3, 2)),
#'   stan_data$alpha
#' )
#' hist(logliks, 150)
posterior_lda_loglik <- function(N, beta_posterior, alpha) {
  logliks <- matrix(0, nrow(N), nrow(beta_posterior))
  for (i in seq_len(nrow(N))) {
    logliks[i, ] <- posterior_lda_loglik_sample(N[i, ], beta_posterior, alpha)
  }
  logliks
}
