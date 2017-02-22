
rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), alpha)
  x / sum(x)
}

#' Get LDA predictive distribution, when theta[i] is unknown
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
#' @examples
#' lda_pred(samples$beta, sum(X[test_ix[1], ]), stan_data$alpha)
lda_pred <- function(beta_samples, n, alpha) {
  n_samples <- nrow(beta_samples)
  V <- dim(beta_samples)[3]
  preds <- matrix(0, n_samples, V)

  for (i in seq_len(n_samples)) {
    theta <- rdirichlet(alpha)
    preds[i, ] <- rmultinom(1, n, prob = t(beta_samples[i,, ]) %*% theta)
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
lda_loglik_data <- function(N, beta, theta) {
  logliks <- seq_len(nrow(N))
  for (i in seq_len(nrow(N))) {
    logliks[i] <- lda_loglik_sample(N[i, ], beta, theta)
  }
  logliks
}

