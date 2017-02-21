data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> D_test; // number of test docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma;
}

parameters {
  simplex[K] theta[D + D_test]; // topic mixtures
  simplex[V] beta[K]; // word dist for k^th topic
  int<lower=0> n_test[D_test, V]; // held out samples
}

model {
  for (d in 1:(D + D_test)) {
    theta[d] ~ dirichlet(alpha);
  }

  for (k in 1:K) {
    beta[k] ~ dirichlet(gamma);
  }

  for (d in 1:(D + D_test)) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];
    for (k in 1:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    if (d <= D) {
      n[d] ~ multinomial(eta);
    } else {
      n_test[d - D] ~ multinomial(eta);
    }
  }

  for (d in 1:D_test) {
    vector[V] eta_test;
    eta_test = beta[1] * theta_test[d, 1];
    for (k in 1:K) {
      eta_test = eta_test + beta[k] * theta[d, k];
    }
    n_test[d] ~ multinomial(eta_test);
  }
}

generated quantities {
  int<lower=0> n_sim[D + D_test, V]; // simulated word counts, for posterior checking

  for (d in 1:(D + D_test)) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    n_sim[d] = multinomial_rng(eta, sum(n[d]));
  }
}
