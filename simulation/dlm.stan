data {
  int<lower=1> T;          // number of time points
  int<lower=1> K;          // number of series (variables)
  vector[T] y;             // vector of observations
  matrix[T, K] x;          // matrix of predictors
}

parameters {
  matrix[T, K] state;      // state matrix
  cov_matrix[K] W;         // covariance matrix for system evolution
  real<lower=0> sigma;     // observation noise standard deviation
}

model {
  // Priors
  for (k in 1:K) {
    state[1, k] ~ normal(0, 10); // prior on initial state
  }
  
  W ~ inv_wishart(10, diag_matrix(rep_vector(10, K))); // prior on W
  sigma ~ inv_gamma(1, 1); // prior on observation noise 

  // State evolution
  for (t in 2:T) {
    state[t,] ~ multi_normal(state[t-1,], W); // state evolution
  }

  // Observations
  for (t in 1:T) {
    y[t] ~ normal(dot_product(state[t,], x[t,]), sigma); // observation model
  }
}
