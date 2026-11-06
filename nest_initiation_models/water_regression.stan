data {
  int<lower=0> N;             // number of years
  vector[N] y;                // nest initiation (day of water year)
  int<lower=0> K;             // number of covariates
  matrix[N, K] X;             // design matrix (standardized covariates)
}

parameters {
  real alpha;                 // intercept
  vector[K] beta;             // covariate coefficients
  real<lower=0> sigma;        // residual SD
}

model {
  // Priors
  alpha ~ normal(306, 20);    // same as baseline: ~mean nesting date in WY days
  beta  ~ normal(0, 10);      // weakly informative; shrinks toward no effect
  sigma ~ exponential(1);

  // Likelihood
  y ~ normal(alpha + X * beta, sigma);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;          // for LOO-IC comparison with baseline
  for (i in 1:N) {
    y_rep[i]   = normal_rng(alpha + X[i] * beta, sigma);
    log_lik[i] = normal_lpdf(y[i] | alpha + X[i] * beta, sigma);
  }
}
