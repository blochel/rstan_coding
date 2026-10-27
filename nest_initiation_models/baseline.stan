//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> N;       // number of years
  vector[N] y;          // nest initiation day of year
}

parameters {
  real alpha;           // intercept (mean nest date)
  real<lower=0> sigma;  // standard deviation
}

model {
  // Priors
  alpha ~ normal(306, 20);  // ~March 1 = day 60, +/- 20 days                      (NOTE: WY starts May 1 - will need to be altered for water dynamics)
  sigma ~ exponential(1);    // weakly informative, forces positive

  // Likelihood
  y ~ normal(alpha, sigma);
}

generated quantities {
  vector[N] y_rep;           // posterior predictive samples
  for (i in 1:N)
    y_rep[i] = normal_rng(alpha, sigma);
}