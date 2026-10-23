library(rstan)

# Your data — replace with your actual nest initiation dates (e.g., day of year)
y <- c(120, 115, 118, 122, 110)  # example: day of year











# baseline model  ---------------------------------------------------------

stan_data <- list(
  N = length(y),
  y = y
)

fit <- stan(
  file = "nest_initiation_models/baseline.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 42
)

print(fit)



# posterior checks  -------------------------------------------------------


# Summary of posterior distributions
print(fit, pars = c("alpha", "sigma"))

# Trace plots — check chains are mixing well
traceplot(fit, pars = c("alpha", "sigma"))

# Extract posterior samples and replicates 
samples <- extract(fit)
y_rep <- extract(fit)$y_rep


# Plot posterior for alpha
hist(samples$alpha, 
     main  = "Posterior: alpha (mean nest date)", 
     xlab  = "Day of Year",
     col   = "steelblue",
     border = "white")

# Plot posterior for sigma
hist(samples$sigma,
     main  = "Posterior: sigma",
     xlab  = "SD",
     col   = "coral",
     border = "white")

# Posterior predictive check — does y_rep look like your real data?
# Dark blue line (y) — the kernel density of your actual observed nest initiation dates
# Light blue lines (y_rep) — 50 simulated datasets drawn from your fitted posterior
bayesplot::ppc_dens_overlay(y, y_rep[1:50, ])


