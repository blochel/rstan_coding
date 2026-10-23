


# packages ----------------------------------------------------------------

library(dplyr)
library(edenR)
library(lubridate)
library(readr)
library(rstan)




# nesting WOST initiation data ---------------------------------------------


url <- "https://raw.githubusercontent.com/weecology/EvergladesWadingBird/main/Nesting/woodstork_dates.csv"
bird_raw <- read_csv(url)



# A south Florida water year begins in May of the preceding year designation
# Robert Kadlec, John Juston,2026, An Everglades water quality recovery model,
# Ecological Engineering, Volume 224, https://doi.org/10.1016/j.ecoleng.2025.107880.



bird_data <- bird_raw |>
  mutate(julian = yday(flight_first_nest_obs),
         # Water year starts in May of preceding year
         # example: water year 1998 = May 1997 to April 1998
         water_year = ifelse(month(flight_first_nest_obs) >= 5, 
                             year(flight_first_nest_obs) + 1, 
                             year(flight_first_nest_obs))) 


# Take earliest nest initiation date per year, convert to Julian day
bird_summary <- bird_data |>
  #only the first observation (less data)
  group_by(year) |>
  summarise(y = min(julian, na.rm = TRUE))

#watchout for tricky water year issues later! 
bird_data |> 
  mutate(correct_Q = year == water_year) |> 
  filter(correct_Q != TRUE)


# baseline model  ---------------------------------------------------------



y <- bird_summary$y

stan_data <- list(
  N = length(y),
  y = y
)

fit <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 42
)

print(fit, pars = c("alpha", "sigma"))



# posterior checks  - baseline --------------------------------------------


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



# Forecasting    -  baseline ----------------------------------------------

# Hold out last year for forecasting
train <- bird_summary %>% filter(water_year < max(water_year))
test  <- bird_summary %>% filter(water_year == max(water_year))

y_train <- train$y

stan_data <- list(
  N = length(y_train),
  y = y_train
)

fit <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 42
)

# Extract posterior predictive for one new observation
samples    <- extract(fit)
y_pred     <- rnorm(4000, mean = samples$alpha, sd = samples$sigma)

# Convert everything back to dates
pred_dates  <- as.Date(paste0(test$water_year, "-01-01")) + y_pred - 1
actual_date <- as.Date(paste0(test$water_year, "-01-01")) + test$y - 1

# Days off — median prediction vs actual
days_off <- round(median(y_pred) - test$y)
sign_label <- ifelse(days_off > 0, "days late", "days early")

# Plot
hist(pred_dates,
     breaks = 30,
     col    = "steelblue",
     border = "white",
     main   = paste("Forecast vs Actual \n Baseline model", test$water_year),
     xlab   = "Nest Initiation Date")
abline(v = as.numeric(actual_date), col = "red", lwd = 2)
legend("topright",
       legend = c("Predicted", 
                  paste("Actual —", format(actual_date, "%b %d")),
                  paste(abs(days_off), sign_label)),
       col    = c("steelblue", "red", NA),
       lwd    = c(NA, 2, NA))



# forecast for next year     -   baseline ---------------------------------

# 2027 forecast — no actual date yet
y_pred_2027 <- rnorm(4000, mean = samples$alpha, sd = samples$sigma)

# Compute quantiles on numeric values first, then convert to dates
ci_95       <- as.Date("2027-01-01") + quantile(y_pred_2027, probs = c(0.025, 0.975)) - 1
ci_50       <- as.Date("2027-01-01") + quantile(y_pred_2027, probs = c(0.25,  0.75))  - 1
median_pred <- as.Date("2027-01-01") + median(y_pred_2027) - 1

pred_dates_2027 <- as.Date("2027-01-01") + y_pred_2027 - 1

# Plot
hist(pred_dates_2027,
     breaks = 30,
     col    = "steelblue",
     border = "white",
     main   = "2027 WOST Nest Initiation Forecast \n Baseline Model",
     xlab   = "Predicted Nest Initiation Date")
abline(v = as.numeric(median_pred), col = "darkred", lwd = 3)
abline(v = as.numeric(ci_50),       col = "darkred", lwd = 2, lty = 2)
abline(v = as.numeric(ci_95),       col = "darkred", lwd = 2, lty = 3)
legend("topright",
       legend = c(paste("Median —", format(median_pred, "%b %d")),
                  paste("50% CI:", format(ci_50[1], "%b %d"), "to", format(ci_50[2], "%b %d")),
                  paste("95% CI:", format(ci_95[1], "%b %d"), "to", format(ci_95[2], "%b %d"))),
       col = "darkred", lwd = c(2, 1, 1), lty = c(1, 2, 3))


