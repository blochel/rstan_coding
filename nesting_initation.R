


# packages ----------------------------------------------------------------

library(dplyr)
library(edenR)
library(ggplot2)
library(lubridate)
library(purrr)
library(readr)
library(rstan)
library(tibble)




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

# hydro data --------------------------------------------------------------

get_data_water <- function(eden_path = "WaterData", update = FALSE, max_age_days = 7) {
  
  water_file <- file.path(eden_path, "eden_covariates.csv")
  
  # Auto-trigger update if file is missing or older than max_age_days
  needs_update <- update ||
    !file.exists(water_file) || {
      stamp <- file.info(water_file)$mtime
      is.na(stamp) || difftime(Sys.time(), stamp, units = "days") > max_age_days
    }
  
  if (needs_update) {
    cat("📥 Downloading fresh EDEN water data...\n")
    edenR::update_water(eden_path)
    water <- get_eden_covariates(eden_path = eden_path, level = "subregions") |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "all")) |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "wcas")) |>
      dplyr::select(year, region = Name, variable, value) |>
      as.data.frame() |>
      dplyr::select(-geometry) |>
      pivot_wider(names_from = "variable", values_from = "value") |>
      mutate(year = as.integer(year)) |>
      arrange(year, region)
    
    write_csv(water, water_file)
    cat(glue::glue(
      "✓ Water data updated and saved to {water_file}\n",
      "  Years: {min(water$year)}-{max(water$year)}\n",
      "  Regions: {paste(unique(water$region), collapse = ', ')}\n"
    ))
    
  } else {
    water <- read_csv(water_file, show_col_types = FALSE)
    # 
    age_days <- round(as.numeric(difftime(Sys.time(),
                                          file.info(water_file)$mtime,
                                          units = "days")), 1)
    cat(glue::glue(
      "✓ Water data loaded from: {water_file}\n",
      "  Years: {min(water$year)}-{max(water$year)}\n",
      "  Regions: {paste(unique(water$region), collapse = ', ')}\n",
      "  Last modified: {format(file.info(water_file)$mtime, '%Y-%m-%d')} ({age_days} days old)\n"
    ))
  }
  
  return(water)
}


water <- get_data_water(eden_path = "WaterData", update = FALSE)

# Inspect what regions and variables you have
water |> distinct(region)
water |> names()


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
train <- bird_summary %>% filter(year < max(year))
test  <- bird_summary %>% filter(year == max(year))

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
pred_dates  <- as.Date(paste0(test$year, "-01-01")) + y_pred - 1
actual_date <- as.Date(paste0(test$year, "-01-01")) + test$y - 1

pred_df <- tibble(date = pred_dates)

ggplot(pred_df, aes(x = date)) +
  geom_histogram(aes(y = after_stat(density), fill = "Predicted"),
                 bins   = 30,
                 color  = "white") +
  geom_vline(aes(xintercept = as.numeric(actual_date), color = "Actual"),
             linewidth = 1.2) +
  annotate("text",
           x     = actual_date + 3,
           y     = Inf,
           label = paste(abs(days_off), sign_label),
           hjust = 0, vjust = 1.5,
           color = "red", size = 3.5) +
  scale_fill_manual(name   = NULL, values = c("Predicted" = "steelblue")) +
  scale_color_manual(name  = NULL, values = c("Actual"    = "red")) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
  labs(title = paste("Forecast vs Actual\nBaseline Model", test$year),
       x     = "Nest Initiation Date",
       y     = "Density") +
  theme_minimal() +
  theme(legend.position   = "top",
        legend.direction  = "horizontal",
        plot.title        = element_text(hjust = 0.5, face = "bold"))


# forecast for next year     -   baseline ---------------------------------

# 2027 forecast — no actual date yet
y_pred_2027 <- rnorm(4000, mean = samples$alpha, sd = samples$sigma)

# Compute quantiles on numeric values first, then convert to dates
ci_95       <- as.Date("2027-01-01") + quantile(y_pred_2027, probs = c(0.025, 0.975)) - 1
ci_50       <- as.Date("2027-01-01") + quantile(y_pred_2027, probs = c(0.25,  0.75))  - 1
median_pred <- as.Date("2027-01-01") + median(y_pred_2027) - 1

pred_dates_2027 <- as.Date("2027-01-01") + y_pred_2027 - 1
pred_df_2027 <- tibble(date = pred_dates_2027)

ggplot(pred_df_2027, aes(x = date)) +
  geom_histogram(aes(y = after_stat(density), fill = "Predicted"),
                 bins  = 30,
                 color = "white") +
  geom_vline(aes(xintercept = as.numeric(median_pred),  color = "Median"),    linewidth = 1.2, linetype = "solid") +
  geom_vline(aes(xintercept = as.numeric(ci_50["25%"]), color = "50% CI"),    linewidth = 0.9, linetype = "dashed") +
  geom_vline(aes(xintercept = as.numeric(ci_50["75%"]), color = "50% CI"),    linewidth = 0.9, linetype = "dashed") +
  geom_vline(aes(xintercept = as.numeric(ci_95["2.5%"]),  color = "95% CI"),  linewidth = 0.9, linetype = "dotted") +
  geom_vline(aes(xintercept = as.numeric(ci_95["97.5%"]), color = "95% CI"),  linewidth = 0.9, linetype = "dotted") +
  scale_fill_manual(name   = NULL, values = c("Predicted" = "steelblue")) +
  scale_color_manual(name  = NULL,
                     values = c("Median" = "darkred", "50% CI" = "darkred", "95% CI" = "darkred"),
                     labels = c(
                       "Median" = paste("Median —",   format(median_pred,   "%b %d")),
                       "50% CI" = paste("50% CI:",    format(ci_50[1], "%b %d"), "to", format(ci_50[2], "%b %d")),
                       "95% CI" = paste("95% CI:",    format(ci_95[1], "%b %d"), "to", format(ci_95[2], "%b %d"))
                     )) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
  labs(title = "2027 WOST Nest Initiation Forecast\nBaseline Model",
       x     = "Predicted Nest Initiation Date",
       y     = "Density") +
  theme_minimal() +
  theme(legend.position  = "top",
        legend.direction = "horizontal",
        plot.title       = element_text(hjust = 0.5, face = "bold"))











# testing where model fails -----------------------------------------------




# ── LOO expanding window cross-validation ────────────────────────────────────

# Minimum years needed to fit a model before we start forecasting
min_train <- 5

years <- sort(unique(bird_summary$year))
loo_years <- years[years >= years[min_train + 1]]  # first forecastable year

run_loo_fold <- function(hold_out_year) {
  train <- bird_summary |> filter(year < hold_out_year)
  test  <- bird_summary |> filter(year == hold_out_year)
  
  stan_data <- list(N = nrow(train), y = train$y)
  
  fit <- stan(
    file   = "nest_initiation_models/baseline.stan",
    data   = stan_data,
    chains = 4,
    iter   = 2000,
    warmup = 1000,
    seed   = 42,
    refresh = 0  # suppress per-fold output
  )
  
  samp   <- extract(fit)
  y_pred <- rnorm(4000, mean = samp$alpha, sd = samp$sigma)
  
  tibble(
    year        = hold_out_year,
    actual      = test$y,
    pred_median = median(y_pred),
    pred_mean   = mean(y_pred),
    pred_sd     = sd(y_pred),
    ci50_lo     = quantile(y_pred, 0.25),
    ci50_hi     = quantile(y_pred, 0.75),
    ci95_lo     = quantile(y_pred, 0.025),
    ci95_hi     = quantile(y_pred, 0.975),
    residual    = actual - pred_median  # positive = nesting earlier than predicted
  )
}

# Run all folds — this will take a few minutes
loo_results <- map_dfr(loo_years, run_loo_fold)

# Flag off years
loo_results <- loo_results |>
  mutate(off_year = abs(residual) > 14,
         direction = case_when(
           residual >  14 ~ "Earlier than predicted",
           residual < -14 ~ "Later than predicted",
           TRUE           ~ "Within 14 days"
         ))

# ── LOO residual plot ─────────────────────────────────────────────────────────

ggplot(loo_results, aes(x = year, y = residual)) +
  geom_hline(yintercept = 0,   color = "grey50", linetype = "dashed") +
  geom_hline(yintercept =  14, color = "coral",  linetype = "dotted", linewidth = 0.8) +
  geom_hline(yintercept = -14, color = "coral",  linetype = "dotted", linewidth = 0.8) +
  geom_ribbon(aes(ymin = ci95_lo - pred_median, ymax = ci95_hi - pred_median),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = ci50_lo - pred_median, ymax = ci50_hi - pred_median),
              fill = "steelblue", alpha = 0.30) +
  geom_point(aes(color = direction), size = 3) +
  geom_line(color = "grey40", linewidth = 0.5) +
  ggrepel::geom_label_repel(
    data  = filter(loo_results, off_year),
    aes(label = year, color = direction),
    size  = 3, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Earlier than predicted" = "steelblue",
    "Later than predicted"   = "coral",
    "Within 14 days"         = "grey40"
  )) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(title    = "LOO Expanding Window — Baseline Model Residuals",
       subtitle = "Residual = Actual − Predicted median  |  Dotted lines = ±14 day threshold",
       x        = "Year",
       y        = "Residual (days)",
       color    = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"))

# ── Pull off years for hydro investigation ────────────────────────────────────

off_years <- loo_results |> filter(off_year) |> pull(year)
cat("Off years (|residual| > 14 days):", paste(off_years, collapse = ", "), "\n")

# ── Pick a region and join to LOO results ─────────────────────────────────────
# "all" = system-wide; swap for a subregion (e.g. "ENP", "WCA3A") as needed
water_system <- water |>
  filter(region == "all") |>
  mutate(year = as.integer(year))

# Join to LOO residuals — water year lags one year behind nesting year
# e.g. nesting year 2005 is driven by water year Oct 2004 – Apr 2005
# The EDEN covariates are already annual summaries keyed to year,
# so join on year directly and revisit lag structure once you see the variables
loo_hydro <- loo_results |>
  left_join(water_system, by = "year")

# ── Quick look at what came in ────────────────────────────────────────────────
glimpse(loo_hydro)

# ── Hydro pattern plot for off years ─────────────────────────────────────────
# Pivot all numeric hydro columns long so we can facet by variable
hydro_vars <- water_system |>
  select(where(is.numeric), -year) |>
  names()

loo_hydro_long <- loo_hydro |>
  select(year, residual, off_year, direction, all_of(hydro_vars)) |>
  pivot_longer(cols = all_of(hydro_vars),
               names_to  = "variable",
               values_to = "value")

# Scatter: each hydro variable vs residual, off years labelled
ggplot(loo_hydro_long, aes(x = value, y = residual)) +
  geom_hline(yintercept = c(-14, 0, 14),
             linetype = c("dotted", "dashed", "dotted"),
             color    = "grey50") +
  geom_point(aes(color = direction), size = 2, alpha = 0.8) +
  ggrepel::geom_label_repel(
    data  = filter(loo_hydro_long, off_year),
    aes(label = year, color = direction),
    size  = 2.5, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Earlier than predicted" = "steelblue",
    "Later than predicted"   = "coral",
    "Within 14 days"         = "grey60"
  )) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(title    = "Hydro Covariates vs. LOO Residuals",
       subtitle = "Dotted lines = ±14 day threshold",
       x        = "Covariate value",
       y        = "Residual (days)  —  positive = earlier than predicted",
       color    = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"))