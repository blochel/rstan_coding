# packages ----------------------------------------------------------------
library(dplyr)
library(edenR)
library(ggplot2)
library(lubridate)
library(purrr)
library(readr)
library(rstan)
library(tibble)
library(tidyr)
extract <- rstan::extract

# nesting WOST initiation data ---------------------------------------------
url <- "https://raw.githubusercontent.com/weecology/EvergladesWadingBird/main/Nesting/woodstork_dates.csv"
bird_raw <- read_csv(url)

# A south Florida water year begins in May of the preceding year designation
# Robert Kadlec, John Juston,2026, An Everglades water quality recovery model,
# Ecological Engineering, Volume 224, https://doi.org/10.1016/j.ecoleng.2025.107880.
bird_data <- bird_raw |>
  mutate(
    water_year       = ifelse(month(flight_first_nest_obs) >= 6,
                              year(flight_first_nest_obs) + 1,
                              year(flight_first_nest_obs)),
    water_year_start = as.Date(paste0(water_year - 1, "-05-01")),
    julian_wy        = as.numeric(flight_first_nest_obs - water_year_start) + 1
  )

bird_summary <- bird_data |>
  group_by(water_year) |>
  summarise(y = min(julian_wy, na.rm = TRUE)) |>
  rename(year = water_year)

# hydro data --------------------------------------------------------------
get_data_water <- function(eden_path = "WaterData", update = FALSE, max_age_days = 7) {
  water_file  <- file.path(eden_path, "eden_covariates.csv")
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
  } else {
    water <- read_csv(water_file, show_col_types = FALSE)
  }
  return(water)
}

water <- get_data_water(eden_path = "WaterData", update = FALSE)

# baseline model  ---------------------------------------------------------
y         <- bird_summary$y
stan_data <- list(N = length(y), y = y)

fit <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 42
)

# posterior checks  - baseline --------------------------------------------
print(fit, pars = c("alpha", "sigma"))
traceplot(fit, pars = c("alpha", "sigma"))

samples <- extract(fit)
y_rep   <- extract(fit)$y_rep

hist(samples$alpha,
     main = "Posterior: alpha (mean nest date)",
     xlab = "Day of Water Year", col = "steelblue", border = "white")

hist(samples$sigma,
     main = "Posterior: sigma",
     xlab = "SD", col = "coral", border = "white")

bayesplot::ppc_dens_overlay(y, y_rep[1:50, ])

# forecasting - baseline --------------------------------------------------
train   <- bird_summary |> filter(year < max(year))
test    <- bird_summary |> filter(year == max(year))
y_train <- train$y

stan_data <- list(N = length(y_train), y = y_train)

fit <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 42
)

samples <- extract(fit)
y_pred  <- rnorm(4000, mean = samples$alpha, sd = samples$sigma)

# water year anchor: May 1 of (test year - 1)
wy_anchor   <- as.Date(paste0(test$year - 1, "-05-01"))
pred_dates  <- wy_anchor + y_pred - 1
actual_date <- wy_anchor + test$y  - 1
median_date <- wy_anchor + median(y_pred) - 1
pred_df     <- tibble(date = pred_dates)

days_off   <- round(median(y_pred) - test$y)
sign_label <- ifelse(days_off > 0, "days late", "days early")

ggplot(pred_df, aes(x = date)) +
  geom_histogram(aes(y = after_stat(density), fill = "Predicted"),
                 bins = 30, color = "white") +
  geom_vline(xintercept = actual_date, color = "red",     linewidth = 1.2) +
  geom_vline(xintercept = median_date, color = "darkblue", linewidth = 1.2,
             linetype = "dashed") +
  annotate("text",
           x = actual_date + 3, y = Inf,
           label = format(actual_date, "%b %d"),
           hjust = 0, vjust = 3.2, color = "red", size = 3.5) +
  annotate("text",
           x = actual_date + 3, y = Inf,
           label = paste(abs(days_off), sign_label),
           hjust = 0, vjust = 1.5, color = "red", size = 3.5) +
  scale_fill_manual(name  = NULL, values = c("Predicted" = "steelblue")) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
  labs(title = paste("Forecast vs Actual\nBaseline Model", test$year),
       x = "Nest Initiation Date", y = "Density") +
  theme_minimal() +
  theme(legend.position  = "top",
        legend.direction = "horizontal",
        plot.title       = element_text(hjust = 0.5, face = "bold"))

# 2027 forecast -----------------------------------------------------------
y_pred_2027     <- rnorm(4000, mean = samples$alpha, sd = samples$sigma)
wy_anchor_2027  <- as.Date("2026-05-01")

ci_95           <- wy_anchor_2027 + quantile(y_pred_2027, probs = c(0.025, 0.975)) - 1
ci_50           <- wy_anchor_2027 + quantile(y_pred_2027, probs = c(0.25,  0.75))  - 1
median_pred     <- wy_anchor_2027 + median(y_pred_2027) - 1
pred_dates_2027 <- wy_anchor_2027 + y_pred_2027 - 1
pred_df_2027    <- tibble(date = pred_dates_2027)

ggplot(pred_df_2027, aes(x = date)) +
  geom_histogram(aes(y = after_stat(density), fill = "Predicted"),
                 bins = 30, color = "white") +
  geom_vline(aes(xintercept = as.numeric(median_pred),    color = "Median"),  linewidth = 1.2, linetype = "solid") +
  geom_vline(aes(xintercept = as.numeric(ci_50["25%"]),   color = "50% CI"), linewidth = 0.9, linetype = "dashed") +
  geom_vline(aes(xintercept = as.numeric(ci_50["75%"]),   color = "50% CI"), linewidth = 0.9, linetype = "dashed") +
  geom_vline(aes(xintercept = as.numeric(ci_95["2.5%"]),  color = "95% CI"), linewidth = 0.9, linetype = "dotted") +
  geom_vline(aes(xintercept = as.numeric(ci_95["97.5%"]), color = "95% CI"), linewidth = 0.9, linetype = "dotted") +
  scale_fill_manual(name  = NULL, values = c("Predicted" = "steelblue")) +
  scale_color_manual(name = NULL,
                     values = c("Median" = "darkred", "50% CI" = "darkred", "95% CI" = "darkred"),
                     labels = c(
                       "Median" = paste("Median —", format(median_pred, "%b %d")),
                       "50% CI" = paste("50% CI:", format(ci_50[1], "%b %d"),
                                        "to", format(ci_50[2], "%b %d")),
                       "95% CI" = paste("95% CI:", format(ci_95[1], "%b %d"),
                                        "to", format(ci_95[2], "%b %d"))
                     )) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
  labs(title = "2027 WOST Nest Initiation Forecast\nBaseline Model",
       x = "Predicted Nest Initiation Date", y = "Density") +
  theme_minimal() +
  theme(legend.position  = "top",
        legend.direction = "horizontal",
        plot.title       = element_text(hjust = 0.5, face = "bold"))

# LOO expanding window ----------------------------------------------------
min_train <- 5
years     <- sort(unique(bird_summary$year))
loo_years <- years[years >= years[min_train + 1]]

run_loo_fold <- function(hold_out_year) {
  train <- bird_summary |> filter(year < hold_out_year)
  test  <- bird_summary |> filter(year == hold_out_year)
  
  fit <- stan(
    file    = "nest_initiation_models/baseline.stan",
    data    = list(N = nrow(train), y = train$y),
    chains  = 4,
    iter    = 2000,
    warmup  = 1000,
    seed    = 42,
    refresh = 0
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
    residual    = actual - pred_median  # positive = nesting LATER than predicted
  )
}

loo_results <- map_dfr(loo_years, run_loo_fold)

loo_results <- loo_results |>
  mutate(
    off_year  = abs(residual) > 14,
    direction = case_when(
      residual >  14 ~ "Later than predicted",
      residual < -14 ~ "Earlier than predicted",
      TRUE           ~ "Within 14 days"
    )
  )

# LOO residual plot -------------------------------------------------------
ggplot(loo_results, aes(x = year, y = residual)) +
  geom_hline(yintercept =   0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept =  14, linetype = "dotted", color = "coral") +
  geom_hline(yintercept = -14, linetype = "dotted", color = "coral") +
  geom_ribbon(aes(ymin = ci95_lo - pred_median, ymax = ci95_hi - pred_median),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = ci50_lo - pred_median, ymax = ci50_hi - pred_median),
              fill = "steelblue", alpha = 0.30) +
  geom_line(color = "grey40", linewidth = 0.5) +
  geom_point(aes(color = direction), size = 3) +
  ggrepel::geom_label_repel(
    data = filter(loo_results, off_year),
    aes(label = year, color = direction),
    size = 3, show.legend = FALSE
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
       y        = "Residual (days)  —  positive = later than predicted",
       color    = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"))

# hydro covariate exploration ---------------------------------------------
off_years <- loo_results |> filter(off_year) |> pull(year)
cat("Off years (|residual| > 14 days):", paste(off_years, collapse = ", "), "\n")

water_system <- water |>
  filter(region == "all") |>
  mutate(year = as.integer(year))

loo_hydro <- loo_results |>
  left_join(water_system, by = "year")

hydro_vars <- water_system |>
  select(where(is.numeric), -year) |>
  names()

loo_hydro_long <- loo_hydro |>
  select(year, residual, off_year, direction, all_of(hydro_vars)) |>
  pivot_longer(cols      = all_of(hydro_vars),
               names_to  = "variable",
               values_to = "value")

ggplot(loo_hydro_long, aes(x = value, y = residual)) +
  geom_hline(yintercept =   0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept =  14, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = -14, linetype = "dotted", color = "grey50") +
  geom_smooth(method = "lm", color = "grey30", fill = "grey80", alpha = 0.4) +
  geom_point(aes(color = direction), size = 2, alpha = 0.8) +
  ggrepel::geom_label_repel(
    data = filter(loo_hydro_long, off_year),
    aes(label = year, color = direction),
    size = 2.5, show.legend = FALSE
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
       y        = "Residual (days)  —  positive = later than predicted",
       color    = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"))