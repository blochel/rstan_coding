# fit_models.R ------------------------------------------------------------
# Loads data, fits baseline and water regression Stan models,
# and saves all fit objects and processed data to outputs/

# packages ----------------------------------------------------------------
library(dplyr)
library(edenR)
library(loo)
library(lubridate)
library(purrr)
library(readr)
library(rstan)
library(tibble)
library(tidyr)
extract <- rstan::extract

dir.create("outputs", showWarnings = FALSE)

# nesting WOST initiation data --------------------------------------------
# A south Florida water year begins in May of the preceding year designation
# Robert Kadlec, John Juston,2026, An Everglades water quality recovery model,
# Ecological Engineering, Volume 224, https://doi.org/10.1016/j.ecoleng.2025.107880.
url      <- "https://raw.githubusercontent.com/weecology/EvergladesWadingBird/main/Nesting/woodstork_dates.csv"
bird_raw <- read_csv(url, show_col_types = FALSE)

bird_data <- bird_raw |>
  mutate(
    water_year       = ifelse(month(flight_first_nest_obs) >= 6,
                              year(flight_first_nest_obs) + 1,
                              year(flight_first_nest_obs)),
    water_year_start = as.Date(paste0(water_year - 1, "-05-01")),
    julian_wy        = as.numeric(flight_first_nest_obs - water_year_start) + 1
  )

bird_summary <- bird_data |>
  # filter(colony != 'corkscrew_swamp', colony != 'butcher') |>
  group_by(water_year) |>
  summarise(y = min(julian_wy, na.rm = TRUE)) |>
  rename(year = water_year)

# hydro data --------------------------------------------------------------
get_data_water <- function(eden_path = "WaterData", update = FALSE, max_age_days = 7) {
  water_file   <- file.path(eden_path, "eden_covariates.csv")
  needs_update <- update ||
    !file.exists(water_file) || {
      stamp <- file.info(water_file)$mtime
      is.na(stamp) || difftime(Sys.time(), stamp, units = "days") > max_age_days
    }
  
  if (needs_update) {
    cat("📥 Downloading fresh EDEN water data...\n")
    edenR::update_water(eden_path)
    water <- get_eden_covariates(eden_path = eden_path, level = "subregions") |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "all"))    |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "wcas"))   |>
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

water        <- get_data_water(eden_path = "WaterData", update = FALSE)
water_system <- water |> filter(region == "all") |> mutate(year = as.integer(year))

# covariate setup ---------------------------------------------------------
covariate_names <- c("pre_recession", "init_depth")

model_data <- bird_summary |>
  left_join(water_system, by = "year") |>
  filter(if_all(all_of(covariate_names), ~ !is.na(.)))

# baseline — full data ----------------------------------------------------
stan_data_baseline <- list(N = nrow(model_data), y = model_data$y)

fit_baseline_full <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = stan_data_baseline,
  chains = 4, iter = 2000, warmup = 1000, seed = 42
)

# baseline — train/test split ---------------------------------------------
train_base <- bird_summary |> filter(year < max(year))
test_base  <- bird_summary |> filter(year == max(year))

fit_baseline_forecast <- stan(
  file   = "nest_initiation_models/baseline.stan",
  data   = list(N = nrow(train_base), y = train_base$y),
  chains = 4, iter = 2000, warmup = 1000, seed = 42
)

# regression — full data --------------------------------------------------
X_full   <- scale(model_data |> select(all_of(covariate_names)))
scale_center_full <- attr(X_full, "scaled:center")
scale_sd_full     <- attr(X_full, "scaled:scale")

stan_data_reg_full <- list(
  N = nrow(model_data),
  y = model_data$y,
  K = length(covariate_names),
  X = X_full
)

fit_reg_full <- stan(
  file   = "nest_initiation_models/water_regression.stan",
  data   = stan_data_reg_full,
  chains = 4, iter = 2000, warmup = 1000, seed = 42
)

# regression — train/test split -------------------------------------------
train_reg <- model_data |> filter(year < max(year))
test_reg  <- model_data |> filter(year == max(year))

X_train_raw  <- train_reg |> select(all_of(covariate_names))
scale_center <- colMeans(X_train_raw)
scale_sd     <- apply(X_train_raw, 2, sd)
X_train      <- scale(X_train_raw, center = scale_center, scale = scale_sd)

fit_reg_forecast <- stan(
  file   = "nest_initiation_models/water_regression.stan",
  data   = list(N = nrow(train_reg), y = train_reg$y,
                K = length(covariate_names), X = X_train),
  chains = 4, iter = 2000, warmup = 1000, seed = 42
)

# LOO expanding window — baseline -----------------------------------------
min_train <- 10
years     <- sort(unique(bird_summary$year))
loo_years <- years[years >= years[min_train + 1]]

run_loo_fold_baseline <- function(hold_out_year) {
  tr  <- bird_summary |> filter(year < hold_out_year)
  te  <- bird_summary |> filter(year == hold_out_year)
  fit <- stan(
    file    = "nest_initiation_models/baseline.stan",
    data    = list(N = nrow(tr), y = tr$y),
    chains  = 4, iter = 2000, warmup = 1000, seed = 42, refresh = 0
  )
  samp   <- extract(fit)
  y_pred <- rnorm(4000, mean = samp$alpha, sd = samp$sigma)
  tibble(
    year        = hold_out_year,
    actual      = te$y,
    pred_median = median(y_pred),
    pred_mean   = mean(y_pred),
    pred_sd     = sd(y_pred),
    ci50_lo     = quantile(y_pred, 0.25),
    ci50_hi     = quantile(y_pred, 0.75),
    ci95_lo     = quantile(y_pred, 0.025),
    ci95_hi     = quantile(y_pred, 0.975),
    residual    = actual - pred_median
  )
}

loo_results_baseline <- map_dfr(loo_years, run_loo_fold_baseline) |>
  mutate(
    off_year  = abs(residual) > 14,
    direction = case_when(
      residual >  14 ~ "Later than predicted",
      residual < -14 ~ "Earlier than predicted",
      TRUE           ~ "Within 14 days"
    )
  )

# save outputs ------------------------------------------------------------
saveRDS(bird_summary,           "outputs/bird_summary.rds")
saveRDS(water,                  "outputs/water.rds")
saveRDS(water_system,           "outputs/water_system.rds")
saveRDS(model_data,             "outputs/model_data.rds")
saveRDS(covariate_names,        "outputs/covariate_names.rds")
saveRDS(fit_baseline_full,      "outputs/fit_baseline_full.rds")
saveRDS(fit_baseline_forecast,  "outputs/fit_baseline_forecast.rds")
saveRDS(train_base,             "outputs/train_base.rds")
saveRDS(test_base,              "outputs/test_base.rds")
saveRDS(fit_reg_full,           "outputs/fit_reg_full.rds")
saveRDS(fit_reg_forecast,       "outputs/fit_reg_forecast.rds")
saveRDS(train_reg,              "outputs/train_reg.rds")
saveRDS(test_reg,               "outputs/test_reg.rds")
saveRDS(scale_center,           "outputs/scale_center.rds")
saveRDS(scale_sd,               "outputs/scale_sd.rds")
saveRDS(scale_center_full,      "outputs/scale_center_full.rds")
saveRDS(scale_sd_full,          "outputs/scale_sd_full.rds")
saveRDS(loo_results_baseline,   "outputs/loo_results_baseline.rds")

cat("✅ fit_models.R complete\n")