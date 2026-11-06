# evaluate_models.R -------------------------------------------------------
# Evaluate existing fitted models and run expanding-window water forecasts.
#
# Required: run fit_models.R first.
# Full baseline and regression fits must use model_data in the same row order.
# Saved scale_center/scale_sd must correspond to fit_reg_forecast.
#
# Expanding-window evaluation and PSIS-LOO are saved separately.
# Water covariates must be available at the intended forecast date.

library(dplyr)
library(loo)
library(rstan)
library(tibble)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

read_output <- function(name) {
  path <- file.path("outputs", paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run fit_models.R first.")
  }
  readRDS(path)
}

save_output <- function(object, name) {
  saveRDS(object, file.path("outputs", paste0(name, ".rds")))
}

# Load -------------------------------------------------------------------
bird_summary          <- read_output("bird_summary")
water_system          <- read_output("water_system")
model_data            <- read_output("model_data")
covariate_names       <- read_output("covariate_names")
fit_baseline_full     <- read_output("fit_baseline_full")
fit_baseline_forecast <- read_output("fit_baseline_forecast")
fit_reg_full          <- read_output("fit_reg_full")
fit_reg_forecast      <- read_output("fit_reg_forecast")
train_base            <- read_output("train_base")
test_base             <- read_output("test_base")
train_reg             <- read_output("train_reg")
test_reg              <- read_output("test_reg")
scale_center          <- read_output("scale_center")
scale_sd              <- read_output("scale_sd")
loo_results_baseline  <- read_output("loo_results_baseline")

# Validation -------------------------------------------------------------
required_columns <- c("year", "y", covariate_names)

if (!all(required_columns %in% names(model_data))) {
  stop("model_data is missing required columns.")
}

if (anyDuplicated(model_data$year)) {
  stop("model_data must have exactly one row per water year.")
}

if (anyDuplicated(water_system$year)) {
  stop("water_system must have exactly one row per water year.")
}

if (anyDuplicated(loo_results_baseline$year)) {
  stop("loo_results_baseline has duplicated evaluation years.")
}

if (nrow(test_base) != 1L || nrow(test_reg) != 1L) {
  stop("The saved single-year forecast test sets must each have one row.")
}

numeric_columns <- c("y", covariate_names)
if (!all(vapply(model_data[numeric_columns], is.numeric, logical(1)))) {
  stop("Response and covariates must be numeric.")
}

water_eval_data <- model_data |>
  filter(if_all(all_of(numeric_columns), ~ is.finite(.x))) |>
  arrange(year)

if (nrow(water_eval_data) <= 10L) {
  stop("Need more than 10 complete years for expanding-window evaluation.")
}

# Helpers ----------------------------------------------------------------
classify_residuals <- function(data) {
  data |>
    mutate(
      off_year = abs(residual) > 14,
      direction = case_when(
        residual > 14  ~ "Later than predicted",
        residual < -14 ~ "Earlier than predicted",
        TRUE           ~ "Within 14 days"
      )
    ) |>
    arrange(year)
}

standardize_new <- function(data, center, spread) {
  # Preserve the model's covariate order when vectors are named.
  if (!is.null(names(center))) {
    center <- center[covariate_names]
  }
  if (!is.null(names(spread))) {
    spread <- spread[covariate_names]
  }
  
  if (length(center) != length(covariate_names) ||
      length(spread) != length(covariate_names) ||
      any(!is.finite(center)) ||
      any(!is.finite(spread)) ||
      any(spread <= 0)) {
    stop("Invalid covariate scaling parameters.")
  }
  
  X <- as.matrix(data[, covariate_names, drop = FALSE])
  
  if (any(!is.finite(X))) {
    stop("Prediction covariates contain missing or non-finite values.")
  }
  
  scale(X, center = center, scale = spread)
}

predict_baseline <- function(fit, seed) {
  set.seed(seed)
  draws <- rstan::extract(fit, pars = c("alpha", "sigma"))
  rnorm(length(draws$alpha), draws$alpha, draws$sigma)
}

predict_water <- function(fit, X, seed) {
  if (nrow(X) != 1L) {
    stop("predict_water expects one prediction row.")
  }
  
  set.seed(seed)
  draws <- rstan::extract(fit, pars = c("alpha", "beta", "sigma"))
  beta <- matrix(
    draws$beta,
    nrow = length(draws$alpha),
    ncol = ncol(X)
  )
  
  mu <- draws$alpha + drop(beta %*% as.numeric(X))
  rnorm(length(mu), mu, draws$sigma)
}

summarize_prediction <- function(year, actual, predictions) {
  tibble(
    year = year,
    actual = actual,
    pred_median = median(predictions),
    pred_mean = mean(predictions),
    pred_sd = sd(predictions),
    ci50_lo = unname(quantile(predictions, 0.25)),
    ci50_hi = unname(quantile(predictions, 0.75)),
    ci95_lo = unname(quantile(predictions, 0.025)),
    ci95_hi = unname(quantile(predictions, 0.975)),
    residual = actual - median(predictions)
  )
}

report_forecast <- function(test, predictions, label) {
  anchor <- as.Date(paste0(test$year - 1, "-05-01"))
  actual_date <- anchor + test$y - 1
  predicted_date <- anchor + round(median(predictions)) - 1
  difference <- round(median(predictions) - test$y)
  
  timing <- if (difference > 0) {
    "days late"
  } else if (difference < 0) {
    "days early"
  } else {
    "days off"
  }
  
  cat(
    "\n--- ", label, " forecast ---\n",
    "Water year: ", test$year,
    " | Actual: ", format(actual_date, "%Y-%m-%d"),
    " | Predicted: ", format(predicted_date, "%Y-%m-%d"),
    " | Off by: ", abs(difference), " ", timing, "\n",
    sep = ""
  )
}

# Posterior summaries ----------------------------------------------------
cat("\n--- Baseline: full fit ---\n")
print(fit_baseline_full, pars = c("alpha", "sigma"))

cat("\n--- Water regression: full fit ---\n")
print(fit_reg_full, pars = c("alpha", "beta", "sigma"))

# PSIS-LOO ---------------------------------------------------------------
# Full fits must share observations AND observation order.
# This is distinct from the chronological expanding-window evaluation.
#
# If an older Stan model lacks log_lik, derive it from posterior draws.
# model_data$y is assumed to be the response used by both full fits.

get_log_lik <- function(fit, regression = FALSE) {
  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep
  
  if (ncol(y_rep) != nrow(model_data)) {
    stop(
      "Full fit response length differs from model_data. ",
      "Refit both full models on the same observations."
    )
  }
  
  ll <- tryCatch(
    loo::extract_log_lik(
      fit,
      parameter_name = "log_lik",
      merge_chains = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(ll)) {
    if (dim(ll)[3] != nrow(model_data)) {
      stop("log_lik does not match model_data.")
    }
    return(ll)
  }
  
  message("log_lik absent; deriving it from posterior draws.")
  
  parameter_names <- if (regression) {
    c("alpha", "beta", "sigma")
  } else {
    c("alpha", "sigma")
  }
  
  draws <- as.array(fit, pars = parameter_names)
  iterations <- dim(draws)[1]
  chains <- dim(draws)[2]
  
  # Derive likelihood from saved posterior predictive conditional means.
  # For the regression, use Stan's lp-independent mean if present;
  # otherwise require log_lik in the model to avoid assuming its scaling.
  if (regression) {
    stop(
      "Water model lacks log_lik. Add log_lik to generated quantities ",
      "and rerun fit_models.R before evaluating."
    )
  }
  
  ll <- array(
    NA_real_,
    dim = c(iterations, chains, nrow(model_data))
  )
  
  for (i in seq_len(nrow(model_data))) {
    ll[, , i] <- matrix(
      dnorm(
        model_data$y[i],
        mean = as.vector(draws[, , "alpha"]),
        sd = as.vector(draws[, , "sigma"]),
        log = TRUE
      ),
      nrow = iterations,
      ncol = chains
    )
  }
  
  ll
}

compute_loo <- function(fit, regression = FALSE) {
  ll <- get_log_lik(fit, regression)
  
  # Shift each observation's log likelihood for numerical stability.
  shifted <- ll
  for (i in seq_len(dim(ll)[3])) {
    shifted[, , i] <- ll[, , i] - max(ll[, , i])
  }
  
  r_eff <- loo::relative_eff(exp(shifted))
  loo::loo(ll, r_eff = r_eff)
}

# Do not let unavailable PSIS-LOO prevent chronological evaluation.
loo_baseline <- tryCatch(
  compute_loo(fit_baseline_full),
  error = function(e) {
    warning("Baseline PSIS-LOO skipped: ", conditionMessage(e))
    NULL
  }
)

loo_reg <- tryCatch(
  compute_loo(fit_reg_full, regression = TRUE),
  error = function(e) {
    warning("Water PSIS-LOO skipped: ", conditionMessage(e))
    NULL
  }
)

loo_comparison <- NULL

if (!is.null(loo_baseline)) {
  cat("\n--- Baseline PSIS-LOO ---\n")
  print(loo_baseline)
}

if (!is.null(loo_reg)) {
  cat("\n--- Water PSIS-LOO ---\n")
  print(loo_reg)
}

loo_comparison <- NULL

if (!is.null(loo_baseline) && !is.null(loo_reg)) {
  loo_comparison <- loo::loo_compare(
    x = list(
      baseline = loo_baseline,
      water = loo_reg
    )
  )
  
  cat("\n--- PSIS-LOO comparison ---\n")
  print(loo_comparison)
}

# Saved single-year forecasts -------------------------------------------
y_pred_base <- predict_baseline(fit_baseline_forecast, seed = 101)

X_test <- standardize_new(test_reg, scale_center, scale_sd)
y_pred_reg <- predict_water(fit_reg_forecast, X_test, seed = 102)

report_forecast(test_base, y_pred_base, "Baseline")
report_forecast(test_reg, y_pred_reg, "Water regression")

forecast_accuracy <- bind_rows(
  summarize_prediction(test_base$year, test_base$y, y_pred_base) |>
    mutate(model = "Baseline"),
  summarize_prediction(test_reg$year, test_reg$y, y_pred_reg) |>
    mutate(model = "Water regression")
)

# Expanding-window water evaluation -------------------------------------
# Use the baseline's evaluation years where water data are available.
# Every fold estimates scaling using only its training observations.
# Reuse the compiled Stan model from the saved water fit.

min_train <- 10L
evaluation_years <- sort(unique(loo_results_baseline$year))
water_stan_model <- rstan::get_stanmodel(fit_reg_full)

water_folds <- list()
skipped_folds <- list()

for (hold_out_year in evaluation_years) {
  train <- water_eval_data |> filter(year < hold_out_year)
  test <- water_eval_data |> filter(year == hold_out_year)
  
  if (nrow(test) != 1L || nrow(train) < min_train) {
    skipped_folds[[length(skipped_folds) + 1L]] <- tibble(
      year = hold_out_year,
      reason = if (nrow(test) != 1L) {
        "Missing complete water-model observation"
      } else {
        "Fewer than 10 complete training years"
      }
    )
    next
  }
  
  X_raw <- as.matrix(train[, covariate_names, drop = FALSE])
  fold_center <- colMeans(X_raw)
  fold_sd <- apply(X_raw, 2, sd)
  
  if (any(!is.finite(fold_sd)) || any(fold_sd <= 0)) {
    skipped_folds[[length(skipped_folds) + 1L]] <- tibble(
      year = hold_out_year,
      reason = "A training covariate has zero or invalid SD"
    )
    next
  }
  
  cat(
    "\nWater fold:", hold_out_year,
    "| Training years:", nrow(train), "\n"
  )
  
  X_train <- standardize_new(train, fold_center, fold_sd)
  X_holdout <- standardize_new(test, fold_center, fold_sd)
  
  fit_fold <- rstan::sampling(
    object = water_stan_model,
    data = list(
      N = nrow(train),
      y = train$y,
      K = length(covariate_names),
      X = X_train
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    seed = as.integer(hold_out_year),
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )
  
  # Report diagnostics for each fold.
  rstan::check_hmc_diagnostics(fit_fold)
  
  predictions <- predict_water(
    fit_fold,
    X_holdout,
    seed = as.integer(hold_out_year + 10000)
  )
  
  water_folds[[length(water_folds) + 1L]] <-
    summarize_prediction(hold_out_year, test$y, predictions) |>
    mutate(n_train = nrow(train))
}

if (length(water_folds) == 0L) {
  stop("No eligible expanding-window water folds were available.")
}

loo_results_regression <- bind_rows(water_folds) |>
  classify_residuals()

loo_results_baseline <- classify_residuals(loo_results_baseline)
skipped_water_folds <- bind_rows(skipped_folds)

# Comparisons restricted to shared evaluation years ----------------------
common_years <- intersect(
  loo_results_baseline$year,
  loo_results_regression$year
)

expanding_results <- bind_rows(
  loo_results_baseline |>
    filter(year %in% common_years) |>
    mutate(model = "Baseline"),
  loo_results_regression |>
    filter(year %in% common_years) |>
    mutate(model = "Water regression")
)

expanding_metrics <- expanding_results |>
  group_by(model) |>
  summarise(
    n_years = n(),
    MAE = mean(abs(residual)),
    RMSE = sqrt(mean(residual^2)),
    bias = mean(residual),
    within_14_days = mean(abs(residual) <= 14),
    .groups = "drop"
  )

cat("\n--- Expanding-window accuracy: shared years ---\n")
print(expanding_metrics)

if (nrow(skipped_water_folds) > 0L) {
  cat("\n--- Skipped water folds ---\n")
  print(skipped_water_folds)
}

# Hydro residual tables --------------------------------------------------
loo_hydro <- loo_results_baseline |>
  left_join(water_system, by = "year")

loo_hydro_regression <- loo_results_regression |>
  left_join(water_system, by = "year")

# Future forecasts -------------------------------------------------------
# Baseline uses the full fit, rather than the held-out-year forecast fit.
forecast_year <- 2027L
y_pred_2027_base <- predict_baseline(fit_baseline_full, seed = 202701)

# Do NOT automatically substitute 2026 covariates for a 2027 prediction.
#
# Optional input:
# outputs/future_water_covariates_2027.rds
#
# This must be a one-row data frame containing the named covariates,
# constructed for the intended forecast date. It may represent an explicit
# scenario; predictions are conditional on these supplied values.
#
# If provided, refit on all complete historical model_data with known scaling.

future_path <- "outputs/future_water_covariates_2027.rds"
y_pred_2027_reg <- NULL

if (file.exists(future_path)) {
  future_water <- readRDS(future_path)
  
  if (!is.data.frame(future_water) ||
      nrow(future_water) != 1L ||
      !all(covariate_names %in% names(future_water))) {
    stop("Future water covariates must be a one-row data frame.")
  }
  
  future_train <- water_eval_data |>
    filter(year < forecast_year)
  
  X_future_train_raw <- as.matrix(
    future_train[, covariate_names, drop = FALSE]
  )
  
  future_center <- colMeans(X_future_train_raw)
  future_sd <- apply(X_future_train_raw, 2, sd)
  
  X_future_train <- standardize_new(
    future_train, future_center, future_sd
  )
  X_future <- standardize_new(
    future_water, future_center, future_sd
  )
  
  fit_reg_future <- rstan::sampling(
    object = water_stan_model,
    data = list(
      N = nrow(future_train),
      y = future_train$y,
      K = length(covariate_names),
      X = X_future_train
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    seed = 202703,
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )
  
  rstan::check_hmc_diagnostics(fit_reg_future)
  y_pred_2027_reg <- predict_water(
    fit_reg_future, X_future, seed = 202702
  )
} else {
  message(
    "Skipping 2027 water forecast: no explicitly supplied future covariates."
  )
}

# Save -------------------------------------------------------------------
save_output(loo_baseline, "loo_baseline")
save_output(loo_reg, "loo_reg")
save_output(loo_comparison, "loo_comparison")
save_output(y_pred_base, "y_pred_base")
save_output(y_pred_reg, "y_pred_reg")
save_output(forecast_accuracy, "forecast_accuracy")
save_output(loo_results_baseline, "loo_results_baseline")
save_output(loo_results_regression, "loo_results_regression")
save_output(skipped_water_folds, "skipped_water_folds")
save_output(expanding_results, "expanding_results")
save_output(expanding_metrics, "expanding_metrics")
save_output(loo_hydro, "loo_hydro")
save_output(loo_hydro_regression, "loo_hydro_regression")
save_output(y_pred_2027_base, "y_pred_2027_base")
save_output(y_pred_2027_reg, "y_pred_2027_reg")
save_output(forecast_year, "forecast_year")

write.csv(
  loo_results_regression,
  "outputs/loo_results_regression.csv",
  row.names = FALSE
)
write.csv(
  expanding_metrics,
  "outputs/expanding_metrics.csv",
  row.names = FALSE
)

cat("\nevaluate_models.R complete.\n")