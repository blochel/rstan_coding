# plot_models.R -----------------------------------------------------------
# Create baseline and water-model figures from saved evaluation outputs.
#
# Required: run evaluate_models.R first.
# No predictions are redrawn here.

library(bayesplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rstan)
library(scales)
library(tibble)
library(tidyr)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

read_output <- function(name) {
  path <- file.path("outputs", paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run evaluate_models.R first.")
  }
  readRDS(path)
}

save_figure <- function(plot, filename, width = 8, height = 5) {
  ggsave(
    filename = file.path("figures", filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# Load -------------------------------------------------------------------
model_data            <- read_output("model_data")
water_system          <- read_output("water_system")
fit_baseline_full     <- read_output("fit_baseline_full")
fit_reg_full          <- read_output("fit_reg_full")
test_base             <- read_output("test_base")
test_reg              <- read_output("test_reg")
y_pred_base           <- read_output("y_pred_base")
y_pred_reg            <- read_output("y_pred_reg")
loo_results_baseline  <- read_output("loo_results_baseline")
loo_results_regression <- read_output("loo_results_regression")
expanding_results     <- read_output("expanding_results")
loo_hydro             <- read_output("loo_hydro")
loo_hydro_regression  <- read_output("loo_hydro_regression")
forecast_year         <- read_output("forecast_year")
y_pred_2027_base       <- read_output("y_pred_2027_base")
y_pred_2027_reg        <- read_output("y_pred_2027_reg")

direction_colors <- c(
  "Earlier than predicted" = "steelblue",
  "Later than predicted" = "coral",
  "Within 14 days" = "grey40"
)

model_colors <- c(
  "Baseline" = "grey40",
  "Water regression" = "steelblue"
)

figure_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Helpers ----------------------------------------------------------------
make_ppc <- function(fit, y, title) {
  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep
  
  if (ncol(y_rep) != length(y)) {
    stop(
      "Posterior predictive observations do not match model_data. ",
      "Refit full models on the same observations."
    )
  }
  
  rows <- seq_len(min(50L, nrow(y_rep)))
  
  bayesplot::ppc_dens_overlay(
    y = y,
    yrep = y_rep[rows, , drop = FALSE]
  ) +
    labs(title = title) +
    figure_theme()
}

make_forecast_plot <- function(predictions, test, model_label) {
  if (nrow(test) != 1L) {
    stop("Forecast plotting expects a single test year.")
  }
  
  anchor <- as.Date(paste0(test$year - 1, "-05-01"))
  actual_date <- anchor + test$y - 1
  median_date <- anchor + round(median(predictions)) - 1
  days_off <- round(median(predictions) - test$y)
  
  timing <- if (days_off > 0) {
    "days late"
  } else if (days_off < 0) {
    "days early"
  } else {
    "days off"
  }
  
  tibble(date = anchor + predictions - 1) |>
    ggplot(aes(x = date)) +
    geom_histogram(
      aes(y = after_stat(density), fill = "Predicted"),
      bins = 30,
      color = "white"
    ) +
    geom_vline(
      xintercept = actual_date,
      color = "red",
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = median_date,
      color = "darkblue",
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    annotate(
      "text",
      x = actual_date + 3,
      y = Inf,
      label = format(actual_date, "%b %d"),
      hjust = 0,
      vjust = 3.2,
      color = "red",
      size = 3.5
    ) +
    annotate(
      "text",
      x = actual_date + 3,
      y = Inf,
      label = paste(abs(days_off), timing),
      hjust = 0,
      vjust = 1.5,
      color = "red",
      size = 3.5
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("Predicted" = "steelblue")
    ) +
    scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
    labs(
      title = paste0(
        "Forecast vs Actual\n", model_label, " — ", test$year
      ),
      x = "Nest Initiation Date",
      y = "Density"
    ) +
    figure_theme()
}

make_future_plot <- function(predictions, year, model_label) {
  anchor <- as.Date(paste0(year - 1, "-05-01"))
  q <- unname(quantile(
    predictions,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
  ))
  dates <- anchor + round(q) - 1
  
  lines <- tibble(
    date = dates,
    interval = factor(
      c(
        "95% predictive interval",
        "50% predictive interval",
        "Median",
        "50% predictive interval",
        "95% predictive interval"
      ),
      levels = c(
        "Median",
        "50% predictive interval",
        "95% predictive interval"
      )
    )
  )
  
  interval_labels <- c(
    "Median" = paste("Median:", format(dates[3], "%b %d")),
    "50% predictive interval" = paste(
      "50%:", format(dates[2], "%b %d"),
      "to", format(dates[4], "%b %d")
    ),
    "95% predictive interval" = paste(
      "95%:", format(dates[1], "%b %d"),
      "to", format(dates[5], "%b %d")
    )
  )
  
  tibble(date = anchor + predictions - 1) |>
    ggplot(aes(x = date)) +
    geom_histogram(
      aes(y = after_stat(density)),
      fill = "steelblue",
      bins = 30,
      color = "white"
    ) +
    geom_vline(
      data = lines,
      aes(xintercept = date, linetype = interval),
      color = "darkred",
      linewidth = 0.9
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Median" = "solid",
        "50% predictive interval" = "dashed",
        "95% predictive interval" = "dotted"
      ),
      labels = interval_labels
    ) +
    scale_x_date(date_labels = "%b %d", date_breaks = "2 months") +
    labs(
      title = paste0(
        year, " WOST Nest Initiation Forecast\n", model_label
      ),
      x = "Predicted Nest Initiation Date",
      y = "Density"
    ) +
    figure_theme() +
    theme(legend.position = "bottom")
}

make_residual_plot <- function(results, model_label) {
  results <- results |> arrange(year)
  
  ggplot(results, aes(x = year, y = residual)) +
    geom_hline(
      yintercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = c(-14, 14),
      color = "coral",
      linetype = "dashed"
    ) +
    geom_line(color = "grey40", linewidth = 0.5) +
    geom_point(aes(color = direction), size = 3) +
    ggrepel::geom_label_repel(
      data = results |> filter(off_year),
      aes(label = year, color = direction),
      size = 3,
      show.legend = FALSE,
      seed = 42,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      name = NULL,
      values = direction_colors,
      limits = names(direction_colors),
      drop = FALSE
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
      title = paste0(
        "Expanding Window — ", model_label, " Residuals"
      ),
      subtitle = paste(
        "Residual = Actual − Predicted median",
        "| Dashed lines = ±14-day threshold"
      ),
      x = "Year",
      y = "Residual (days) — positive = later than predicted"
    ) +
    figure_theme()
}

make_hydro_plot <- function(results, model_label) {
  hydro_vars <- water_system |>
    select(where(is.numeric), -year) |>
    names()
  
  long_data <- results |>
    select(
      year, residual, off_year, direction,
      all_of(hydro_vars)
    ) |>
    pivot_longer(
      cols = all_of(hydro_vars),
      names_to = "variable",
      values_to = "value"
    ) |>
    filter(is.finite(value), is.finite(residual))
  
  ggplot(long_data, aes(x = value, y = residual)) +
    geom_hline(
      yintercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = c(-14, 14),
      color = "coral",
      linetype = "dashed"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = "grey30",
      se = TRUE
    ) +
    geom_point(aes(color = direction), size = 2, alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = long_data |> filter(off_year),
      aes(label = year, color = direction),
      size = 2.5,
      show.legend = FALSE,
      seed = 42,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      name = NULL,
      values = direction_colors,
      limits = names(direction_colors),
      drop = FALSE
    ) +
    facet_wrap(~ variable, scales = "free_x") +
    labs(
      title = paste0(
        "Hydro Covariates vs. Expanding-Window Residuals\n",
        model_label
      ),
      subtitle = "Dashed lines = ±14-day threshold",
      x = "Covariate value",
      y = "Residual (days) — positive = later than predicted"
    ) +
    figure_theme()
}

# 1–2. Posterior predictive checks ---------------------------------------
p_ppc_base <- make_ppc(
  fit_baseline_full,
  model_data$y,
  "Posterior Predictive Check — Baseline Model"
)

p_ppc_reg <- make_ppc(
  fit_reg_full,
  model_data$y,
  "Posterior Predictive Check — Water Regression Model"
)

save_figure(p_ppc_base, "01_ppc_baseline.png")
save_figure(p_ppc_reg, "02_ppc_regression.png")

# 3–4. Forecast vs actual -------------------------------------------------
p_forecast_base <- make_forecast_plot(
  y_pred_base, test_base, "Baseline Model"
)

p_forecast_reg <- make_forecast_plot(
  y_pred_reg, test_reg, "Water Regression Model"
)

save_figure(p_forecast_base, "03_forecast_baseline.png")
save_figure(p_forecast_reg, "04_forecast_regression.png")

# 5–6. Future forecasts ---------------------------------------------------
p_2027_base <- make_future_plot(
  y_pred_2027_base, forecast_year, "Baseline Model"
)

save_figure(p_2027_base, "05_forecast_2027_baseline.png")

future_reg_filename <- file.path(
  "figures", "06_forecast_2027_regression.png"
)

if (!is.null(y_pred_2027_reg)) {
  p_2027_reg <- make_future_plot(
    y_pred_2027_reg,
    forecast_year,
    "Water Regression — Conditional on Supplied Covariates"
  )
  
  save_figure(p_2027_reg, "06_forecast_2027_regression.png")
} else {
  # Remove stale figures from an earlier run.
  if (file.exists(future_reg_filename)) {
    unlink(future_reg_filename)
  }
  
  message("2027 water forecast figure skipped: no future covariates supplied.")
}

# 7. Baseline residuals ---------------------------------------------------
p_loo_resid <- make_residual_plot(
  loo_results_baseline, "Baseline Model"
)

save_figure(
  p_loo_resid,
  "07_loo_residuals_baseline.png",
  width = 9,
  height = 5
)

# 8. Baseline hydro covariates vs residuals -------------------------------
p_hydro <- make_hydro_plot(loo_hydro, "Baseline Model")

save_figure(
  p_hydro,
  "08_hydro_loo_residuals.png",
  width = 12,
  height = 8
)

# 9. NEW: Water regression residuals --------------------------------------
p_loo_resid_reg <- make_residual_plot(
  loo_results_regression, "Water Regression Model"
)

save_figure(
  p_loo_resid_reg,
  "09_loo_residuals_regression.png",
  width = 9,
  height = 5
)

# 10. NEW: Side-by-side residual comparison -------------------------------
# Use only years evaluated by both models.
comparison_data <- expanding_results |>
  mutate(
    model = factor(
      model,
      levels = c("Baseline", "Water regression")
    )
  ) |>
  arrange(model, year)

p_residual_comparison <- ggplot(
  comparison_data,
  aes(x = year, y = residual)
) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-14, 14),
    color = "coral",
    linetype = "dashed"
  ) +
  geom_line(color = "grey40", linewidth = 0.5) +
  geom_point(aes(color = direction), size = 3) +
  ggrepel::geom_label_repel(
    data = comparison_data |> filter(off_year),
    aes(label = year, color = direction),
    size = 3,
    show.legend = FALSE,
    seed = 42,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    name = NULL,
    values = direction_colors,
    limits = names(direction_colors),
    drop = FALSE
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    title = "Expanding-Window Residuals — Model Comparison",
    subtitle = "Shared evaluation years | Dashed lines = ±14-day threshold",
    x = "Year",
    y = "Residual (days) — positive = later than predicted"
  ) +
  figure_theme()

save_figure(
  p_residual_comparison,
  "10_loo_residuals_comparison.png",
  width = 10,
  height = 8
)

# 11. NEW: Water hydro covariates vs residuals -----------------------------
p_hydro_reg <- make_hydro_plot(
  loo_hydro_regression, "Water Regression Model"
)

save_figure(
  p_hydro_reg,
  "11_hydro_loo_residuals_regression.png",
  width = 12,
  height = 8
)

# 12. NEW: Absolute forecast errors ---------------------------------------
p_absolute_errors <- ggplot(
  comparison_data,
  aes(x = year, y = abs(residual), color = model, group = model)
) +
  geom_hline(
    yintercept = 14,
    color = "coral",
    linetype = "dashed"
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(name = NULL, values = model_colors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    title = "Expanding-Window Absolute Forecast Errors",
    subtitle = "Shared evaluation years | Lower is better",
    x = "Year",
    y = "Absolute error (days)"
  ) +
  figure_theme()

save_figure(
  p_absolute_errors,
  "12_absolute_error_comparison.png",
  width = 9,
  height = 5
)

cat("\nplot_models.R complete — figures saved to figures/.\n")