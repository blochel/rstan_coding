# Wood Stork Nest Initiation Models

## File Structure

```text
project/
├── master.R
├── fit_models.R
├── evaluate_models.R
├── plot_models.R
├── nest_initiation_models/
│   ├── baseline.stan
│   └── water_regression.stan
├── WaterData/
├── outputs/
└── figures/
```

## Baseline (Intercept-Only Normal) Model

**The model:**

$$
y_i \sim \mathrm{Normal}(\mu_i, \sigma)
$$

$$
\mu_i = \alpha
$$

- **$y_i$** = first observed nest initiation date in water year $i$, expressed as day of the water year
- **$\mu_i$** = expected nest initiation date
- **$\alpha$** = intercept
- **$\sigma$** = residual standard deviation, in days

Dates are expressed relative to May 1 of the preceding calendar year.
For example, water year 2027 begins on May 1, 2026.

## Water Regression Model

The water regression model allows expected nesting date to vary with
pre-recession and initial water depth.

**The model:**

$$
y_i \sim \mathrm{Normal}(\mu_i, \sigma)
$$

$$
\mu_i = \alpha + \beta_1 z_{1,i} + \beta_2 z_{2,i}
$$

- **$y_i$** = first observed nest initiation date in water year $i$, expressed as day of the water year
- **$\mu_i$** = expected nest initiation date in water year $i$
- **$\alpha$** = expected nesting date when both covariates are at their training-data means
- **$z_{1,i}$** = standardized `pre_recession`
- **$z_{2,i}$** = standardized `init_depth`
- **$\beta_1$** = change in expected nesting date, in days, for a one-standard-deviation increase in pre-recession, holding initial depth constant
- **$\beta_2$** = change in expected nesting date, in days, for a one-standard-deviation increase in initial depth, holding pre-recession constant
- **$\sigma$** = residual standard deviation, in days

Positive coefficients indicate later nesting; negative coefficients
indicate earlier nesting.

### Standardized Covariates

Each covariate is standardized using its training-data mean and standard
deviation:

$$
z_{k,i} = \frac{x_{k,i} - \bar{x}_{k}}{s_k}
$$

- **$x_{k,i}$** = original value of covariate $k$ in year $i$
- **$\bar{x}_{k}$** = mean of covariate $k$ in the training data
- **$s_k$** = standard deviation of covariate $k$ in the training data

The same training-data means and standard deviations are applied to
prediction data. In expanding-window evaluation, scaling is recalculated
using only the training years in each fold.

For genuine forecasting, covariates must be available at the forecast
date. Covariates calculated using observations after nesting begins
should not be treated as advance predictors.

## Expanding-Window Forecast Residuals

For each evaluation year $i$, the model is fitted using only earlier
years. The fitted model generates posterior predictive nesting dates for
the evaluation year.

The point prediction is the median of these posterior predictive draws:

$$
\widehat{y}_i = \mathrm{median}(y_i^{\mathrm{rep}})
$$

The residual is the observed nesting date minus the predicted median:

$$
r_i = y_i - \widehat{y}_i
$$

- **$r_i$** = forecast residual, in days
- **$y_i$** = observed nest initiation date
- **$y_i^{\mathrm{rep}}$** = posterior predictive draws for year $i$, generated from a model fitted to earlier years
- **$\widehat{y}_i$** = median posterior predictive nesting date

Positive residuals indicate nesting occurred later than predicted.
Negative residuals indicate nesting occurred earlier than predicted.
A residual of zero indicates the predicted median matched the observed date.

### Residual Classification

| Residual | Classification |
|---|---|
| $r_i < -14$ | Earlier than predicted |
| $-14 \leq r_i \leq 14$ | Within 14 days |
| $r_i > 14$ | Later than predicted |

Horizontal dashed lines at $-14$ and $+14$ days mark the forecast-error
thresholds. Years outside these thresholds are labeled in the residual plots.

This calculation applies to both the baseline and water regression models.
Model-comparison figures and accuracy summaries use shared evaluation years.
Water-model folds require complete covariates and at least 10 complete
training years.

## Model Evaluation

### PSIS-LOO

PSIS-LOO evaluates predictive performance by approximating the exclusion
of each observation from the full model fit.

- **ELPD:** higher values indicate better estimated predictive performance.
- **LOOIC:** lower values indicate better estimated predictive performance.
- **Pareto k:** diagnostic values used to assess the reliability of the approximation.

Comparisons require both models to use the same observations in the same order.
The ELPD difference should be interpreted alongside its standard error.

### Expanding-Window Evaluation

Expanding-window evaluation trains on earlier years and predicts a later
year. Unlike PSIS-LOO, it preserves the chronological training–prediction
structure.

Accuracy summaries include mean absolute error, root mean squared error,
mean residual, and the proportion of predictions within 14 days.

## Figures

Figures are saved in `figures/`.

| File | Description |
|---|---|
| `01_ppc_baseline.png` | Baseline posterior predictive check |
| `02_ppc_regression.png` | Water regression posterior predictive check |
| `03_forecast_baseline.png` | Baseline held-out-year forecast |
| `04_forecast_regression.png` | Water regression held-out-year forecast |
| `05_forecast_2027_baseline.png` | Baseline 2027 forecast |
| `06_forecast_2027_regression.png` | Conditional water regression 2027 forecast, if future covariates are supplied |
| `07_loo_residuals_baseline.png` | Baseline expanding-window residuals |
| `08_hydro_loo_residuals.png` | Hydro covariates versus baseline residuals |
| `09_loo_residuals_regression.png` | Water regression expanding-window residuals |
| `10_loo_residuals_comparison.png` | Residual comparison on shared evaluation years |
| `11_hydro_loo_residuals_regression.png` | Hydro covariates versus water regression residuals |
| `12_absolute_error_comparison.png` | Absolute forecast errors on shared evaluation years |

## Future Modeling Ideas

- Temperature as an additional predictor.
- Multiscale habitat selection and water dynamics near nesting locations.
- Additional water data for Corkscrew Swamp and Big Cypress.

Potential monitoring locations:

- https://waterdata.usgs.gov/monitoring-location/USGS-262159081290400
- https://waterdata.usgs.gov/monitoring-location/USGS-260355080541401