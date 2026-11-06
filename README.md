
**File Structure **
```
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

$$y \sim \text{Normal}(\mu, \sigma)$$ $$\mu = \alpha$$

- **y** = nest initiation date (Julian day of year)
- **μ** = mean nest initiation date
- **α** = intercept
- **σ** = standard deviation


## Water Regression Model

**The model:**

$$
y_i \sim \text{Normal}(\mu_i, \sigma)
$$

$$
\mu_i = \alpha
+ \beta_{\text{pre}} z_{\text{pre\_recession},i}
+ \beta_{\text{depth}} z_{\text{init\_depth},i}
$$

**Standardized covariates:**

$$
z_{k,i} = \frac{x_{k,i} - \bar{x}_{k,\text{train}}}{s_{k,\text{train}}}
$$

- **$y_i$** = first observed nest initiation date in water year $i$, expressed as day of the water year
- **$\mu_i$** = expected nest initiation date in water year $i$
- **$\alpha$** = intercept: expected nesting date when both covariates are at their training-data means
- **$z_{\text{pre\_recession},i}$** = standardized pre-recession covariate
- **$z_{\text{init\_depth},i}$** = standardized initial water depth
- **$\beta_{\text{pre}}, \beta_{\text{depth}}$** = changes in expected nesting date, in days, for a one-standard-deviation increase in each covariate, holding the other constant
- **$\sigma$** = residual standard deviation, in days

Positive coefficients indicate later nesting; negative coefficients indicate earlier nesting. Covariates are standardized using training-data means and standard deviations, which are also applied to prediction data.






## Expanding-Window Forecast Residuals

For each evaluation year $i$, the model is fitted using only earlier years. The residual is the observed nesting date minus the median posterior predictive nesting date:

$$
r_i = y_i - \widehat{y}_i
$$

$$
\widehat{y}_i =
\operatorname{median}\left(
y_i^{\mathrm{rep}} \mid \text{data from years before } i
\right)
$$

- **$r_i$** = forecast residual, in days
- **$y_i$** = observed nest initiation date, expressed as day of the water year
- **$\widehat{y}_i$** = median posterior predictive nesting date for year $i$

Positive residuals indicate nesting occurred later than predicted;
negative residuals indicate nesting occurred earlier than predicted.
A residual of zero indicates the prediction matched the observed date.

**Residual classification:**

$$
\text{Category}_i =
\begin{cases}
\text{Earlier than predicted}, & r_i < -14 \\
\text{Within 14 days}, & -14 \leq r_i \leq 14 \\
\text{Later than predicted}, & r_i > 14
\end{cases}
$$

The horizontal lines at $-14$ and $+14$ days mark the forecast-error
thresholds. Years outside these thresholds are labeled in the plot.
This calculation applies to both the baseline and water regression models.




# ideas 
- Temperature
- multiscale habitat selection (water dynamics in area of initiation)
- missing water data for Corkscrew and Big Cypress...
      - https://waterdata.usgs.gov/monitoring-location/USGS-262159081290400#dataTypeId=continuous-63160-0&period=P365D&showMedian=true&showFieldMeasurements=true
      - https://waterdata.usgs.gov/monitoring-location/USGS-260355080541401/#dataTypeId=continuous-00065-0&period=P7D&showFieldMeasurements=true 
