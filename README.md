# EventProjection
Predicts enrollment and event outcomes at the design or interim-analysis stage using parametric enrollment, dropout, and time-to-event models with simulation-based projection. Extending the functionality of the 'eventPred' package, this package incorporates cure rate and delayed treatment effect models adapted from 'EventPredInCure', and further adds support for user-specified fixed cure rates. It provides a unified framework for model fitting, prediction, and visualization across design scenarios and ongoing clinical trials.

## Main features

- Unified workflow centered on `fitEnrollment()`, `fitEvent()`, `fitDropout()`, `predictEnrollment()`, `predictEvent()`, and `getPrediction()`.
- Supports both design-stage projection and interim-analysis projection.
- Preserves standard non-cure projection workflows from `eventPred`.
- Adds cure-rate and delayed-treatment-effect event models adapted from `EventPredInCure`.
- Supports either estimated cure fractions or user-specified fixed cure rates.
- Includes a Shiny app for interactive model fitting and projection.
- Includes example datasets for quick start and reproducible demonstrations.

## Installation

You can install the package from GitHub:
```r
library(devtools)
devtools::install_github("wonderzhm/EventProjection")
```
If you want to run the Shiny app, install the optional app dependencies as well:

```r
install.packages(c(
  "shinyMatrix", "shinyFeedback", "shinyjs", "shinybusy",
  "readxl", "writexl", "DT", "prompter"
))
```
## Quick start

```r
library(EventProjection)

data(interimData2)

set.seed(3000)
pred <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 100,
  showplot = FALSE
)

pred
```

## Example workflow

A typical analysis-stage workflow is:

1. Summarize the observed study status with `summarizeObserved()`.
2. Fit one or more enrollment models with `fitEnrollment()` if enrollment is ongoing.
3. Fit one or more event models with `fitEvent()`.
4. Fit one or more dropout models with `fitDropout()`.
5. Project future enrollment and events with `predictEnrollment()`, `predictEvent()`, or `getPrediction()`.

A typical design-stage workflow is:

1. Specify prior distributions for enrollment, event, and optionally dropout.
2. Pass the prior objects to `predictEnrollment()`, `predictEvent()`, or `getPrediction()`.
3. Use simulated outputs and prediction intervals to summarize uncertainty in study timelines.

## Supported model classes

### Enrollment models

- `"poisson"`
- `"time-decay"`
- `"b-spline"`
- `"piecewise poisson"`
- `"piecewise uniform"`

### Standard event models

- `"exponential"`
- `"weibull"`
- `"log-logistic"`
- `"log-normal"`
- `"piecewise exponential"`
- `"model averaging"`
- `"spline"`
- `"cox"`

### Cure-aware event models

- `"exponential with cured population"`
- `"weibull with cured population"`
- `"log-normal with cured population"`
- `"log-logistic with cured population"`
- `"piecewise exponential with cured population"`

### Dropout models

- `"exponential"`
- `"weibull"`
- `"log-logistic"`
- `"log-normal"`
- `"piecewise exponential"`
- `"model averaging"`
- `"cox"`

## Shiny app

To launch the interactive app:

```r
EventProjection::runShinyApp()
```

For backward compatibility, the legacy launcher name is also available:

```r
EventProjection::runShinyApp_eventPred()
```

The app supports importing study data from `.csv` and `.rda`/`.RData` files.

## Example datasets

The package includes the following example datasets:

- `interimData1`
- `interimData2`
- `finalData`

These datasets can be loaded with `data()` and used to reproduce examples from the documentation and vignette.

## Documentation

A longer technical tutorial is included as a vignette:

```r
vignette("eventprojection_tutorial", package = "EventProjection")
```

## Notes

- Core modeling functions are intended to remain close to `eventPred` for standard non-cure, non-delay workflows.
- Cure-aware and delayed-treatment-effect functionality is intended to remain close to `EventPredInCure` while using the unified `EventProjection` interface.
- The package is currently prepared for GitHub distribution and internal/research use.

## License

GPL (>= 2)

