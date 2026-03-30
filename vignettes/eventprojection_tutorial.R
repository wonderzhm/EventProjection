## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## -----------------------------------------------------------------------------
library(EventProjection)


## -----------------------------------------------------------------------------
data(interimData1)
data(interimData2)
data(finalData)


## -----------------------------------------------------------------------------
obs1 <- summarizeObserved(
  df = interimData1,
  to_predict = "enrollment and event",
  showplot = FALSE,
  by_treatment = FALSE
)

obs2 <- summarizeObserved(
  df = interimData2,
  to_predict = "event only",
  showplot = FALSE,
  by_treatment = FALSE
)


## -----------------------------------------------------------------------------
fit_enroll_bs <- fitEnrollment(
  df = interimData1,
  enroll_model = "b-spline",
  nknots = 1,
  showplot = FALSE
)

fit_enroll_td <- fitEnrollment(
  df = interimData1,
  enroll_model = "time-decay",
  showplot = FALSE
)

fit_enroll_pw <- fitEnrollment(
  df = interimData1,
  enroll_model = "piecewise poisson",
  accrualTime = c(0, 90, 180),
  showplot = FALSE
)


## -----------------------------------------------------------------------------
pred_enroll <- predictEnrollment(
  df = interimData1,
  target_n = 300,
  enroll_fit = fit_enroll_bs$fit,
  lags = 30,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
design_enroll_fit <- list(
  model = "piecewise poisson",
  theta = log(26 / 9 * seq(1, 9) / 30.4375),
  vtheta = diag(9) * 1e-8,
  accrualTime = seq(0, 8) * 30.4375
)

design_enroll_pred <- predictEnrollment(
  target_n = 300,
  enroll_fit = design_enroll_fit,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_weibull <- fitEvent(
  df = interimData2,
  event_model = "weibull",
  showplot = FALSE
)

fit_event_piecewise <- fitEvent(
  df = interimData2,
  event_model = "piecewise exponential",
  piecewiseSurvivalTime = c(0, 140, 352),
  showplot = FALSE
)

fit_event_ma <- fitEvent(
  df = interimData2,
  event_model = "model averaging",
  showplot = FALSE
)

fit_event_cox <- fitEvent(
  df = interimData2,
  event_model = "cox",
  m = 5,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_dropout_exp <- fitDropout(
  df = interimData2,
  dropout_model = "exponential",
  showplot = FALSE
)

fit_dropout_piecewise <- fitDropout(
  df = interimData2,
  dropout_model = "piecewise exponential",
  piecewiseDropoutTime = c(0, 180),
  showplot = FALSE
)


## -----------------------------------------------------------------------------
pred_event_complete <- predictEvent(
  df = interimData2,
  target_d = 200,
  event_fit = fit_event_piecewise$fit,
  dropout_fit = fit_dropout_exp$fit,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_interim <- fitEvent(
  df = interimData1,
  event_model = "weibull",
  showplot = FALSE
)

fit_dropout_interim <- fitDropout(
  df = interimData1,
  dropout_model = "exponential",
  showplot = FALSE
)

pred_event_interim <- predictEvent(
  df = interimData1,
  target_d = 200,
  newSubjects = pred_enroll$newSubjects,
  event_fit = fit_event_interim$fit,
  dropout_fit = fit_dropout_interim$fit,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
pred_full_interim <- getPrediction(
  df = interimData1,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_model = "time-decay",
  event_model = "weibull",
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
pred_full_complete <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "cox",
  m = 5,
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
event_prior <- list(
  model = "piecewise exponential",
  theta = log(c(0.0533, 0.0309) / 30.4375),
  vtheta = diag(2) * 1e-8,
  piecewiseSurvivalTime = c(0, 6) * 30.4375
)

dropout_prior <- list(
  model = "exponential",
  theta = log((-log(1 - 0.05) / 12) / 30.4375),
  vtheta = 1e-8
)

pred_design <- getPrediction(
  df = NULL,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_prior = design_enroll_fit,
  event_prior = event_prior,
  dropout_prior = dropout_prior,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
set.seed(101)
interim_cov <- as.data.frame(interimData2)
interim_cov$age <- round(rnorm(nrow(interim_cov), mean = 62, sd = 8))
interim_cov$sex <- factor(sample(c("Female", "Male"), nrow(interim_cov), replace = TRUE))

pred_cov <- getPrediction(
  df = interim_cov,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull",
  dropout_model = "exponential",
  covariates_event = c("age", "sex"),
  covariates_dropout = c("age"),
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_cure_est <- fitEvent(
  df = interimData2,
  event_model = "weibull with cured population",
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_cure_fixed <- fitEvent(
  df = interimData2,
  event_model = "weibull with cured population",
  cure_rate = 0.25,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_cure_fixed_by_trt <- fitEvent(
  df = interimData2,
  event_model = "weibull with cured population",
  by_treatment = TRUE,
  cure_rate = c(0.20, 0.35),
  showplot = FALSE
)


## -----------------------------------------------------------------------------
fit_event_cure_piecewise_fixed <- fitEvent(
  df = interimData2,
  event_model = "piecewise exponential with cured population",
  cure_rate = 0.30,
  piecewiseSurvivalTime = c(0, 120, 240),
  showplot = FALSE
)


## -----------------------------------------------------------------------------
pred_event_cure_fixed <- getPrediction(
  df = interimData2,
  to_predict = "event only",
  target_d = 200,
  event_model = "weibull with cured population",
  cure_rate = 0.25,
  dropout_model = "exponential",
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
event_prior_cure_fixed <- list(
  model = "weibull with cured population",
  theta = c(log(1.2), log(300)),
  vtheta = diag(c(1e-8, 1e-8))
)

pred_design_cure_fixed <- getPrediction(
  df = NULL,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_prior = design_enroll_fit,
  event_prior = event_prior_cure_fixed,
  dropout_prior = dropout_prior,
  cure_rate = 0.30,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
event_prior_cure_fixed_pw <- list(
  model = "piecewise exponential with cured population",
  theta = log(c(0.0040, 0.0025, 0.0015)),
  vtheta = diag(3) * 1e-8,
  piecewiseSurvivalTime = c(0, 120, 240)
)

pred_design_cure_fixed_pw <- getPrediction(
  df = NULL,
  to_predict = "enrollment and event",
  target_n = 300,
  target_d = 200,
  enroll_prior = design_enroll_fit,
  event_prior = event_prior_cure_fixed_pw,
  dropout_prior = dropout_prior,
  cure_rate = 0.30,
  pilevel = 0.90,
  nreps = 200,
  showsummary = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
event_prior_delay <- list()
event_prior_delay[[1]] <- list(
  model = "weibull with cured population and delayed treatment",
  theta = c(qlogis(0.20), log(1.1), log(320), 0, 1),
  vtheta = matrix(0, 5, 5)
)
event_prior_delay[[2]] <- list(
  model = "weibull with cured population and delayed treatment",
  theta = c(qlogis(0.35), log(1.1), log(320), 60, 0.70),
  vtheta = matrix(0, 5, 5)
)


## -----------------------------------------------------------------------------
fit_event_by_trt <- fitEvent(
  df = interimData2,
  event_model = "weibull",
  by_treatment = TRUE,
  showplot = FALSE
)

fit_dropout_by_trt <- fitDropout(
  df = interimData2,
  dropout_model = "exponential",
  by_treatment = TRUE,
  showplot = FALSE
)


## -----------------------------------------------------------------------------
sessionInfo()

