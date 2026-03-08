testthat::test_that("event prediction at design stage", {
  set.seed(1000)

  enroll_prior <- list(
    model = "piecewise poisson",
    accrualTime = seq(0, 8)*30.4375,
    theta = log(26/9*seq(1, 9)/30.4375),
    vtheta = diag(9)*1e-8)

  event_prior <- list(
    model = "piecewise exponential",
    theta = log(c(0.0533, 0.0309)/30.4375),
    vtheta = diag(2)*1e-8,
    piecewiseSurvivalTime = c(0,6)*30.4375)

  dropout_prior <- list(
    model = "exponential",
    theta = log((-log(1-0.05)/12)/30.4375),
    vtheta = 1e-8)

  pred1 <- getPrediction(
    df = NULL,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 200,
    enroll_prior = enroll_prior,
    event_prior = event_prior,
    dropout_prior = dropout_prior,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1352, 1180, 1527)

  testthat::expect_equal(a, b)
})


testthat::test_that("event prediction during enrollment phase", {
  set.seed(2000)

  pred1 <- getPrediction(
    df = interimData1,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 200,
    enroll_model = "time-decay",
    event_model = "weibull",
    dropout_model = "exponential",
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1231, 835, 2018)

  testthat::expect_equal(a, b)
})


testthat::test_that("event prediction after enrollment completion", {
  set.seed(3000)

  pred1 <- getPrediction(
    df = interimData2, to_predict = "event only",
    target_d = 200,
    event_model = "cox", m = 5,
    dropout_model = "exponential",
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$event_pred$event_pred_day)
  b <- c(1027, 1003, 1099)

  testthat::expect_equal(a, b)
})
