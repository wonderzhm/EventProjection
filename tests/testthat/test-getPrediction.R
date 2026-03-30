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

  testthat::expect_true("event_pred_time" %in% names(pred1$event_pred))
  testthat::expect_length(as.numeric(pred1$event_pred$event_pred_time), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred$event_pred_time))))
  testthat::expect_equal(
    as.numeric(pred1$event_pred$event_pred_day),
    EventProjection:::.ep_study_time_to_day(pred1$event_pred$event_pred_time)
  )
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

  testthat::expect_true("event_pred_time" %in% names(pred1$event_pred))
  testthat::expect_length(as.numeric(pred1$event_pred$event_pred_time), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred$event_pred_time))))
  testthat::expect_equal(
    as.numeric(pred1$event_pred$event_pred_day),
    EventProjection:::.ep_study_time_to_day(pred1$event_pred$event_pred_time)
  )
  testthat::expect_equal(
    as.numeric(pred1$event_pred$event_pred_date - as.Date(interimData1$trialsdt[1]) + 1),
    as.numeric(pred1$event_pred$event_pred_day)
  )
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

  testthat::expect_true("event_pred_time" %in% names(pred1$event_pred))
  testthat::expect_length(as.numeric(pred1$event_pred$event_pred_time), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred$event_pred_time))))
  testthat::expect_equal(
    as.numeric(pred1$event_pred$event_pred_day),
    EventProjection:::.ep_study_time_to_day(pred1$event_pred$event_pred_time)
  )
  testthat::expect_equal(
    as.numeric(pred1$event_pred$event_pred_date - as.Date(interimData2$trialsdt[1]) + 1),
    as.numeric(pred1$event_pred$event_pred_day)
  )
})


testthat::test_that("getPrediction can omit duplicated large outputs", {
  set.seed(2000)

  pred1 <- getPrediction(
    df = interimData1,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 200,
    enroll_model = "time-decay",
    event_model = "weibull",
    dropout_model = "exponential",
    pilevel = 0.90, nreps = 20,
    showsummary = FALSE, showplot = FALSE,
    return_subject_data = FALSE,
    return_simulation_data = FALSE)

  testthat::expect_false("subject_data" %in% names(pred1))
  testthat::expect_false("newSubjects" %in% names(pred1$enroll_pred))
  testthat::expect_true("enroll_pred_df" %in% names(pred1$enroll_pred))
  testthat::expect_false("newEvents" %in% names(pred1$event_pred))
  testthat::expect_true("enroll_pred_df" %in% names(pred1$event_pred))
  testthat::expect_true("event_pred_df" %in% names(pred1$event_pred))
  testthat::expect_true("dropout_pred_df" %in% names(pred1$event_pred))
  testthat::expect_true("ongoing_pred_df" %in% names(pred1$event_pred))
  testthat::expect_true("event_pred_day" %in% names(pred1$event_pred))
  testthat::expect_true("event_pred_time" %in% names(pred1$event_pred))
})


testthat::test_that("getPrediction returns a prior-informed event fit", {
  set.seed(3000)

  event_prior <- list(
    model = "weibull",
    theta = c(log(500), 0),
    vtheta = diag(c(0.05, 0.05))
  )

  pred1 <- getPrediction(
    df = interimData2,
    to_predict = "event only",
    target_d = 200,
    event_model = "weibull",
    event_prior = event_prior,
    dropout_model = "none",
    pilevel = 0.90,
    nreps = 20,
    showsummary = FALSE,
    showplot = FALSE
  )

  testthat::expect_true("event_fit_posterior" %in% names(pred1))
  testthat::expect_true(all(c("fit", "kmdf", "dffit", "text") %in%
                              names(pred1$event_fit_posterior)))
  testthat::expect_equal(
    pred1$event_fit_posterior$fit$model,
    pred1$event_fit$fit$model
  )
  testthat::expect_gt(
    sum(abs(pred1$event_fit$dffit$surv - pred1$event_fit_posterior$dffit$surv)),
    0
  )
})


testthat::test_that("getPrediction returns a prior-informed dropout fit", {
  set.seed(3000)

  dropout_prior <- list(
    model = "exponential",
    theta = log(1 / 800),
    vtheta = 0.05
  )

  pred1 <- getPrediction(
    df = interimData2,
    to_predict = "event only",
    target_d = 200,
    event_model = "weibull",
    dropout_model = "exponential",
    dropout_prior = dropout_prior,
    pilevel = 0.90,
    nreps = 20,
    showsummary = FALSE,
    showplot = FALSE
  )

  testthat::expect_true("dropout_fit_posterior" %in% names(pred1))
  testthat::expect_true(all(c("fit", "kmdf", "dffit", "text") %in%
                              names(pred1$dropout_fit_posterior)))
  testthat::expect_equal(
    pred1$dropout_fit_posterior$fit$model,
    pred1$dropout_fit$fit$model
  )
  testthat::expect_gt(
    sum(abs(pred1$dropout_fit$dffit$surv - pred1$dropout_fit_posterior$dffit$surv)),
    0
  )
})


testthat::test_that("getPrediction rejects design-stage event-only requests", {
  testthat::expect_error(
    getPrediction(
      df = NULL,
      to_predict = "event only",
      target_d = 200,
      showsummary = FALSE,
      showplot = FALSE
    ),
    "does not support design-stage event-only prediction"
  )
})
