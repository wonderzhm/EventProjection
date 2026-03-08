testthat::test_that("cured-population event workflow runs on observed data", {
  fit1 <- fitEvent(
    df = interimData2,
    event_model = "weibull with cured population",
    showplot = FALSE,
    by_treatment = FALSE
  )

  testthat::expect_equal(fit1$fit$model, "weibull with cured population")
  testthat::expect_length(fit1$fit$theta, 3)
  testthat::expect_true(all(is.finite(fit1$fit$theta)))
  testthat::expect_true(all(c("time", "surv") %in% names(fit1$dffit)))
  testthat::expect_gt(nrow(fit1$dffit), 0)

  dropout_fits <- fitDropout(
    df = interimData2,
    dropout_model = "exponential",
    showplot = FALSE
  )

  set.seed(4101)
  pred1 <- predictEvent(
    df = interimData2,
    target_d = 200,
    event_fit = fit1$fit,
    dropout_fit = dropout_fits$fit,
    pilevel = 0.90,
    nreps = 20,
    showsummary = FALSE,
    showplot = FALSE
  )

  testthat::expect_length(as.numeric(pred1$event_pred_day), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred_day))))
  testthat::expect_gt(nrow(pred1$newEvents), 0)
})


testthat::test_that("fixed cure_rate is accepted in fitEvent and design-stage prediction", {
  fixed_cure <- 0.30

  fit1 <- fitEvent(
    df = interimData2,
    event_model = "weibull with cured population",
    cure_rate = fixed_cure,
    showplot = FALSE,
    by_treatment = FALSE
  )

  testthat::expect_equal(fit1$fit$model, "weibull with cured population")
  testthat::expect_equal(unname(fit1$fit$theta[1]), stats::qlogis(fixed_cure), tolerance = 1e-10)
  testthat::expect_equal(as.numeric(fit1$fit$vtheta[1, ]), c(0, 0, 0), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(fit1$fit$vtheta[, 1]), c(0, 0, 0), tolerance = 1e-12)

  enroll_prior <- list(
    model = "piecewise poisson",
    theta = log(26 / 9 * seq(1, 9) / 30.4375),
    vtheta = diag(9) * 1e-8,
    accrualTime = seq(0, 8) * 30.4375
  )

  event_prior <- list(
    model = "weibull with cured population",
    theta = c(log(1.2), log(300)),
    vtheta = diag(2) * 1e-8
  )

  dropout_prior <- list(
    model = "exponential",
    theta = log((-log(1 - 0.05) / 12) / 30.4375),
    vtheta = 1e-8
  )

  set.seed(4102)
  pred1 <- getPrediction(
    df = NULL,
    to_predict = "enrollment and event",
    target_n = 300,
    target_d = 150,
    enroll_prior = enroll_prior,
    event_prior = event_prior,
    dropout_prior = dropout_prior,
    cure_rate = fixed_cure,
    pilevel = 0.90,
    nreps = 20,
    showsummary = FALSE,
    showplot = FALSE
  )

  testthat::expect_equal(pred1$stage, "Design stage")
  testthat::expect_true("event_pred" %in% names(pred1))
  testthat::expect_length(as.numeric(pred1$event_pred$event_pred_day), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred$event_pred_day))))
  testthat::expect_gt(nrow(pred1$subject_data), 0)
})


testthat::test_that("piecewise uniform enrollment runs at the design stage", {
  enroll_prior <- list(
    model = "piecewise uniform",
    theta = log(c(25, 45) / 30.4375),
    vtheta = diag(2) * 1e-8,
    accrualTime = c(0, 120)
  )

  pred1 <- predictEnrollment(
    target_n = 80,
    enroll_fit = enroll_prior,
    pilevel = 0.90,
    nreps = 5,
    fix_parameter = TRUE,
    showsummary = FALSE,
    showplot = FALSE
  )

  testthat::expect_length(as.numeric(pred1$enroll_pred_day), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$enroll_pred_day))))
  testthat::expect_gt(nrow(pred1$newSubjects), 0)

  draw1 <- pred1$newSubjects[pred1$newSubjects$draw == 1, ][["arrivalTime"]]
  draw2 <- pred1$newSubjects[pred1$newSubjects$draw == 2, ][["arrivalTime"]]

  testthat::expect_true(all(diff(draw1) >= 0))
  testthat::expect_equal(draw1, draw2)
})


testthat::test_that("delayed-treatment design-stage branch runs by treatment", {
  enroll_prior <- list(
    model = "piecewise poisson",
    theta = log(26 / 9 * seq(1, 9) / 30.4375),
    vtheta = diag(9) * 1e-8,
    accrualTime = seq(0, 8) * 30.4375
  )

  event_prior_delay <- list(
    list(
      model = "weibull with cured population and delayed treatment",
      theta = c(stats::qlogis(0.20), log(1.1), log(320), 0, 1),
      vtheta = matrix(0, 5, 5)
    ),
    list(
      model = "weibull with cured population and delayed treatment",
      theta = c(stats::qlogis(0.35), log(1.1), log(320), 60, 0.70),
      vtheta = matrix(0, 5, 5)
    )
  )

  dropout_prior_one <- list(
    model = "exponential",
    theta = log((-log(1 - 0.05) / 12) / 30.4375),
    vtheta = 1e-8
  )

  dropout_prior <- list(
    dropout_prior_one,
    dropout_prior_one
  )

  set.seed(4103)
  pred1 <- getPrediction(
    df = NULL,
    to_predict = "enrollment and event",
    target_n = 120,
    target_d = 50,
    enroll_prior = enroll_prior,
    event_prior = event_prior_delay,
    dropout_prior = dropout_prior,
    pilevel = 0.90,
    nreps = 20,
    showsummary = FALSE,
    showplot = FALSE,
    by_treatment = TRUE,
    ngroups = 2,
    alloc = c(1, 1),
    treatment_label = c("Control", "Treatment")
  )

  testthat::expect_equal(pred1$stage, "Design stage")
  testthat::expect_true(all(c("treatment", "treatment_description") %in% names(pred1$subject_data)))
  testthat::expect_true(all(c("Control", "Treatment") %in% pred1$subject_data$treatment_description))
  testthat::expect_length(as.numeric(pred1$event_pred$event_pred_day), 3)
  testthat::expect_true(all(is.finite(as.numeric(pred1$event_pred$event_pred_day))))
})


testthat::test_that("app launcher builds a shiny app object when app dependencies are installed", {
  app_pkgs <- c(
    "shinyMatrix",
    "shinyFeedback",
    "shinyjs",
    "shinybusy",
    "readxl",
    "writexl",
    "DT",
    "prompter"
  )

  for (pkg in app_pkgs) {
    testthat::skip_if_not_installed(pkg)
  }

  app_runner <- if (exists("runShinyApp", mode = "function")) {
    runShinyApp
  } else {
    runShinyApp_eventPred
  }

  app <- app_runner()

  testthat::expect_s3_class(app, "shiny.appobj")
})
