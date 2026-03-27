testthat::test_that("enrollment prediction at design stage", {
  set.seed(1000)

  pred1 <- predictEnrollment(
    target_n = 300,
    enroll_fit = list(
      model = "piecewise poisson",
      theta = log(26/9*seq(1, 9)/30.4375),
      vtheta = diag(9)*1e-8,
      accrualTime = seq(0, 8)*30.4375),
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$enroll_pred_day)
  b <- c(469, 443, 509)

  testthat::expect_equal(a, b)
})


testthat::test_that("enrollment prediction during enrollment stage", {
  set.seed(2000)

  enroll_fits <- fitEnrollment(
    df = interimData1, enroll_model = "b-spline",
    nknots = 1, showplot = FALSE)

  pred1 <- predictEnrollment(
    df = interimData1,
    target_n = 300,
    enroll_fit = enroll_fits$fit,
    lags = 30,
    pilevel = 0.90, nreps = 100,
    showsummary = FALSE, showplot = FALSE)

  a <- as.numeric(pred1$enroll_pred_day)
  b <- c(439, 420, 471)

  testthat::expect_equal(a, b)
})


testthat::test_that("enrollment prediction can omit large subject-level output", {
  set.seed(2000)

  enroll_fits <- fitEnrollment(
    df = interimData1, enroll_model = "b-spline",
    nknots = 1, showplot = FALSE)

  pred1 <- predictEnrollment(
    df = interimData1,
    target_n = 300,
    enroll_fit = enroll_fits$fit,
    lags = 30,
    pilevel = 0.90, nreps = 20,
    showsummary = FALSE, showplot = FALSE,
    return_new_subjects = FALSE)

  testthat::expect_false("newSubjects" %in% names(pred1))
  testthat::expect_true("enroll_pred_df" %in% names(pred1))
  testthat::expect_true("enroll_pred_day" %in% names(pred1))
})


testthat::test_that("b-spline enrollment prediction handles lags beyond observed history", {
  set.seed(2000)

  enroll_fits <- fitEnrollment(
    df = interimData1, enroll_model = "b-spline",
    nknots = 1, showplot = FALSE)

  testthat::expect_no_error(
    predictEnrollment(
      df = interimData1,
      target_n = 300,
      enroll_fit = enroll_fits$fit,
      lags = 1000,
      pilevel = 0.90, nreps = 20,
      showsummary = FALSE, showplot = FALSE
    )
  )
})
