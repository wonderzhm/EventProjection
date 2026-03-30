testthat::test_that("continuous study-time dates are mapped consistently", {
  origin_date <- as.Date("2020-01-01")

  testthat::expect_equal(
    EventProjection:::.ep_study_time_to_date(c(1, 1.9, 2, 2.75), origin_date),
    as.Date(c("2020-01-01", "2020-01-01", "2020-01-02", "2020-01-02"))
  )

  testthat::expect_equal(
    EventProjection:::.ep_followup_time_to_date(
      c(1, 1.25, 2, 3.6),
      as.Date("2020-06-15")
    ),
    as.Date(c("2020-06-15", "2020-06-15", "2020-06-16", "2020-06-17"))
  )
})


testthat::test_that("study-time day summaries preserve continuous-time internals", {
  testthat::expect_equal(
    EventProjection:::.ep_study_time_to_day(c(1, 1.2, 2.8, 10.01)),
    c(1, 1, 2, 10)
  )
})


testthat::test_that("same-date prediction rows keep the latest time point", {
  dt <- data.table::data.table(
    t = c(31, 31.2, 32),
    n = c(10, 11, 12),
    lower = c(9, 10, 11),
    upper = c(11, 12, 13),
    pilevel = 0.9,
    mean = c(10, 11, 12),
    var = c(0, 0, 0),
    parameter = "Enrollment",
    date = EventProjection:::.ep_study_time_to_date(
      c(31, 31.2, 32),
      as.Date("2020-01-01")
    )
  )

  out <- EventProjection:::.ep_collapse_same_date_rows(dt)

  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(out$t, c(31.2, 32))
  testthat::expect_equal(out$n, c(11, 12))
})
