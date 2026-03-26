#' @title Fit enrollment model
#' @description Fits a specified enrollment model to subject-level enrollment data.
#'
#' @details
#' This function is adapted from the original implementations in the
#' \pkg{eventPred} and \pkg{EventPredInCure} packages (original author: Kaifeng Lu).
#'
#' **Merged/updated features**
#' \itemize{
#'   \item Keeps \pkg{eventPred} plotting controls: \code{showplot}, \code{generate_plot},
#'         \code{interactive_plot}.
#'   \item Adds \code{criterion} from \pkg{EventPredInCure} to control whether the plot
#'         annotates AIC, BIC, or both.
#'   \item Makes column-name handling more robust by lowercasing input column names
#'         (requires \code{trialsdt}, \code{randdt}, \code{cutoffdt} after lowercasing).
#'   \item Removes aliases \code{enroll_fit} / \code{enroll_fit_plot} to keep the
#'         return object minimal; downstream callers should use \code{$fit} and
#'         \code{$fit_plot}.
#' }
#'
#' For the time-decay model, the mean function is
#' \deqn{\mu(t) = (\mu/\delta)\left(t - (1/\delta)(1 - \exp(-\delta t))\right)}
#' and the rate function is
#' \deqn{\lambda(t) = (\mu/\delta)(1 - \exp(-\delta t)).}
#'
#' For the B-spline model, the daily enrollment rate is
#' \deqn{\lambda(t) = \exp(B(t)^\top \theta),}
#' where \eqn{B(t)} represents the B-spline basis functions.
#'
#' @param df A \code{data.frame} containing subject-level enrollment information.
#'   Must include (case-insensitive) columns \code{trialsdt}, \code{randdt}, and \code{cutoffdt}.
#' @param enroll_model Enrollment model: one of \code{"Poisson"}, \code{"Time-decay"},
#'   \code{"B-spline"}, or \code{"Piecewise Poisson"} (case-insensitive).
#' @param nknots Number of inner knots for the B-spline enrollment model (non-negative integer).
#' @param accrualTime Accrual time cutpoints for piecewise Poisson (must start at 0 and be increasing).
#'   Example: \code{c(0, 30)} defines intervals [0,30) and [30,Inf).
#' @param showplot Logical; if \code{TRUE} and \code{generate_plot=TRUE}, prints the plot.
#' @param generate_plot Logical; if \code{TRUE}, create a fitted enrollment curve plot.
#' @param interactive_plot Logical; if \code{TRUE}, uses \pkg{plotly}; otherwise uses \pkg{ggplot2}.
#' @param criterion Character; which information criterion to display in the plot annotations:
#'   \code{"aic"}, \code{"bic"}, or \code{"both"}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{fit}: list containing \code{model}, \code{theta}, \code{vtheta}, \code{aic}, \code{bic},
#'         and optionally \code{x} (B-spline) or \code{accrualTime} (piecewise Poisson).
#'   \item \code{fit_plot}: (if \code{generate_plot=TRUE}) plot object (\pkg{plotly} or \pkg{ggplot2}).
#'   \item \code{enrolldf}: de-duplicated observed cumulative enrollment by day.
#'   \item \code{dffit}: fitted cumulative enrollment by day.
#'   \item \code{text}: character vector used in plot annotation (model + AIC/BIC depending on \code{criterion}).
#' }
#'
#' @references
#' Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for accrual in clinical trials.
#' Statistics in Medicine. 2010;29:649-658.
#'
#' @examples
#' \dontrun{
#' enroll_fit <- fitEnrollment(df = interimData1, 
#' enroll_model = "b-spline", nknots = 1, showplot = FALSE)
#' }
#'
#' @export
#'
fitEnrollment <- function(df,
                          enroll_model = "b-spline",
                          nknots = 0,
                          accrualTime = 0,
                          showplot = TRUE,
                          generate_plot = TRUE,
                          interactive_plot = TRUE,
                          criterion = "both") {
  
  # ---- Input checks ----
  erify::check_class(df, "data.frame")
  
  enroll_model_lc <- tolower(enroll_model)
  erify::check_content(enroll_model_lc,
                       c("poisson", "time-decay", "b-spline", "piecewise poisson"))
  
  criterion_lc <- tolower(criterion)
  erify::check_content(criterion_lc, c("aic", "bic", "both"))
  
  erify::check_n(nknots, zero = TRUE)
  erify::check_bool(showplot)
  erify::check_bool(generate_plot)
  erify::check_bool(interactive_plot)
  
  if (accrualTime[1] != 0) stop("accrualTime must start with 0")
  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }
  
  # UPDATED: make column handling case-insensitive
  df <- as.data.frame(df)
  names(df) <- tolower(names(df))
  
  required_cols <- c("trialsdt", "randdt", "cutoffdt")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  df$trialsdt <- as.Date(df$trialsdt)
  df$randdt   <- as.Date(df$randdt)
  df$cutoffdt <- as.Date(df$cutoffdt)
  
  trialsdt <- df$trialsdt[1]
  cutoffdt <- df$cutoffdt[1]
  n0 <- nrow(df)
  
  erify::check_positive(
    n0,
    supplement = "The number of subjects must be positive to fit an enrollment model."
  )
  
  if (any(df$randdt < trialsdt)) stop("randdt must be greater than or equal to trialsdt")
  if (any(df$randdt > cutoffdt)) stop("randdt must be less than or equal to cutoffdt")
  
  # up to the last randomization date to account for enrollment completion
  t0 <- as.numeric(max(df$randdt) - trialsdt + 1)
  
  # ---- Build observed cumulative enrollment (deduplicated by day) ----
  o <- order(df$randdt)
  df1 <- df[o, , drop = FALSE]
  df1$t <- as.numeric(df1$randdt - trialsdt + 1)
  df1$n <- seq_len(nrow(df1))
  
  # keep last subject for each randomization date
  idx_last <- tapply(seq_len(nrow(df1)), df1$randdt, max)
  df1u <- df1[as.integer(idx_last), c("t", "n"), drop = FALSE]
  df1u <- df1u[order(df1u$t), , drop = FALSE]
  
  # add day 1
  df1u <- rbind(data.frame(t = 1, n = 0), df1u)
  
  # ---- Fit enrollment model ----
  if (enroll_model_lc == "poisson") {
    fit1 <- list(
      model = "Poisson",
      theta = log(n0 / t0),
      vtheta = 1 / n0,
      aic = -2 * (-n0 + n0 * log(n0 / t0)) + 2,
      bic = -2 * (-n0 + n0 * log(n0 / t0)) + log(n0)
    )
    dffit1 <- data.frame(t = seq(1, t0))
    dffit1$n <- exp(fit1$theta) * dffit1$t
    
  } else if (enroll_model_lc == "time-decay") {
    
    fmu_td <- function(t, theta) {
      mu <- exp(theta[1])
      delta <- exp(theta[2])
      mu / delta * (t - 1 / delta * (1 - exp(-delta * t)))
    }
    
    llik_td <- function(theta, t, df) {
      mu <- exp(theta[1])
      delta <- exp(theta[2])
      a1 <- -mu / delta * (t - 1 / delta * (1 - exp(-delta * t)))
      a2 <- sum(log(mu / delta) + log(1 - exp(-delta * df$t)))
      a1 + a2
    }
    
    # initialization follows eventPred implementation
    tail_n <- df1u$n[df1u$t >= 0.75 * t0]
    beta <- (n0 - tail_n[1]) / (t0 / 4)
    mu0 <- 2 * n0 / t0^2
    delta0 <- mu0 / beta
    theta0 <- c(log(mu0), log(delta0))
    
    opt1 <- stats::optim(
      theta0, llik_td, gr = NULL, t = t0, df = df1,
      control = c(fnscale = -1),
      hessian = TRUE
    )
    
    fit1 <- list(
      model = "Time-decay",
      theta = opt1$par,
      vtheta = solve(-opt1$hessian),
      aic = -2 * opt1$value + 4,
      bic = -2 * opt1$value + 2 * log(n0)
    )
    
    dffit1 <- data.frame(t = seq(1, t0))
    dffit1$n <- fmu_td(dffit1$t, fit1$theta)
    
  } else if (enroll_model_lc == "b-spline") {
    
    K <- nknots
    days <- seq(1, t0)
    n_by_day <- as.numeric(table(factor(df1$t, levels = days)))
    
    x <- splines::bs(days, df = K + 4, intercept = TRUE)
    
    llik_bs <- function(theta, n, x) {
      lambda <- exp(as.vector(x %*% theta))
      -sum(lambda) + sum(n * log(lambda))
    }
    
    theta0 <- as.vector(solve(t(x) %*% x, t(x) %*% log(pmax(n_by_day, 0.1))))
    
    opt1 <- stats::optim(
      theta0, llik_bs, gr = NULL, n = n_by_day, x = x,
      control = c(fnscale = -1),
      hessian = TRUE
    )
    
    fit1 <- list(
      model = "B-spline",
      theta = opt1$par,
      vtheta = solve(-opt1$hessian),
      aic = -2 * opt1$value + 2 * (K + 4),
      bic = -2 * opt1$value + (K + 4) * log(n0),
      x = x
    )
    
    fmu_bs <- function(t, theta, x) {
      lambda <- exp(as.vector(x %*% theta))
      cumsum(lambda)[t]
    }
    
    dffit1 <- data.frame(t = seq(1, t0))
    dffit1$n <- fmu_bs(dffit1$t, fit1$theta, x)
    
  } else if (enroll_model_lc == "piecewise poisson") {
    
    u <- accrualTime[accrualTime < t0]
    u2 <- c(u, t0)
    
    factors <- cut(df1$t, breaks = u2)
    n_int <- as.numeric(table(factors))
    t_int <- diff(u2)
    
    vtheta <- if (length(u) > 1) diag(1 / n_int) else (1 / n_int) * diag(1)
    
    fit1 <- list(
      model = "Piecewise Poisson",
      theta = log(n_int / t_int),
      vtheta = vtheta,
      aic = -2 * sum(-n_int + n_int * log(n_int / t_int)) + 2 * length(u),
      bic = -2 * sum(-n_int + n_int * log(n_int / t_int)) + length(u) * log(n0),
      accrualTime = u
    )
    
    lambda <- n_int / t_int
    psum <- c(0, cumsum(n_int))
    time <- seq(1, t0)
    j <- findInterval(time, u)
    m <- psum[j] + lambda[j] * (time - u[j])
    
    dffit1 <- data.frame(t = time, n = m)
  }
  
  text_vec <- .ep_build_enrollment_fit_text(
    fit = fit1,
    accrual_time = accrualTime,
    nknots = nknots,
    criterion_lc = criterion_lc
  )
  
  # ---- Plot (optional) ----
  if (generate_plot) {
    
    fittedEnroll <- .ep_build_enrollment_fit_plot(
      enrolldf = df1u,
      dffit = dffit1,
      text_vec = text_vec,
      interactive_plot = interactive_plot
    )
    
    if (showplot) print(fittedEnroll)
    return(.ep_build_enrollment_fit_output(
      fit = fit1,
      fit_plot = fittedEnroll,
      enrolldf = df1u,
      dffit = dffit1,
      text_vec = text_vec
    ))
      
      # add 1–3 annotation lines depending on criterion
  }
  
  # no plot
  .ep_build_enrollment_fit_output(
    fit = fit1,
    enrolldf = df1u,
    dffit = dffit1,
    text_vec = text_vec
  )
}
