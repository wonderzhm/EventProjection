#' @title Fit time-to-event model
#' @description Fits a specified time-to-event model to the event data.
#'
#' @details
#' This merged implementation is based on the original \code{fitEvent()} from the
#' \pkg{eventPred} package by Kaifeng Lu, with additional features adapted from
#' \pkg{EventPredInCure}:
#' \itemize{
#'   \item Added \code{criterion} ("aic", "bic", or "both") to control which
#'   information criteria are shown in the figure annotation.
#'   \item Added support for Chen (2016) style mixture-cure event models:
#'   "exponential with cured population", "weibull with cured population",
#'   "log-normal with cured population", "log-logistic with cured population",
#'   and "piecewise exponential with cured population".
#'   \item Added an optional fixed cure-rate feature via \code{cure_rate}: 
#'   when provided, the cure fraction is fixed rather than estimated for mixture-cure models; 
#'   when \code{by_treatment=TRUE}, \code{cure_rate} may be a scalar (same for all arms) 
#'   or a vector with one value per treatment arm.
#'   \item Added censoring markers to the fitted survival curve plots.
#'   \item Standardized column names in \code{df} to lower case to make the
#'   function more robust to input case.
#' }
#'
#' The core fitting logic, covariate support, Cox option, and return structure
#' follow \pkg{eventPred} to preserve downstream compatibility.
#'
#' @param df The subject-level event data, including \code{time} and \code{event}.
#'   The data should also include \code{treatment} coded as 1, 2, and so on, and
#'   optionally \code{treatment_description} for fitting by treatment.
#' @param event_model The event model used to analyze the event data.
#'   One of: "exponential", "weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", "spline", "cox",
#'   "exponential with cured population", "weibull with cured population",
#'   "log-normal with cured population", "log-logistic with cured population",
#'   "piecewise exponential with cured population".
#' @param cure_rate Optional fixed cured fraction for mixture-cure event models.
#'   If \code{NULL} (default), the cure rate is estimated from the data.
#'   If provided, it must be numeric in (0, 1). When \code{by_treatment=FALSE},
#'   supply a single value. When \code{by_treatment=TRUE}, \code{cure_rate} may be
#'   a single value (applied to all arms) or a vector with one value per treatment arm.
#' @param piecewiseSurvivalTime A vector that specifies the time intervals for
#'   the piecewise exponential survival distribution (with or without cure).
#'   Must start with 0, e.g., \code{c(0, 60)} breaks the time axis into
#'   intervals [0, 60) and [60, Inf). Default is 0.
#' @param k Number of inner knots for the spline model. Default 0.
#' @param scale Spline scale: "hazard", "odds", or "normal". Default "hazard".
#' @param m Number of event-time intervals to extrapolate the hazard beyond the
#'   last observed event time when \code{event_model="cox"}. Default 5.
#' @param showplot Logical; if \code{TRUE} and \code{generate_plot=TRUE}, prints plots.
#' @param by_treatment Logical; fit by treatment groups if \code{TRUE}. Default \code{FALSE}.
#' @param covariates Optional character vector of baseline covariate names to include.
#'   (Not supported for cured-population models in this merged version.)
#' @param generate_plot Logical; whether to generate plots. Default \code{TRUE}.
#' @param interactive_plot Logical; if \code{TRUE}, uses \pkg{plotly}; otherwise \pkg{ggplot2}.
#' @param criterion Character; which information criterion to display in plot annotations:
#'   "aic", "bic", or "both". Default "both".
#'
#' @return
#' If \code{by_treatment=FALSE}, a list with components:
#' \itemize{
#'   \item \code{fit}: model fit object (model, theta, vtheta, aic, bic, and extras).
#'   \item \code{fit_plot}: plot object if \code{generate_plot=TRUE}.
#'   \item \code{kmdf}: Kaplan-Meier curve data (includes \code{n.censor}).
#'   \item \code{dffit}: fitted survival curve data.
#'   \item \code{text}: character vector used in plot annotation.
#' }
#'
#' If \code{by_treatment=TRUE}, returns a list of such lists, one per treatment.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' \itemize{
#' \item Royston, P. and Parmar, M.K.B. (2002). Statistics in Medicine 21:2175-2197.
#' \item Chen, T.-T. (2016). BMC Medical Research Methodology 16:1-10.
#' }
#'
#' @examples
#'
#' event_fit <- fitEvent(
#'   df = interimData2,
#'   event_model = "piecewise exponential",
#'   piecewiseSurvivalTime = c(0, 180), showplot = FALSE)
#'
#' @export
fitEvent <- function(df, event_model = "model averaging",
                     cure_rate = NULL,
                     piecewiseSurvivalTime = 0,
                     k = 0, scale = "hazard", m = 5,
                     showplot = TRUE, by_treatment = FALSE,
                     covariates = NULL,
                     generate_plot = TRUE,
                     interactive_plot = TRUE,
                     criterion = "both") {
  
  # ---- Compatibility / validation ----
  erify::check_class(df, "data.frame")
  
  event_model_lc <- tolower(event_model)
  criterion_lc <- tolower(criterion)
  
  erify::check_content(
    event_model_lc,
    c("exponential", "weibull", "log-logistic",
      "log-normal", "piecewise exponential",
      "model averaging", "spline", "cox",
      "exponential with cured population",
      "weibull with cured population",
      "log-normal with cured population",
      "log-logistic with cured population",
      "piecewise exponential with cured population")
  )
  
  erify::check_content(criterion_lc, c("aic", "bic", "both"))
  
  if (piecewiseSurvivalTime[1] != 0) stop("piecewiseSurvivalTime must start with 0")
  if (length(piecewiseSurvivalTime) > 1 && any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }
  
  erify::check_n(k, zero = TRUE)
  erify::check_content(tolower(scale), c("hazard", "odds", "normal"))
  erify::check_n(m)
  
  erify::check_bool(showplot)
  erify::check_bool(by_treatment)
  erify::check_bool(generate_plot)
  erify::check_bool(interactive_plot)
  
  # ---- Standardize column names (merged change) ----
  dt <- data.table::as.data.table(data.table::copy(df))
  data.table::setnames(dt, tolower(names(dt)))
  
  if (!all(c("time", "event") %in% names(dt))) {
    stop("df must contain columns 'time' and 'event' (case-insensitive).")
  }
  if (by_treatment && !("treatment" %in% names(dt))) {
    stop("df must contain column 'treatment' (case-insensitive) when by_treatment=TRUE.")
  }
  # treat empty covariates vector as NULL (merged robustness)
  if (!is.null(covariates) && length(covariates) == 0) covariates <- NULL
  
  # ---- Cure models: currently no covariate support (merged design choice) ----
  is_cure_model <- grepl("with cured population$", event_model_lc)
  
  if (!is.null(cure_rate) && !is_cure_model) {
    stop("cure_rate is only applicable when event_model is a cured-population model (i.e., ends with 'with cured population').")
  }
  if (!is.null(cure_rate)) {
    if (!is.numeric(cure_rate) || any(is.na(cure_rate))) {
      stop("cure_rate must be numeric with no missing values, and each value must be in (0, 1).")
    }
    if (any(cure_rate <= 0) || any(cure_rate >= 1)) {
      stop("Each cure_rate value must be strictly between 0 and 1 (exclusive).")
    }
  }
  
  if (is_cure_model && !is.null(covariates)) {
    stop("Cured-population event models do not support covariates in this merged version.")
  }
  
  # construct the formula for survival analysis (eventPred behavior)
  if (!is.null(covariates)) {
    covariates <- tolower(covariates)
    if (!all(covariates %in% names(dt))) stop("All covariates must exist in df (case-insensitive).")
    xnames <- paste(covariates, collapse = "+")
    formula <- stats::as.formula(paste("survival::Surv(time, event) ~", xnames))
  } else {
    formula <- survival::Surv(time, event) ~ 1
  }
  
  if (by_treatment) {
    ngroups <- dt[, data.table::uniqueN(get("treatment"))]
    if (!("treatment_description" %in% names(dt))) {
      dt[, `:=`(treatment_description = paste("Treatment", get("treatment")))]
    }
  } else {
    ngroups <- 1
    dt[, `:=`(treatment = 1)]
  }
  
  if (ngroups == 1) by_treatment <- FALSE
  
  # ---- Fixed cure-rate handling (merged update) ----
  # cure_rate can be scalar (applied to all arms) or a vector (one per arm) when by_treatment=TRUE.
  cure_rate_vec <- NULL
  if (is_cure_model && !is.null(cure_rate)) {
    if (by_treatment) {
      if (length(cure_rate) == 1L) {
        cure_rate_vec <- rep(as.numeric(cure_rate), ngroups)
      } else if (length(cure_rate) == ngroups) {
        cure_rate_vec <- as.numeric(cure_rate)
      } else {
        stop("When by_treatment=TRUE, cure_rate must be a single value or a vector with length equal to the number of treatment groups (", ngroups, ").")
      }
    } else {
      if (length(cure_rate) != 1L) {
        stop("When by_treatment=FALSE, cure_rate must be a single value.")
      }
      cure_rate_vec <- rep(as.numeric(cure_rate), ngroups)
    }
  }
  
  
  # ---- Fit by treatment group ----
  event_fit <- list()
  
  for (i in 1:ngroups) {
    
    # Use arm-specific fixed cure rate if provided; otherwise estimate cure rate.
    if (!is.null(cure_rate_vec)) {
      cure_rate <- cure_rate_vec[i]
    } else {
      cure_rate <- NULL
    }
    
    df1 <- dt[get("treatment") == i]
    
    n0 <- nrow(df1)
    d0 <- df1[, sum(get("event"))]
    ex0 <- df1[, sum(get("time"))]
    
    erify::check_positive(d0, supplement = "The number of events must be positive to fit an event model.")
    
    x <- model.matrix(formula, df1)
    q <- ncol(x) - 1
    
    kmfit <- survival::survfit(survival::Surv(time, event) ~ 1, data = df1)
    kmdf <- data.table::data.table(time = kmfit$time, surv = kmfit$surv, n.censor = kmfit$n.censor)
    kmdf <- data.table::rbindlist(
      list(data.table::data.table(time = 0, surv = 1, n.censor = 0), kmdf),
      use.names = TRUE
    )
    kmdf_censor <- kmdf[get("n.censor") > 0]
    
    # ---- Fit model ----
    p_cure <- NULL  # cured fraction (only for cure models)
    
    if (event_model_lc == "exponential") {
      erify::check_positive(d0 - q, supplement = paste(
        "The number of events must be >=", q + 1, "to fit an exponential model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "exponential")
      
      fit2 <- list(model = "Exponential",
                   theta = -as.numeric(reg$coefficients),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 1),
                   bic = -2 * reg$loglik[1] + (q + 1) * log(n0))
      
      rate <- exp(as.numeric(x %*% fit2$theta))
      dffit2 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::pexp(t, rate, lower.tail = FALSE))))]
      
    } else if (event_model_lc == "weibull") {
      erify::check_positive(d0 - q - 1, supplement = paste(
        "The number of events must be >=", q + 2, "to fit a Weibull model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "weibull")
      
      fit2 <- list(model = "Weibull",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      shape <- exp(-fit2$theta[q + 2])
      scale2 <- exp(as.numeric(x %*% fit2$theta[1:(q + 1)]))
      
      dffit2 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::pweibull(t, shape, scale2, lower.tail = FALSE))))]
      
    } else if (event_model_lc == "log-logistic") {
      erify::check_positive(d0 - q - 1, supplement = paste(
        "The number of events must be >=", q + 2, "to fit a log-logistic model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "loglogistic")
      
      fit2 <- list(model = "Log-logistic",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      location <- as.numeric(x %*% fit2$theta[1:(q + 1)])
      scale2 <- exp(fit2$theta[q + 2])
      
      dffit2 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::plogis(log(t), location, scale2, lower.tail = FALSE))))]
      
    } else if (event_model_lc == "log-normal") {
      erify::check_positive(d0 - q - 1, supplement = paste(
        "The number of events must be >=", q + 2, "to fit a log-normal model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "lognormal")
      
      fit2 <- list(model = "Log-normal",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      meanlog <- as.numeric(x %*% fit2$theta[1:(q + 1)])
      sdlog <- exp(fit2$theta[q + 2])
      
      dffit2 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::plnorm(t, meanlog, sdlog, lower.tail = FALSE))))]
      
    } else if (event_model_lc == "piecewise exponential") {
      J <- length(piecewiseSurvivalTime)
      
      erify::check_positive(d0 - J - q + 1, supplement = paste(
        "The number of events must be >=", J + q, "to fit a piecewise exponential model."))
      
      fit2 <- pwexpreg(df1$time, df1$event, J, piecewiseSurvivalTime, q, x)
      
      time <- seq(0, max(df1$time))
      surv <- purrr::map(1:n0, function(l)
        ppwexp(time, fit2$theta, J, fit2$piecewiseSurvivalTime,
               q, x[l,], lower.tail = FALSE))
      surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      
      dffit2 <- data.table::data.table(time, surv)
      
    } else if (event_model_lc == "model averaging") {
      erify::check_positive(d0 - q - 1, supplement = paste(
        "The number of events must be >=", q + 2, "to fit a model averaging model."))
      
      reg1 <- survival::survreg(formula, data = df1, dist = "weibull")
      reg2 <- survival::survreg(formula, data = df1, dist = "lognormal")
      aic1 <- -2 * reg1$loglik[1] + 2 * (q + 2)
      aic2 <- -2 * reg2$loglik[1] + 2 * (q + 2)
      bic1 <- -2 * reg1$loglik[1] + (q + 2) * log(n0)
      bic2 <- -2 * reg2$loglik[1] + (q + 2) * log(n0)
      
      w1 <- 1 / (1 + exp(-0.5 * (bic2 - bic1)))
      
      theta <- c(as.numeric(reg1$coefficients), log(reg1$scale),
                 as.numeric(reg2$coefficients), log(reg2$scale))
      vtheta <- as.matrix(Matrix::bdiag(reg1$var, reg2$var))
      
      fit2 <- list(model = "Model averaging",
                   theta = theta,
                   vtheta = vtheta,
                   aic = w1 * aic1 + (1 - w1) * aic2,
                   bic = w1 * bic1 + (1 - w1) * bic2,
                   w1 = w1)
      
      time <- seq(0, max(df1$time))
      surv <- purrr::map(1:n0, function(l)
        pmodavg(time, fit2$theta, w1, q, x[l,], lower.tail = FALSE))
      surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      
      dffit2 <- data.table::data.table(time, surv)
      
    } else if (event_model_lc == "spline") {
      erify::check_positive(d0 - k - q - 1, supplement = paste(
        "The number of events must be >=", k + q + 2, "to fit a spline model."))
      
      spl <- flexsurv::flexsurvspline(
        formula, data = df1, k = k, scale = scale, method = "Nelder-Mead")
      
      fit2 <- list(model = "Spline",
                   theta = as.numeric(spl$coefficients),
                   vtheta = spl$cov,
                   aic = -2 * spl$loglik + 2 * (k + q + 2),
                   bic = -2 * spl$loglik + (k + q + 2) * log(n0),
                   knots = spl$knots,
                   scale = spl$scale)
      
      time <- seq(0, max(df1$time))
      
      if (q > 0) {
        xbeta <- as.numeric(as.matrix(x[, -1]) %*%
                              fit2$theta[(k + 3):(k + q + 2)])
        
        surv <- purrr::map(1:n0, function(l)
          flexsurv::psurvspline(
            time, gamma = fit2$theta[1:(k + 2)], knots = fit2$knots,
            scale = fit2$scale, offset = xbeta[l], lower.tail = FALSE))
        surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      } else {
        surv <- flexsurv::psurvspline(
          time, gamma = fit2$theta, knots = fit2$knots, scale = fit2$scale,
          lower.tail = FALSE)
      }
      
      dffit2 <- data.table::data.table(time, surv)
      
    } else if (event_model_lc == "cox") {
      erify::check_positive(d0 - q, supplement = paste(
        "The number of events must be >=", q + 1, "to fit a Cox model."))
      erify::check_positive(d0 - m + 1, supplement = paste("m must be <= the observed number of events", d0))
      
      reg <- phregr(df1, time = "time", event = "event", covariates = covariates)
      
      bh <- data.table::setDT(reg$basehaz)[get("nevent") > 0]
      haz <- bh$haz
      vhaz <- bh$varhaz
      
      if (q > 0) {
        if (q == 1) {
          ghaz <- matrix(bh$gradhaz, ncol = 1)
        } else {
          ghaz <- do.call(cbind, lapply(1:q, function(ii) bh[[paste0("gradhaz.", ii)]]))
        }
      }
      
      M <- nrow(bh)
      tcut <- c(0, bh$time)
      dtcut <- diff(tcut)
      lambda1 <- haz / dtcut
      
      d <- bh$nevent
      llik <- reg$sumstat$loglik1 + sum(d * (log(d / dtcut) - 1))
      
      if (q > 0) {
        theta <- c(log(lambda1), as.numeric(reg$beta))
        vbeta <- reg$vbeta
        dimnames(vbeta) <- NULL
        vtheta <- matrix(0, M + q, M + q)
        vtheta[(M + 1):(M + q), (M + 1):(M + q)] <- vbeta
        vtheta[1:M, 1:M] <- diag(vhaz / (haz * haz)) + ghaz %*% vbeta %*% t(ghaz)
        vtheta[1:M, (M + 1):(M + q)] <- ghaz %*% vbeta
        vtheta[(M + 1):(M + q), 1:M] <- t(vtheta[1:M, (M + 1):(M + q)])
      } else {
        theta <- log(lambda1)
        vtheta <- diag(vhaz / (haz * haz))
      }
      
      fit2 <- list(model = "Cox",
                   theta = theta,
                   vtheta = vtheta,
                   aic = -2 * llik + 2 * (M + q),
                   bic = -2 * llik + (M + q) * log(n0),
                   piecewiseSurvivalTime = tcut)
      
      time <- seq(0, max(df1$time))
      lambda2 <- sum(bh$haz[(M - m + 1):M]) / (bh$time[M] - bh$time[M - m])
      lambda <- c(lambda1, lambda2)
      
      s1 <- sapply(time, function(t)
        ppwexp(t, log(lambda), M + 1, tcut, lower.tail = FALSE))
      
      if (q > 0) {
        xbeta <- as.numeric(as.matrix(x[, -1]) %*% reg$beta)
        surv <- apply(outer(s1, exp(xbeta), `^`), 1, mean)
      } else {
        surv <- s1
      }
      
      dffit2 <- data.table::data.table(time, surv)
    } else if (event_model_lc == "exponential with cured population") {
      # Chen (2016) mixture-cure exponential model (EventPredInCure logic)
      if (is.null(cure_rate)) {
        temp <- stats::nlminb(
          start = c(-2, log(d0 / ex0)),
          objective = loglik_Chen_exponential,
          df = as.data.frame(df1)
        )
        theta <- temp$par
        
        fit2 <- list(
          model = "exponential with cured population",
          theta = theta,
          vtheta = MASS::ginv(numDeriv::hessian(loglik_Chen_exponential, theta, df = as.data.frame(df1))),
          bic = 2 * temp$objective + log(n0) * 2,
          aic = 2 * temp$objective + 2 * 2
        )
        
        p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
      } else {
        temp <- stats::nlminb(
          start = c(log(d0 / ex0)),
          objective = loglik_Chen_exponential2,
          df = as.data.frame(df1),
          p = cure_rate
        )
        
        theta <- c(stats::qlogis(cure_rate), temp$par)
        
        vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
        vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
          loglik_Chen_exponential2,
          temp$par,
          df = as.data.frame(df1),
          p = cure_rate
        ))
        
        fit2 <- list(
          model = "exponential with cured population",
          theta = theta,
          vtheta = vtheta,
          bic = 2 * temp$objective + log(n0) * 1,
          aic = 2 * temp$objective + 2 * 1
        )
        
        p_cure <- cure_rate
      }
      
      tt <- seq(0, max(df1$time))
      dffit2 <- data.table::data.table(time = tt, surv = SP_Chen_exponential(theta, tt))
      
    } else if (event_model_lc == "weibull with cured population") {
      # Chen (2016) mixture-cure Weibull model (EventPredInCure logic)
      reg0 <- survival::survreg(survival::Surv(time, event) ~ 1, data = df1, dist = "weibull")
      
      # Note: loglik_Chen_weibull* parameterizes Weibull by (log(shape), log(scale)).
      # A good starting value for log(shape) is log(1/reg0$scale).
      if (is.null(cure_rate)) {
        temp <- stats::nlminb(
          start = c(-2, log(1 / reg0$scale), as.numeric(reg0$coeff)),
          objective = loglik_Chen_weibull,
          df = as.data.frame(df1)
        )
        theta <- temp$par
        
        fit2 <- list(
          model = "weibull with cured population",
          theta = theta,
          vtheta = MASS::ginv(numDeriv::hessian(loglik_Chen_weibull, theta, df = as.data.frame(df1))),
          bic = 2 * temp$objective + log(n0) * 3,
          aic = 2 * temp$objective + 2 * 3
        )
        
        p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
      } else {
        temp <- stats::nlminb(
          start = c(log(1 / reg0$scale), as.numeric(reg0$coeff)),
          objective = loglik_Chen_weibull2,
          df = as.data.frame(df1),
          p = cure_rate
        )
        
        theta <- c(stats::qlogis(cure_rate), temp$par)
        
        vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
        vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
          loglik_Chen_weibull2,
          temp$par,
          df = as.data.frame(df1),
          p = cure_rate
        ))
        
        fit2 <- list(
          model = "weibull with cured population",
          theta = theta,
          vtheta = vtheta,
          bic = 2 * temp$objective + log(n0) * 2,
          aic = 2 * temp$objective + 2 * 2
        )
        
        p_cure <- cure_rate
      }
      
      tt <- seq(0, max(df1$time))
      dffit2 <- data.table::data.table(time = tt, surv = SP_Chen_weibull(theta, tt))
      
    } else if (event_model_lc == "log-logistic with cured population") {
      # Chen (2016) mixture-cure log-logistic model (EventPredInCure logic)
      reg0 <- survival::survreg(survival::Surv(time, event) ~ 1, data = df1, dist = "loglogistic")
      
      if (is.null(cure_rate)) {
        temp <- stats::nlminb(
          start = c(-2, log(1 / reg0$scale), as.numeric(reg0$coeff)),
          objective = loglik_Chen_log_logistic,
          df = as.data.frame(df1)
        )
        theta <- temp$par
        
        fit2 <- list(
          model = "log-logistic with cured population",
          theta = theta,
          vtheta = MASS::ginv(numDeriv::hessian(loglik_Chen_log_logistic, theta, df = as.data.frame(df1))),
          bic = 2 * temp$objective + log(n0) * 3,
          aic = 2 * temp$objective + 2 * 3
        )
        
        p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
      } else {
        temp <- stats::nlminb(
          start = c(log(1 / reg0$scale), as.numeric(reg0$coeff)),
          objective = loglik_Chen_log_logistic2,
          df = as.data.frame(df1),
          p = cure_rate
        )
        
        theta <- c(stats::qlogis(cure_rate), temp$par)
        
        vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
        vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
          loglik_Chen_log_logistic2,
          temp$par,
          df = as.data.frame(df1),
          p = cure_rate
        ))
        
        fit2 <- list(
          model = "log-logistic with cured population",
          theta = theta,
          vtheta = vtheta,
          bic = 2 * temp$objective + log(n0) * 2,
          aic = 2 * temp$objective + 2 * 2
        )
        
        p_cure <- cure_rate
      }
      
      tt <- seq(0, max(df1$time))
      dffit2 <- data.table::data.table(time = tt, surv = SP_Chen_log_logistic(theta, tt))
      
    } else if (event_model_lc == "log-normal with cured population") {
      # Chen (2016) mixture-cure log-normal model (EventPredInCure logic)
      reg0 <- survival::survreg(survival::Surv(time, event) ~ 1, data = df1, dist = "lognormal")
      
      if (is.null(cure_rate)) {
        temp <- stats::nlminb(
          start = c(-2, as.numeric(reg0$coeff), log(reg0$scale)),
          objective = loglik_Chen_log_normal,
          df = as.data.frame(df1)
        )
        theta <- temp$par
        
        fit2 <- list(
          model = "log-normal with cured population",
          theta = theta,
          vtheta = MASS::ginv(numDeriv::hessian(loglik_Chen_log_normal, theta, df = as.data.frame(df1))),
          bic = 2 * temp$objective + log(n0) * 3,
          aic = 2 * temp$objective + 2 * 3
        )
        
        p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
      } else {
        temp <- stats::nlminb(
          start = c(as.numeric(reg0$coeff), log(reg0$scale)),
          objective = loglik_Chen_log_normal2,
          df = as.data.frame(df1),
          p = cure_rate
        )
        
        theta <- c(stats::qlogis(cure_rate), temp$par)
        
        vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
        vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
          loglik_Chen_log_normal2,
          temp$par,
          df = as.data.frame(df1),
          p = cure_rate
        ))
        
        fit2 <- list(
          model = "log-normal with cured population",
          theta = theta,
          vtheta = vtheta,
          bic = 2 * temp$objective + log(n0) * 2,
          aic = 2 * temp$objective + 2 * 2
        )
        
        p_cure <- cure_rate
      }
      
      tt <- seq(0, max(df1$time))
      dffit2 <- data.table::data.table(time = tt, surv = SP_Chen_log_normal(theta, tt))
      
    } else if (event_model_lc == "piecewise exponential with cured population") {
      # Chen (2016) mixture-cure piecewise exponential model (EventPredInCure logic)
      u <- piecewiseSurvivalTime[piecewiseSurvivalTime < max(df1$time)]
      ucut <- c(u, max(df1$time))
      J <- length(u)
      
      # Use J hazards corresponding to the J intervals defined by ucut.
      haz_start <- rep(log(0.0030), J)
      
      if (is.null(cure_rate)) {
        temp <- stats::nlminb(
          start = c(-2, haz_start),
          objective = loglik_Chen_piecewise_exponential,
          df = as.data.frame(df1),
          piecewiseSurvivalTime = piecewiseSurvivalTime
        )
        theta <- temp$par
        
        fit2 <- list(
          model = "piecewise exponential with cured population",
          theta = theta,
          vtheta = MASS::ginv(numDeriv::hessian(
            loglik_Chen_piecewise_exponential,
            theta,
            df = as.data.frame(df1),
            piecewiseSurvivalTime = piecewiseSurvivalTime
          )),
          bic = 2 * temp$objective + (J + 1) * log(n0),
          aic = 2 * temp$objective + (J + 1) * 2,
          piecewiseSurvivalTime = u
        )
        
        p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
      } else {
        temp <- stats::nlminb(
          start = haz_start,
          objective = loglik_Chen_piecewise_exponential2,
          df = as.data.frame(df1),
          piecewiseSurvivalTime = piecewiseSurvivalTime,
          p = cure_rate
        )
        
        theta <- c(stats::qlogis(cure_rate), temp$par)
        
        vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
        vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
          loglik_Chen_piecewise_exponential2,
          temp$par,
          df = as.data.frame(df1),
          piecewiseSurvivalTime = piecewiseSurvivalTime,
          p = cure_rate
        ))
        
        fit2 <- list(
          model = "piecewise exponential with cured population",
          theta = theta,
          vtheta = vtheta,
          bic = 2 * temp$objective + J * log(n0),
          aic = 2 * temp$objective + J * 2,
          piecewiseSurvivalTime = u
        )
        
        p_cure <- cure_rate
      }
      
      tt <- seq(0, max(df1$time))
      dffit2 <- data.table::data.table(
        time = tt,
        surv = SP_Chen_piecewise_exponential(theta, tt, piecewiseSurvivalTime = piecewiseSurvivalTime)
      )
      
    } else {
      stop("Unsupported event_model.")
    }
    
    # ---- Build plot annotation text (merged change: criterion + cure fraction) ----
    if (grepl("^piecewise exponential", tolower(fit2$model))) {
      ps <- if (!is.null(fit2$piecewiseSurvivalTime)) fit2$piecewiseSurvivalTime else piecewiseSurvivalTime
      modeltext <- paste0(fit2$model, "(", paste(ps, collapse = ","), ")")
      
    } else if (tolower(fit2$model) == "spline") {
      modeltext <- paste0(fit2$model, "(k = ", k, ", scale = '", scale, "')")
    } else if (tolower(fit2$model) == "cox") {
      modeltext <- paste0(fit2$model, "(m = ", m, ")")
    } else {
      modeltext <- fit2$model
    }
    
    if (event_model_lc == "model averaging") {
      aictext <- paste("Weighted AIC:", formatC(fit2$aic, format = "f", digits = 2))
      bictext <- paste("Weighted BIC:", formatC(fit2$bic, format = "f", digits = 2))
    } else {
      aictext <- paste("AIC:", formatC(fit2$aic, format = "f", digits = 2))
      bictext <- paste("BIC:", formatC(fit2$bic, format = "f", digits = 2))
    }
    
    text_vec <- c(modeltext)
    if (criterion_lc %in% c("aic", "both")) text_vec <- c(text_vec, aictext)
    if (criterion_lc %in% c("bic", "both")) text_vec <- c(text_vec, bictext)
    if (!is.null(p_cure)) text_vec <- c(text_vec, paste("Prop of Cured:", formatC(p_cure, format = "f", digits = 2)))
    
    # ---- Plot (optional) ----
    if (generate_plot) {
      if (interactive_plot) {
        y_seq <- 0.95 - 0.05 * (0:(length(text_vec) - 1))
        fittedEvent <- plotly::plot_ly() %>%
          plotly::add_lines(
            data = kmdf, x = ~time, y = ~surv, name = "Kaplan-Meier",
            line = list(shape = "hv")) %>%
          plotly::add_markers(
            data = kmdf_censor, x = ~time, y = ~surv, name = "censoring",
            showlegend = FALSE,
            marker = list(symbol = "cross-thin", size = 6,
                          line = list(color = "black", width = 2))) %>%
          plotly::add_lines(
            data = dffit2, x = ~time, y = ~surv, name = "fitted") %>%
          plotly::layout(
            xaxis = list(title = "Days since randomization", zeroline = FALSE),
            yaxis = list(title = "Survival probability", zeroline = FALSE),
            title = list(text = "Fitted time to event survival curve"),
            annotations = list(
              x = rep(0.7, length(text_vec)),
              y = y_seq, xref = "paper", yref = "paper",
              text = paste0("<i>", text_vec, "</i>"),
              xanchor = "left", font = list(size = 14, color = "red"),
              showarrow = FALSE)) %>%
          plotly::hide_legend()
        
      } else {
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
          stop("Package 'ggplot2' is required for interactive_plot=FALSE.")
        }
        
        fittedEvent <- ggplot2::ggplot() +
          ggplot2::geom_step(data = kmdf[-1, ], ggplot2::aes(x = .data$time, y = .data$surv)) +
          ggplot2::geom_point(data = kmdf_censor, ggplot2::aes(x = .data$time, y = .data$surv),
                              shape = 4, inherit.aes = FALSE) +
          ggplot2::geom_line(data = dffit2[-1, ], ggplot2::aes(x = .data$time, y = .data$surv),
                             colour = "red") +
          ggplot2::labs(
            x = "Days since randomization",
            y = "Survival probability",
            title = "Fitted time to event survival curve") +
          ggplot2::theme(legend.position = "none")
        
        x_pos <- max(kmdf$time, na.rm = TRUE)
        y_pos <- max(kmdf$surv, na.rm = TRUE)
        
        for (j in seq_along(text_vec)) {
          fittedEvent <- fittedEvent +
            ggplot2::annotate("text", x = x_pos, y = y_pos, label = text_vec[j],
                              hjust = 1.1, vjust = 1.5 + 2 * (j - 1),
                              size = 5, colour = "red", fontface = "italic")
        }
      }
    }
    
    # by-treatment labels (eventPred style)
    if (by_treatment) {
      if (generate_plot) {
        if (interactive_plot) {
          fittedEvent <- fittedEvent %>%
            plotly::layout(annotations = c(
              list(
                x = 0.5, y = 1,
                text = paste0("<b>", df1$treatment_description[1], "</b>"),
                xanchor = "center", yanchor = "middle", showarrow = FALSE,
                xref = "paper", yref = "paper"
              ),
              fittedEvent$x$layout$annotations
            ))
        } else {
          fittedEvent <- fittedEvent + ggplot2::labs(subtitle = df1$treatment_description[1])
        }
      }
      fit2$treatment <- df1$treatment[1]
      fit2$treatment_description <- df1$treatment_description[1]
    }
    
    # assemble return element
    if (generate_plot) {
      event_fit[[i]] <- list(fit = fit2, fit_plot = fittedEvent,
                             kmdf = kmdf, dffit = dffit2,
                             text = text_vec)
    } else {
      event_fit[[i]] <- list(fit = fit2,
                             kmdf = kmdf, dffit = dffit2,
                             text = text_vec)
    }
  }
  
  # ensure that the sub plots share the same x axis range
  if (by_treatment) {
    if (generate_plot) {
      x_range <- range(dt$time, na.rm = TRUE)
      for (i in 1:ngroups) {
        
        # Use arm-specific fixed cure rate if provided; otherwise estimate cure rate.
        if (!is.null(cure_rate_vec)) {
          cure_rate <- cure_rate_vec[i]
        } else {
          cure_rate <- NULL
        }
        
        if (interactive_plot) {
          event_fit[[i]]$fit_plot <- event_fit[[i]]$fit_plot %>%
            plotly::layout(xaxis = list(range = x_range))
        } else {
          event_fit[[i]]$fit_plot <- event_fit[[i]]$fit_plot +
            ggplot2::scale_x_continuous(limits = x_range)
        }
      }
    }
  } else {
    # unwrap singleton to match eventPred behavior
    event_fit <- event_fit[[1]]
  }
  
  if (generate_plot && showplot) {
    if (by_treatment) {
      for (i in 1:length(event_fit)) print(event_fit[[i]]$fit_plot)
    } else {
      print(event_fit$fit_plot)
    }
  }
  
  event_fit
}