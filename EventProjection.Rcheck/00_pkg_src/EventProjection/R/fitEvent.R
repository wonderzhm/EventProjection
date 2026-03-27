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
    
    event_model_fit <- .ep_fit_event_model_branch(
      df1 = df1,
      formula = formula,
      event_model_lc = event_model_lc,
      n0 = n0,
      d0 = d0,
      ex0 = ex0,
      q = q,
      x = x,
      piecewise_survival_time = piecewiseSurvivalTime,
      k = k,
      scale = scale,
      m = m,
      covariates = covariates,
      cure_rate = cure_rate
    )
    fit2 <- event_model_fit$fit
    dffit2 <- event_model_fit$dffit
    p_cure <- event_model_fit$p_cure
    
    text_vec <- .ep_build_event_fit_text(
      fit = fit2,
      event_model_lc = event_model_lc,
      criterion_lc = criterion_lc,
      piecewise_survival_time = piecewiseSurvivalTime,
      k = k,
      scale = scale,
      m = m,
      p_cure = p_cure
    )
    
    # ---- Plot (optional) ----
    if (generate_plot) {
      fittedEvent <- .ep_build_event_fit_plot(
        kmdf = kmdf,
        kmdf_censor = kmdf_censor,
        dffit = dffit2,
        text_vec = text_vec,
        interactive_plot = interactive_plot
      )
    }
    
    # by-treatment labels (eventPred style)
    if (by_treatment) {
      if (generate_plot) {
        fittedEvent <- .ep_add_event_fit_treatment_label(
          plot_obj = fittedEvent,
          label = df1$treatment_description[1],
          interactive_plot = interactive_plot
        )
      }
      fit2$treatment <- df1$treatment[1]
      fit2$treatment_description <- df1$treatment_description[1]
    }
    
    # assemble return element
    event_fit[[i]] <- .ep_build_event_fit_output(
      fit = fit2,
      fit_plot = if (generate_plot) fittedEvent else NULL,
      kmdf = kmdf,
      dffit = dffit2,
      text_vec = text_vec
    )
  }
  
  # ensure that the sub plots share the same x axis range
  if (by_treatment) {
    if (generate_plot) {
      x_range <- range(dt$time, na.rm = TRUE)
      event_fit <- .ep_sync_event_fit_plot_ranges(
        event_fit = event_fit,
        x_range = x_range,
        interactive_plot = interactive_plot
      )
    }
  } else {
    # unwrap singleton to match eventPred behavior
    event_fit <- event_fit[[1]]
  }
  
  if (generate_plot && showplot) {
    if (by_treatment) {
      for (i in seq_along(event_fit)) print(event_fit[[i]]$fit_plot)
    } else {
      print(event_fit$fit_plot)
    }
  }
  
  event_fit
}
