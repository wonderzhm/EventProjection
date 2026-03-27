#' @title Enrollment and event prediction
#' @description Performs enrollment and event prediction by utilizing
#'   observed data and specified enrollment and event models.
#'
#' @param df The subject-level enrollment and event data, including
#'   \code{trialsdt}, \code{usubjid}, \code{randdt}, and \code{cutoffdt} for
#'   enrollment prediction, and, additionally, \code{time}, \code{event},
#'   and \code{dropout} for event prediction. The data should also include
#'   \code{treatment} coded as 1, 2, and so on, and
#'   \code{treatment_description} for enrollment and
#'   event prediction by treatment. By default, it is set to
#'   \code{NULL} for enrollment and event prediction at the design stage.
#' @param to_predict Specifies what to predict: "enrollment only", "event
#'   only", or "enrollment and event". By default, it is set to
#'   "enrollment and event". When \code{df = NULL}, design-stage
#'   \code{"event only"} prediction is not supported through
#'   \code{getPrediction()}.
#' @param target_n The target number of subjects to enroll in the study.
#' @param target_d The target number of events to reach in the study.
#' @param enroll_model The enrollment model which can be specified as
#'   "Poisson", "Time-decay", "B-spline",
#'   "Piecewise Poisson", or "Piecewise Uniform". UPDATED (merge):
#'   add "Piecewise Uniform" (from EventPredInCure). By default, it is set to "B-spline".
#' @param nknots The number of inner knots for the B-spline enrollment
#'   model. By default, it is set to 0.
#' @param lags The day lags to compute the average enrollment rate to
#'   carry forward for the B-spline enrollment model. By default,
#'   it is set to 30.
#' @param accrualTime The accrual time intervals for the piecewise
#'   Poisson model. Must start with 0, e.g., c(0, 30) breaks the
#'   time axis into 2 accrual intervals: [0, 30) and [30, Inf).
#'   By default, it is set to 0.
#' @param enroll_prior The prior of enrollment model parameters.
#' @param event_model The event model used to analyze the event data
#'   which can be set to one of the following options:
#'   "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", "spline", or "cox".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. By default, it is set to
#'   "model averaging". UPDATED (merge): also supports Chen-style mixture-cure event models
#'   ("exponential with cured population","weibull with cured population",
#'   "log-normal with cured population","log-logistic with cured population", or
#'   "piecewise exponential with cured population") and delayed-treatment variants
#'   ("exponential with cured population and delayed treatment",
#'   "weibull with cured population and delayed treatment",
#'   "log-normal with cured population and delayed treatment",
#'   "log-logistic with cured population and delayed treatment") for design-stage
#'   simulation.
#' @param piecewiseSurvivalTime A vector that specifies the time
#'   intervals for the piecewise exponential survival distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k The number of inner knots of the spline event model of
#'   Royston and Parmar (2002). The default
#'   \code{k=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale If "hazard", the log cumulative hazard is modeled
#'   as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param m The number of event time intervals to extrapolate the hazard
#'   function beyond the last observed event time.
#' @param event_prior The prior of event model parameters.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options:
#'   "none", "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", "spline", or "cox".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. By default, it is set to
#'   "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k_dropout The number of inner knots of the spline dropout model of
#'   Royston and Parmar (2002). The default
#'   \code{k_dropout=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale_dropout} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale_dropout If "hazard", the log cumulative hazard for dropout
#'   is modeled as a spline function. If "odds", the log cumulative odds is
#'   modeled as a spline function. If "normal", -qnorm(S(t)) is
#'   modeled as a spline function.
#' @param m_dropout The number of dropout time intervals to extrapolate
#'   the hazard function beyond the last observed dropout time.
#' @param dropout_prior The prior of dropout model parameters.
#' @param fixedFollowup A Boolean variable indicating whether a fixed
#'   follow-up design is used. By default, it is set to \code{FALSE}
#'   for a variable follow-up design.
#' @param followupTime The follow-up time for a fixed
#'   follow-up design, in days. By default, it is set to 365.
#' @param pilevel The prediction interval level. By default,
#'   it is set to 0.90.
#' @param nyears The number of years after the data cut for prediction.
#'   By default, it is set to 4.
#' @param target_t The target number of days after the data cutoff
#'   used to predict both the number of events and the probability
#'   of achieving the target event count.
#' @param nreps The number of replications for simulation. By default,
#'   it is set to 500.
#' @param showEnrollment A Boolean variable to control whether or not to
#'   show the number of enrolled subjects. By default, it is set to
#'   \code{TRUE}.
#' @param showEvent A Boolean variable to control whether or not to
#'   show the number of events. By default, it is set to
#'   \code{TRUE}.
#' @param showDropout A Boolean variable to control whether or not to
#'   show the number of dropouts. By default, it is set to
#'   \code{FALSE}.
#' @param showOngoing A Boolean variable to control whether or not to
#'   show the number of ongoing subjects. By default, it is set to
#'   \code{FALSE}.
#' @param showsummary A Boolean variable to control whether or not to
#'   show the prediction summary. By default, it is set to \code{TRUE}.
#' @param showplot A Boolean variable to control whether or not to
#'   show the plots. By default, it is set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   predict by treatment group. By default, it is set to \code{FALSE}.
#' @param ngroups The number of treatment groups for enrollment prediction
#'   at the design stage. By default, it is set to 1.
#'   It is replaced with the actual number of
#'   treatment groups in the observed data if \code{df} is not \code{NULL}.
#' @param alloc The treatment allocation in a randomization block.
#'   By default, it is set to \code{NULL}, which yields equal allocation
#'   among the treatment groups.
#' @param treatment_label The treatment labels for treatments in a
#'   randomization block for design stage prediction.
#'   It is replaced with the treatment_description
#'   in the observed data if \code{df} is not \code{NULL}.
#' @param covariates_event The names of baseline covariates from the input
#'   data frame to include in the event model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param event_prior_with_covariates The prior of event model parameters
#'   in the presence of covariates.
#' @param covariates_dropout The names of baseline covariates from the input
#'   data frame to include in the dropout model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param dropout_prior_with_covariates The prior of dropout model
#'   parameters in the presence of covariates.
#' @param fix_parameter Whether to fix parameters at the maximum
#'   likelihood estimates when generating new data for prediction.
#'   Defaults to \code{FALSE}, in which case, parameters will be drawn from
#'   their approximate posterior distribution.
#' @param generate_plot Whether to generate plots.
#' @param interactive_plot Whether to produce interactive plots using
#'   plotly or static plots using ggplot2.
#' @param criterion Information criterion(s) used in figure annotations when
#'   model comparison is involved: "aic", "bic", or "both". Default is "both".
#'   UPDATED (merge): adapted from EventPredInCure.
#' @param seed.num Optional random seed for simulation-based prediction
#'   (enrollment and/or event). If NULL, no seed is set. UPDATED (merge):
#'   adapted from EventPredInCure.
#' @param cure_rate Optional fixed cure rate for Chen-style mixture-cure event
#'   models (e.g., "exponential with cured population"). If provided, the cure
#'   fraction is fixed at this value rather than estimated. When
#'   \code{by_treatment=TRUE}, \code{cure_rate} may be a numeric vector of length
#'   \code{ngroups}, with each element corresponding to one treatment arm.
#'   UPDATED (merge): new feature to support fixed cure rates.
#' @param return_subject_data Logical; if \code{TRUE} (default), keep the
#'   combined subject-level simulation output in the returned object. Setting
#'   this to \code{FALSE} can materially reduce memory use for large
#'   \code{nreps}.
#' @param return_simulation_data Logical; if \code{TRUE} (default), keep nested
#'   subject-level simulated data such as \code{newSubjects} and
#'   \code{newEvents} in \code{enroll_pred} and \code{event_pred}. Setting this
#'   to \code{FALSE} reduces the size of the returned object.
#'
#' @details
#' For the time-decay model, the mean function is
#' \eqn{\mu(t) = (\mu/\delta)(t - (1/\delta)(1 - \exp(-\delta t)))}
#' and the rate function is
#' \eqn{\lambda(t) = (\mu/\delta)(1 - \exp(-\delta t))}.
#' For the B-spline model, the daily enrollment rate is approximated as
#' \eqn{\lambda(t) = \exp(B(t)' \theta)},
#' where \code{B(t)} represents the B-spline basis functions.
#'
#' The \code{enroll_prior} variable should be a list that
#' includes \code{model} to specify the enrollment model
#' (poisson, time-decay, or piecewise poisson),
#' \code{theta} and \code{vtheta} to indicate the parameter
#' values and the covariance matrix. One can use a very small
#' value of \code{vtheta} to fix the parameter values.
#' For the piecewise Poisson enrollment model, the list
#' should also include \code{accrualTime}. It should be noted
#' that the B-spline model is not appropriate for use as prior.
#'
#' For event prediction by treatment with prior information,
#' the \code{event_prior} (\code{dropout_prior}) variable should be
#' a list with one element per treatment. For each treatment, the
#' element should include \code{model} to specify the event (dropout)
#' model (exponential, weibull, log-logistic, log-normal,
#' or piecewise exponential), and \code{theta} and \code{vtheta} to
#' indicate the parameter values and the covariance matrix.
#' For the piecewise exponential event (dropout) model, the list
#' should also include \code{piecewiseSurvivalTime}
#' (\code{piecewiseDropoutTime}) to indicate the location of knots.
#' It should be noted that the model averaging, spline, and
#' cox options are not appropriate for use as prior.
#'
#' If the event prediction is not by treatment while the prior
#' information is given by treatment, then each element of
#' \code{event_prior} (\code{dropout_prior}) should also include
#' \code{w} to specify the weight of the treatment in a
#' randomization block. If the prediction is not by treatment and
#' the prior is given for the overall study, then \code{event_prior}
#' (\code{dropout_prior}) is a flat list with \code{model},
#' \code{theta}, and \code{vtheta}. For the piecewise exponential
#' event (dropout) model, it should also include
#' \code{piecewiseSurvivalTime} (\code{piecewiseDropoutTime}) to
#' indicate the location of knots.
#'
#' For analysis-stage enrollment and event prediction, the
#' \code{enroll_prior}, \code{event_prior}, and
#' \code{dropout_prior} are either set to \code{NULL} to
#' use the observed data only, or specify the prior distribution
#' of model parameters to be combined with observed data likelihood
#' for enhanced modeling flexibility.
#'
#' @note Origin (merge): merged from eventPred::getPrediction and EventPredInCure::getPrediction.
#' Please retain original package authorship/credit per upstream licenses.
#' @return A list that includes the fits of observed data models,
#' as well as simulated enrollment data for new subjects and
#' simulated event data for ongoing and new subjects. Large returned simulation
#' objects can be omitted with \code{return_subject_data} and
#' \code{return_simulation_data}.
#'
#'
#' This merged implementation is based on the original \code{getPrediction()}
#' from the \pkg{eventPred} package by Kaifeng Lu, with additional features
#' adapted from \pkg{EventPredInCure}. Key merged updates include:
#' \itemize{
#'   \item Added \code{criterion} ("aic", "bic", or "both") to control which
#'   information criteria are shown in figure annotations (passed through to
#'   \code{fitEnrollment()}, \code{fitEvent()}, \code{fitDropout()} when supported).
#'   \item Added support for Chen (2016) style mixture-cure and delayed-treatment
#'   event models in input validation and design-stage simulation.
#'   \item Added optional \code{seed.num} for reproducible simulation-based predictions.
#'   \item Standardized \code{df} column names to lower case for robustness.
#'   \item UPDATED (merge): added an optional fixed \code{cure_rate} interface
#'   (passed to \code{fitEvent()}) for cured models.
#' }
#'
#' @examples
#' # Event prediction after enrollment completion
#' set.seed(3000)
#'
#' pred <- getPrediction(
#'   df = interimData2, to_predict = "event only",
#'   target_d = 200,
#'   event_model = "weibull",
#'   dropout_model = "exponential",
#'   pilevel = 0.90, nreps = 100, showplot = FALSE)
#'
#' @export
#'
getPrediction <- function(
    df = NULL, to_predict = "enrollment and event",
    target_n = NA, target_d = NA,
    enroll_model = "b-spline", nknots = 0, lags = 30,
    accrualTime = 0,
    enroll_prior = NULL,
    event_model = "model averaging", piecewiseSurvivalTime = 0,
    k = 0, scale = "hazard", m = 5,
    event_prior = NULL,
    dropout_model = "exponential", piecewiseDropoutTime = 0,
    k_dropout = 0, scale_dropout = "hazard", m_dropout = 5,
    dropout_prior = NULL,
    fixedFollowup = FALSE, followupTime = 365,
    pilevel = 0.90, nyears = 4,
    target_t = NA, nreps = 500,
    showEnrollment = TRUE, showEvent = TRUE,
    showDropout = FALSE, showOngoing = FALSE,
    showsummary = TRUE, showplot = TRUE,
    by_treatment = FALSE, ngroups = 1, alloc = NULL,
    treatment_label = NULL,
    covariates_event = NULL,
    event_prior_with_covariates = NULL,
    covariates_dropout = NULL,
    dropout_prior_with_covariates = NULL,
    fix_parameter = FALSE,
    generate_plot = TRUE,
    interactive_plot = TRUE,
    criterion = "both",
    seed.num = NULL,
    cure_rate = NULL,
    return_subject_data = TRUE,
    return_simulation_data = TRUE
) {
  .old_options <- options(datatable.showProgress = FALSE)
  on.exit(options(.old_options), add = TRUE)
  to_predict_lc <- tolower(to_predict)
  dt <- NULL

  # ---- Input validation ----
  if (!is.null(df)) erify::check_class(df, "data.frame")

  erify::check_content(to_predict_lc,
                       c("enrollment only", "event only",
                         "enrollment and event"))

  if (is.null(df) && to_predict_lc == "event only") {
    stop(paste(
      "getPrediction() does not support design-stage event-only prediction",
      "when df is NULL. Use predictEvent() directly when subject-level",
      "enrollment simulation input is already available."
    ))
  }

  if (!is.na(target_n)) erify::check_n(target_n)
  if (!is.na(target_d)) erify::check_n(target_d)
  if (is.na(target_n) && is.na(target_d))
    stop("At least one of target_n and target_d must be specified")
  if (!is.na(target_n) && !is.na(target_d) && target_d > target_n)
    stop("target_d cannot exceed target_n")
  erify::check_bool(return_subject_data)
  erify::check_bool(return_simulation_data)


  # check by_treatment, ngroups, and alloc
  erify::check_bool(by_treatment)

  if (is.null(df) && by_treatment) {
    erify::check_n(ngroups)
  }

  if (is.null(df)) by_treatment = TRUE

  if (!is.null(df)) {
    # ---- Preprocess observed data ----
    # UPDATED (merge): normalize column names to lower-case for robustness
    dt <- data.table::setDT(data.table::copy(df))
    data.table::setnames(dt, tolower(names(dt)))
  }


  if (by_treatment) {
    if (!is.null(df)) {
      ngroups = dt[, data.table::uniqueN(get("treatment"))]
    }

    if (is.null(alloc)) {
      alloc = rep(1, ngroups)
    } else {
      if (length(alloc) != ngroups) {
        stop("length of alloc must be equal to the number of treatments")
      }

      if (any(alloc <= 0 | alloc != round(alloc))) {
        stop("elements of alloc must be positive integers")
      }
    }
  } else {
    ngroups = 1
  }

  if (ngroups == 1) {
    by_treatment = FALSE
  }

  # -----------------------------------------------------------------------
  # UPDATED (merge): additional input validation for merged functionality
  # -----------------------------------------------------------------------
  if (!is.null(seed.num)) {
    erify::check_n(seed.num, zero = TRUE)
    seed.num <- as.integer(seed.num)
  }

  erify::check_content(tolower(criterion), c("aic", "bic", "both"))

  # Fixed cure fraction for mixture-cure event models
  if (!is.null(cure_rate)) {
    if (any(!is.finite(cure_rate))) stop("cure_rate must be finite")
    if (any(cure_rate <= 0 | cure_rate >= 1)) stop("cure_rate must be in (0, 1)")
    if (by_treatment) {
      if (length(cure_rate) == 1) {
        cure_rate <- rep(cure_rate, ngroups)
      } else if (length(cure_rate) != ngroups) {
        stop("length of cure_rate must be 1 or ngroups when by_treatment=TRUE")
      }
    } else {
      if (length(cure_rate) != 1) stop("length of cure_rate must be 1 when by_treatment=FALSE")
    }
  }

  if (!is.null(treatment_label) && length(treatment_label) != ngroups) {
    stop(paste("length of treatment_label must be equal to",
               "the number of treatments"))
  }
  
  needs_event_prediction <- grepl("event", to_predict, ignore.case = TRUE)
  needs_enrollment_subjects <- needs_event_prediction ||
    return_simulation_data || return_subject_data
  needs_event_subjects <- return_simulation_data || return_subject_data
  observed <- NULL
  enroll_fit <- NULL
  event_fit <- NULL
  event_fit_w_x <- NULL
  dropout_fit <- NULL
  dropout_fit_w_x <- NULL
  enroll_pred <- NULL
  event_pred <- NULL
  subject_data <- NULL

  erify::check_content(tolower(enroll_model),
                       c("poisson", "time-decay", "b-spline",
                         "piecewise poisson", "piecewise uniform"))  # UPDATED (merge)

  erify::check_n(nknots, zero = TRUE)
  erify::check_n(lags, zero = TRUE)

  if (accrualTime[1] != 0) {
    stop("accrualTime must start with 0")
  }
  if (length(accrualTime) > 1 && any(diff(accrualTime) <= 0)) {
    stop("accrualTime should be increasing")
  }

  # check enrollment model prior
  if (!is.null(enroll_prior)) {
    erify::check_class(enroll_prior, "list")
    erify::check_content(tolower(enroll_prior$model), c(
      "poisson", "time-decay", "piecewise poisson", "piecewise uniform"))  # UPDATED (merge)

    model = tolower(enroll_prior$model)
    p = length(enroll_prior$theta)
    vtheta = enroll_prior$vtheta

    if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                   ncol(vtheta) != p)) ||
        (p == 1 && length(c(vtheta)) != 1)) {
      stop(paste("Dimensions of vtheta must be compatible with the length",
                 "of theta in enroll_prior"))
    }

    if ((model == "poisson" && p != 1) ||
        (model == "time-decay" && p != 2) ||
        (model == "piecewise poisson" &&
         p != length(enroll_prior$accrualTime)) ||
        (model == "piecewise uniform" &&
         p != length(enroll_prior$accrualTime))) {  # UPDATED (merge)
      stop(paste("Length of theta must be compatible with model",
                 "in enroll_prior"))
    }

    if (model == "piecewise poisson" || model == "piecewise uniform") {  # UPDATED (merge)
      if (enroll_prior$accrualTime[1] != 0) {
        stop("accrualTime must start with 0 in enroll_prior")
      }
      if (length(enroll_prior$accrualTime) > 1 &&
          any(diff(enroll_prior$accrualTime) <= 0)) {
        stop("accrualTime should be increasing in enroll_prior")
      }
    }

    if (!is.null(df)) {
      if (tolower(enroll_prior$model) != tolower(enroll_model)) {
        stop("Prior and likelihood must use the same enrollment model")
      }

      if ((tolower(enroll_prior$model) == "piecewise poisson" ||
           tolower(enroll_prior$model) == "piecewise uniform") &&  # UPDATED (merge)
          (length(enroll_prior$accrualTime) < length(accrualTime) ||
           !all.equal(enroll_prior$accrualTime[1:length(accrualTime)],
                      accrualTime))) {
        stop(paste("accrualTime of piecewise Poisson must be a subset of",
                   "that in enroll_prior"))
      }
    }
  }


  erify::check_content(
    tolower(event_model),
    c("exponential", "weibull", "log-logistic", "log-normal",
      "piecewise exponential", "model averaging", "spline", "cox",
      # UPDATED (merge): Chen (2016) mixture-cure models
      "exponential with cured population", "weibull with cured population",
      "log-normal with cured population", "log-logistic with cured population",
      "piecewise exponential with cured population",
      # UPDATED (merge): delayed-treatment variants (design-stage simulation)
      "exponential with cured population and delayed treatment",
      "weibull with cured population and delayed treatment",
      "log-normal with cured population and delayed treatment",
      "log-logistic with cured population and delayed treatment",
      "piecewise exponential with cured population and delayed treatment")
  )


  # UPDATED (merge): delayed-treatment event models are for design-stage simulation.
  if (!is.null(df) && grepl("delayed treatment", event_model, ignore.case = TRUE)) {
    stop("Delayed-treatment event models are only supported when df is NULL (design-stage simulation).")
  }

  if (piecewiseSurvivalTime[1] != 0) {
    stop("piecewiseSurvivalTime must start with 0")
  }
  if (length(piecewiseSurvivalTime) > 1 &&
      any(diff(piecewiseSurvivalTime) <= 0)) {
    stop("piecewiseSurvivalTime should be increasing")
  }

  erify::check_n(k, zero = TRUE)
  erify::check_content(tolower(scale), c("hazard", "odds", "normal"))


  # check event model prior
  if (!is.null(event_prior)) {
    erify::check_class(event_prior, "list")

    if (by_treatment) {
      if (length(event_prior) != ngroups) {
        stop("event_prior must be a list with one element per treatment")
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(event_prior)) {
      event_prior2 <- list()
      event_prior2[[1]] <- event_prior
    } else {
      event_prior2 <- event_prior
    }

    for (j in seq_along(event_prior2)) {
      erify::check_content(tolower(event_prior2[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential",
                             # UPDATED (merge): mixture-cure + delayed variants
                             "exponential with cured population",
                             "weibull with cured population",
                             "log-normal with cured population",
                             "log-logistic with cured population",
                             "piecewise exponential with cured population",
                             "exponential with cured population and delayed treatment",
                             "weibull with cured population and delayed treatment",
                             "log-normal with cured population and delayed treatment",
                             "log-logistic with cured population and delayed treatment",
                             "piecewise exponential with cured population and delayed treatment"))

      model <- tolower(event_prior2[[j]]$model)
      theta <- event_prior2[[j]]$theta
      vtheta <- event_prior2[[j]]$vtheta
      p <- length(theta)

      is_cure <- grepl("with cured population", model, fixed = TRUE)
      is_piecewise <- grepl("piecewise exponential", model, fixed = TRUE)
      is_delayed <- grepl("delayed treatment", model, fixed = TRUE)

      # UPDATED (merge): compute expected parameter dimension using the base model
      # (strip the delayed-treatment suffix when present)
      base_model <- if (is_delayed) {
        sub(" and delayed treatment$", "", model)
      } else {
        model
      }

      # Expected parameter dimension by model (eventPred/EventPredInCure parameterization)
      if (is_piecewise) {
        if (is.null(event_prior2[[j]]$piecewiseSurvivalTime)) {
          stop("piecewiseSurvivalTime must be provided in event_prior for piecewise exponential models.")
        }
        p_expected <- length(event_prior2[[j]]$piecewiseSurvivalTime) + ifelse(is_cure, 1, 0)
      } else if (base_model == "exponential") {
        p_expected <- 1
      } else if (base_model == "weibull" || base_model == "log-normal" || base_model == "log-logistic") {
        p_expected <- 2
      } else if (base_model == "exponential with cured population") {
        p_expected <- 2
      } else if (base_model == "weibull with cured population" ||
                 base_model == "log-normal with cured population" ||
                 base_model == "log-logistic with cured population") {
        p_expected <- 3
      } else {
        stop("Unknown event_prior model: ", model)
      }
      if (is_delayed) {
        p_expected <- p_expected + 2
      }

      # ---------------------------------------------------------------------
      # UPDATED (merge): allow cured-model priors to omit the cure parameter
      # when cure_rate is provided, and pad theta/vtheta accordingly.
      # ---------------------------------------------------------------------
      if (is_cure && !is.null(cure_rate)) {
        cr_j <- if (by_treatment) cure_rate[j] else cure_rate
        if (p == p_expected - 1) {
          # Pad cure parameter as the first element: logit(cure_rate)
          theta <- c(stats::qlogis(cr_j), theta)

          vmat <- if (is.matrix(vtheta)) vtheta else matrix(vtheta, nrow = p, ncol = p)
          vtheta <- rbind(c(0, rep(0, p)), cbind(rep(0, p), vmat))

          event_prior2[[j]]$theta <- theta
          event_prior2[[j]]$vtheta <- vtheta
          p <- length(theta)
        }
      }

      # Dimension check for vtheta (after any padding)
      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p || ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in event_prior"))
      }

      # Final parameter-length check
      if (p != p_expected) {
        stop(paste("Length of theta in event_prior is not compatible with model:",
                   model))
      }

      # If df is provided, ensure piecewiseSurvivalTime of prior is compatible with the analysis model
      if (!is.null(df) &&
          grepl("piecewise exponential", event_model, ignore.case = TRUE) &&
          grepl("piecewise exponential", model, ignore.case = TRUE)) {
        if (!all(event_prior2[[j]]$piecewiseSurvivalTime %in% piecewiseSurvivalTime)) {
          stop(paste("piecewiseSurvivalTime of piecewise exponential model",
                     "must be a subset of that in event_prior"))
        }
      }
    }

    # Ensure downstream design-stage prediction uses the validated/padded
    # prior object (e.g., when cure_rate is supplied and the cure parameter
    # is omitted from event_prior$theta).
    if ("model" %in% names(event_prior)) {
      event_prior <- event_prior2[[1]]
    } else {
      event_prior <- event_prior2
    }
  }


  erify::check_content(tolower(dropout_model),
                       c("none", "exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline", "cox"))

  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0")
  }
  if (length(piecewiseDropoutTime) > 1 &&
      any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }

  erify::check_n(k_dropout, zero = TRUE)
  erify::check_content(tolower(scale_dropout), c("hazard", "odds", "normal"))

  # check dropout model prior
  if (!is.null(dropout_prior)) {
    erify::check_class(dropout_prior, "list")

    if (by_treatment) {
      if (length(dropout_prior) != ngroups) {
        stop("dropout_prior must be a list with one element per treatment")
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(dropout_prior)) {
      dropout_prior2 <- list()
      dropout_prior2[[1]] <- dropout_prior
    } else {
      dropout_prior2 <- dropout_prior
    }

    for (j in 1:length(dropout_prior2)) {
      erify::check_content(tolower(dropout_prior2[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(dropout_prior2[[j]]$model)
      p = length(dropout_prior2[[j]]$theta)
      vtheta = dropout_prior2[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in dropout_prior"))
      }

      if ((model == "exponential" && p != 1) ||
          (model == "weibull" && p != 2) ||
          (model == "log-logistic" && p != 2) ||
          (model == "log-normal" && p != 2) ||
          (model == "piecewise exponential" &&
           p != length(dropout_prior2[[j]]$piecewiseDropoutTime))) {
        stop(paste("Length of theta must be compatible with model",
                   "in dropout_prior"))
      }

      if (model == "piecewise exponential") {
        if (dropout_prior2[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_prior"))
        }
        if (length(dropout_prior2[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_prior2[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_prior"))
        }
      }


      if (!is.null(df)) {
        if (tolower(dropout_prior2[[j]]$model) != tolower(dropout_model)) {
          stop("Prior and likelihood must use the same dropout model")
        }

        if (tolower(dropout_prior2[[j]]$model) == "piecewise exponential" &&
            (length(dropout_prior2[[j]]$piecewiseDropoutTime) <
             length(piecewiseDropoutTime) ||
             !all.equal(dropout_prior2[[j]]$piecewiseDropoutTime[
               1:length(piecewiseDropoutTime)], piecewiseDropoutTime))) {
          stop(paste("piecewiseDropoutTime of piecewise exponential model",
                     "must be a subset of that in dropout_prior"))
        }
      }


      if (!is.null(event_prior) && "w" %in% names(event_prior2[[j]])) {
        if (!("w" %in% names(dropout_prior2[[j]])) ||
            event_prior2[[j]]$w != dropout_prior2[[j]]$w) {
          stop("w must be equal between event_prior and dropout_prior")
        }
      }
    }
  }


  erify::check_n(m)
  erify::check_n(m_dropout)
  erify::check_bool(fixedFollowup)
  erify::check_positive(followupTime)
  erify::check_positive(pilevel)
  erify::check_positive(1-pilevel)
  erify::check_positive(nyears)
  erify::check_n(nreps)
  erify::check_bool(showEnrollment)
  erify::check_bool(showEvent)
  erify::check_bool(showDropout)
  erify::check_bool(showOngoing)
  erify::check_bool(showsummary)
  erify::check_bool(showplot)
  erify::check_bool(fix_parameter)

  if (!all(is.na(target_t))) {
    target_t <- target_t[!is.na(target_t)]
    if (any(target_t <= 0 | target_t > nyears*365)) {
      stop("target_t must be positive and less than nyears*365")
    }
  }


  if (!is.null(covariates_event)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_event")
    }

    if (!all(covariates_event %in% colnames(dt))) {
      stop("All covariates_event must exist in df")
    }

    xnames = paste(covariates_event, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_event = model.matrix(formula, dt)  # design matrix with intercept
    q_event = ncol(x_event) - 1
  }


  # check event_prior_with_covariates
  event_prior_w_x <- event_prior_with_covariates
  if (!is.null(event_prior_w_x)) {
    # UPDATED (merge): covariate priors not supported for cured/delayed event models
    if (grepl("with cured population", event_model, ignore.case = TRUE) ||
        grepl("delayed treatment", event_model, ignore.case = TRUE)) {
      stop("event_prior_with_covariates is not supported for cured/delayed event models.")
    }

    if (is.null(covariates_event)) {
      stop(paste("covariates_event must be provided",
                 "for event_prior_with_covariates"))
    }

    erify::check_class(event_prior_w_x, "list")

    if (by_treatment) {
      if (length(event_prior_w_x) != ngroups) {
        stop(paste("event_prior_with_covariates must be a list with",
                   "one element per treatment"))
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(event_prior_w_x)) {
      event_prior2_w_x <- list()
      event_prior2_w_x[[1]] <- event_prior_w_x
    } else {
      event_prior2_w_x <- event_prior_w_x
    }

    for (j in 1:length(event_prior2_w_x)) {
      erify::check_content(tolower(event_prior2_w_x[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(event_prior2_w_x[[j]]$model)
      p = length(event_prior2_w_x[[j]]$theta)
      vtheta = event_prior2_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in event_prior_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_event) ||
          (model == "weibull" && p != 2 + q_event) ||
          (model == "log-logistic" && p != 2 + q_event) ||
          (model == "log-normal" && p != 2 + q_event) ||
          (model == "piecewise exponential" &&
           p != length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) +
           q_event)) {
        stop(paste("Length of theta must be compatible with model",
                   "in event_prior_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (event_prior2_w_x[[j]]$piecewiseSurvivalTime[1] != 0) {
          stop(paste("piecewiseSurvivalTime must start with 0",
                     "in event_prior_with_covariates"))
        }
        if (length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) > 1 &&
            any(diff(event_prior2_w_x[[j]]$piecewiseSurvivalTime) <= 0)) {
          stop(paste("piecewiseSurvivalTime should be increasing",
                     "in event_prior_with_covariates"))
        }
      }


      if (!is.null(df)) {
        if (tolower(event_prior2_w_x[[j]]$model) != tolower(event_model)) {
          stop("Prior and likelihood must use the same event model")
        }

        if (tolower(event_prior2_w_x[[j]]$model) == "piecewise exponential" &&
            (length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) <
             length(piecewiseSurvivalTime) ||
             !all.equal(event_prior2_w_x[[j]]$piecewiseSurvivalTime[
               1:length(piecewiseSurvivalTime)], piecewiseSurvivalTime))) {
          stop(paste("piecewiseSurvivalTime of piecewise exponential model",
                     "must be a subset of that in",
                     "event_prior_with_covariates"))
        }
      }
    }
  }


  if (!is.null(covariates_dropout)) {
    if (is.null(df)) {
      stop("df must be provided in the presence of covariates_dropout")
    }

    if (!all(covariates_dropout %in% colnames(dt))) {
      stop("All covariates_dropout must exist in df")
    }

    xnames = paste(covariates_dropout, collapse = "+")
    formula = as.formula(paste("~", xnames))
    x_dropout = model.matrix(formula, dt)  # design matrix with intercept
    q_dropout = ncol(x_dropout) - 1
  }


  # check event_prior_with_covariates
  dropout_prior_w_x <- dropout_prior_with_covariates
  if (!is.null(dropout_prior_w_x)) {
    if (is.null(covariates_dropout)) {
      stop(paste("covariates_dropout must be provided for",
                 "dropout_prior_with_covariates"))
    }

    erify::check_class(dropout_prior_w_x, "list")

    if (by_treatment) {
      if (length(dropout_prior_w_x) != ngroups) {
        stop(paste("dropout_prior_with_covariates must be a list with",
                   "one element per treatment"))
      }
    }

    # convert to a list with one element per treatment
    if ("model" %in% names(dropout_prior_w_x)) {
      dropout_prior2_w_x <- list()
      dropout_prior2_w_x[[1]] <- dropout_prior_w_x
    } else {
      dropout_prior2_w_x <- dropout_prior_w_x
    }

    for (j in 1:length(dropout_prior2_w_x)) {
      erify::check_content(tolower(dropout_prior2_w_x[[j]]$model),
                           c("exponential", "weibull", "log-logistic",
                             "log-normal", "piecewise exponential"))

      model = tolower(dropout_prior2_w_x[[j]]$model)
      p = length(dropout_prior2_w_x[[j]]$theta)
      vtheta = dropout_prior2_w_x[[j]]$vtheta

      if ((p > 1 && (!is.matrix(vtheta) || nrow(vtheta) != p ||
                     ncol(vtheta) != p)) ||
          (p == 1 && length(c(vtheta)) != 1)) {
        stop(paste("Dimensions of vtheta must be compatible with",
                   "the length of theta in dropout_prior_with_covariates"))
      }

      if ((model == "exponential" && p != 1 + q_dropout) ||
          (model == "weibull" && p != 2 + q_dropout) ||
          (model == "log-logistic" && p != 2 + q_dropout) ||
          (model == "log-normal" && p != 2 + q_dropout) ||
          (model == "piecewise exponential" &&
           p != length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) +
           q_dropout)) {
        stop(paste("Length of theta must be compatible with model",
                   "in dropout_prior_with_covariates"))
      }

      if (model == "piecewise exponential") {
        if (dropout_prior2_w_x[[j]]$piecewiseDropoutTime[1] != 0) {
          stop(paste("piecewiseDropoutTime must start with 0",
                     "in dropout_prior_with_covariates"))
        }
        if (length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) > 1 &&
            any(diff(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) <= 0)) {
          stop(paste("piecewiseDropoutTime should be increasing",
                     "in dropout_prior_with_covariates"))
        }
      }


      if (!is.null(df)) {
        if (tolower(dropout_prior2_w_x[[j]]$model) !=
            tolower(dropout_model)) {
          stop("Prior and likelihood must use the same dropout model")
        }

        if (tolower(dropout_prior2_w_x[[j]]$model) == "piecewise exponential"
            && (length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) <
                length(piecewiseDropoutTime) ||
                !all.equal(dropout_prior2_w_x[[j]]$piecewiseDropoutTime[
                  1:length(piecewiseDropoutTime)], piecewiseDropoutTime))) {
          stop(paste("piecewiseDropoutTime of piecewise exponential model",
                     "must be a subset of that in",
                     "dropout_prior_with_covariates"))
        }
      }
    }
  }


  # check input data set to ensure it has all the required columns
  if (!is.null(df)) {
    cols = colnames(dt)

    if (tolower(to_predict) == "enrollment only") {
      req_cols = c("trialsdt", "usubjid", "randdt", "cutoffdt")
    } else {
      req_cols = c("trialsdt", "usubjid", "randdt", "time", "event",
                   "dropout", "cutoffdt")
    }

    if (by_treatment) {
      req_cols <- c(req_cols, "treatment")
    }

    if (!all(req_cols %in% cols)) {
      stop(paste("The following columns are missing from df:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
    }

    if (any(is.na(dt[, mget(req_cols)]))) {
      stop(paste("The following columns of df have missing values:",
                 paste(req_cols[sapply(dt, function(x) any(is.na(x)))],
                       collapse = ", ")))
    }

    if ("treatment" %in% cols && !("treatment_description" %in% cols)) {
      dt[, `:=`(treatment_description = paste("Treatment", get("treatment")))]
    }
  }


  if (!is.null(df)) {
    dt$trialsdt <- as.Date(dt$trialsdt)
    dt$randdt <- as.Date(dt$randdt)
    dt$cutoffdt <- as.Date(dt$cutoffdt)

    trialsdt = dt[1, get("trialsdt")]
    cutoffdt = dt[1, get("cutoffdt")]

    # summarize observed data (counts + KM)
    observed <- summarizeObserved(dt, to_predict, showplot,
                                  by_treatment,
                                  generate_plot, interactive_plot)
  }

  if (!is.null(covariates_event)) {
    covariates_event <- tolower(covariates_event)
  }

  # UPDATED (merge): cure/delayed event models currently do not support covariates
  if (!is.null(covariates_event) &&
      grepl("with cured population", event_model, ignore.case = TRUE)) {
    stop("covariates_event is not supported for cured event models.")
  }
  if (!is.null(covariates_event) &&
      grepl("delayed treatment", event_model, ignore.case = TRUE)) {
    stop("covariates_event is not supported for delayed-treatment event models.")
  }


  if (!is.null(covariates_dropout)) {
    covariates_dropout <- tolower(covariates_dropout)
  }


  # ---------------------------------------------------------------------
  # Enrollment model: fit/prior update + prediction
  # ---------------------------------------------------------------------
  # ---------------------------------------------------------------------
  # Enrollment prediction
  #   - Fit enrollment model (if df provided)
  #   - Combine prior + likelihood (optional)
  #   - Predict future enrollment via predictEnrollment()
  # ---------------------------------------------------------------------
  # fit and predict enrollment
  if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) {
      erify::check_n(target_n - observed$n0,
                     supplement = "Enrollment target reached.")

      enroll_fit <- fitEnrollment(
        df = dt,
        enroll_model = enroll_model,
        nknots = nknots,
        accrualTime = accrualTime,
        showplot = showplot,
        generate_plot = generate_plot,
        interactive_plot = interactive_plot,
        criterion = criterion  # UPDATED (merge)
      )
      enroll_fit1 <- enroll_fit$fit

      # combine prior and likelihood to yield posterior
      if (!is.null(enroll_prior)) {
        if (tolower(enroll_model) == "piecewise poisson" &&
            length(enroll_prior$accrualTime) > length(accrualTime)) {
          l1 <- length(accrualTime)

          # information from prior
          info = solve(enroll_prior$vtheta)
          # add information from data
          info[1:l1,1:l1] = info[1:l1,1:l1] +
            solve(enroll_fit1$vtheta)

          mu = solve(enroll_prior$vtheta,
                     enroll_prior$theta)
          mu[1:l1] = mu[1:l1] + solve(enroll_fit1$vtheta,
                                      enroll_fit1$theta)
        } else {
          info = solve(enroll_prior$vtheta) +
            solve(enroll_fit1$vtheta)
          mu = solve(enroll_prior$vtheta,
                     enroll_prior$theta) +
            solve(enroll_fit1$vtheta,
                  enroll_fit1$theta)
        }

        enroll_fit1 <- enroll_prior
        enroll_fit1$vtheta = solve(info)
        enroll_fit1$theta = as.numeric(enroll_fit1$vtheta %*% mu)
      }

      # enrollment prediction at the analysis stage
      enroll_pred <- predictEnrollment(
        df = dt,
        target_n = target_n,
        enroll_fit = enroll_fit1,
        lags = lags,
        pilevel = pilevel,
        nyears = nyears,
        nreps = nreps,
        showsummary = showsummary,
        showplot = FALSE,
        by_treatment = by_treatment,
        ngroups = ngroups,
        alloc = alloc,
        treatment_label = treatment_label,
        fix_parameter = fix_parameter,
        generate_plot = generate_plot,
        interactive_plot = interactive_plot,
        seed.num = seed.num,  # UPDATED (merge)
        return_new_subjects = needs_enrollment_subjects
      )
    } else {
      # enrollment prediction at the design stage
      enroll_pred <- predictEnrollment(
        df = NULL,
        target_n = target_n,
        enroll_fit = enroll_prior,
        lags = lags,
        pilevel = pilevel,
        nyears = nyears,
        nreps = nreps,
        showsummary = showsummary,
        showplot = FALSE,
        by_treatment = by_treatment,
        ngroups = ngroups,
        alloc = alloc,
        treatment_label = treatment_label,
        fix_parameter = fix_parameter,
        generate_plot = generate_plot,
        interactive_plot = interactive_plot,
        seed.num = seed.num,  # UPDATED (merge)
        return_new_subjects = needs_enrollment_subjects
      )
    }
  }


  # ---------------------------------------------------------------------
  # Event model fitting + event prediction
  #   - Optionally fit event model (if df provided)
  #   - Optionally combine prior + likelihood
  #   - Predict future events via predictEvent()
  # ---------------------------------------------------------------------
  # fit and predict event
  if (grepl("event", to_predict, ignore.case = TRUE)) {
    if (!is.null(df)) { # event prediction at analysis stage
      erify::check_n(target_d - observed$d0,
                     supplement = "Event target reached.")

      if (by_treatment) {
        sum_by_trt <- dt[, list(
          n0 = .N, d0 = sum(get("event")), c0 = sum(get("dropout")),
          r0 = sum(!(get("event") | get("dropout")))),
          by = "treatment"]
      }


      # -----------------------------------------------------------------------
      # UPDATED (merge): convert prior-by-treatment to an overall prior (blinded)
      # Using moment-matching pooling adapted from EventPredInCure.
      # NOTE: For Chen-style mixture-cure models, pooling is not supported; users
      # should provide an overall prior directly when by_treatment=FALSE.
      # -----------------------------------------------------------------------
      if (!is.null(event_prior) && !by_treatment && !"model" %in% names(event_prior)) {
        if (length(event_prior) != ngroups) {
          stop("event_prior must have one element per treatment.",
               call. = FALSE)
        }
        w <- vapply(event_prior, function(x) x$w, numeric(1))
        if (any(w <= 0) || any(!is.finite(w))) stop("w must be positive and finite in event_prior")
        w <- w / sum(w)

        model0 <- tolower(event_prior[[1]]$model)
        if (any(vapply(event_prior, function(x) tolower(x$model) != model0, logical(1)))) {
          stop("All elements of event_prior must have the same model when pooling.")
        }

        if (grepl("with cured population", model0, fixed = TRUE)) {
          stop("Pooling event_prior across treatments is not supported for cured models. ",
               "Please provide an overall event_prior when by_treatment=FALSE.")
        }

        # Helper: mean/var for parametric distributions under eventPred parameterization
        .weibull_mean_var <- function(theta) {
          sc <- exp(theta[1]); sh <- exp(-theta[2])
          mu <- sc * gamma(1 + 1/sh)
          va <- sc^2 * (gamma(1 + 2/sh) - gamma(1 + 1/sh)^2)
          c(mu = mu, var = va)
        }
        .lnorm_mean_var <- function(theta) {
          mu <- exp(theta[1] + exp(theta[2])^2 / 2)
          va <- (exp(exp(theta[2])^2) - 1) * exp(2 * theta[1] + exp(theta[2])^2)
          c(mu = mu, var = va)
        }

        if (model0 == "exponential") {
          lam <- exp(vapply(event_prior, function(x) x$theta, numeric(1)))
          lam1 <- 1 / sum(w / lam)
          vtheta <- vapply(event_prior, function(x) x$vtheta, numeric(1))
          vtheta1 <- sum((w/lam)^2 * vtheta) * lam1^2
          event_prior1 <- list(model = "exponential", theta = log(lam1), vtheta = vtheta1)

        } else if (model0 == "weibull") {
          mm <- t(vapply(event_prior, function(x) .weibull_mean_var(x$theta), numeric(2)))
          m1 <- mm[, "mu"]; v1 <- mm[, "var"]
          mu1 <- sum(w * m1)
          sig2 <- sum(w * v1) + sum(w * (m1 - mu1)^2)

          froot <- function(th2) {
            sh <- exp(-th2)
            mu_w <- gamma(1 + 1/sh)
            va_w <- gamma(1 + 2/sh) - mu_w^2
            va_w - sig2 / mu1^2 * mu_w^2
          }
          th2 <- stats::uniroot(froot, interval = c(-10, 10))$root
          sh <- exp(-th2)
          th1 <- log(mu1 / gamma(1 + 1/sh))
          theta1 <- c(th1, th2)

          gmean <- function(th) {
            sc <- exp(th[1]); sh <- exp(-th[2])
            mu <- sc * gamma(1 + 1/sh)
            va <- sc^2 * (gamma(1 + 2/sh) - gamma(1 + 1/sh)^2)
            c(mu = mu, var = va)
          }
          gm <- function(th) numDeriv::jacobian(gmean, th)

          vtheta1 <- matrix(0, 2, 2)
          for (i in seq_len(ngroups)) {
            gi <- gm(event_prior[[i]]$theta)
            vmi <- gi %*% event_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge): matrix mult
            vtheta1 <- vtheta1 + w[i]^2 * vmi
          }
          ginv <- solve(gm(theta1))
          vtheta1 <- ginv %*% vtheta1 %*% t(ginv)
          event_prior1 <- list(model = "weibull", theta = theta1, vtheta = vtheta1)

        } else if (model0 == "log-normal") {
          mm <- t(vapply(event_prior, function(x) .lnorm_mean_var(x$theta), numeric(2)))
          m1 <- mm[, "mu"]; v1 <- mm[, "var"]
          mu1 <- sum(w * m1)
          sig2 <- sum(w * v1) + sum(w * (m1 - mu1)^2)

          sig2log <- log(sig2 / mu1^2 + 1)
          theta1 <- c(log(mu1) - sig2log/2, log(sqrt(sig2log)))

          gm <- function(th) numDeriv::jacobian(.lnorm_mean_var, th)

          vtheta1 <- matrix(0, 2, 2)
          for (i in seq_len(ngroups)) {
            gi <- gm(event_prior[[i]]$theta)
            vmi <- gi %*% event_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge)
            vtheta1 <- vtheta1 + w[i]^2 * vmi
          }
          ginv <- solve(gm(theta1))
          vtheta1 <- ginv %*% vtheta1 %*% t(ginv)
          event_prior1 <- list(model = "log-normal", theta = theta1, vtheta = vtheta1)

        } else if (model0 == "log-logistic") {
          # Quantile-matching approach adapted from EventPredInCure (with bug fixes).
          theta_all <- unlist(lapply(event_prior, function(x) x$theta))  # FIXED (merge)
          k2 <- length(theta_all) / 2
          # Map (theta1,theta2) -> q-th quantile on original time scale
          fq <- function(th, q) exp(th[1] + qlogis(q)/exp(th[2]))

          qbar <- function(th_all) {
            mat <- matrix(th_all, nrow = 2, ncol = k2)
            q25 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.25), numeric(1))
            q50 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.50), numeric(1))
            q75 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.75), numeric(1))
            c(sum(w * log(q25)), sum(w * log(q50)), sum(w * log(q75)))
          }
          qlog <- qbar(theta_all)

          obj <- function(th) {
            q <- vapply(c(0.25, 0.50, 0.75), function(p) fq(th, p), numeric(1))
            sum((log(q) - qlog)^2)
          }
          theta1 <- stats::optim(c(qlog[2], 0), obj)$par

          # Approximate vtheta via numerical Jacobian of pooled estimate wrt per-arm thetas
          g <- numDeriv::jacobian(function(th_all) {
            qlog2 <- qbar(th_all)
            obj2 <- function(th) {
              q <- vapply(c(0.25,0.50,0.75), function(p) fq(th, p), numeric(1))
              sum((log(q) - qlog2)^2)
            }
            stats::optim(c(qlog2[2], 0), obj2)$par
          }, theta_all)

          vtheta1 <- matrix(0, 2, 2)
          for (i in seq_len(ngroups)) {
            gi <- g[, (2*i-1):(2*i), drop = FALSE]
            vtheta1 <- vtheta1 + gi %*% event_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge)
          }
          event_prior1 <- list(model = "log-logistic", theta = theta1, vtheta = vtheta1)

        } else if (model0 == "piecewise exponential") {
          np <- length(event_prior[[1]]$theta)
          lam <- matrix(NA_real_, nrow = ngroups, ncol = np)
          for (i in seq_len(ngroups)) lam[i, ] <- exp(event_prior[[i]]$theta)
          lam1 <- 1 / colSums((w / lam))

          vtheta1 <- matrix(0, np, np)
          for (j in seq_len(np)) {
            vj <- vapply(event_prior, function(x) x$vtheta[j, j], numeric(1))
            vtheta1[j, j] <- sum((w/lam[, j])^2 * vj) * lam1[j]^2
          }
          event_prior1 <- list(model = "piecewise exponential",
                               piecewiseSurvivalTime = event_prior[[1]]$piecewiseSurvivalTime,
                               theta = log(lam1), vtheta = vtheta1)
        } else {
          stop("Pooling for model ", model0, " is not implemented.")
        }

      } else if (!is.null(event_prior)) {
        event_prior1 <- event_prior
      }

      if (!is.null(event_prior_w_x) && !by_treatment &&
          !("model" %in% names(event_prior_w_x))) {

        m0 <- length(event_prior_w_x)
        w <- sapply(event_prior_w_x, function(sub_list) sub_list$w)
        if (any(w <= 0)) {
          stop("w must be positive in event_prior_with_covariates")
        }
        w <- w/sum(w)

        # check model consistency across treatments
        model <- tolower(event_prior_w_x[[1]]$model)
        if (m0 > 1) {
          for (j in 2:m0) {
            if (tolower(event_prior_w_x[[j]]$model) != model) {
              stop(paste("Event model must be equal across treatments in",
                         "event_prior_with_covariates"))
            }
          }
        }
        # check piecewiseSurvivalTime consistency (if applicable)
        if (model == "piecewise exponential") {
          if (m0 > 1) {
            for (j in 2:m0) {
              if (length(event_prior_w_x[[j]]$piecewiseSurvivalTime) !=
                  length(event_prior_w_x[[1]]$piecewiseSurvivalTime)) {
                stop(paste("Length of piecewiseSurvivalTime must be equal across",
                           "treatments in event_prior_with_covariates"))
              }
              if (any(event_prior_w_x[[j]]$piecewiseSurvivalTime !=
                      event_prior_w_x[[1]]$piecewiseSurvivalTime)) {
                stop(paste("piecewiseSurvivalTime must be equal across treatments in",
                           "event_prior_with_covariates"))
              }
            }
          }
        }

        # UPDATED (merge): pooled prior via moment matching helper
        pooled <- .ep_pool_priors_moment_match(event_prior_w_x, w)
        theta <- pooled$theta
        vtheta <- pooled$vtheta

        if (model %in% c("exponential", "weibull", "log-logistic",
                         "log-normal")) {
          event_prior1_w_x <- list(
            model = model, theta = theta, vtheta = vtheta)
        } else if (model == "piecewise exponential") {
          event_prior1_w_x <- list(
            model = model, theta = theta, vtheta = vtheta,
            piecewiseSurvivalTime =
              event_prior_w_x[[1]]$piecewiseSurvivalTime)
        }
      } else if (!is.null(event_prior_w_x)) {
        event_prior1_w_x <- event_prior_w_x
      }


      # ----- Event model fitting (optional) -----


      # Fit event model and optionally combine with prior


      #


      # fit the event model without covariates
      if (!(to_predict == "event only" && !is.null(covariates_event))) {
        event_fit <- fitEvent(
          df = dt,
          event_model = event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime,
          k = k,
          scale = scale,
          m = m,
          showplot = showplot,
          by_treatment = by_treatment,
          covariates = NULL,
          generate_plot = generate_plot,
          interactive_plot = interactive_plot,
          criterion = criterion,  # UPDATED (merge)
          cure_rate = cure_rate   # UPDATED (merge)
        )
        if (!by_treatment) {
          event_fit1 <- list(event_fit$fit)
        } else {
          event_fit1 <- purrr::map(event_fit, "fit")  # UPDATED (merge)
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(event_prior)) {
          if (!by_treatment) {
            event_prior2 <- list()
            event_prior2[[1]] <- event_prior1
          } else {
            event_prior2 <- event_prior1
          }

          for (j in 1:ngroups) {
            if (tolower(event_model) == "piecewise exponential" &&
                length(event_prior2[[j]]$piecewiseSurvivalTime) >
                length(piecewiseSurvivalTime)) {
              l1 <- length(piecewiseSurvivalTime)

              # information from prior
              info = solve(event_prior2[[j]]$vtheta)
              # add information from data
              info[1:l1,1:l1] = info[1:l1,1:l1] +
                solve(event_fit1[[j]]$vtheta)

              mu = solve(event_prior2[[j]]$vtheta,
                         event_prior2[[j]]$theta)
              mu[1:l1] = mu[1:l1] + solve(event_fit1[[j]]$vtheta,
                                          event_fit1[[j]]$theta)
            } else {
              info = solve(event_prior2[[j]]$vtheta) +
                solve(event_fit1[[j]]$vtheta)
              mu = solve(event_prior2[[j]]$vtheta,
                         event_prior2[[j]]$theta) +
                solve(event_fit1[[j]]$vtheta,
                      event_fit1[[j]]$theta)
            }

            event_fit2 <- event_prior2[[j]]
            event_fit2$vtheta = solve(info)
            event_fit2$theta = as.numeric(event_fit2$vtheta %*% mu)

            event_fit1[[j]] <- event_fit2
          }
        }

        if (!by_treatment) {
          event_fit1 <- event_fit1[[1]]
        }
      } else {
        event_fit1 <- NULL
      }


      # fit the event model with covariates
      if (!is.null(covariates_event)) {
        event_fit_w_x <- fitEvent(
          df = dt,
          event_model = event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime,
          k = k,
          scale = scale,
          m = m,
          showplot = showplot,
          by_treatment = by_treatment,
          covariates = covariates_event,
          generate_plot = generate_plot,
          interactive_plot = interactive_plot,
          criterion = criterion,  # UPDATED (merge)
          cure_rate = cure_rate   # UPDATED (merge)
        )
        if (!by_treatment) {
          event_fit1_w_x <- list()
          event_fit1_w_x[[1]] <- event_fit_w_x$fit
        } else {
          event_fit1_w_x <- purrr::map(event_fit_w_x,
                                       function(fit) fit$fit)
        }

        # combine prior and likelihood to yield posterior
        if (!is.null(event_prior_w_x)) {
          if (!by_treatment) {
            event_prior2_w_x <- list()
            event_prior2_w_x[[1]] <- event_prior1_w_x
          } else {
            event_prior2_w_x <- event_prior1_w_x
          }

          for (j in 1:ngroups) {
            if (tolower(event_model) == "piecewise exponential" &&
                length(event_prior2_w_x[[j]]$piecewiseSurvivalTime) >
                length(piecewiseSurvivalTime)) {
              l1 <- length(piecewiseSurvivalTime)

              # information from prior
              info = solve(event_prior2_w_x[[j]]$vtheta)
              # add information from data
              info[1:l1,1:l1] = info[1:l1,1:l1] +
                solve(event_fit1_w_x[[j]]$vtheta)

              mu = solve(event_prior2_w_x[[j]]$vtheta,
                         event_prior2_w_x[[j]]$theta)
              mu[1:l1] = mu[1:l1] + solve(event_fit1_w_x[[j]]$vtheta,
                                          event_fit1_w_x[[j]]$theta)
            } else {
              info = solve(event_prior2_w_x[[j]]$vtheta) +
                solve(event_fit1_w_x[[j]]$vtheta)
              mu = solve(event_prior2_w_x[[j]]$vtheta,
                         event_prior2_w_x[[j]]$theta) +
                solve(event_fit1_w_x[[j]]$vtheta,
                      event_fit1_w_x[[j]]$theta)
            }

            event_fit2_w_x <- event_prior2_w_x[[j]]
            event_fit2_w_x$vtheta = solve(info)
            event_fit2_w_x$theta = as.numeric(event_fit2_w_x$vtheta
                                              %*% mu)

            event_fit1_w_x[[j]] <- event_fit2_w_x
          }
        }

        if (!by_treatment) {
          event_fit1_w_x <- event_fit1_w_x[[1]]
        }
      } else {
        event_fit1_w_x <- NULL
      }



      # ---------------------------------------------------------------------
      # Dropout model: fit/prior update + prediction (optional)
      # ---------------------------------------------------------------------
      # whether to include dropout model
      if (tolower(dropout_model) != "none") {


        # ---------------------------------------------------------------------
        # UPDATED (merge): convert dropout prior-by-treatment to an overall prior
        # Using moment-matching pooling adapted from EventPredInCure (bug-fixed).
        # ---------------------------------------------------------------------
        if (!is.null(dropout_prior) && !by_treatment &&
            !("model" %in% names(dropout_prior))) {
          if (length(dropout_prior) != ngroups) {
            stop("dropout_prior must have one element per treatment.",
                 call. = FALSE)
          }
          w <- vapply(dropout_prior, function(x) x$w, numeric(1))
          if (any(w <= 0) || any(!is.finite(w))) stop("w must be positive and finite in dropout_prior")
          w <- w / sum(w)

          model0 <- tolower(dropout_prior[[1]]$model)
          if (any(vapply(dropout_prior, function(x) tolower(x$model) != model0, logical(1)))) {
            stop("All elements of dropout_prior must have the same model when pooling.")
          }

          .weibull_mean_var <- function(theta) {
            sc <- exp(theta[1]); sh <- exp(-theta[2])
            mu <- sc * gamma(1 + 1/sh)
            va <- sc^2 * (gamma(1 + 2/sh) - gamma(1 + 1/sh)^2)
            c(mu = mu, var = va)
          }
          .lnorm_mean_var <- function(theta) {
            mu <- exp(theta[1] + exp(theta[2])^2 / 2)
            va <- (exp(exp(theta[2])^2) - 1) * exp(2 * theta[1] + exp(theta[2])^2)
            c(mu = mu, var = va)
          }

          if (model0 == "exponential") {
            lam <- exp(vapply(dropout_prior, function(x) x$theta, numeric(1)))
            lam1 <- 1 / sum(w / lam)
            vtheta <- vapply(dropout_prior, function(x) x$vtheta, numeric(1))
            vtheta1 <- sum((w/lam)^2 * vtheta) * lam1^2
            dropout_prior1 <- list(model = "exponential", theta = log(lam1), vtheta = vtheta1)

          } else if (model0 == "weibull") {
            mm <- t(vapply(dropout_prior, function(x) .weibull_mean_var(x$theta), numeric(2)))
            m1 <- mm[, "mu"]; v1 <- mm[, "var"]
            mu1 <- sum(w * m1)
            sig2 <- sum(w * v1) + sum(w * (m1 - mu1)^2)

            froot <- function(th2) {
              sh <- exp(-th2)
              mu_w <- gamma(1 + 1/sh)
              va_w <- gamma(1 + 2/sh) - mu_w^2
              va_w - sig2 / mu1^2 * mu_w^2
            }
            th2 <- stats::uniroot(froot, interval = c(-10, 10))$root
            sh <- exp(-th2)
            th1 <- log(mu1 / gamma(1 + 1/sh))
            theta1 <- c(th1, th2)

            gmean <- function(th) {
              sc <- exp(th[1]); sh <- exp(-th[2])
              mu <- sc * gamma(1 + 1/sh)
              va <- sc^2 * (gamma(1 + 2/sh) - gamma(1 + 1/sh)^2)
              c(mu = mu, var = va)
            }
            gm <- function(th) numDeriv::jacobian(gmean, th)

            vtheta1 <- matrix(0, 2, 2)
            for (i in seq_len(ngroups)) {
              gi <- gm(dropout_prior[[i]]$theta)
              vmi <- gi %*% dropout_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge)
              vtheta1 <- vtheta1 + w[i]^2 * vmi
            }
            ginv <- solve(gm(theta1))
            vtheta1 <- ginv %*% vtheta1 %*% t(ginv)
            dropout_prior1 <- list(model = "weibull", theta = theta1, vtheta = vtheta1)

          } else if (model0 == "log-normal") {
            mm <- t(vapply(dropout_prior, function(x) .lnorm_mean_var(x$theta), numeric(2)))
            m1 <- mm[, "mu"]; v1 <- mm[, "var"]
            mu1 <- sum(w * m1)
            sig2 <- sum(w * v1) + sum(w * (m1 - mu1)^2)

            sig2log <- log(sig2 / mu1^2 + 1)
            theta1 <- c(log(mu1) - sig2log/2, log(sqrt(sig2log)))

            gm <- function(th) numDeriv::jacobian(.lnorm_mean_var, th)

            vtheta1 <- matrix(0, 2, 2)
            for (i in seq_len(ngroups)) {
              gi <- gm(dropout_prior[[i]]$theta)
              vmi <- gi %*% dropout_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge)
              vtheta1 <- vtheta1 + w[i]^2 * vmi
            }
            ginv <- solve(gm(theta1))
            vtheta1 <- ginv %*% vtheta1 %*% t(ginv)
            dropout_prior1 <- list(model = "log-normal", theta = theta1, vtheta = vtheta1)

          } else if (model0 == "log-logistic") {
            theta_all <- unlist(lapply(dropout_prior, function(x) x$theta))  # FIXED (merge)
            k2 <- length(theta_all) / 2
            fq <- function(th, q) exp(th[1] + qlogis(q)/exp(th[2]))

            qbar <- function(th_all) {
              mat <- matrix(th_all, nrow = 2, ncol = k2)
              q25 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.25), numeric(1))
              q50 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.50), numeric(1))
              q75 <- vapply(seq_len(k2), function(i) fq(mat[, i], 0.75), numeric(1))
              c(sum(w * log(q25)), sum(w * log(q50)), sum(w * log(q75)))
            }
            qlog <- qbar(theta_all)

            obj <- function(th) {
              q <- vapply(c(0.25, 0.50, 0.75), function(p) fq(th, p), numeric(1))
              sum((log(q) - qlog)^2)
            }
            theta1 <- stats::optim(c(qlog[2], 0), obj)$par

            g <- numDeriv::jacobian(function(th_all) {
              qlog2 <- qbar(th_all)
              obj2 <- function(th) {
                q <- vapply(c(0.25,0.50,0.75), function(p) fq(th, p), numeric(1))
                sum((log(q) - qlog2)^2)
              }
              stats::optim(c(qlog2[2], 0), obj2)$par
            }, theta_all)

            vtheta1 <- matrix(0, 2, 2)
            for (i in seq_len(ngroups)) {
              gi <- g[, (2*i-1):(2*i), drop = FALSE]
              vtheta1 <- vtheta1 + gi %*% dropout_prior[[i]]$vtheta %*% t(gi)  # FIXED (merge)
            }
            dropout_prior1 <- list(model = "log-logistic", theta = theta1, vtheta = vtheta1)

          } else if (model0 == "piecewise exponential") {
            np <- length(dropout_prior[[1]]$theta)
            lam <- matrix(NA_real_, nrow = ngroups, ncol = np)
            for (i in seq_len(ngroups)) lam[i, ] <- exp(dropout_prior[[i]]$theta)
            lam1 <- 1 / colSums((w / lam))

            vtheta1 <- matrix(0, np, np)
            for (j in seq_len(np)) {
              vj <- vapply(dropout_prior, function(x) x$vtheta[j, j], numeric(1))
              vtheta1[j, j] <- sum((w/lam[, j])^2 * vj) * lam1[j]^2
            }
            dropout_prior1 <- list(model = "piecewise exponential",
                                   piecewiseDropoutTime = dropout_prior[[1]]$piecewiseDropoutTime,
                                   theta = log(lam1), vtheta = vtheta1)
          } else {
            stop("Pooling for model ", model0, " is not implemented.")
          }

        } else if (!is.null(dropout_prior)) {
          dropout_prior1 <- dropout_prior
        }

        if (!is.null(dropout_prior_w_x) && !by_treatment &&
            !("model" %in% names(dropout_prior_w_x))) {

          m0 <- length(dropout_prior_w_x)
          w <- sapply(dropout_prior_w_x, function(sub_list) sub_list$w)
          if (any(w <= 0)) {
            stop("w must be positive in dropout_prior_with_covariates")
          }
          w <- w/sum(w)

          # check model consistency across treatments
          model <- tolower(dropout_prior_w_x[[1]]$model)
          if (m0 > 1) {
            for (j in 2:m0) {
              if (tolower(dropout_prior_w_x[[j]]$model) != model) {
                stop(paste("dropout model must be equal across treatments in",
                           "dropout_prior_with_covariates"))
              }
            }
          }
          # check piecewiseDropoutTime consistency (if applicable)
          if (model == "piecewise exponential") {
            if (m0 > 1) {
              for (j in 2:m0) {
                if (length(dropout_prior_w_x[[j]]$piecewiseDropoutTime) !=
                    length(dropout_prior_w_x[[1]]$piecewiseDropoutTime)) {
                  stop(paste("Length of piecewiseDropoutTime must be equal across",
                             "treatments in dropout_prior_with_covariates"))
                }
                if (any(dropout_prior_w_x[[j]]$piecewiseDropoutTime !=
                        dropout_prior_w_x[[1]]$piecewiseDropoutTime)) {
                  stop(paste("piecewiseDropoutTime must be equal across treatments in",
                             "dropout_prior_with_covariates"))
                }
              }
            }
          }

          # UPDATED (merge): pooled prior via moment matching helper
          pooled <- .ep_pool_priors_moment_match(dropout_prior_w_x, w)
          theta <- pooled$theta
          vtheta <- pooled$vtheta

          if (model %in% c("exponential", "weibull", "log-logistic",
                           "log-normal")) {
            dropout_prior1_w_x <- list(
              model = model, theta = theta, vtheta = vtheta)
          } else if (model == "piecewise exponential") {
            dropout_prior1_w_x <- list(
              model = model, theta = theta, vtheta = vtheta,
              piecewiseDropoutTime =
                dropout_prior_w_x[[1]]$piecewiseDropoutTime)
          }
        } else if (!is.null(dropout_prior_w_x)) {
          dropout_prior1_w_x <- dropout_prior_w_x
        }


        # ----- Dropout model fitting (optional) -----


        # Fit dropout model and optionally combine with prior; used for event projection


        #


        # fit the dropout model without covariates
        if (!(to_predict == "event only" && !is.null(covariates_dropout))) {
          dropout_fit <- fitDropout(
            df = dt,
            dropout_model = dropout_model,
            piecewiseDropoutTime = piecewiseDropoutTime,
            k_dropout = k_dropout,
            scale_dropout = scale_dropout,
            m_dropout = m_dropout,
            showplot = showplot,
            by_treatment = by_treatment,
            covariates = NULL,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            criterion = criterion  # UPDATED (merge)
          )
          if (!by_treatment) {
            dropout_fit1 <- list(dropout_fit$fit)
          } else {
            dropout_fit1 <- purrr::map(dropout_fit, "fit")  # UPDATED (merge)
          }

          # combine prior and likelihood to yield posterior
          if (!is.null(dropout_prior)) {
            if (!by_treatment) {
              dropout_prior2 <- list()
              dropout_prior2[[1]] <- dropout_prior1
            } else {
              dropout_prior2 <- dropout_prior1
            }

            for (j in 1:ngroups) {
              if (tolower(dropout_model) == "piecewise exponential" &&
                  length(dropout_prior2[[j]]$piecewiseDropoutTime) >
                  length(piecewiseDropoutTime)) {
                l1 <- length(piecewiseDropoutTime)

                # information from prior
                info = solve(dropout_prior2[[j]]$vtheta)
                # add information from data
                info[1:l1,1:l1] = info[1:l1,1:l1] +
                  solve(dropout_fit1[[j]]$vtheta)

                mu = solve(dropout_prior2[[j]]$vtheta,
                           dropout_prior2[[j]]$theta)
                mu[1:l1] = mu[1:l1] + solve(dropout_fit1[[j]]$vtheta,
                                            dropout_fit1[[j]]$theta)
              } else {
                info = solve(dropout_prior2[[j]]$vtheta) +
                  solve(dropout_fit1[[j]]$vtheta)
                mu = solve(dropout_prior2[[j]]$vtheta,
                           dropout_prior2[[j]]$theta) +
                  solve(dropout_fit1[[j]]$vtheta,
                        dropout_fit1[[j]]$theta)
              }

              dropout_fit2 <- dropout_prior2[[j]]
              dropout_fit2$vtheta = solve(info)
              dropout_fit2$theta = as.numeric(dropout_fit2$vtheta %*% mu)

              dropout_fit1[[j]] <- dropout_fit2
            }
          }

          if (!by_treatment) {
            dropout_fit1 <- dropout_fit1[[1]]
          }
        } else {
          dropout_fit1 <- NULL
        }


        # fit the dropout model with covariates
        if (!is.null(covariates_dropout)) {
          dropout_fit_w_x <- fitDropout(
            df = dt,
            dropout_model = dropout_model,
            piecewiseDropoutTime = piecewiseDropoutTime,
            k_dropout = k_dropout,
            scale_dropout = scale_dropout,
            m_dropout = m_dropout,
            showplot = showplot,
            by_treatment = by_treatment,
            covariates = covariates_dropout,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            criterion = criterion  # UPDATED (merge)
          )
          if (!by_treatment) {
            dropout_fit1_w_x <- list()
            dropout_fit1_w_x[[1]] <- dropout_fit_w_x$fit
          } else {
            dropout_fit1_w_x <- purrr::map(dropout_fit_w_x,
                                           function(fit) fit$fit)
          }

          # combine prior and likelihood to yield posterior
          if (!is.null(dropout_prior_w_x)) {
            if (!by_treatment) {
              dropout_prior2_w_x <- list()
              dropout_prior2_w_x[[1]] <- dropout_prior1_w_x
            } else {
              dropout_prior2_w_x = dropout_prior1_w_x
            }

            for (j in 1:ngroups) {
              if (tolower(dropout_model) == "piecewise exponential" &&
                  length(dropout_prior2_w_x[[j]]$piecewiseDropoutTime) >
                  length(piecewiseDropoutTime)) {
                l1 <- length(piecewiseDropoutTime)

                # information from prior
                info = solve(dropout_prior2_w_x[[j]]$vtheta)
                # add information from data
                info[1:l1,1:l1] = info[1:l1,1:l1] +
                  solve(dropout_fit1_w_x[[j]]$vtheta)

                mu = solve(dropout_prior2_w_x[[j]]$vtheta,
                           dropout_prior2_w_x[[j]]$theta)
                mu[1:l1] = mu[1:l1] + solve(dropout_fit1_w_x[[j]]$vtheta,
                                            dropout_fit1_w_x[[j]]$theta)
              } else {
                info = solve(dropout_prior2_w_x[[j]]$vtheta) +
                  solve(dropout_fit1_w_x[[j]]$vtheta)
                mu = solve(dropout_prior2_w_x[[j]]$vtheta,
                           dropout_prior2_w_x[[j]]$theta) +
                  solve(dropout_fit1_w_x[[j]]$vtheta,
                        dropout_fit1_w_x[[j]]$theta)
              }

              dropout_fit2_w_x <- dropout_prior2_w_x[[j]]
              dropout_fit2_w_x$vtheta = solve(info)
              dropout_fit2_w_x$theta = as.numeric(dropout_fit2_w_x$vtheta
                                                  %*% mu)

              dropout_fit1_w_x[[j]] <- dropout_fit2_w_x
            }
          }

          if (!by_treatment) {
            dropout_fit1_w_x <- dropout_fit1_w_x[[1]]
          }
        } else {
          dropout_fit1_w_x <- NULL
        }

        # ----- Event projection -----

        # Project time to target #events (or events by time) using fitted models

        #

        # event prediction with a dropout model
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = dt,
            target_d = target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            m = m,
            dropout_fit = dropout_fit1,
            m_dropout = m_dropout,
            fixedFollowup = fixedFollowup,
            followupTime = followupTime,
            pilevel = pilevel,
            nyears = nyears,
            target_t = target_t,
            nreps = nreps,
            showEnrollment = showEnrollment,
            showEvent = showEvent,
            showDropout = showDropout,
            showOngoing = showOngoing,
            showsummary = showsummary,
            showplot = FALSE,
            by_treatment = by_treatment,
            covariates_event = covariates_event,
            event_fit_with_covariates = event_fit1_w_x,
            covariates_dropout = covariates_dropout,
            dropout_fit_with_covariates = dropout_fit1_w_x,
            fix_parameter = fix_parameter,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            seed.num = seed.num,  # UPDATED (merge)
            return_new_events = needs_event_subjects
          )
        } else {
          event_pred <- predictEvent(
            df = dt,
            target_d = target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            m = m,
            dropout_fit = dropout_fit1,
            m_dropout = m_dropout,
            fixedFollowup = fixedFollowup,
            followupTime = followupTime,
            pilevel = pilevel,
            nyears = nyears,
            target_t = target_t,
            nreps = nreps,
            showEnrollment = showEnrollment,
            showEvent = showEvent,
            showDropout = showDropout,
            showOngoing = showOngoing,
            showsummary = showsummary,
            showplot = FALSE,
            by_treatment = by_treatment,
            covariates_event = covariates_event,
            event_fit_with_covariates = event_fit1_w_x,
            covariates_dropout = covariates_dropout,
            dropout_fit_with_covariates = dropout_fit1_w_x,
            fix_parameter = fix_parameter,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            seed.num = seed.num,  # UPDATED (merge)
            return_new_events = needs_event_subjects
          )
        }
      } else {  # no dropout model
        if (grepl("enrollment", to_predict, ignore.case = TRUE)) {
          event_pred <- predictEvent(
            df = dt,
            target_d = target_d,
            newSubjects = enroll_pred$newSubjects,
            event_fit = event_fit1,
            m = m,
            dropout_fit = NULL,
            m_dropout = m_dropout,
            fixedFollowup = fixedFollowup,
            followupTime = followupTime,
            pilevel = pilevel,
            nyears = nyears,
            target_t = target_t,
            nreps = nreps,
            showEnrollment = showEnrollment,
            showEvent = showEvent,
            showDropout = showDropout,
            showOngoing = showOngoing,
            showsummary = showsummary,
            showplot = FALSE,
            by_treatment = by_treatment,
            covariates_event = covariates_event,
            event_fit_with_covariates = event_fit1_w_x,
            covariates_dropout = covariates_dropout,
            dropout_fit_with_covariates = NULL,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            seed.num = seed.num,  # UPDATED (merge)
            return_new_events = needs_event_subjects
          )
        } else {
          event_pred <- predictEvent(
            df = dt,
            target_d = target_d,
            newSubjects = NULL,
            event_fit = event_fit1,
            m = m,
            dropout_fit = NULL,
            m_dropout = m_dropout,
            fixedFollowup = fixedFollowup,
            followupTime = followupTime,
            pilevel = pilevel,
            nyears = nyears,
            target_t = target_t,
            nreps = nreps,
            showEnrollment = showEnrollment,
            showEvent = showEvent,
            showDropout = showDropout,
            showOngoing = showOngoing,
            showsummary = showsummary,
            showplot = FALSE,
            by_treatment = by_treatment,
            covariates_event = covariates_event,
            event_fit_with_covariates = event_fit1_w_x,
            covariates_dropout = covariates_dropout,
            dropout_fit_with_covariates = NULL,
            generate_plot = generate_plot,
            interactive_plot = interactive_plot,
            seed.num = seed.num,  # UPDATED (merge)
            return_new_events = needs_event_subjects
          )
        }

      }
    } else { # event prediction at design stage
      if (!is.null(dropout_prior)) {
        event_pred <- predictEvent(
          df = NULL,
          target_d = target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior,
          m = m,
          dropout_fit = dropout_prior,
          m_dropout = m_dropout,
          fixedFollowup = fixedFollowup,
          followupTime = followupTime,
          pilevel = pilevel,
          nyears = nyears,
          target_t = target_t,
          nreps = nreps,
          showEnrollment = showEnrollment,
          showEvent = showEvent,
          showDropout = showDropout,
          showOngoing = showOngoing,
          showsummary = showsummary,
          showplot = FALSE,
          by_treatment = by_treatment,
          covariates_event = covariates_event,
          event_fit_with_covariates = event_prior_w_x,
          covariates_dropout = covariates_dropout,
          dropout_fit_with_covariates = dropout_prior_w_x,
          fix_parameter = fix_parameter,
          generate_plot = generate_plot,
          interactive_plot = interactive_plot,
          seed.num = seed.num,  # UPDATED (merge)
          return_new_events = needs_event_subjects
        )
      } else {
        event_pred <- predictEvent(
          df = NULL,
          target_d = target_d,
          newSubjects = enroll_pred$newSubjects,
          event_fit = event_prior,
          m = m,
          dropout_fit = NULL,
          m_dropout = m_dropout,
          fixedFollowup = fixedFollowup,
          followupTime = followupTime,
          pilevel = pilevel,
          nyears = nyears,
          target_t = target_t,
          nreps = nreps,
          showEnrollment = showEnrollment,
          showEvent = showEvent,
          showDropout = showDropout,
          showOngoing = showOngoing,
          showsummary = showsummary,
          showplot = FALSE,
          by_treatment = by_treatment,
          covariates_event = covariates_event,
          event_fit_with_covariates = event_prior_w_x,
          covariates_dropout = covariates_dropout,
          dropout_fit_with_covariates = dropout_prior_w_x,
          fix_parameter = fix_parameter,
          generate_plot = generate_plot,
          interactive_plot = interactive_plot,
          seed.num = seed.num,  # UPDATED (merge)
          return_new_events = needs_event_subjects
        )

      }
    }
  }


  # obtain subject-level data from all subjects
  if (return_subject_data) {
    subject_data <- .ep_build_prediction_subject_data(
      dt = if (!is.null(df)) dt else NULL,
      to_predict = to_predict,
      enroll_pred = enroll_pred,
      event_pred = event_pred,
      by_treatment = by_treatment
    )
  }
  
  enroll_pred <- .ep_trim_prediction_output(enroll_pred, return_simulation_data)
  event_pred <- .ep_trim_prediction_output(event_pred, return_simulation_data)

  observed_out <- observed
  enroll_fit_out <- if (!is.null(df)) enroll_fit else enroll_prior
  event_fit_out <- if (!is.null(event_fit)) event_fit else event_prior
  event_fit_w_x_out <- event_fit_w_x
  dropout_fit_out <- if (!is.null(dropout_fit)) dropout_fit else dropout_prior
  dropout_fit_w_x_out <- dropout_fit_w_x

  .ep_maybe_print_prediction_plot(
    to_predict_lc = to_predict_lc,
    generate_plot = generate_plot,
    showplot = showplot,
    enroll_pred = enroll_pred,
    event_pred = event_pred
  )

  .ep_build_get_prediction_result(
    has_observed_data = !is.null(df),
    to_predict_lc = to_predict_lc,
    observed = observed_out,
    enroll_fit = enroll_fit_out,
    enroll_pred = enroll_pred,
    event_fit = event_fit_out,
    event_fit_with_covariates = event_fit_w_x_out,
    dropout_fit = dropout_fit_out,
    dropout_fit_with_covariates = dropout_fit_w_x_out,
    event_pred = event_pred,
    dropout_model = dropout_model,
    covariates_event = covariates_event,
    covariates_dropout = covariates_dropout,
    subject_data = subject_data,
    return_subject_data = return_subject_data
  )
}
