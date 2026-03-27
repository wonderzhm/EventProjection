#' @title Fit time-to-dropout model
#' @description Fits a specified time-to-dropout model to the dropout data.
#'
#' @details
#' This merged implementation is based on the original \code{fitDropout()} from the
#' \pkg{eventPred} package by Kaifeng Lu, with additional user-interface features
#' adapted from \pkg{EventPredInCure}:
#' \itemize{
#'   \item Added \code{criterion} ("aic", "bic", or "both") to control which
#'   information criteria are shown in the figure annotation.
#'   \item Standardized column names in \code{df} to lower case to make the
#'   function more robust to input case.
#'   \item Added a small positional-argument compatibility shim for legacy
#'   \pkg{EventPredInCure} calls of the form
#'   \code{fitDropout(df, dropout_model, piecewiseDropoutTime, by_treatment, criterion)}.
#' }
#'
#' The statistical model fitting logic and returned object structure follow
#' \pkg{eventPred} to retain full functionality (covariates, model averaging,
#' spline, Cox, interactive/static plots).
#'
#' @param df The subject-level dropout data, including \code{time} and
#'   \code{dropout}. The data should also include \code{treatment}
#'   coded as 1, 2, and so on, and \code{treatment_description}
#'   for fitting the dropout model by treatment.
#' @param dropout_model The dropout model used to analyze the dropout data
#'   which can be set to one of the following options:
#'   "exponential", "Weibull", "log-logistic", "log-normal",
#'   "piecewise exponential", "model averaging", "spline", or "cox".
#'   The model averaging uses the \code{exp(-bic/2)} weighting and
#'   combines Weibull and log-normal models. The spline model of
#'   Royston and Parmar (2002) assumes that a transformation of
#'   the survival function is modeled as a natural cubic spline
#'   function of log time. By default, it is set to "exponential".
#' @param piecewiseDropoutTime A vector that specifies the time
#'   intervals for the piecewise exponential dropout distribution.
#'   Must start with 0, e.g., c(0, 60) breaks the time axis into 2
#'   event intervals: [0, 60) and [60, Inf). By default, it is set to 0.
#' @param k_dropout The number of inner knots of the spline. The default
#'   \code{k_dropout=0} gives a Weibull, log-logistic or log-normal model,
#'   if \code{scale_dropout} is "hazard", "odds", or "normal", respectively.
#'   The knots are chosen as equally-spaced quantiles of the log
#'   uncensored survival times. The boundary knots are chosen as the
#'   minimum and maximum log uncensored survival times.
#' @param scale_dropout  The scale of the spline. The default is "hazard",
#'   in which case the log cumulative hazard is modeled as a spline
#'   function. If \code{scale = "odds"}, the log cumulative odds is
#'   modeled as a spline function. If \code{scale = "normal"},
#'   \code{-qnorm(S(t))} is modeled as a spline function.
#' @param m_dropout The number of dropout time intervals to extrapolate
#'   the hazard function beyond the last observed dropout time when
#'   \code{dropout_model = "cox"}.
#' @param showplot A Boolean variable to control whether or not to
#'   show the fitted time-to-dropout survival curve. By default, it is
#'   set to \code{TRUE}.
#' @param by_treatment A Boolean variable to control whether or not to
#'   fit the time-to-dropout data by treatment group. By default,
#'   it is set to \code{FALSE}.
#' @param covariates The names of baseline covariates from the input
#'   data frame to include in the dropout model, e.g., c("age", "sex").
#'   Factor variables need to be declared in the input data frame.
#' @param generate_plot Whether to generate plots.
#' @param interactive_plot Whether to produce interactive plots using
#'   plotly or static plots using ggplot2.
#' @param criterion A character variable to denote which information criterion
#'   is shown in the plot annotation, which can be set to one of the following
#'   options: "aic", "bic", or "both". By default, it is set to \code{"both"}.
#'   This parameter only affects plot annotations; model fitting is unchanged.
#'
#' @return
#' A list of results from the model fit including key information
#' such as the dropout model (\code{model}), the estimated model parameters
#' (\code{theta}), covariance matrix (\code{vtheta}), and \code{aic}/\code{bic}.
#'
#' If \code{generate_plot=TRUE}, the fitted time-to-dropout survival curve is
#' also returned (interactive plotly object or static ggplot object).
#'
#' When fitting the dropout model by treatment, the output is a list of lists,
#' with one element per treatment group.
#'
#' @references
#' Patrick Royston and Mahesh K. B. Parmar. Flexible parametric
#' proportional-hazards and proportional-odds models for censored
#' survival data, with application to prognostic modelling and
#' estimation of treatment effects. Stat in Med. 2002; 21:2175-2197.
#'
#' @examples
#' dropout_fit <- fitDropout(
#'   df = interimData2, showplot = FALSE,
#'   dropout_model = "exponential")
#'
#' @export
#'
fitDropout <- function(df, dropout_model = "exponential",
                       piecewiseDropoutTime = 0,
                       k_dropout = 0, scale_dropout = "hazard",
                       m_dropout = 5,
                       showplot = TRUE, by_treatment = FALSE,
                       covariates = NULL,
                       generate_plot = TRUE,
                       interactive_plot = TRUE,
                       criterion = "both") {
  
  # ---- Compatibility shim (EventPredInCure legacy positional usage) ----
  # Legacy call pattern:
  #   fitDropout(df, dropout_model, piecewiseDropoutTime, by_treatment, criterion)
  # In that case:
  #   k_dropout receives a logical (by_treatment),
  #   scale_dropout receives a string in c("aic","bic","both") (criterion).
  if (is.logical(k_dropout) && length(k_dropout) == 1L && missing(by_treatment)) {
    by_treatment <- k_dropout
    k_dropout <- 0
  }
  if (is.character(scale_dropout) &&
      length(scale_dropout) == 1L &&
      tolower(scale_dropout) %in% c("aic", "bic", "both") &&
      missing(criterion)) {
    criterion <- tolower(scale_dropout)
    scale_dropout <- "hazard"
  }
  
  erify::check_class(df, "data.frame")
  
  erify::check_content(tolower(criterion), c("aic", "bic", "both"))
  
  erify::check_content(tolower(dropout_model),
                       c("exponential", "weibull", "log-logistic",
                         "log-normal", "piecewise exponential",
                         "model averaging", "spline", "cox"))
  
  if (piecewiseDropoutTime[1] != 0) {
    stop("piecewiseDropoutTime must start with 0")
  }
  if (length(piecewiseDropoutTime) > 1 && any(diff(piecewiseDropoutTime) <= 0)) {
    stop("piecewiseDropoutTime should be increasing")
  }
  
  erify::check_n(k_dropout, zero = TRUE)
  erify::check_content(tolower(scale_dropout), c("hazard", "odds", "normal"))
  erify::check_n(m_dropout)
  
  erify::check_bool(showplot)
  erify::check_bool(by_treatment)
  
  # ---- Standardize input column names (merged change) ----
  dt <- data.table::as.data.table(data.table::copy(df))
  data.table::setnames(dt, tolower(names(dt)))
  
  # minimal required columns
  if (!all(c("time", "dropout") %in% names(dt))) {
    stop("df must contain columns 'time' and 'dropout' (case-insensitive).")
  }
  if (by_treatment && !("treatment" %in% names(dt))) {
    stop("df must contain column 'treatment' (case-insensitive) when by_treatment=TRUE.")
  }
  
  # construct the formula for survival analysis
  if (!is.null(covariates)) {
    covariates <- tolower(covariates)
    if (!all(covariates %in% names(dt))) {
      stop("All covariates must exist in df (case-insensitive).")
    }
    xnames <- paste(covariates, collapse = "+")
    formula <- stats::as.formula(paste("survival::Surv(time, dropout) ~", xnames))
  } else {
    formula <- survival::Surv(time, dropout) ~ 1
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
  
  if (ngroups == 1) {
    by_treatment <- FALSE
  }
  
  # ---- Fit by treatment group ----
  dropout_fit <- list()
  
  for (i in 1:ngroups) {
    df1 <- dt[get("treatment") == i]
    
    n0 <- nrow(df1)
    c0 <- df1[, sum(get("dropout"))]
    ex0 <- df1[, sum(get("time"))]
    
    x <- model.matrix(formula, df1)
    q <- ncol(x) - 1
    
    kmfit <- survival::survfit(survival::Surv(time, dropout) ~ 1, data = df1)
    kmdf <- data.table::data.table(time = kmfit$time, surv = kmfit$surv)
    df0 <- data.table::data.table(time = 0, surv = 1)
    kmdf <- data.table::rbindlist(list(df0, kmdf), use.names = TRUE)
    
    if (tolower(dropout_model) == "exponential") {
      erify::check_positive(c0 - q, supplement = paste(
        "The number of dropouts must be >=", q + 1,
        "to fit an exponential model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "exponential")
      
      fit3 <- list(model = "Exponential",
                   theta = -as.numeric(reg$coefficients),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 1),
                   bic = -2 * reg$loglik[1] + (q + 1) * log(n0))
      
      rate <- exp(as.numeric(x %*% fit3$theta))
      
      dffit3 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::pexp(t, rate, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "weibull") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a Weibull model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "weibull")
      
      fit3 <- list(model = "Weibull",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      shape <- exp(-fit3$theta[q + 2])
      scale <- exp(as.numeric(x %*% fit3$theta[1:(q + 1)]))
      
      dffit3 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::pweibull(t, shape, scale, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "log-logistic") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a log-logistic model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "loglogistic")
      
      fit3 <- list(model = "Log-logistic",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      location <- as.numeric(x %*% fit3$theta[1:(q + 1)])
      scale <- exp(fit3$theta[q + 2])
      
      dffit3 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::plogis(log(t), location, scale, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "log-normal") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a log-normal model."))
      
      reg <- survival::survreg(formula, data = df1, dist = "lognormal")
      
      fit3 <- list(model = "Log-normal",
                   theta = c(as.numeric(reg$coefficients), log(reg$scale)),
                   vtheta = reg$var,
                   aic = -2 * reg$loglik[1] + 2 * (q + 2),
                   bic = -2 * reg$loglik[1] + (q + 2) * log(n0))
      
      meanlog <- as.numeric(x %*% fit3$theta[1:(q + 1)])
      sdlog <- exp(fit3$theta[q + 2])
      
      dffit3 <- data.table::data.table(time = seq(0, max(df1$time)))[
        , `:=`(surv = sapply(get("time"), function(t)
          mean(stats::plnorm(t, meanlog, sdlog, lower.tail = FALSE))))]
    } else if (tolower(dropout_model) == "piecewise exponential") {
      J <- length(piecewiseDropoutTime)
      
      erify::check_positive(c0 - J - q + 1, supplement = paste(
        "The number of dropouts must be >=", J + q,
        "to fit a piecewise exponential model."))
      
      fit3 <- pwexpreg(df1$time, df1$dropout, J, piecewiseDropoutTime, q, x)
      
      time <- seq(0, max(df1$time))
      
      surv <- purrr::map(1:n0, function(l)
        ppwexp(time, fit3$theta, J, fit3$piecewiseDropoutTime, q, x[l,],
               lower.tail = FALSE))
      surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      
      dffit3 <- data.table::data.table(time, surv)
    } else if (tolower(dropout_model) == "model averaging") {
      erify::check_positive(c0 - q - 1, supplement = paste(
        "The number of dropouts must be >=", q + 2,
        "to fit a model averaging model."))
      
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
      
      fit3 <- list(model = "Model averaging",
                   theta = theta,
                   vtheta = vtheta,
                   aic = w1 * aic1 + (1 - w1) * aic2,
                   bic = w1 * bic1 + (1 - w1) * bic2,
                   w1 = w1)
      
      time <- seq(0, max(df1$time))
      
      surv <- purrr::map(1:n0, function(l)
        pmodavg(time, fit3$theta, w1, q, x[l,], lower.tail = FALSE))
      surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      
      dffit3 <- data.table::data.table(time, surv)
    } else if (tolower(dropout_model) == "spline") {
      erify::check_positive(c0 - k_dropout - q - 1, supplement = paste(
        "The number of dropouts must be >=", k_dropout + q + 2,
        "to fit a spline model."))
      
      spl <- flexsurv::flexsurvspline(
        formula, data = df1, k = k_dropout, scale = scale_dropout,
        method = "Nelder-Mead")
      
      fit3 <- list(model = "Spline",
                   theta = as.numeric(spl$coefficients),
                   vtheta = spl$cov,
                   aic = -2 * spl$loglik + 2 * (k_dropout + q + 2),
                   bic = -2 * spl$loglik + (k_dropout + q + 2) * log(n0),
                   knots = spl$knots,
                   scale = spl$scale)
      
      time <- seq(0, max(df1$time))
      
      if (q > 0) {
        xbeta <- as.numeric(as.matrix(x[, -1]) %*%
                              fit3$theta[(k_dropout + 3):(k_dropout + q + 2)])
        
        surv <- purrr::map(1:n0, function(l)
          flexsurv::psurvspline(
            time, gamma = fit3$theta[1:(k_dropout + 2)], knots = fit3$knots,
            scale = fit3$scale, offset = xbeta[l], lower.tail = FALSE))
        surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
      } else {
        surv <- flexsurv::psurvspline(
          time, gamma = fit3$theta, knots = fit3$knots, scale = fit3$scale,
          lower.tail = FALSE)
      }
      
      dffit3 <- data.table::data.table(time, surv)
    } else if (tolower(dropout_model) == "cox") {
      erify::check_positive(c0 - q, supplement = paste(
        "The number of dropouts must be >=", q + 1,
        "to fit a Cox model."))
      
      erify::check_positive(c0 - m_dropout + 1, supplement = paste(
        "m_dropout must be <= the observed number of dropouts", c0))
      
      reg <- phregr(df1, time = "time", event = "dropout",
                    covariates = covariates)
      
      bh <- data.table::setDT(reg$basehaz)[get("nevent") > 0]
      haz <- bh$haz
      vhaz <- bh$varhaz
      
      if (q > 0) {
        if (q == 1) {
          ghaz <- matrix(bh$gradhaz, ncol = 1)
        } else {
          ghaz <- do.call(cbind, lapply(1:q, function(ii)
            bh[[paste0("gradhaz.", ii)]]))
        }
      }
      
      M <- nrow(bh)
      tcut <- c(0, bh$time)
      dt_int <- diff(tcut)
      lambda1 <- haz / dt_int
      
      d <- bh$nevent
      llik <- reg$sumstat$loglik1 + sum(d * (log(d / dt_int) - 1))
      
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
      
      fit3 <- list(model = "Cox",
                   theta = theta,
                   vtheta = vtheta,
                   aic = -2 * llik + 2 * (M + q),
                   bic = -2 * llik + (M + q) * log(n0),
                   piecewiseDropoutTime = tcut)
      
      time <- seq(0, max(df1$time))
      
      lambda2 <- sum(bh$haz[(M - m_dropout + 1):M]) /
        (bh$time[M] - bh$time[M - m_dropout])
      
      lambda <- c(lambda1, lambda2)
      
      s1 <- sapply(time, function(t)
        ppwexp(t, log(lambda), M + 1, tcut, lower.tail = FALSE))
      
      if (q > 0) {
        xbeta <- as.numeric(as.matrix(x[, -1]) %*% reg$beta)
        surv <- apply(outer(s1, exp(xbeta), `^`), 1, mean)
      } else {
        surv <- s1
      }
      
      dffit3 <- data.table::data.table(time, surv)
    } else {
      stop("Unsupported dropout_model.")
    }
    
    # ---- Build plot annotation text (merged change: criterion) ----
    if (tolower(fit3$model) == "piecewise exponential") {
      modeltext <- paste0(paste0(fit3$model, "("),
                          paste(piecewiseDropoutTime, collapse = ","), ")")
    } else if (tolower(fit3$model) == "spline") {
      modeltext <- paste0(fit3$model, "(k = ", k_dropout, ", ", "scale = '",
                          scale_dropout, "')")
    } else if (tolower(fit3$model) == "cox") {
      modeltext <- paste0(fit3$model, "(m = ", m_dropout, ")")
    } else {
      modeltext <- fit3$model
    }
    
    if (tolower(dropout_model) == "model averaging") {
      aictext <- paste("Weighted AIC:", formatC(fit3$aic, format = "f", digits = 2))
      bictext <- paste("Weighted BIC:", formatC(fit3$bic, format = "f", digits = 2))
    } else {
      aictext <- paste("AIC:", formatC(fit3$aic, format = "f", digits = 2))
      bictext <- paste("BIC:", formatC(fit3$bic, format = "f", digits = 2))
    }
    
    crit <- tolower(criterion)
    text_vec <- c(modeltext)
    if (crit %in% c("aic", "both")) text_vec <- c(text_vec, aictext)
    if (crit %in% c("bic", "both")) text_vec <- c(text_vec, bictext)
    
    # plot the survival curve
    if (generate_plot) {
      if (interactive_plot) {
        y_seq <- 0.95 - 0.05 * (0:(length(text_vec) - 1))
        fittedDropout <- plotly::plot_ly() %>%
          plotly::add_lines(
            data = kmdf, x = ~time, y = ~surv, name = "Kaplan-Meier",
            line = list(shape = "hv")) %>%
          plotly::add_lines(
            data = dffit3, x = ~time, y = ~surv, name = "fitted") %>%
          plotly::layout(
            xaxis = list(title = "Days since randomization",
                         zeroline = FALSE),
            yaxis = list(title = "Survival probability", zeroline = FALSE),
            title = list(text = "Fitted time to dropout survival curve"),
            annotations = list(
              x = rep(0.7, length(text_vec)),
              y = y_seq,
              xref = "paper",
              yref = "paper",
              text = paste0("<i>", text_vec, "</i>"),
              xanchor = "left",
              font = list(size = 14, color = "red"),
              showarrow = FALSE)) %>%
          plotly::hide_legend()
      } else {
        x_pos <- max(kmdf$time, na.rm = TRUE)
        y_pos <- max(kmdf$surv, na.rm = TRUE)
        
        fittedDropout <- ggplot2::ggplot() +
          ggplot2::geom_step(data = kmdf[-1, ], ggplot2::aes(
            x = .data$time, y = .data$surv)) +
          ggplot2::geom_line(data = dffit3[-1, ], ggplot2::aes(
            x = .data$time, y = .data$surv), colour = "red") +
          ggplot2::labs(
            x = "Days since randomization",
            y = "Survival probability",
            title = "Fitted time to dropout survival curve") +
          ggplot2::theme(legend.position = "none")
        
        for (j in seq_along(text_vec)) {
          fittedDropout <- fittedDropout +
            ggplot2::annotate("text", x = x_pos, y = y_pos, label = text_vec[j],
                              hjust = 1.1, vjust = 1.5 + 2 * (j - 1),
                              size = 5, colour = "red", fontface = "italic")
        }
      }
    }
    
    if (by_treatment) {
      if (generate_plot) {
        if (interactive_plot) {
          fittedDropout <- fittedDropout %>%
            plotly::layout(annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", df1$treatment_description[1], "</b>"),
              xanchor = "center", yanchor = "middle", showarrow = FALSE,
              xref = "paper", yref = "paper"))
        } else {
          fittedDropout <- fittedDropout +
            ggplot2::labs(subtitle = df1$treatment_description[1])
        }
      }
      
      fit3$treatment <- df1$treatment[1]
      fit3$treatment_description <- df1$treatment_description[1]
    }
    
    if (generate_plot) {
      dropout_fit[[i]] <- list(fit = fit3, fit_plot = fittedDropout,
                               kmdf = kmdf, dffit = dffit3,
                               text = text_vec)
    } else {
      dropout_fit[[i]] <- list(fit = fit3, kmdf = kmdf, dffit = dffit3,
                               text = text_vec)
    }
  }
  
  # ensure that the sub plots share the same x axis range
  if (by_treatment) {
    if (generate_plot) {
      x_range <- range(dt$time)
      for (i in 1:ngroups) {
        if (interactive_plot) {
          dropout_fit[[i]]$fit_plot <- dropout_fit[[i]]$fit_plot %>%
            plotly::layout(xaxis = list(range = x_range))
        } else {
          dropout_fit[[i]]$fit_plot <- dropout_fit[[i]]$fit_plot +
            ggplot2::scale_x_continuous(limits = x_range)
        }
      }
    }
  } else {
    if (generate_plot) {
      dropout_fit <- list(fit = fit3, fit_plot = fittedDropout,
                          kmdf = kmdf, dffit = dffit3,
                          text = text_vec)
    } else {
      dropout_fit <- list(fit = fit3, kmdf = kmdf, dffit = dffit3,
                          text = text_vec)
    }
  }
  
  if (generate_plot && showplot) {
    if (by_treatment) {
      for (i in 1:ngroups) {
        print(dropout_fit[[i]]$fit_plot)
      }
    } else {
      print(dropout_fit$fit_plot)
    }
  }
  
  dropout_fit
}
