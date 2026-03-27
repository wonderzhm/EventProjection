.ep_build_ic_text <- function(modeltext, aic, bic, criterion_lc,
                              weighted = FALSE, extra_lines = NULL) {
  prefix <- if (weighted) "Weighted " else ""
  text_vec <- c(modeltext)

  if (criterion_lc %in% c("aic", "both")) {
    text_vec <- c(
      text_vec,
      paste0(prefix, "AIC: ", formatC(aic, format = "f", digits = 2))
    )
  }

  if (criterion_lc %in% c("bic", "both")) {
    text_vec <- c(
      text_vec,
      paste0(prefix, "BIC: ", formatC(bic, format = "f", digits = 2))
    )
  }

  if (!is.null(extra_lines)) {
    text_vec <- c(text_vec, extra_lines)
  }

  text_vec
}


.ep_standardize_columns <- function(x, canonical) {
  if (is.null(x)) {
    return(x)
  }

  nms <- names(x)
  low <- tolower(nms)

  for (nm in canonical) {
    j <- match(tolower(nm), low)
    if (!is.na(j) && nms[j] != nm) {
      names(x)[j] <- nm
    }
  }

  x
}


.ep_format_enrollment_fit_label <- function(fit, accrual_time, nknots) {
  model_lc <- tolower(fit$model)

  if (model_lc == "piecewise poisson") {
    return(paste0(fit$model, "(", paste(accrual_time, collapse = " "), ")"))
  }

  if (model_lc == "b-spline") {
    return(paste0(fit$model, "(nknots = ", nknots, ")"))
  }

  fit$model
}


.ep_build_enrollment_fit_text <- function(fit, accrual_time, nknots,
                                          criterion_lc) {
  .ep_build_ic_text(
    modeltext = .ep_format_enrollment_fit_label(
      fit = fit,
      accrual_time = accrual_time,
      nknots = nknots
    ),
    aic = fit$aic,
    bic = fit$bic,
    criterion_lc = criterion_lc
  )
}


.ep_build_enrollment_fit_plot <- function(enrolldf, dffit, text_vec,
                                          interactive_plot) {
  if (interactive_plot) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive_plot=TRUE.")
    }

    y_pos <- seq(0.95, by = -0.05, length.out = length(text_vec))
    x_pos <- rep(0.05, length(text_vec))

    return(
      plotly::plot_ly() %>%
        plotly::add_lines(
          data = enrolldf,
          x = ~t,
          y = ~n,
          name = "observed",
          line = list(shape = "hv")
        ) %>%
        plotly::add_lines(
          data = dffit,
          x = ~t,
          y = ~n,
          name = "fitted"
        ) %>%
        plotly::layout(
          xaxis = list(title = "Days since trial start", zeroline = FALSE),
          yaxis = list(title = "Subjects", zeroline = FALSE),
          title = list(text = "Fitted enrollment curve"),
          annotations = list(
            x = x_pos,
            y = y_pos,
            xref = "paper",
            yref = "paper",
            text = paste0("<i>", text_vec, "</i>"),
            xanchor = "left",
            font = list(size = 14, color = "red"),
            showarrow = FALSE
          )
        ) %>%
        plotly::hide_legend()
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for interactive_plot=FALSE.")
  }

  fitted_enroll <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = enrolldf,
      ggplot2::aes(x = .data$t, y = .data$n)
    ) +
    ggplot2::geom_line(
      data = dffit,
      ggplot2::aes(x = .data$t, y = .data$n),
      colour = "red"
    ) +
    ggplot2::labs(
      x = "Days since trial start",
      y = "Subjects",
      title = "Fitted enrollment curve"
    ) +
    ggplot2::theme(legend.position = "none")

  for (i in seq_along(text_vec)) {
    fitted_enroll <- fitted_enroll +
      ggplot2::annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = text_vec[i],
        hjust = -0.1,
        vjust = 1.5 + 2 * (i - 1),
        colour = "red",
        size = 5,
        fontface = "italic"
      )
  }

  fitted_enroll
}


.ep_build_enrollment_fit_output <- function(fit, enrolldf, dffit, text_vec,
                                            fit_plot = NULL) {
  if (is.null(fit_plot)) {
    return(list(
      fit = fit,
      enrolldf = enrolldf,
      dffit = dffit,
      text = text_vec
    ))
  }

  list(
    fit = fit,
    fit_plot = fit_plot,
    enrolldf = enrolldf,
    dffit = dffit,
    text = text_vec
  )
}


.ep_format_event_fit_label <- function(fit, piecewise_survival_time,
                                       k, scale, m) {
  model_lc <- tolower(fit$model)

  if (grepl("^piecewise exponential", model_lc)) {
    ps <- fit$piecewiseSurvivalTime
    if (is.null(ps)) {
      ps <- piecewise_survival_time
    }
    return(paste0(fit$model, "(", paste(ps, collapse = ","), ")"))
  }

  if (model_lc == "spline") {
    return(paste0(fit$model, "(k = ", k, ", scale = '", scale, "')"))
  }

  if (model_lc == "cox") {
    return(paste0(fit$model, "(m = ", m, ")"))
  }

  fit$model
}


.ep_build_event_fit_text <- function(fit, event_model_lc, criterion_lc,
                                     piecewise_survival_time, k, scale, m,
                                     p_cure = NULL) {
  extra_lines <- NULL
  if (!is.null(p_cure)) {
    extra_lines <- paste(
      "Prop of Cured:",
      formatC(p_cure, format = "f", digits = 2)
    )
  }

  .ep_build_ic_text(
    modeltext = .ep_format_event_fit_label(
      fit = fit,
      piecewise_survival_time = piecewise_survival_time,
      k = k,
      scale = scale,
      m = m
    ),
    aic = fit$aic,
    bic = fit$bic,
    criterion_lc = criterion_lc,
    weighted = identical(event_model_lc, "model averaging"),
    extra_lines = extra_lines
  )
}


.ep_build_event_fit_plot <- function(kmdf, kmdf_censor, dffit, text_vec,
                                     interactive_plot) {
  if (interactive_plot) {
    y_seq <- 0.95 - 0.05 * (0:(length(text_vec) - 1))

    return(
      plotly::plot_ly() %>%
        plotly::add_lines(
          data = kmdf,
          x = ~time,
          y = ~surv,
          name = "Kaplan-Meier",
          line = list(shape = "hv")
        ) %>%
        plotly::add_markers(
          data = kmdf_censor,
          x = ~time,
          y = ~surv,
          name = "censoring",
          showlegend = FALSE,
          marker = list(
            symbol = "cross-thin",
            size = 6,
            line = list(color = "black", width = 2)
          )
        ) %>%
        plotly::add_lines(
          data = dffit,
          x = ~time,
          y = ~surv,
          name = "fitted"
        ) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
          yaxis = list(title = "Survival probability", zeroline = FALSE),
          title = list(text = "Fitted time to event survival curve"),
          annotations = list(
            x = rep(0.7, length(text_vec)),
            y = y_seq,
            xref = "paper",
            yref = "paper",
            text = paste0("<i>", text_vec, "</i>"),
            xanchor = "left",
            font = list(size = 14, color = "red"),
            showarrow = FALSE
          )
        ) %>%
        plotly::hide_legend()
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for interactive_plot=FALSE.")
  }

  fitted_event <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = kmdf[-1, ],
      ggplot2::aes(x = .data$time, y = .data$surv)
    ) +
    ggplot2::geom_point(
      data = kmdf_censor,
      ggplot2::aes(x = .data$time, y = .data$surv),
      shape = 4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = dffit[-1, ],
      ggplot2::aes(x = .data$time, y = .data$surv),
      colour = "red"
    ) +
    ggplot2::labs(
      x = "Days since randomization",
      y = "Survival probability",
      title = "Fitted time to event survival curve"
    ) +
    ggplot2::theme(legend.position = "none")

  x_pos <- max(kmdf$time, na.rm = TRUE)
  y_pos <- max(kmdf$surv, na.rm = TRUE)

  for (j in seq_along(text_vec)) {
    fitted_event <- fitted_event +
      ggplot2::annotate(
        "text",
        x = x_pos,
        y = y_pos,
        label = text_vec[j],
        hjust = 1.1,
        vjust = 1.5 + 2 * (j - 1),
        size = 5,
        colour = "red",
        fontface = "italic"
      )
  }

  fitted_event
}


.ep_add_event_fit_treatment_label <- function(plot_obj, label,
                                              interactive_plot) {
  if (interactive_plot) {
    existing_annotations <- plot_obj$x$layout$annotations

    return(
      plot_obj %>%
        plotly::layout(annotations = c(
          list(list(
            x = 0.5,
            y = 1,
            text = paste0("<b>", label, "</b>"),
            xanchor = "center",
            yanchor = "middle",
            showarrow = FALSE,
            xref = "paper",
            yref = "paper"
          )),
          existing_annotations
        ))
    )
  }

  plot_obj + ggplot2::labs(subtitle = label)
}


.ep_build_event_fit_output <- function(fit, kmdf, dffit, text_vec,
                                       fit_plot = NULL) {
  if (is.null(fit_plot)) {
    return(list(
      fit = fit,
      kmdf = kmdf,
      dffit = dffit,
      text = text_vec
    ))
  }

  list(
    fit = fit,
    fit_plot = fit_plot,
    kmdf = kmdf,
    dffit = dffit,
    text = text_vec
  )
}


.ep_sync_event_fit_plot_ranges <- function(event_fit, x_range,
                                           interactive_plot) {
  for (i in seq_along(event_fit)) {
    if (interactive_plot) {
      event_fit[[i]]$fit_plot <- event_fit[[i]]$fit_plot %>%
        plotly::layout(xaxis = list(range = x_range))
    } else {
      event_fit[[i]]$fit_plot <- event_fit[[i]]$fit_plot +
        ggplot2::scale_x_continuous(limits = x_range)
    }
  }

  event_fit
}


.ep_build_prediction_output <- function(..., subject_data = NULL,
                                        return_subject_data = TRUE) {
  out <- list(...)

  if (return_subject_data) {
    out$subject_data <- subject_data
  }

  out
}


.ep_trim_prediction_output <- function(pred, return_simulation_data) {
  if (is.null(pred)) {
    return(pred)
  }

  if (!return_simulation_data) {
    pred$newSubjects <- NULL
    pred$newEvents <- NULL
  }

  pred
}


.ep_pool_priors_moment_match <- function(prior_list, w) {
  if (length(prior_list) != length(w)) {
    stop("w must have same length as prior_list")
  }

  w <- w / sum(w)
  theta_bar <- 0
  second_moment <- 0

  for (i in seq_along(prior_list)) {
    th <- prior_list[[i]]$theta
    vt <- prior_list[[i]]$vtheta

    if (!is.matrix(vt)) {
      vt <- matrix(vt, nrow = length(th), ncol = length(th))
    }

    theta_bar <- theta_bar + w[i] * th
    second_moment <- second_moment + w[i] * (vt + th %*% t(th))
  }

  vtheta_bar <- second_moment - theta_bar %*% t(theta_bar)
  list(theta = theta_bar, vtheta = vtheta_bar)
}


.ep_build_prediction_subject_data <- function(dt = NULL, to_predict,
                                              enroll_pred = NULL,
                                              event_pred = NULL,
                                              by_treatment = FALSE) {
  to_predict_lc <- tolower(to_predict)

  if (to_predict_lc == "enrollment only") {
    subject_data <- enroll_pred$newSubjects

    if (!is.null(dt)) {
      observed_dt <- data.table::copy(dt)
      observed_dt[, `:=`(
        arrivalTime = as.numeric(get("randdt") - get("trialsdt") + 1),
        draw = 0
      )]

      observed_cols <- c("draw", "usubjid", "arrivalTime")
      if (by_treatment) {
        observed_cols <- c(
          observed_cols,
          "treatment",
          "treatment_description"
        )
      }

      subject_data <- data.table::rbindlist(
        list(
          observed_dt[, mget(observed_cols)],
          subject_data
        ),
        use.names = TRUE
      )
    }
  } else {
    subject_data <- event_pred$newEvents

    if (!is.null(dt)) {
      observed_dt <- data.table::copy(dt)
      observed_dt[, `:=`(
        arrivalTime = as.numeric(get("randdt") - get("trialsdt") + 1),
        totalTime = as.numeric(get("randdt") - get("trialsdt")) + get("time"),
        draw = 0
      )]

      observed_cols <- c(
        "draw",
        "usubjid",
        "arrivalTime",
        "time",
        "event",
        "dropout",
        "totalTime"
      )
      if (by_treatment) {
        observed_cols <- c(
          "draw",
          "usubjid",
          "arrivalTime",
          "treatment",
          "treatment_description",
          "time",
          "event",
          "dropout",
          "totalTime"
        )
      }

      subject_data <- data.table::rbindlist(
        list(
          observed_dt[get("event") | get("dropout"), mget(observed_cols)],
          subject_data
        ),
        use.names = TRUE
      )
    }
  }

  if (!is.null(dt)) {
    varnames <- c(setdiff(names(dt), names(subject_data)), "usubjid")
    subject_data <- merge(
      dt[, mget(varnames)],
      subject_data,
      by = "usubjid",
      all.y = TRUE
    )[order(get("draw"))]
  }

  subject_data
}


.ep_get_prediction_target_label <- function(to_predict_lc) {
  switch(
    to_predict_lc,
    "enrollment only" = "Enrollment only",
    "enrollment and event" = "Enrollment and event",
    "event only" = "Event only",
    stop("Unsupported to_predict value.")
  )
}


.ep_get_prediction_stage_label <- function(has_observed_data, to_predict_lc) {
  if (!has_observed_data) {
    return("Design stage")
  }

  if (to_predict_lc == "event only") {
    return("Real-time after enrollment completion")
  }

  "Real-time before enrollment completion"
}


.ep_maybe_print_prediction_plot <- function(to_predict_lc, generate_plot,
                                            showplot, enroll_pred, event_pred) {
  if (!(generate_plot && showplot)) {
    return(invisible(NULL))
  }

  if (to_predict_lc == "enrollment only") {
    print(enroll_pred$enroll_pred_plot)
  } else {
    print(event_pred$event_pred_plot)
  }

  invisible(NULL)
}


.ep_build_get_prediction_result <- function(has_observed_data, to_predict_lc,
                                            observed = NULL,
                                            enroll_fit = NULL,
                                            enroll_pred = NULL,
                                            event_fit = NULL,
                                            event_fit_posterior = NULL,
                                            event_fit_with_covariates = NULL,
                                            dropout_fit = NULL,
                                            dropout_fit_posterior = NULL,
                                            dropout_fit_with_covariates = NULL,
                                            event_pred = NULL,
                                            dropout_model = "none",
                                            covariates_event = NULL,
                                            covariates_dropout = NULL,
                                            subject_data = NULL,
                                            return_subject_data = TRUE) {
  out <- list(
    stage = .ep_get_prediction_stage_label(
      has_observed_data = has_observed_data,
      to_predict_lc = to_predict_lc
    ),
    to_predict = .ep_get_prediction_target_label(to_predict_lc)
  )

  if (!has_observed_data) {
    if (to_predict_lc == "event only") {
      stop("Internal error: design-stage event-only prediction is not supported.")
    }

    out$enroll_fit <- enroll_fit
    out$enroll_pred <- enroll_pred

    if (to_predict_lc == "enrollment and event") {
      out$event_fit <- event_fit
      if (!is.null(dropout_fit)) {
        out$dropout_fit <- dropout_fit
      }
      out$event_pred <- event_pred
    }

    return(do.call(
      .ep_build_prediction_output,
      c(
        out,
        list(
          subject_data = subject_data,
          return_subject_data = return_subject_data
        )
      )
    ))
  }

  out$observed <- observed

  if (to_predict_lc == "enrollment only") {
    out$enroll_fit <- enroll_fit
    out$enroll_pred <- enroll_pred
  } else if (to_predict_lc == "enrollment and event") {
    out$enroll_fit <- enroll_fit
    out$enroll_pred <- enroll_pred
    out$event_fit <- event_fit
    if (!is.null(event_fit_posterior)) {
      out$event_fit_posterior <- event_fit_posterior
    }

    if (!is.null(covariates_event)) {
      out$event_fit_with_covariates <- event_fit_with_covariates
    }

    if (tolower(dropout_model) != "none") {
      out$dropout_fit <- dropout_fit
      if (!is.null(dropout_fit_posterior)) {
        out$dropout_fit_posterior <- dropout_fit_posterior
      }

      if (!is.null(covariates_dropout)) {
        out$dropout_fit_with_covariates <- dropout_fit_with_covariates
      }
    }

    out$event_pred <- event_pred
  } else if (to_predict_lc == "event only") {
    if (!is.null(covariates_event)) {
      out$event_fit_with_covariates <- event_fit_with_covariates
    } else {
      out$event_fit <- event_fit
    }
    if (!is.null(event_fit_posterior)) {
      out$event_fit_posterior <- event_fit_posterior
    }

    if (tolower(dropout_model) != "none") {
      if (!is.null(covariates_dropout)) {
        out$dropout_fit_with_covariates <- dropout_fit_with_covariates
      } else {
        out$dropout_fit <- dropout_fit
      }
      if (!is.null(dropout_fit_posterior)) {
        out$dropout_fit_posterior <- dropout_fit_posterior
      }
    }

    out$event_pred <- event_pred
  } else {
    stop("Unsupported to_predict value.")
  }

  do.call(
    .ep_build_prediction_output,
    c(
      out,
      list(
        subject_data = subject_data,
        return_subject_data = return_subject_data
      )
    )
  )
}
