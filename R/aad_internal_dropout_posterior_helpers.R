.ep_format_dropout_fit_label <- function(fit, piecewise_dropout_time,
                                         k_dropout, scale_dropout,
                                         m_dropout) {
  model_lc <- tolower(fit$model)

  if (grepl("^piecewise exponential", model_lc)) {
    ps <- fit$piecewiseDropoutTime
    if (is.null(ps)) {
      ps <- piecewise_dropout_time
    }
    return(paste0(fit$model, "(", paste(ps, collapse = ","), ")"))
  }

  if (model_lc == "spline") {
    return(paste0(
      fit$model,
      "(k = ",
      k_dropout,
      ", scale = '",
      scale_dropout,
      "')"
    ))
  }

  if (model_lc == "cox") {
    return(paste0(fit$model, "(m = ", m_dropout, ")"))
  }

  fit$model
}


.ep_get_dropout_fit_time_grid <- function(dropout_fit_output) {
  if (!is.null(dropout_fit_output$dffit) &&
      "time" %in% names(dropout_fit_output$dffit)) {
    return(dropout_fit_output$dffit$time)
  }

  if (!is.null(dropout_fit_output$kmdf) &&
      "time" %in% names(dropout_fit_output$kmdf)) {
    return(sort(unique(dropout_fit_output$kmdf$time)))
  }

  0
}


.ep_compute_dropout_survival_from_fit <- function(fit, time_grid,
                                                  piecewise_dropout_time = 0) {
  model_lc <- tolower(fit$model)

  surv <- switch(
    model_lc,
    "exponential" = {
      stats::pexp(time_grid, exp(fit$theta[1]), lower.tail = FALSE)
    },
    "weibull" = {
      stats::pweibull(
        time_grid,
        shape = exp(-fit$theta[2]),
        scale = exp(fit$theta[1]),
        lower.tail = FALSE
      )
    },
    "log-logistic" = {
      stats::plogis(
        log(time_grid),
        location = fit$theta[1],
        scale = exp(fit$theta[2]),
        lower.tail = FALSE
      )
    },
    "log-normal" = {
      stats::plnorm(
        time_grid,
        meanlog = fit$theta[1],
        sdlog = exp(fit$theta[2]),
        lower.tail = FALSE
      )
    },
    "piecewise exponential" = {
      ps <- fit$piecewiseDropoutTime
      if (is.null(ps)) {
        ps <- piecewise_dropout_time
      }
      ppwexp(time_grid, fit$theta, length(ps), ps, lower.tail = FALSE)
    },
    NULL
  )

  if (is.null(surv)) {
    return(NULL)
  }

  data.table::data.table(time = time_grid, surv = surv)
}


.ep_build_dropout_fit_plot <- function(kmdf, dffit, text_vec,
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
        plotly::add_lines(
          data = dffit,
          x = ~time,
          y = ~surv,
          name = "fitted"
        ) %>%
        plotly::layout(
          xaxis = list(title = "Days since randomization", zeroline = FALSE),
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
            showarrow = FALSE
          )
        ) %>%
        plotly::hide_legend()
    )
  }

  fitted_dropout <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = kmdf[-1, ],
      ggplot2::aes(x = .data$time, y = .data$surv)
    ) +
    ggplot2::geom_line(
      data = dffit[-1, ],
      ggplot2::aes(x = .data$time, y = .data$surv),
      colour = "red"
    ) +
    ggplot2::labs(
      x = "Days since randomization",
      y = "Survival probability",
      title = "Fitted time to dropout survival curve"
    ) +
    ggplot2::theme(legend.position = "none")

  x_pos <- max(kmdf$time, na.rm = TRUE)
  y_pos <- max(kmdf$surv, na.rm = TRUE)

  for (j in seq_along(text_vec)) {
    fitted_dropout <- fitted_dropout +
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

  fitted_dropout
}


.ep_add_dropout_fit_treatment_label <- function(plot_obj, label,
                                                interactive_plot) {
  if (interactive_plot) {
    return(
      plot_obj %>%
        plotly::layout(
          annotations = list(
            x = 0.5,
            y = 1,
            text = paste0("<b>", label, "</b>"),
            xanchor = "center",
            yanchor = "middle",
            showarrow = FALSE,
            xref = "paper",
            yref = "paper"
          )
        )
    )
  }

  plot_obj + ggplot2::labs(subtitle = label)
}


.ep_build_dropout_fit_output <- function(fit, kmdf, dffit, text_vec,
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


.ep_posterior_dropout_fit_text <- function(fit, piecewise_dropout_time,
                                           k_dropout, scale_dropout,
                                           m_dropout) {
  c(
    .ep_format_dropout_fit_label(
      fit = fit,
      piecewise_dropout_time = piecewise_dropout_time,
      k_dropout = k_dropout,
      scale_dropout = scale_dropout,
      m_dropout = m_dropout
    ),
    "Posterior mean fit: prior + likelihood"
  )
}


.ep_build_single_posterior_dropout_fit <- function(dropout_fit_output,
                                                   posterior_fit,
                                                   piecewise_dropout_time,
                                                   k_dropout,
                                                   scale_dropout,
                                                   m_dropout,
                                                   generate_plot,
                                                   interactive_plot) {
  if (is.null(dropout_fit_output) || is.null(posterior_fit)) {
    return(NULL)
  }

  time_grid <- .ep_get_dropout_fit_time_grid(dropout_fit_output)
  dffit <- .ep_compute_dropout_survival_from_fit(
    fit = posterior_fit,
    time_grid = time_grid,
    piecewise_dropout_time = piecewise_dropout_time
  )

  if (is.null(dffit)) {
    return(NULL)
  }

  fit_for_output <- posterior_fit
  fit_for_output$model <- dropout_fit_output$fit$model

  if (!is.null(dropout_fit_output$fit$treatment)) {
    fit_for_output$treatment <- dropout_fit_output$fit$treatment
    fit_for_output$treatment_description <-
      dropout_fit_output$fit$treatment_description
  }

  kmdf <- dropout_fit_output$kmdf
  text_vec <- .ep_posterior_dropout_fit_text(
    fit = fit_for_output,
    piecewise_dropout_time = piecewise_dropout_time,
    k_dropout = k_dropout,
    scale_dropout = scale_dropout,
    m_dropout = m_dropout
  )

  fit_plot <- NULL
  if (generate_plot) {
    fit_plot <- .ep_build_dropout_fit_plot(
      kmdf = kmdf,
      dffit = dffit,
      text_vec = text_vec,
      interactive_plot = interactive_plot
    )

    if (!is.null(fit_for_output$treatment_description)) {
      fit_plot <- .ep_add_dropout_fit_treatment_label(
        plot_obj = fit_plot,
        label = fit_for_output$treatment_description,
        interactive_plot = interactive_plot
      )
    }
  }

  .ep_build_dropout_fit_output(
    fit = fit_for_output,
    fit_plot = fit_plot,
    kmdf = kmdf,
    dffit = dffit,
    text_vec = text_vec
  )
}


.ep_build_posterior_dropout_fit_output <- function(dropout_fit,
                                                   posterior_fit,
                                                   piecewise_dropout_time,
                                                   k_dropout,
                                                   scale_dropout,
                                                   m_dropout,
                                                   generate_plot,
                                                   interactive_plot) {
  if (is.null(dropout_fit) || is.null(posterior_fit)) {
    return(NULL)
  }

  if (!.ep_is_event_fit_output_collection(dropout_fit)) {
    return(
      .ep_build_single_posterior_dropout_fit(
        dropout_fit_output = dropout_fit,
        posterior_fit = posterior_fit,
        piecewise_dropout_time = piecewise_dropout_time,
        k_dropout = k_dropout,
        scale_dropout = scale_dropout,
        m_dropout = m_dropout,
        generate_plot = generate_plot,
        interactive_plot = interactive_plot
      )
    )
  }

  posterior_dropout_fit <- purrr::map2(
    dropout_fit,
    posterior_fit,
    ~ .ep_build_single_posterior_dropout_fit(
      dropout_fit_output = .x,
      posterior_fit = .y,
      piecewise_dropout_time = piecewise_dropout_time,
      k_dropout = k_dropout,
      scale_dropout = scale_dropout,
      m_dropout = m_dropout,
      generate_plot = generate_plot,
      interactive_plot = interactive_plot
    )
  )

  if (generate_plot && length(posterior_dropout_fit) > 1) {
    x_range <- range(unlist(lapply(posterior_dropout_fit, function(x) {
      x$dffit$time
    })), na.rm = TRUE)

    posterior_dropout_fit <- .ep_sync_event_fit_plot_ranges(
      event_fit = posterior_dropout_fit,
      x_range = x_range,
      interactive_plot = interactive_plot
    )
  }

  posterior_dropout_fit
}
