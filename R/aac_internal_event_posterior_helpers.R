.ep_is_event_fit_output_collection <- function(event_fit) {
  is.list(event_fit) &&
    length(event_fit) > 0 &&
    is.list(event_fit[[1]]) &&
    "fit" %in% names(event_fit[[1]])
}


.ep_get_event_fit_time_grid <- function(event_fit_output) {
  if (!is.null(event_fit_output$dffit) &&
      "time" %in% names(event_fit_output$dffit)) {
    return(event_fit_output$dffit$time)
  }

  if (!is.null(event_fit_output$kmdf) &&
      "time" %in% names(event_fit_output$kmdf)) {
    return(sort(unique(event_fit_output$kmdf$time)))
  }

  0
}


.ep_compute_event_survival_from_fit <- function(fit, time_grid,
                                                piecewise_survival_time = 0) {
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
      ps <- fit$piecewiseSurvivalTime
      if (is.null(ps)) {
        ps <- piecewise_survival_time
      }
      ppwexp(time_grid, fit$theta, length(ps), ps, lower.tail = FALSE)
    },
    "exponential with cured population" = {
      SP_Chen_exponential(fit$theta, time_grid)
    },
    "weibull with cured population" = {
      SP_Chen_weibull(fit$theta, time_grid)
    },
    "log-logistic with cured population" = {
      SP_Chen_log_logistic(fit$theta, time_grid)
    },
    "log-normal with cured population" = {
      SP_Chen_log_normal(fit$theta, time_grid)
    },
    "piecewise exponential with cured population" = {
      ps <- fit$piecewiseSurvivalTime
      if (is.null(ps)) {
        ps <- piecewise_survival_time
      }
      SP_Chen_piecewise_exponential(fit$theta, time_grid, ps)
    },
    NULL
  )

  if (is.null(surv)) {
    return(NULL)
  }

  data.table::data.table(time = time_grid, surv = surv)
}


.ep_posterior_event_fit_text <- function(fit, piecewise_survival_time, k,
                                         scale, m) {
  text_vec <- c(
    .ep_format_event_fit_label(
      fit = fit,
      piecewise_survival_time = piecewise_survival_time,
      k = k,
      scale = scale,
      m = m
    ),
    "Posterior mean fit: prior + likelihood"
  )

  if (grepl("with cured population", tolower(fit$model), fixed = TRUE)) {
    text_vec <- c(
      text_vec,
      paste(
        "Prop of Cured:",
        formatC(stats::plogis(fit$theta[1]), format = "f", digits = 2)
      )
    )
  }

  text_vec
}


.ep_build_single_posterior_event_fit <- function(event_fit_output,
                                                 posterior_fit,
                                                 piecewise_survival_time,
                                                 k, scale, m, generate_plot,
                                                 interactive_plot) {
  if (is.null(event_fit_output) || is.null(posterior_fit)) {
    return(NULL)
  }

  time_grid <- .ep_get_event_fit_time_grid(event_fit_output)
  dffit <- .ep_compute_event_survival_from_fit(
    fit = posterior_fit,
    time_grid = time_grid,
    piecewise_survival_time = piecewise_survival_time
  )

  if (is.null(dffit)) {
    return(NULL)
  }

  fit_for_output <- posterior_fit
  fit_for_output$model <- event_fit_output$fit$model

  if (!is.null(event_fit_output$fit$treatment)) {
    fit_for_output$treatment <- event_fit_output$fit$treatment
    fit_for_output$treatment_description <-
      event_fit_output$fit$treatment_description
  }

  kmdf <- event_fit_output$kmdf
  kmdf_censor <- kmdf[get("n.censor") > 0]
  text_vec <- .ep_posterior_event_fit_text(
    fit = fit_for_output,
    piecewise_survival_time = piecewise_survival_time,
    k = k,
    scale = scale,
    m = m
  )

  fit_plot <- NULL
  if (generate_plot) {
    fit_plot <- .ep_build_event_fit_plot(
      kmdf = kmdf,
      kmdf_censor = kmdf_censor,
      dffit = dffit,
      text_vec = text_vec,
      interactive_plot = interactive_plot
    )

    if (!is.null(fit_for_output$treatment_description)) {
      fit_plot <- .ep_add_event_fit_treatment_label(
        plot_obj = fit_plot,
        label = fit_for_output$treatment_description,
        interactive_plot = interactive_plot
      )
    }
  }

  .ep_build_event_fit_output(
    fit = fit_for_output,
    fit_plot = fit_plot,
    kmdf = kmdf,
    dffit = dffit,
    text_vec = text_vec
  )
}


.ep_build_posterior_event_fit_output <- function(event_fit, posterior_fit,
                                                 piecewise_survival_time, k,
                                                 scale, m, generate_plot,
                                                 interactive_plot) {
  if (is.null(event_fit) || is.null(posterior_fit)) {
    return(NULL)
  }

  if (!.ep_is_event_fit_output_collection(event_fit)) {
    return(
      .ep_build_single_posterior_event_fit(
        event_fit_output = event_fit,
        posterior_fit = posterior_fit,
        piecewise_survival_time = piecewise_survival_time,
        k = k,
        scale = scale,
        m = m,
        generate_plot = generate_plot,
        interactive_plot = interactive_plot
      )
    )
  }

  posterior_event_fit <- purrr::map2(
    event_fit,
    posterior_fit,
    ~ .ep_build_single_posterior_event_fit(
      event_fit_output = .x,
      posterior_fit = .y,
      piecewise_survival_time = piecewise_survival_time,
      k = k,
      scale = scale,
      m = m,
      generate_plot = generate_plot,
      interactive_plot = interactive_plot
    )
  )

  if (generate_plot && length(posterior_event_fit) > 1) {
    x_range <- range(unlist(lapply(posterior_event_fit, function(x) {
      x$dffit$time
    })), na.rm = TRUE)

    posterior_event_fit <- .ep_sync_event_fit_plot_ranges(
      event_fit = posterior_event_fit,
      x_range = x_range,
      interactive_plot = interactive_plot
    )
  }

  posterior_event_fit
}
