.ep_fit_survreg_event_model <- function(df1, formula, event_model_lc,
                                        n0, d0, q, x) {
  time_grid <- seq(0, max(df1$time))

  if (event_model_lc == "exponential") {
    erify::check_positive(d0 - q, supplement = paste(
      "The number of events must be >=", q + 1, "to fit an exponential model."))

    reg <- survival::survreg(formula, data = df1, dist = "exponential")

    fit2 <- list(
      model = "Exponential",
      theta = -as.numeric(reg$coefficients),
      vtheta = reg$var,
      aic = -2 * reg$loglik[1] + 2 * (q + 1),
      bic = -2 * reg$loglik[1] + (q + 1) * log(n0)
    )

    rate <- exp(as.numeric(x %*% fit2$theta))
    surv <- sapply(time_grid, function(t) {
      mean(stats::pexp(t, rate, lower.tail = FALSE))
    })

  } else if (event_model_lc == "weibull") {
    erify::check_positive(d0 - q - 1, supplement = paste(
      "The number of events must be >=", q + 2, "to fit a Weibull model."))

    reg <- survival::survreg(formula, data = df1, dist = "weibull")

    fit2 <- list(
      model = "Weibull",
      theta = c(as.numeric(reg$coefficients), log(reg$scale)),
      vtheta = reg$var,
      aic = -2 * reg$loglik[1] + 2 * (q + 2),
      bic = -2 * reg$loglik[1] + (q + 2) * log(n0)
    )

    shape <- exp(-fit2$theta[q + 2])
    scale2 <- exp(as.numeric(x %*% fit2$theta[1:(q + 1)]))
    surv <- sapply(time_grid, function(t) {
      mean(stats::pweibull(t, shape, scale2, lower.tail = FALSE))
    })

  } else if (event_model_lc == "log-logistic") {
    erify::check_positive(d0 - q - 1, supplement = paste(
      "The number of events must be >=", q + 2, "to fit a log-logistic model."))

    reg <- survival::survreg(formula, data = df1, dist = "loglogistic")

    fit2 <- list(
      model = "Log-logistic",
      theta = c(as.numeric(reg$coefficients), log(reg$scale)),
      vtheta = reg$var,
      aic = -2 * reg$loglik[1] + 2 * (q + 2),
      bic = -2 * reg$loglik[1] + (q + 2) * log(n0)
    )

    location <- as.numeric(x %*% fit2$theta[1:(q + 1)])
    scale2 <- exp(fit2$theta[q + 2])
    surv <- sapply(time_grid, function(t) {
      mean(stats::plogis(log(t), location, scale2, lower.tail = FALSE))
    })

  } else if (event_model_lc == "log-normal") {
    erify::check_positive(d0 - q - 1, supplement = paste(
      "The number of events must be >=", q + 2, "to fit a log-normal model."))

    reg <- survival::survreg(formula, data = df1, dist = "lognormal")

    fit2 <- list(
      model = "Log-normal",
      theta = c(as.numeric(reg$coefficients), log(reg$scale)),
      vtheta = reg$var,
      aic = -2 * reg$loglik[1] + 2 * (q + 2),
      bic = -2 * reg$loglik[1] + (q + 2) * log(n0)
    )

    meanlog <- as.numeric(x %*% fit2$theta[1:(q + 1)])
    sdlog <- exp(fit2$theta[q + 2])
    surv <- sapply(time_grid, function(t) {
      mean(stats::plnorm(t, meanlog, sdlog, lower.tail = FALSE))
    })

  } else {
    stop("Unsupported survreg event model.")
  }

  list(
    fit = fit2,
    dffit = data.table::data.table(time = time_grid, surv = surv),
    p_cure = NULL
  )
}


.ep_fit_piecewise_event_model <- function(df1, n0, d0, q, x,
                                          piecewise_survival_time) {
  j <- length(piecewise_survival_time)

  erify::check_positive(d0 - j - q + 1, supplement = paste(
    "The number of events must be >=", j + q,
    "to fit a piecewise exponential model."))

  fit2 <- pwexpreg(
    df1$time,
    df1$event,
    j,
    piecewise_survival_time,
    q,
    x
  )

  time_grid <- seq(0, max(df1$time))
  surv <- purrr::map(seq_len(n0), function(idx) {
    ppwexp(
      time_grid,
      fit2$theta,
      j,
      fit2$piecewiseSurvivalTime,
      q,
      x[idx, ],
      lower.tail = FALSE
    )
  })
  surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)

  list(
    fit = fit2,
    dffit = data.table::data.table(time = time_grid, surv = surv),
    p_cure = NULL
  )
}


.ep_fit_model_averaging_event_model <- function(df1, formula, n0, d0, q, x) {
  erify::check_positive(d0 - q - 1, supplement = paste(
    "The number of events must be >=", q + 2,
    "to fit a model averaging model."))

  reg1 <- survival::survreg(formula, data = df1, dist = "weibull")
  reg2 <- survival::survreg(formula, data = df1, dist = "lognormal")
  aic1 <- -2 * reg1$loglik[1] + 2 * (q + 2)
  aic2 <- -2 * reg2$loglik[1] + 2 * (q + 2)
  bic1 <- -2 * reg1$loglik[1] + (q + 2) * log(n0)
  bic2 <- -2 * reg2$loglik[1] + (q + 2) * log(n0)

  w1 <- 1 / (1 + exp(-0.5 * (bic2 - bic1)))

  theta <- c(
    as.numeric(reg1$coefficients), log(reg1$scale),
    as.numeric(reg2$coefficients), log(reg2$scale)
  )
  vtheta <- as.matrix(Matrix::bdiag(reg1$var, reg2$var))

  fit2 <- list(
    model = "Model averaging",
    theta = theta,
    vtheta = vtheta,
    aic = w1 * aic1 + (1 - w1) * aic2,
    bic = w1 * bic1 + (1 - w1) * bic2,
    w1 = w1
  )

  time_grid <- seq(0, max(df1$time))
  surv <- purrr::map(seq_len(n0), function(idx) {
    pmodavg(time_grid, fit2$theta, w1, q, x[idx, ], lower.tail = FALSE)
  })
  surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)

  list(
    fit = fit2,
    dffit = data.table::data.table(time = time_grid, surv = surv),
    p_cure = NULL
  )
}


.ep_fit_spline_event_model <- function(df1, formula, n0, d0, q, x, k, scale) {
  erify::check_positive(d0 - k - q - 1, supplement = paste(
    "The number of events must be >=", k + q + 2,
    "to fit a spline model."))

  spl <- flexsurv::flexsurvspline(
    formula,
    data = df1,
    k = k,
    scale = scale,
    method = "Nelder-Mead"
  )

  fit2 <- list(
    model = "Spline",
    theta = as.numeric(spl$coefficients),
    vtheta = spl$cov,
    aic = -2 * spl$loglik + 2 * (k + q + 2),
    bic = -2 * spl$loglik + (k + q + 2) * log(n0),
    knots = spl$knots,
    scale = spl$scale
  )

  time_grid <- seq(0, max(df1$time))

  if (q > 0) {
    xbeta <- as.numeric(as.matrix(x[, -1]) %*% fit2$theta[(k + 3):(k + q + 2)])

    surv <- purrr::map(seq_len(n0), function(idx) {
      flexsurv::psurvspline(
        time_grid,
        gamma = fit2$theta[1:(k + 2)],
        knots = fit2$knots,
        scale = fit2$scale,
        offset = xbeta[idx],
        lower.tail = FALSE
      )
    })
    surv <- apply(matrix(purrr::list_c(surv), ncol = n0), 1, mean)
  } else {
    surv <- flexsurv::psurvspline(
      time_grid,
      gamma = fit2$theta,
      knots = fit2$knots,
      scale = fit2$scale,
      lower.tail = FALSE
    )
  }

  list(
    fit = fit2,
    dffit = data.table::data.table(time = time_grid, surv = surv),
    p_cure = NULL
  )
}


.ep_fit_cox_event_model <- function(df1, n0, d0, q, x, m, covariates) {
  erify::check_positive(d0 - q, supplement = paste(
    "The number of events must be >=", q + 1, "to fit a Cox model."))
  erify::check_positive(d0 - m + 1, supplement = paste(
    "m must be <= the observed number of events", d0))

  reg <- phregr(df1, time = "time", event = "event", covariates = covariates)

  bh <- data.table::setDT(reg$basehaz)[get("nevent") > 0]
  haz <- bh$haz
  vhaz <- bh$varhaz

  if (q > 0) {
    if (q == 1) {
      ghaz <- matrix(bh$gradhaz, ncol = 1)
    } else {
      ghaz <- do.call(cbind, lapply(seq_len(q), function(ii) {
        bh[[paste0("gradhaz.", ii)]]
      }))
    }
  }

  m_events <- nrow(bh)
  tcut <- c(0, bh$time)
  dtcut <- diff(tcut)
  lambda1 <- haz / dtcut

  d <- bh$nevent
  llik <- reg$sumstat$loglik1 + sum(d * (log(d / dtcut) - 1))

  if (q > 0) {
    theta <- c(log(lambda1), as.numeric(reg$beta))
    vbeta <- reg$vbeta
    dimnames(vbeta) <- NULL
    vtheta <- matrix(0, m_events + q, m_events + q)
    vtheta[(m_events + 1):(m_events + q), (m_events + 1):(m_events + q)] <- vbeta
    vtheta[1:m_events, 1:m_events] <-
      diag(vhaz / (haz * haz)) + ghaz %*% vbeta %*% t(ghaz)
    vtheta[1:m_events, (m_events + 1):(m_events + q)] <- ghaz %*% vbeta
    vtheta[(m_events + 1):(m_events + q), 1:m_events] <-
      t(vtheta[1:m_events, (m_events + 1):(m_events + q)])
  } else {
    theta <- log(lambda1)
    vtheta <- diag(vhaz / (haz * haz))
  }

  fit2 <- list(
    model = "Cox",
    theta = theta,
    vtheta = vtheta,
    aic = -2 * llik + 2 * (m_events + q),
    bic = -2 * llik + (m_events + q) * log(n0),
    piecewiseSurvivalTime = tcut
  )

  time_grid <- seq(0, max(df1$time))
  lambda2 <- sum(bh$haz[(m_events - m + 1):m_events]) /
    (bh$time[m_events] - bh$time[m_events - m])
  lambda <- c(lambda1, lambda2)

  s1 <- sapply(time_grid, function(t) {
    ppwexp(t, log(lambda), m_events + 1, tcut, lower.tail = FALSE)
  })

  if (q > 0) {
    xbeta <- as.numeric(as.matrix(x[, -1]) %*% reg$beta)
    surv <- apply(outer(s1, exp(xbeta), `^`), 1, mean)
  } else {
    surv <- s1
  }

  list(
    fit = fit2,
    dffit = data.table::data.table(time = time_grid, surv = surv),
    p_cure = NULL
  )
}


.ep_fit_exponential_cure_event_model <- function(df1, n0, d0, ex0, cure_rate) {
  df1_df <- as.data.frame(df1)

  if (is.null(cure_rate)) {
    temp <- stats::nlminb(
      start = c(-2, log(d0 / ex0)),
      objective = loglik_Chen_exponential,
      df = df1_df
    )
    theta <- temp$par

    fit2 <- list(
      model = "exponential with cured population",
      theta = theta,
      vtheta = MASS::ginv(numDeriv::hessian(
        loglik_Chen_exponential,
        theta,
        df = df1_df
      )),
      bic = 2 * temp$objective + log(n0) * 2,
      aic = 2 * temp$objective + 2 * 2
    )

    p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
  } else {
    temp <- stats::nlminb(
      start = c(log(d0 / ex0)),
      objective = loglik_Chen_exponential2,
      df = df1_df,
      p = cure_rate
    )

    theta <- c(stats::qlogis(cure_rate), temp$par)
    vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
    vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
      loglik_Chen_exponential2,
      temp$par,
      df = df1_df,
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

  time_grid <- seq(0, max(df1$time))

  list(
    fit = fit2,
    dffit = data.table::data.table(
      time = time_grid,
      surv = SP_Chen_exponential(theta, time_grid)
    ),
    p_cure = p_cure
  )
}


.ep_fit_weibull_cure_event_model <- function(df1, n0, cure_rate) {
  df1_df <- as.data.frame(df1)
  reg0 <- survival::survreg(survival::Surv(time, event) ~ 1, data = df1, dist = "weibull")

  if (is.null(cure_rate)) {
    temp <- stats::nlminb(
      start = c(-2, log(1 / reg0$scale), as.numeric(reg0$coeff)),
      objective = loglik_Chen_weibull,
      df = df1_df
    )
    theta <- temp$par

    fit2 <- list(
      model = "weibull with cured population",
      theta = theta,
      vtheta = MASS::ginv(numDeriv::hessian(
        loglik_Chen_weibull,
        theta,
        df = df1_df
      )),
      bic = 2 * temp$objective + log(n0) * 3,
      aic = 2 * temp$objective + 2 * 3
    )

    p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
  } else {
    temp <- stats::nlminb(
      start = c(log(1 / reg0$scale), as.numeric(reg0$coeff)),
      objective = loglik_Chen_weibull2,
      df = df1_df,
      p = cure_rate
    )

    theta <- c(stats::qlogis(cure_rate), temp$par)
    vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
    vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
      loglik_Chen_weibull2,
      temp$par,
      df = df1_df,
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

  time_grid <- seq(0, max(df1$time))

  list(
    fit = fit2,
    dffit = data.table::data.table(
      time = time_grid,
      surv = SP_Chen_weibull(theta, time_grid)
    ),
    p_cure = p_cure
  )
}


.ep_fit_log_logistic_cure_event_model <- function(df1, n0, cure_rate) {
  df1_df <- as.data.frame(df1)
  reg0 <- survival::survreg(
    survival::Surv(time, event) ~ 1,
    data = df1,
    dist = "loglogistic"
  )

  if (is.null(cure_rate)) {
    temp <- stats::nlminb(
      start = c(-2, log(1 / reg0$scale), as.numeric(reg0$coeff)),
      objective = loglik_Chen_log_logistic,
      df = df1_df
    )
    theta <- temp$par

    fit2 <- list(
      model = "log-logistic with cured population",
      theta = theta,
      vtheta = MASS::ginv(numDeriv::hessian(
        loglik_Chen_log_logistic,
        theta,
        df = df1_df
      )),
      bic = 2 * temp$objective + log(n0) * 3,
      aic = 2 * temp$objective + 2 * 3
    )

    p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
  } else {
    temp <- stats::nlminb(
      start = c(log(1 / reg0$scale), as.numeric(reg0$coeff)),
      objective = loglik_Chen_log_logistic2,
      df = df1_df,
      p = cure_rate
    )

    theta <- c(stats::qlogis(cure_rate), temp$par)
    vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
    vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
      loglik_Chen_log_logistic2,
      temp$par,
      df = df1_df,
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

  time_grid <- seq(0, max(df1$time))

  list(
    fit = fit2,
    dffit = data.table::data.table(
      time = time_grid,
      surv = SP_Chen_log_logistic(theta, time_grid)
    ),
    p_cure = p_cure
  )
}


.ep_fit_log_normal_cure_event_model <- function(df1, n0, cure_rate) {
  df1_df <- as.data.frame(df1)
  reg0 <- survival::survreg(survival::Surv(time, event) ~ 1, data = df1, dist = "lognormal")

  if (is.null(cure_rate)) {
    temp <- stats::nlminb(
      start = c(-2, as.numeric(reg0$coeff), log(reg0$scale)),
      objective = loglik_Chen_log_normal,
      df = df1_df
    )
    theta <- temp$par

    fit2 <- list(
      model = "log-normal with cured population",
      theta = theta,
      vtheta = MASS::ginv(numDeriv::hessian(
        loglik_Chen_log_normal,
        theta,
        df = df1_df
      )),
      bic = 2 * temp$objective + log(n0) * 3,
      aic = 2 * temp$objective + 2 * 3
    )

    p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
  } else {
    temp <- stats::nlminb(
      start = c(as.numeric(reg0$coeff), log(reg0$scale)),
      objective = loglik_Chen_log_normal2,
      df = df1_df,
      p = cure_rate
    )

    theta <- c(stats::qlogis(cure_rate), temp$par)
    vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
    vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
      loglik_Chen_log_normal2,
      temp$par,
      df = df1_df,
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

  time_grid <- seq(0, max(df1$time))

  list(
    fit = fit2,
    dffit = data.table::data.table(
      time = time_grid,
      surv = SP_Chen_log_normal(theta, time_grid)
    ),
    p_cure = p_cure
  )
}


.ep_fit_piecewise_exponential_cure_event_model <- function(df1, n0,
                                                           piecewise_survival_time,
                                                           cure_rate) {
  df1_df <- as.data.frame(df1)
  u <- piecewise_survival_time[piecewise_survival_time < max(df1$time)]
  j <- length(u)
  haz_start <- rep(log(0.0030), j)

  if (is.null(cure_rate)) {
    temp <- stats::nlminb(
      start = c(-2, haz_start),
      objective = loglik_Chen_piecewise_exponential,
      df = df1_df,
      piecewiseSurvivalTime = piecewise_survival_time
    )
    theta <- temp$par

    fit2 <- list(
      model = "piecewise exponential with cured population",
      theta = theta,
      vtheta = MASS::ginv(numDeriv::hessian(
        loglik_Chen_piecewise_exponential,
        theta,
        df = df1_df,
        piecewiseSurvivalTime = piecewise_survival_time
      )),
      bic = 2 * temp$objective + (j + 1) * log(n0),
      aic = 2 * temp$objective + (j + 1) * 2,
      piecewiseSurvivalTime = u
    )

    p_cure <- exp(theta[1]) / (1 + exp(theta[1]))
  } else {
    temp <- stats::nlminb(
      start = haz_start,
      objective = loglik_Chen_piecewise_exponential2,
      df = df1_df,
      piecewiseSurvivalTime = piecewise_survival_time,
      p = cure_rate
    )

    theta <- c(stats::qlogis(cure_rate), temp$par)
    vtheta <- matrix(0, nrow = length(theta), ncol = length(theta))
    vtheta[-1, -1] <- MASS::ginv(numDeriv::hessian(
      loglik_Chen_piecewise_exponential2,
      temp$par,
      df = df1_df,
      piecewiseSurvivalTime = piecewise_survival_time,
      p = cure_rate
    ))

    fit2 <- list(
      model = "piecewise exponential with cured population",
      theta = theta,
      vtheta = vtheta,
      bic = 2 * temp$objective + j * log(n0),
      aic = 2 * temp$objective + j * 2,
      piecewiseSurvivalTime = u
    )

    p_cure <- cure_rate
  }

  time_grid <- seq(0, max(df1$time))

  list(
    fit = fit2,
    dffit = data.table::data.table(
      time = time_grid,
      surv = SP_Chen_piecewise_exponential(
        theta,
        time_grid,
        piecewiseSurvivalTime = piecewise_survival_time
      )
    ),
    p_cure = p_cure
  )
}


.ep_fit_event_model_branch <- function(df1, formula, event_model_lc,
                                       n0, d0, ex0, q, x,
                                       piecewise_survival_time,
                                       k, scale, m,
                                       covariates, cure_rate) {
  if (event_model_lc %in% c("exponential", "weibull", "log-logistic", "log-normal")) {
    return(.ep_fit_survreg_event_model(
      df1 = df1,
      formula = formula,
      event_model_lc = event_model_lc,
      n0 = n0,
      d0 = d0,
      q = q,
      x = x
    ))
  }

  switch(
    event_model_lc,
    "piecewise exponential" = .ep_fit_piecewise_event_model(
      df1 = df1,
      n0 = n0,
      d0 = d0,
      q = q,
      x = x,
      piecewise_survival_time = piecewise_survival_time
    ),
    "model averaging" = .ep_fit_model_averaging_event_model(
      df1 = df1,
      formula = formula,
      n0 = n0,
      d0 = d0,
      q = q,
      x = x
    ),
    "spline" = .ep_fit_spline_event_model(
      df1 = df1,
      formula = formula,
      n0 = n0,
      d0 = d0,
      q = q,
      x = x,
      k = k,
      scale = scale
    ),
    "cox" = .ep_fit_cox_event_model(
      df1 = df1,
      n0 = n0,
      d0 = d0,
      q = q,
      x = x,
      m = m,
      covariates = covariates
    ),
    "exponential with cured population" = .ep_fit_exponential_cure_event_model(
      df1 = df1,
      n0 = n0,
      d0 = d0,
      ex0 = ex0,
      cure_rate = cure_rate
    ),
    "weibull with cured population" = .ep_fit_weibull_cure_event_model(
      df1 = df1,
      n0 = n0,
      cure_rate = cure_rate
    ),
    "log-logistic with cured population" = .ep_fit_log_logistic_cure_event_model(
      df1 = df1,
      n0 = n0,
      cure_rate = cure_rate
    ),
    "log-normal with cured population" = .ep_fit_log_normal_cure_event_model(
      df1 = df1,
      n0 = n0,
      cure_rate = cure_rate
    ),
    "piecewise exponential with cured population" = .ep_fit_piecewise_exponential_cure_event_model(
      df1 = df1,
      n0 = n0,
      piecewise_survival_time = piecewise_survival_time,
      cure_rate = cure_rate
    ),
    stop("Unsupported event_model.")
  )
}
