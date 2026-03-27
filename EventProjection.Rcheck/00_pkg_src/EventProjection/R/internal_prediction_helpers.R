.ep_split_by_draw_sorted <- function(dt, value_col, nreps) {
  if (nreps < 1) {
    stop("nreps must be at least 1.")
  }

  if (is.null(dt) || nrow(dt) == 0) {
    return(rep(list(numeric(0)), nreps))
  }

  split_vals <- split(dt[[value_col]], factor(dt$draw, levels = seq_len(nreps)))

  lapply(split_vals, function(x) {
    if (length(x) == 0) {
      numeric(0)
    } else {
      sort.int(as.numeric(x), method = "auto")
    }
  })
}


.ep_count_matrix_from_splits <- function(times, value_splits) {
  times <- as.numeric(times)
  count_mat <- matrix(0L, nrow = length(times), ncol = length(value_splits))

  for (i in seq_along(value_splits)) {
    vals <- value_splits[[i]]
    if (length(vals) > 0) {
      count_mat[, i] <- findInterval(times, vals)
    }
  }

  count_mat
}


.ep_count_matrix_from_values <- function(dt, value_col, times, nreps) {
  .ep_count_matrix_from_splits(
    times = times,
    value_splits = .ep_split_by_draw_sorted(dt, value_col, nreps)
  )
}


.ep_summarize_count_matrix <- function(count_mat, times, pilevel) {
  plower <- (1 - pilevel) / 2
  pupper <- 1 - plower

  qmat <- t(apply(
    count_mat,
    1,
    stats::quantile,
    probs = c(0.5, plower, pupper),
    names = FALSE
  ))

  if (ncol(count_mat) > 1) {
    var_vec <- apply(count_mat, 1, stats::var)
  } else {
    var_vec <- rep(0, nrow(count_mat))
  }

  data.table::data.table(
    t = as.numeric(times),
    n = qmat[, 1],
    pilevel = pilevel,
    lower = qmat[, 2],
    upper = qmat[, 3],
    mean = rowMeans(count_mat),
    var = var_vec
  )
}


.ep_add_summary_offset <- function(summary_dt, offset) {
  if (offset != 0) {
    summary_dt[, `:=`(
      n = get("n") + offset,
      lower = get("lower") + offset,
      upper = get("upper") + offset,
      mean = get("mean") + offset
    )]
  }

  summary_dt
}


.ep_observed_summary_from_counts <- function(times, counts, pilevel) {
  counts <- as.numeric(counts)

  data.table::data.table(
    t = as.numeric(times),
    n = counts,
    pilevel = pilevel,
    lower = NA_real_,
    upper = NA_real_,
    mean = counts,
    var = 0
  )
}


.ep_observed_ongoing_summary <- function(times, arrival_values,
                                         completed_values, pilevel) {
  arrival_values <- sort.int(as.numeric(arrival_values), method = "auto")
  completed_values <- sort.int(as.numeric(completed_values), method = "auto")

  arrival_counts <- findInterval(as.numeric(times), arrival_values)
  completed_counts <- findInterval(as.numeric(times), completed_values)

  .ep_observed_summary_from_counts(
    times = times,
    counts = arrival_counts - completed_counts,
    pilevel = pilevel
  )
}


.ep_observed_grouped_ongoing_summary <- function(groups_dt, data_dt,
                                                 group_cols, times,
                                                 pilevel) {
  out_list <- lapply(seq_len(nrow(groups_dt)), function(i) {
    group_row <- groups_dt[i]
    dt_i <- data_dt

    for (col in group_cols) {
      dt_i <- dt_i[get(col) == group_row[[col]][1]]
    }

    summary_i <- .ep_observed_ongoing_summary(
      times = times,
      arrival_values = dt_i$arrivalTime,
      completed_values = dt_i[get("event") | get("dropout")]$totalTime,
      pilevel = pilevel
    )

    for (col in group_cols) {
      summary_i[, (col) := group_row[[col]][1]]
    }

    summary_i
  })

  data.table::rbindlist(out_list, use.names = TRUE)
}


.ep_append_summary_timepoint <- function(summary_dt, tp, group_cols = NULL) {
  if (is.null(group_cols)) {
    return(data.table::rbindlist(
      list(summary_dt, data.table::copy(summary_dt)[.N][, `:=`(t = tp)]),
      use.names = TRUE
    )[, .SD[.N], by = "t"])
  }

  data.table::rbindlist(
    list(
      summary_dt,
      data.table::copy(summary_dt)[, .SD[.N], by = group_cols][, `:=`(t = tp)]
    ),
    use.names = TRUE
  )[, .SD[.N], by = c(group_cols, "t")]
}


.ep_summarize_target_times <- function(count_mat, times, target_times,
                                       pilevel, offset, target_d) {
  idx <- match(as.numeric(target_times), as.numeric(times))
  if (anyNA(idx)) {
    stop("All target_times must be present in times.")
  }

  plower <- (1 - pilevel) / 2
  pupper <- 1 - plower

  dt_list <- lapply(seq_along(idx), function(i) {
    counts <- count_mat[idx[i], ]
    q <- stats::quantile(
      counts,
      probs = c(0.5, plower, pupper),
      names = FALSE
    ) + offset

    data.table::data.table(
      t = as.numeric(target_times[i]),
      n = q[1],
      pilevel = pilevel,
      lower = q[2],
      upper = q[3],
      mean = mean(counts) + offset,
      var = if (length(counts) > 1) stats::var(counts) else 0,
      target_d = target_d,
      prob_gt_target_d = mean(counts >= target_d - offset)
    )
  })

  data.table::rbindlist(dt_list, use.names = TRUE)
}
