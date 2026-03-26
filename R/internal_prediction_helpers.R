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
