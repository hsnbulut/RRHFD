# Core functions for adaptive robust rolling Hotelling T-squared fault detection.

as_matrix <- function(x, arg = "x") {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!all(is.finite(x))) {
    stop(arg, " must contain only finite numeric values.", call. = FALSE)
  }
  x
}

safe_cov <- function(x, ridge = 1e-8) {
  x <- as_matrix(x)
  p <- ncol(x)
  if (nrow(x) <= 1L) {
    return(diag(ridge, p))
  }
  s <- stats::cov(x)
  if (!all(is.finite(s))) {
    s <- diag(p)
  }
  s + diag(ridge, p)
}

safe_inv_quad <- function(cov, diff, ridge = 1e-8) {
  cov <- as_matrix(cov, "cov")
  diff <- as.numeric(diff)
  p <- ncol(cov)
  if (length(diff) != p) {
    stop("diff length must match covariance dimension.", call. = FALSE)
  }
  cov2 <- cov + diag(ridge, p)
  as.numeric(crossprod(diff, solve(cov2, diff)))
}

logdet_safe <- function(cov, ridge = 1e-8) {
  cov <- as_matrix(cov, "cov")
  p <- ncol(cov)
  as.numeric(determinant(cov + diag(ridge, p), logarithm = TRUE)$modulus[1L])
}

#' Adaptive MCD Alpha
#'
#' Compute the adaptive MCD alpha used by the RRH-FD detector.
#'
#' @param w_roll Rolling window length.
#' @param n_spikes Expected number of pre-fault spikes inside the window.
#' @param c_factor Inflation factor applied to the estimated contamination rate.
#' @param lower Lower bound for alpha.
#' @param upper Upper bound for alpha.
#'
#' @return A numeric alpha value between `lower` and `upper`.
#' @examples
#' adaptive_mcd_alpha(50, n_spikes = 5)
adaptive_mcd_alpha <- function(w_roll,
                               n_spikes = 5,
                               c_factor = 1.25,
                               lower = 0.50,
                               upper = 0.95) {
  if (length(w_roll) != 1L || !is.finite(w_roll) || w_roll <= 1) {
    stop("w_roll must be a finite value greater than 1.", call. = FALSE)
  }
  if (lower <= 0 || upper > 1 || lower > upper) {
    stop("lower and upper must satisfy 0 < lower <= upper <= 1.", call. = FALSE)
  }
  eps_hat <- n_spikes / w_roll
  alpha <- 1 - c_factor * eps_hat
  max(lower, min(upper, alpha))
}

#' RRH-FD Control Settings
#'
#' Build a small list of detector settings.
#'
#' @param w_roll Rolling window length.
#' @param alpha_level Target false alarm level.
#' @param n_spikes Expected number of pre-fault spikes inside the window.
#' @param c_factor Inflation factor for adaptive alpha.
#' @param alpha_mcd Optional fixed MCD alpha. If `NULL`, adaptive alpha is used.
#' @param alpha_lower Lower bound for adaptive alpha.
#' @param alpha_upper Upper bound for adaptive alpha.
#' @param ridge Ridge added to covariance matrices.
#' @param use_robustbase Use `robustbase::covMcd()` when available.
#'
#' @return A list of control settings.
rrh_control <- function(w_roll = 50,
                        alpha_level = 0.05,
                        n_spikes = 5,
                        c_factor = 1.25,
                        alpha_mcd = NULL,
                        alpha_lower = 0.50,
                        alpha_upper = 0.95,
                        ridge = 1e-8,
                        use_robustbase = NULL) {
  if (is.null(alpha_mcd)) {
    alpha_mcd <- adaptive_mcd_alpha(
      w_roll = w_roll,
      n_spikes = n_spikes,
      c_factor = c_factor,
      lower = alpha_lower,
      upper = alpha_upper
    )
  }
  list(
    w_roll = as.integer(w_roll),
    alpha_level = alpha_level,
    n_spikes = n_spikes,
    c_factor = c_factor,
    alpha_mcd = alpha_mcd,
    ridge = ridge,
    use_robustbase = use_robustbase
  )
}

trimmed_cov <- function(x, alpha = 0.75, ridge = 1e-8) {
  x <- as_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  h <- max(p + 1L, ceiling(alpha * n))
  h <- min(h, n)

  center0 <- stats::median(x[, 1L])
  center0 <- apply(x, 2L, stats::median)
  scale0 <- apply(x, 2L, stats::mad, constant = 1.4826)
  scale0[!is.finite(scale0) | scale0 <= 0] <- 1
  z <- sweep(sweep(x, 2L, center0, "-"), 2L, scale0, "/")
  d2 <- rowSums(z^2)
  keep <- order(d2)[seq_len(h)]
  x_keep <- x[keep, , drop = FALSE]
  center <- colMeans(x_keep)
  cov <- safe_cov(x_keep, ridge = ridge)

  list(
    center = center,
    cov = cov,
    alpha = alpha,
    h = h,
    subset = keep,
    method = "trimmed"
  )
}

#' Robust MCD Covariance Estimate
#'
#' Estimate a robust center and covariance for a rolling residual window.
#'
#' @param x Numeric matrix whose rows are residual observations.
#' @param alpha MCD alpha value.
#' @param ridge Ridge added to covariance matrices.
#' @param use_robustbase If `TRUE`, require `robustbase::covMcd()`. If `FALSE`,
#'   always use the trimmed fallback. If `NULL`, use `robustbase` when installed.
#'
#' @return A list with `center`, `cov`, `alpha`, `h`, `subset`, and `method`.
#' @examples
#' x <- simulate_residuals(80, 2, seed = 1)
#' rrh_cov_mcd(x, alpha = 0.85)
rrh_cov_mcd <- function(x,
                        alpha = 0.75,
                        ridge = 1e-8,
                        use_robustbase = NULL) {
  x <- as_matrix(x)
  if (nrow(x) <= ncol(x)) {
    stop("x must have more rows than columns.", call. = FALSE)
  }
  if (alpha <= 0 || alpha > 1) {
    stop("alpha must be in (0, 1].", call. = FALSE)
  }

  has_robustbase <- requireNamespace("robustbase", quietly = TRUE)
  if (is.null(use_robustbase)) {
    use_robustbase <- has_robustbase
  }
  if (isTRUE(use_robustbase) && !has_robustbase) {
    stop("robustbase is not installed. Set use_robustbase = FALSE or install robustbase.",
         call. = FALSE)
  }

  if (isTRUE(use_robustbase)) {
    cov_mcd <- getExportedValue("robustbase", "covMcd")
    fit <- cov_mcd(x, alpha = alpha)
    cov <- as.matrix(fit$cov) + diag(ridge, ncol(x))
    return(list(
      center = as.numeric(fit$center),
      cov = cov,
      alpha = alpha,
      h = fit$quan,
      subset = fit$best,
      method = "robustbase::covMcd"
    ))
  }

  trimmed_cov(x, alpha = alpha, ridge = ridge)
}

#' One-Step RRH-FD Statistic
#'
#' Compute the robust rolling Hotelling statistic for one residual vector.
#'
#' @param rk Current residual vector.
#' @param window Historical residual window, one row per time point.
#' @param alpha_mcd MCD alpha.
#' @param ridge Ridge added to covariance matrices.
#' @param use_robustbase Passed to [rrh_cov_mcd()].
#'
#' @return A list with statistic, center, covariance, alpha, and estimator method.
rrh_stat <- function(rk,
                     window,
                     alpha_mcd = 0.75,
                     ridge = 1e-8,
                     use_robustbase = NULL) {
  window <- as_matrix(window, "window")
  rk <- as.numeric(rk)
  if (length(rk) != ncol(window)) {
    stop("rk length must equal ncol(window).", call. = FALSE)
  }
  est <- rrh_cov_mcd(
    x = window,
    alpha = alpha_mcd,
    ridge = ridge,
    use_robustbase = use_robustbase
  )
  diff <- rk - est$center
  stat <- safe_inv_quad(est$cov, diff, ridge = ridge)
  list(
    statistic = stat,
    center = est$center,
    cov = est$cov,
    alpha = est$alpha,
    method = est$method
  )
}

#' Empirical RRH-FD Threshold Calibration
#'
#' Calibrate an RRH-FD threshold under a clean/null residual stream.
#'
#' @param null_residuals Clean residual matrix. If supplied, random time points
#'   are sampled and their preceding windows are used.
#' @param simulate_fn Optional function taking one integer `n` and returning an
#'   `n` by `p` residual matrix.
#' @param B Number of calibration replicates.
#' @param w_roll Rolling window length.
#' @param alpha_level Target false alarm level.
#' @param n_spikes Expected number of spikes for adaptive alpha.
#' @param alpha_mcd Optional fixed MCD alpha. If `NULL`, adaptive alpha is used.
#' @param ridge Ridge added to covariance matrices.
#' @param seed Optional random seed.
#' @param use_robustbase Passed to [rrh_cov_mcd()].
#'
#' @return A list with threshold, null statistics, and settings.
#' @examples
#' clean <- simulate_residuals(200, 2, seed = 1)
#' rrh_calibrate_threshold(clean, B = 50, w_roll = 30)
rrh_calibrate_threshold <- function(null_residuals = NULL,
                                    simulate_fn = NULL,
                                    B = 1000,
                                    w_roll = 50,
                                    alpha_level = 0.05,
                                    n_spikes = 5,
                                    alpha_mcd = NULL,
                                    ridge = 1e-8,
                                    seed = NULL,
                                    use_robustbase = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(alpha_mcd)) {
    alpha_mcd <- adaptive_mcd_alpha(w_roll = w_roll, n_spikes = n_spikes)
  }
  if (is.null(null_residuals) && is.null(simulate_fn)) {
    stop("Provide null_residuals or simulate_fn.", call. = FALSE)
  }

  vals <- numeric(B)
  if (!is.null(null_residuals)) {
    null_residuals <- as_matrix(null_residuals, "null_residuals")
    n <- nrow(null_residuals)
    if (n <= w_roll) {
      stop("null_residuals must have more rows than w_roll.", call. = FALSE)
    }
    idx <- sample(seq.int(w_roll + 1L, n), size = B, replace = TRUE)
    for (b in seq_len(B)) {
      k <- idx[b]
      window <- null_residuals[(k - w_roll):(k - 1L), , drop = FALSE]
      vals[b] <- rrh_stat(
        rk = null_residuals[k, ],
        window = window,
        alpha_mcd = alpha_mcd,
        ridge = ridge,
        use_robustbase = use_robustbase
      )$statistic
    }
  } else {
    for (b in seq_len(B)) {
      sim <- as_matrix(simulate_fn(w_roll + 1L), "simulate_fn output")
      if (nrow(sim) < w_roll + 1L) {
        stop("simulate_fn must return at least w_roll + 1 rows.", call. = FALSE)
      }
      window <- sim[seq_len(w_roll), , drop = FALSE]
      vals[b] <- rrh_stat(
        rk = sim[w_roll + 1L, ],
        window = window,
        alpha_mcd = alpha_mcd,
        ridge = ridge,
        use_robustbase = use_robustbase
      )$statistic
    }
  }

  threshold <- as.numeric(stats::quantile(vals, probs = 1 - alpha_level,
                                          type = 8, na.rm = TRUE))
  list(
    threshold = threshold,
    null_stats = vals,
    alpha_level = alpha_level,
    alpha_mcd = alpha_mcd,
    w_roll = w_roll,
    B = B
  )
}

#' Run Online RRH-FD Detection
#'
#' Compute rolling robust Hotelling statistics and alarms for a residual matrix.
#'
#' @param residuals Numeric residual matrix, one row per time point.
#' @param w_roll Rolling window length.
#' @param threshold Alarm threshold.
#' @param n_spikes Expected number of spikes for adaptive alpha.
#' @param alpha_mcd Optional fixed MCD alpha. If `NULL`, adaptive alpha is used.
#' @param ridge Ridge added to covariance matrices.
#' @param use_robustbase Passed to [rrh_cov_mcd()].
#'
#' @return A data frame with `time`, `stat_rrh`, `alarm`, `alpha_mcd`,
#'   `logBIR`, and estimator `method`.
#' @examples
#' clean <- simulate_residuals(200, 2, seed = 1)
#' cal <- rrh_calibrate_threshold(clean, B = 50, w_roll = 30)
#' det <- rrh_detect(clean, w_roll = 30, threshold = cal$threshold)
#' head(det)
rrh_detect <- function(residuals,
                       w_roll = 50,
                       threshold,
                       n_spikes = 5,
                       alpha_mcd = NULL,
                       ridge = 1e-8,
                       use_robustbase = NULL) {
  residuals <- as_matrix(residuals, "residuals")
  n <- nrow(residuals)
  p <- ncol(residuals)
  if (n <= w_roll) {
    stop("residuals must have more rows than w_roll.", call. = FALSE)
  }
  if (missing(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be a finite scalar.", call. = FALSE)
  }
  if (is.null(alpha_mcd)) {
    alpha_mcd <- adaptive_mcd_alpha(w_roll = w_roll, n_spikes = n_spikes)
  }

  times <- seq.int(w_roll + 1L, n)
  stat <- rep(NA_real_, length(times))
  alarm <- rep(0L, length(times))
  log_bir <- rep(NA_real_, length(times))
  method <- rep(NA_character_, length(times))

  for (i in seq_along(times)) {
    k <- times[i]
    window <- residuals[(k - w_roll):(k - 1L), , drop = FALSE]
    s_classic <- safe_cov(window, ridge = ridge)
    out <- rrh_stat(
      rk = residuals[k, ],
      window = window,
      alpha_mcd = alpha_mcd,
      ridge = ridge,
      use_robustbase = use_robustbase
    )
    stat[i] <- out$statistic
    alarm[i] <- as.integer(stat[i] > threshold)
    log_bir[i] <- logdet_safe(s_classic, ridge = ridge) -
      logdet_safe(out$cov, ridge = ridge)
    method[i] <- out$method
  }

  data.frame(
    time = times,
    stat_rrh = stat,
    alarm = alarm,
    threshold = threshold,
    alpha_mcd = alpha_mcd,
    logBIR = log_bir,
    method = method,
    stringsAsFactors = FALSE
  )
}

#' Fault Detection Metrics
#'
#' Compute simple pre/post fault metrics from an RRH-FD detection result.
#'
#' @param detections Data frame returned by [rrh_detect()].
#' @param fault_time Fault start time.
#' @param horizon Optional horizon for early detection rate.
#'
#' @return A one-row data frame with false alarm rate, missed detection rate,
#'   detection delay, detection flag, and optional horizon detection flag.
rrh_metrics <- function(detections, fault_time, horizon = 50) {
  if (!all(c("time", "alarm") %in% names(detections))) {
    stop("detections must contain time and alarm columns.", call. = FALSE)
  }
  pre <- detections$time < fault_time
  post <- detections$time >= fault_time
  pf <- if (any(pre)) mean(detections$alarm[pre] == 1L, na.rm = TRUE) else NA_real_
  pm <- if (any(post)) mean(detections$alarm[post] == 0L, na.rm = TRUE) else NA_real_
  post_alarm_times <- detections$time[post & detections$alarm == 1L]
  detected <- as.integer(length(post_alarm_times) > 0L)
  delay <- if (detected) min(post_alarm_times) - fault_time else NA_real_
  horizon_end <- fault_time + horizon
  in_horizon <- detections$time >= fault_time & detections$time <= horizon_end
  det_h <- as.integer(any(detections$alarm[in_horizon] == 1L, na.rm = TRUE))
  data.frame(
    pf = pf,
    pm = pm,
    delay = delay,
    detected = detected,
    horizon = horizon,
    detected_in_horizon = det_h
  )
}

#' Simulate Residuals
#'
#' Generate a simple multivariate residual stream for examples and calibration.
#'
#' @param n Number of observations.
#' @param p Residual dimension.
#' @param sigma Covariance matrix. Defaults to identity.
#' @param fault_time Optional fault start time.
#' @param fault_mean Mean shift after `fault_time`.
#' @param contamination Optional list with `index` and `shift` entries.
#' @param seed Optional random seed.
#'
#' @return An `n` by `p` residual matrix.
simulate_residuals <- function(n,
                               p = 2,
                               sigma = diag(p),
                               fault_time = NULL,
                               fault_mean = rep(0, p),
                               contamination = NULL,
                               seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sigma <- as_matrix(sigma, "sigma")
  if (ncol(sigma) != p || nrow(sigma) != p) {
    stop("sigma must be a p by p matrix.", call. = FALSE)
  }
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x <- z %*% chol(sigma)

  if (!is.null(fault_time)) {
    if (fault_time < 1 || fault_time > n) {
      stop("fault_time must be between 1 and n.", call. = FALSE)
    }
    x[fault_time:n, ] <- sweep(
      x[fault_time:n, , drop = FALSE],
      2L,
      fault_mean,
      "+"
    )
  }

  if (!is.null(contamination)) {
    idx <- contamination$index
    shift <- contamination$shift
    if (is.null(idx) || is.null(shift)) {
      stop("contamination must contain index and shift.", call. = FALSE)
    }
    idx <- idx[idx >= 1 & idx <= n]
    x[idx, ] <- sweep(x[idx, , drop = FALSE], 2L, shift, "+")
  }

  colnames(x) <- paste0("r", seq_len(p))
  x
}
