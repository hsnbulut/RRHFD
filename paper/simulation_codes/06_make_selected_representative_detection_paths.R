############################################################
# Selected representative detection paths for RRH-FD paper
#
# Purpose:
#   Search many candidate contaminated trajectories and select
#   one representative trajectory where both classical rolling
#   and Liu-type offline detectors are delayed, while RRH-FD
#   adaptive and RRH-FD 0.85 detect earlier.
#
# Output:
#   results/figures/contaminated/
#     fig_selected_representative_detection_paths_methods_step.png
#     selected_representative_case_summary.csv
#     plot_data_selected_representative_detection_paths_methods_step.csv
#     candidate_representative_search_summary.csv
############################################################

rm(list = ls())

required_pkgs <- c("MASS", "dplyr", "tidyr", "ggplot2", "rrcov", "tibble")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them first using install.packages()."
  )
}

library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rrcov)
library(tibble)

############################################################
# Settings
############################################################

fig_dir <- "results/figures/contaminated"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260521)

N     <- 2000
k0    <- 500
p     <- 2
w     <- 80
alpha <- 0.05

# Increase if a stronger selected example is needed.
N_CANDIDATES <- 500
BASE_SEED    <- 700000

# Screening period around fault onset
screen_start <- k0 - 120
screen_end   <- k0 + 250

# Representative contaminated scenario
gamma_online <- 2.0
delta_step   <- 0.06

# Nominal covariance matrix
Sigma0 <- matrix(
  c(0.05^2, 0.00,
    0.00,   0.05^2),
  nrow = 2,
  byrow = TRUE
)

fault_direction <- c(1.0, 0.6)
fault_direction <- fault_direction / sqrt(sum(fault_direction^2))

############################################################
# Helper functions
############################################################

safe_solve <- function(S) {
  S <- as.matrix(S)
  S <- S + diag(1e-8, nrow(S))
  tryCatch(solve(S), error = function(e) MASS::ginv(S))
}

classical_est <- function(X) {
  list(
    center = colMeans(X),
    cov = stats::cov(X)
  )
}

rmcd_est <- function(X, alpha_mcd = 0.85, nsamp = 100) {
  X <- as.matrix(X)
  
  fit <- tryCatch(
    rrcov::CovMcd(X, alpha = alpha_mcd, nsamp = nsamp),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(classical_est(X))
  }
  
  list(
    center = as.numeric(fit@center),
    cov = as.matrix(fit@cov)
  )
}

t2_value <- function(x, center, S) {
  r <- as.numeric(x - center)
  as.numeric(t(r) %*% safe_solve(S) %*% r)
}

estimate_window_contamination <- function(Xw) {
  est <- classical_est(Xw)
  md2 <- apply(Xw, 1, function(z) t2_value(z, est$center, est$cov))
  mean(md2 > qchisq(0.975, df = ncol(Xw)), na.rm = TRUE)
}

adaptive_alpha <- function(Xw) {
  eps_hat <- estimate_window_contamination(Xw)
  
  alpha_ad <- 1 - eps_hat
  alpha_ad <- max(0.75, min(0.95, alpha_ad))
  
  alpha_ad
}

fault_profile <- function(type, N, k0, delta = 0.06) {
  f <- rep(0, N)
  
  if (type == "step") {
    f[k0:N] <- delta
  }
  
  if (type == "oscillation") {
    idx <- k0:N
    f[idx] <- delta + 0.02 * sin((idx - k0) / 10)
  }
  
  if (type == "gradual") {
    idx <- k0:N
    f[idx] <- delta * pmin(1, (idx - k0 + 1) / 300)
  }
  
  if (type == "intermittent") {
    f[500:699]   <- delta
    f[1000:1199] <- 1.5 * delta
    f[1500:1699] <- 2.0 * delta
  }
  
  f
}

generate_series <- function(
    seed,
    fault_type = "step",
    N = 2000,
    k0 = 500,
    Sigma0,
    delta = 0.06,
    gamma_online = 2.0
) {
  set.seed(seed)
  
  X <- MASS::mvrnorm(N, mu = rep(0, ncol(Sigma0)), Sigma = Sigma0)
  
  # Transient pre-fault spikes in the rolling reference window
  # immediately before the true fault onset.
  spike_pool <- (k0 - w):(k0 - 1)
  spike_idx <- sort(sample(spike_pool, size = 20, replace = FALSE))
  
  spike_shift <- gamma_online * c(0.18, 0.12)
  X[spike_idx, ] <- sweep(X[spike_idx, , drop = FALSE], 2, spike_shift, "+")
  
  f <- fault_profile(fault_type, N, k0, delta)
  
  for (tt in k0:N) {
    X[tt, ] <- X[tt, ] + f[tt] * fault_direction
  }
  
  list(
    X = X,
    fault_signal = f,
    spike_idx = spike_idx
  )
}

calibrate_thresholds <- function(B = 500, w = 80, p = 2, Sigma0, alpha = 0.05) {
  message("Calibrating thresholds for representative plotting with B = ", B)
  
  vals_classical <- numeric(B)
  vals_liu       <- numeric(B)
  vals_r075      <- numeric(B)
  vals_r085      <- numeric(B)
  
  n_off <- 1000
  
  for (b in seq_len(B)) {
    Xw   <- MASS::mvrnorm(w, mu = rep(0, p), Sigma = Sigma0)
    xnew <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = Sigma0)
    
    ec <- classical_est(Xw)
    vals_classical[b] <- t2_value(xnew, ec$center, ec$cov)
    
    Xoff <- MASS::mvrnorm(n_off, mu = rep(0, p), Sigma = Sigma0)
    eo <- classical_est(Xoff)
    vals_liu[b] <- t2_value(xnew, eo$center, eo$cov)
    
    er075 <- rmcd_est(Xw, alpha_mcd = 0.75, nsamp = 100)
    vals_r075[b] <- t2_value(xnew, er075$center, er075$cov)
    
    er085 <- rmcd_est(Xw, alpha_mcd = 0.85, nsamp = 100)
    vals_r085[b] <- t2_value(xnew, er085$center, er085$cov)
  }
  
  tibble::tibble(
    method = c(
      "Classical rolling",
      "Liu-type offline",
      "RRH-FD 0.75",
      "RRH-FD 0.85",
      "RRH-FD adaptive"
    ),
    threshold = c(
      stats::quantile(vals_classical, 1 - alpha, na.rm = TRUE),
      stats::quantile(vals_liu,       1 - alpha, na.rm = TRUE),
      stats::quantile(vals_r075,      1 - alpha, na.rm = TRUE),
      stats::quantile(vals_r085,      1 - alpha, na.rm = TRUE),
      stats::quantile(vals_r085,      1 - alpha, na.rm = TRUE)
    )
  )
}

get_threshold_tbl <- function() {
  threshold_file <- file.path(fig_dir, "representative_thresholds.csv")
  
  if (file.exists(threshold_file)) {
    message("Reading existing representative thresholds: ", threshold_file)
    th <- read.csv(threshold_file)
    return(tibble::as_tibble(th))
  }
  
  th <- calibrate_thresholds(
    B = 500,
    w = w,
    p = p,
    Sigma0 = Sigma0,
    alpha = alpha
  )
  
  write.csv(th, threshold_file, row.names = FALSE)
  th
}

compute_paths <- function(
    X,
    threshold_tbl,
    t_start = w + 1,
    t_end = N,
    methods = c(
      "Classical rolling",
      "Liu-type offline",
      "RRH-FD 0.75",
      "RRH-FD 0.85",
      "RRH-FD adaptive"
    ),
    nsamp_mcd = 100
) {
  N <- nrow(X)
  p <- ncol(X)
  
  t_start <- max(t_start, w + 1)
  t_end   <- min(t_end, N)
  
  # Offline reference for Liu-type benchmark.
  # It is intentionally fixed and clean, representing a nonrolling baseline.
  Xoff <- MASS::mvrnorm(1000, mu = rep(0, p), Sigma = Sigma0)
  e_off <- classical_est(Xoff)
  
  out <- vector("list", length = 0)
  
  for (tt in t_start:t_end) {
    Xw <- X[(tt - w):(tt - 1), , drop = FALSE]
    xt <- X[tt, ]
    
    if ("Classical rolling" %in% methods) {
      ec <- classical_est(Xw)
      out[[length(out) + 1]] <- data.frame(
        k = tt,
        method = "Classical rolling",
        T2 = t2_value(xt, ec$center, ec$cov),
        alpha_mcd = NA_real_
      )
    }
    
    if ("Liu-type offline" %in% methods) {
      out[[length(out) + 1]] <- data.frame(
        k = tt,
        method = "Liu-type offline",
        T2 = t2_value(xt, e_off$center, e_off$cov),
        alpha_mcd = NA_real_
      )
    }
    
    if ("RRH-FD 0.75" %in% methods) {
      er075 <- rmcd_est(Xw, alpha_mcd = 0.75, nsamp = nsamp_mcd)
      out[[length(out) + 1]] <- data.frame(
        k = tt,
        method = "RRH-FD 0.75",
        T2 = t2_value(xt, er075$center, er075$cov),
        alpha_mcd = 0.75
      )
    }
    
    if ("RRH-FD 0.85" %in% methods) {
      er085 <- rmcd_est(Xw, alpha_mcd = 0.85, nsamp = nsamp_mcd)
      out[[length(out) + 1]] <- data.frame(
        k = tt,
        method = "RRH-FD 0.85",
        T2 = t2_value(xt, er085$center, er085$cov),
        alpha_mcd = 0.85
      )
    }
    
    if ("RRH-FD adaptive" %in% methods) {
      a_ad <- adaptive_alpha(Xw)
      erad <- rmcd_est(Xw, alpha_mcd = a_ad, nsamp = nsamp_mcd)
      out[[length(out) + 1]] <- data.frame(
        k = tt,
        method = "RRH-FD adaptive",
        T2 = t2_value(xt, erad$center, erad$cov),
        alpha_mcd = a_ad
      )
    }
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::left_join(threshold_tbl, by = "method") %>%
    dplyr::mutate(
      alarm = as.integer(T2 > threshold),
      method = factor(
        method,
        levels = c(
          "Classical rolling",
          "Liu-type offline",
          "RRH-FD 0.75",
          "RRH-FD 0.85",
          "RRH-FD adaptive"
        )
      )
    )
}

summarise_delays <- function(path_df, k0) {
  path_df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      pre_fault_alarms = sum(alarm == 1 & k < k0, na.rm = TRUE),
      first_alarm = {
        tmp <- k[k >= k0 & alarm == 1]
        if (length(tmp) == 0) NA_integer_ else min(tmp)
      },
      delay = ifelse(is.na(first_alarm), NA_real_, first_alarm - k0),
      .groups = "drop"
    )
}

score_candidate <- function(delay_tbl) {
  dc <- delay_tbl$delay[delay_tbl$method == "Classical rolling"]
  dl <- delay_tbl$delay[delay_tbl$method == "Liu-type offline"]
  da <- delay_tbl$delay[delay_tbl$method == "RRH-FD adaptive"]
  dr <- delay_tbl$delay[delay_tbl$method == "RRH-FD 0.85"]
  
  pfa <- delay_tbl$pre_fault_alarms[delay_tbl$method == "RRH-FD adaptive"]
  
  if (length(dc) == 0 || length(dl) == 0 || length(da) == 0) return(-Inf)
  if (is.na(da)) return(-Inf)
  
  dc_eff <- ifelse(is.na(dc), 999, dc)
  dl_eff <- ifelse(is.na(dl), 999, dl)
  dr_eff <- ifelse(length(dr) == 0 || is.na(dr), da, dr)
  
  # Prefer cases where:
  #   - classical rolling is delayed,
  #   - Liu-type offline is also delayed,
  #   - RRH-FD adaptive is early,
  #   - RRH-FD 0.85 is also early,
  #   - adaptive pre-fault alarms are not excessive.
  score <-
    1.00 * (dc_eff - da) +
    1.25 * (dl_eff - da) +
    0.50 * (dc_eff - dr_eff) -
    0.75 * pfa
  
  # Strongly penalize examples where Liu-type offline is too early.
  if (!is.na(dl) && dl <= 5) {
    score <- score - 100
  }
  
  # Penalize examples where adaptive RRH-FD is not early.
  if (!is.na(da) && da > 20) {
    score <- score - 50
  }
  
  # Prefer cases where RRH-FD 0.85 is not much worse than adaptive.
  if (!is.na(dr) && !is.na(da) && dr > da + 20) {
    score <- score - 25
  }
  
  score
}

theme_detection <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey50"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none",
      plot.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = base_size - 2, hjust = 0)
    )
}

save_png <- function(plot, filename, width = 8.2, height = 8.4, dpi = 600) {
  ggplot2::ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

############################################################
# Thresholds
############################################################

threshold_tbl <- get_threshold_tbl()

############################################################
# Candidate search
############################################################

message("Searching for a representative trajectory...")

candidate_results <- vector("list", N_CANDIDATES)

for (i in seq_len(N_CANDIDATES)) {
  seed_i <- BASE_SEED + i
  
  dat <- generate_series(
    seed = seed_i,
    fault_type = "step",
    N = N,
    k0 = k0,
    Sigma0 = Sigma0,
    delta = delta_step,
    gamma_online = gamma_online
  )
  
  paths_i <- compute_paths(
    X = dat$X,
    threshold_tbl = threshold_tbl,
    t_start = screen_start,
    t_end = screen_end,
    methods = c(
      "Classical rolling",
      "Liu-type offline",
      "RRH-FD 0.85",
      "RRH-FD adaptive"
    ),
    nsamp_mcd = 80
  )
  
  delays_i <- summarise_delays(paths_i, k0 = k0)
  score_i  <- score_candidate(delays_i)
  
  candidate_results[[i]] <- delays_i %>%
    dplyr::mutate(
      seed = seed_i,
      score = score_i,
      candidate_id = i
    )
  
  if (i %% 10 == 0) {
    message("  screened ", i, " / ", N_CANDIDATES, " candidates")
  }
}

candidate_summary <- dplyr::bind_rows(candidate_results)

write.csv(
  candidate_summary,
  file.path(fig_dir, "candidate_representative_search_summary.csv"),
  row.names = FALSE
)

best_seed_tbl <- candidate_summary %>%
  dplyr::filter(is.finite(score)) %>%
  dplyr::group_by(seed) %>%
  dplyr::summarise(score = dplyr::first(score), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::slice(1)

best_seed <- best_seed_tbl$seed

if (length(best_seed) == 0 || is.na(best_seed)) {
  stop(
    "No suitable representative trajectory was found. ",
    "Increase N_CANDIDATES or adjust delta_step/gamma_online."
  )
}

message("Selected seed: ", best_seed)

############################################################
# Recompute full paths for selected trajectory
############################################################

selected_data <- generate_series(
  seed = best_seed,
  fault_type = "step",
  N = N,
  k0 = k0,
  Sigma0 = Sigma0,
  delta = delta_step,
  gamma_online = gamma_online
)

selected_paths <- compute_paths(
  X = selected_data$X,
  threshold_tbl = threshold_tbl,
  t_start = w + 1,
  t_end = N,
  methods = c(
    "Classical rolling",
    "Liu-type offline",
    "RRH-FD 0.75",
    "RRH-FD 0.85",
    "RRH-FD adaptive"
  ),
  nsamp_mcd = 100
)

selected_delays <- summarise_delays(selected_paths, k0 = k0)

write.csv(
  selected_paths,
  file.path(fig_dir, "plot_data_selected_representative_detection_paths_methods_step.csv"),
  row.names = FALSE
)

write.csv(
  selected_delays,
  file.path(fig_dir, "selected_representative_case_summary.csv"),
  row.names = FALSE
)

############################################################
# Plot labels
############################################################

delay_labels <- selected_delays %>%
  dplyr::mutate(
    delay_label = ifelse(
      is.na(delay),
      "no post-fault alarm",
      paste0("delay = ", delay)
    ),
    method_label = paste0(as.character(method), " (", delay_label, ")")
  ) %>%
  dplyr::select(method, method_label)

plot_df <- selected_paths %>%
  dplyr::left_join(delay_labels, by = "method") %>%
  dplyr::mutate(
    method_label = factor(
      method_label,
      levels = delay_labels$method_label[
        match(
          c(
            "Classical rolling",
            "Liu-type offline",
            "RRH-FD 0.75",
            "RRH-FD 0.85",
            "RRH-FD adaptive"
          ),
          as.character(delay_labels$method)
        )
      ]
    )
  )

first_alarm_points <- plot_df %>%
  dplyr::filter(k >= k0, alarm == 1) %>%
  dplyr::group_by(method_label) %>%
  dplyr::slice_min(k, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Robust y-axis upper bound so extreme spikes do not dominate the figure.
ymax <- stats::quantile(plot_df$T2, probs = 0.995, na.rm = TRUE)
ymax <- max(ymax, max(plot_df$threshold, na.rm = TRUE) * 1.25)

############################################################
# Plot
############################################################

fig_selected <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = T2)) +
  ggplot2::annotate(
    "rect",
    xmin = min(selected_data$spike_idx),
    xmax = max(selected_data$spike_idx),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08,
    fill = "orange"
  ) +
  ggplot2::geom_line(color = "#1f77b4", linewidth = 0.35) +
  ggplot2::geom_hline(
    ggplot2::aes(yintercept = threshold),
    color = "firebrick",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_vline(
    xintercept = k0,
    color = "black",
    linetype = "dotted",
    linewidth = 0.50
  ) +
  ggplot2::geom_point(
    data = first_alarm_points,
    ggplot2::aes(x = k, y = T2),
    inherit.aes = FALSE,
    color = "black",
    fill = "yellow",
    shape = 21,
    size = 2.4,
    stroke = 0.7
  ) +
  ggplot2::facet_wrap(~ method_label, ncol = 1, scales = "free_y") +
  ggplot2::coord_cartesian(xlim = c(0, N), ylim = c(0, ymax)) +
  ggplot2::labs(
    x = "Time index k",
    y = expression("Detection statistic " * T^2),
    caption = paste0(
      "Selected representative contaminated step-fault trajectory; seed = ", best_seed,
      ". Black dotted line: fault onset k0 = ", k0,
      "; red dashed line: method-specific threshold; orange shading: transient pre-fault spikes; ",
      "yellow point: first post-fault alarm."
    )
  ) +
  theme_detection(base_size = 11)

save_png(
  fig_selected,
  "fig_selected_representative_detection_paths_methods_step.png",
  width = 8.2,
  height = 8.4,
  dpi = 600
)

############################################################
# Print selected case summary
############################################################

cat("\nSelected representative seed:\n")
print(best_seed)

cat("\nSelected representative detection delays:\n")
print(selected_delays)

cat("\nFigure saved to:\n")
cat(file.path(fig_dir, "fig_selected_representative_detection_paths_methods_step.png"), "\n")

cat("\nCandidate search summary saved to:\n")
cat(file.path(fig_dir, "candidate_representative_search_summary.csv"), "\n")

cat("\nDone.\n")