############################################################
## prepare_rrhfd_main_grid_and_thresholds.R
## Creates scenario grid, job grid, and calibrated thresholds
############################################################

rm(list = ls())

req_pkgs <- c("MASS", "robustbase", "dplyr", "readr")

for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(MASS)
library(robustbase)
library(dplyr)
library(readr)

set.seed(20260516)

############################################################
## Global settings
############################################################

k_final <- 2000
k_fault <- 500
alpha_level <- 0.05

tau <- 0.01
x0_mean <- c(0, 1)
P0 <- (0.01^2) * diag(2)

Q <- (0.01^2) * diag(2)
R <- 0.1 * diag(2)

p <- 2
omega <- 10

N_SPIKES <- 5

## Main simulation size
R_REPS_MAIN <- as.integer(Sys.getenv("R_REPS_MAIN", "1000"))
R_BLOCK <- as.integer(Sys.getenv("R_BLOCK", "100"))

## Threshold calibration size
B_CAL <- as.integer(Sys.getenv("B_CAL", "10000"))

out_table_dir <- "results/tables"
dir.create(out_table_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## RSDTS model
############################################################

y3_fun <- function(k) {
  0.1 * pi + 0.005 * pi * sin(k / 100)
}

A_fun <- function(k) {
  y3 <- y3_fun(k)
  matrix(
    c(
      1,          y3 * tau,
      -y3 * tau,  1
    ),
    nrow = 2,
    byrow = TRUE
  )
}

C_mat <- matrix(
  c(
    1,  0,
    0, -1
  ),
  nrow = 2,
  byrow = TRUE
)

F_vec <- c(1, 0)

generate_ybar <- function(k_final = 2000) {
  xbar <- matrix(0, nrow = k_final, ncol = 2)
  ybar <- matrix(0, nrow = k_final, ncol = 2)
  
  xbar[1, ] <- x0_mean
  ybar[1, ] <- as.numeric(C_mat %*% xbar[1, ])
  
  for (k in 2:k_final) {
    A_prev <- A_fun(k - 1)
    xbar[k, ] <- as.numeric(A_prev %*% xbar[k - 1, ])
    ybar[k, ] <- as.numeric(C_mat %*% xbar[k, ])
  }
  
  list(xbar = xbar, ybar = ybar)
}

YBAR <- generate_ybar(k_final)$ybar

fault_signal <- function(k, type = "step", step_delta = 0.06) {
  if (k < k_fault) return(0)
  
  if (type == "step") {
    return(step_delta)
  }
  
  if (type == "oscillation") {
    return(0.06 + 0.02 * sin((k - k_fault) / 10))
  }
  
  if (type == "gradual") {
    return(0.06 + (0.005^2) * (k - k_fault))
  }
  
  if (type == "intermittent") {
    if (k >= 500 && k < 700) return(0.10)
    if (k >= 1000 && k < 1200) return(0.14)
    if (k >= 1500 && k < 1700) return(0.18)
    return(0)
  }
  
  stop("Unknown fault type.")
}

generate_rsdt_series <- function(
    fault_type = "step",
    step_delta = 0.06,
    online_contamination = c("none", "pre_spikes"),
    gamma_online = 0,
    spike_index = c(480, 485, 490, 495, 498),
    seed = NULL
) {
  online_contamination <- match.arg(online_contamination)
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(0, nrow = k_final, ncol = 2)
  Y <- matrix(0, nrow = k_final, ncol = 2)
  f <- numeric(k_final)
  
  X[1, ] <- MASS::mvrnorm(1, mu = x0_mean, Sigma = P0)
  
  for (k in 1:k_final) {
    f[k] <- fault_signal(k, type = fault_type, step_delta = step_delta)
    
    vk <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = R)
    Y[k, ] <- as.numeric(C_mat %*% X[k, ] + vk + F_vec * f[k])
    
    if (k < k_final) {
      wk <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = Q)
      X[k + 1, ] <- as.numeric(A_fun(k) %*% X[k, ] + wk)
    }
  }
  
  if (online_contamination == "pre_spikes") {
    spike_index <- spike_index[spike_index > 1 & spike_index < k_fault]
    
    for (kk in spike_index) {
      Y[kk, ] <- Y[kk, ] + gamma_online * c(1, 0)
    }
  }
  
  residual <- Y - YBAR
  
  list(
    X = X,
    Y = Y,
    residual = residual,
    fault = f
  )
}

############################################################
## Estimators and thresholds
############################################################

safe_cov <- function(X) {
  S <- stats::cov(X)
  if (any(!is.finite(S))) {
    S <- diag(ncol(X))
  }
  S + 1e-8 * diag(ncol(X))
}

safe_inv_quad <- function(S, e, ridge = 1e-8) {
  S2 <- S + ridge * diag(ncol(S))
  as.numeric(t(e) %*% solve(S2, e))
}

rmcd_est <- function(W, alpha_mcd = 0.75) {
  fit <- tryCatch(
    robustbase::covMcd(W, alpha = alpha_mcd),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(center = colMeans(W), cov = safe_cov(W)))
  }
  
  list(
    center = fit$center,
    cov = fit$cov + 1e-8 * diag(ncol(W))
  )
}

adaptive_mcd_alpha <- function(
    w_roll,
    n_spikes = 5,
    c_factor = 1.25,
    lower = 0.50,
    upper = 0.95
) {
  eps_hat <- n_spikes / w_roll
  alpha_mcd <- 1 - c_factor * eps_hat
  alpha_mcd <- max(lower, min(upper, alpha_mcd))
  alpha_mcd
}

liu_T2_threshold <- function(p = 2, omega = 10, alpha = 0.05) {
  (p * omega / (omega - p + 1)) *
    qf(1 - alpha, df1 = p, df2 = omega - p + 1)
}

calibrate_classical_rolling_threshold <- function(
    B = 10000,
    w_roll = 50,
    alpha = 0.05,
    seed = 1234
) {
  set.seed(seed)
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    dat <- generate_rsdt_series(
      fault_type = "step",
      step_delta = 0.06,
      online_contamination = "none",
      seed = seed + b
    )
    
    Rclean <- dat$residual[1:(k_fault - 1), , drop = FALSE]
    t <- sample((w_roll + 1):(k_fault - 1), 1)
    
    W <- Rclean[(t - w_roll):(t - 1), , drop = FALSE]
    x_t <- Rclean[t, ]
    
    vals[b] <- safe_inv_quad(safe_cov(W), x_t - colMeans(W))
  }
  
  as.numeric(quantile(vals, probs = 1 - alpha, na.rm = TRUE))
}

calibrate_rrh_threshold <- function(
    B = 10000,
    w_roll = 50,
    alpha = 0.05,
    alpha_mcd_value = 0.75,
    seed = 5678
) {
  set.seed(seed)
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    dat <- generate_rsdt_series(
      fault_type = "step",
      step_delta = 0.06,
      online_contamination = "none",
      seed = seed + b
    )
    
    Rclean <- dat$residual[1:(k_fault - 1), , drop = FALSE]
    t <- sample((w_roll + 1):(k_fault - 1), 1)
    
    W <- Rclean[(t - w_roll):(t - 1), , drop = FALSE]
    x_t <- Rclean[t, ]
    
    est <- rmcd_est(W, alpha_mcd = alpha_mcd_value)
    vals[b] <- safe_inv_quad(est$cov, x_t - est$center)
  }
  
  as.numeric(quantile(vals, probs = 1 - alpha, na.rm = TRUE))
}

############################################################
## Scenario grid
############################################################

## Clean benchmark
clean_grid <- expand.grid(
  study = "clean",
  fault_type = c("step", "oscillation", "gradual", "intermittent"),
  step_delta = c(0.06),
  w_roll = c(50, 80),
  gamma_online = c(0),
  online_contamination = c("none"),
  offline_contamination = c("none"),
  eps_offline = c(0),
  stringsAsFactors = FALSE
)

## Weak step faults under contamination
step_cont_grid <- expand.grid(
  study = "contaminated_rolling",
  fault_type = c("step"),
  step_delta = c(0.03, 0.04, 0.06),
  w_roll = c(50, 80),
  gamma_online = c(1.5, 2.0),
  online_contamination = c("pre_spikes"),
  offline_contamination = c("none", "random"),
  eps_offline = c(0.05, 0.10),
  stringsAsFactors = FALSE
)

step_cont_grid <- step_cont_grid %>%
  filter(
    offline_contamination == "random" |
      (offline_contamination == "none" & eps_offline == 0.05)
  ) %>%
  mutate(
    eps_offline = ifelse(offline_contamination == "none", 0, eps_offline)
  )

## Gradual and intermittent under contamination
other_cont_grid <- expand.grid(
  study = "contaminated_rolling",
  fault_type = c("gradual", "intermittent", "oscillation"),
  step_delta = c(0.06),
  w_roll = c(50, 80),
  gamma_online = c(1.5, 2.0),
  online_contamination = c("pre_spikes"),
  offline_contamination = c("none", "random"),
  eps_offline = c(0.05, 0.10),
  stringsAsFactors = FALSE
)

other_cont_grid <- other_cont_grid %>%
  filter(
    offline_contamination == "random" |
      (offline_contamination == "none" & eps_offline == 0.05)
  ) %>%
  mutate(
    eps_offline = ifelse(offline_contamination == "none", 0, eps_offline)
  )

scenario_grid <- bind_rows(
  clean_grid,
  step_cont_grid,
  other_cont_grid
) %>%
  mutate(
    scenario_id = row_number()
  ) %>%
  select(
    scenario_id,
    everything()
  )

write_csv(
  scenario_grid,
  file.path(out_table_dir, "scenario_grid_rrhfd_main.csv")
)

cat("Number of scenarios:", nrow(scenario_grid), "\n")

############################################################
## Job grid: scenario x block
############################################################

n_blocks <- ceiling(R_REPS_MAIN / R_BLOCK)

job_grid <- expand.grid(
  scenario_id = scenario_grid$scenario_id,
  block_id = seq_len(n_blocks)
) %>%
  left_join(scenario_grid, by = "scenario_id") %>%
  arrange(scenario_id, block_id) %>%
  mutate(
    job_id = row_number(),
    rep_start = (block_id - 1) * R_BLOCK + 1,
    rep_end = pmin(block_id * R_BLOCK, R_REPS_MAIN)
  ) %>%
  select(job_id, scenario_id, block_id, rep_start, rep_end, everything())

write_csv(
  job_grid,
  file.path(out_table_dir, "job_grid_rrhfd_main.csv")
)

cat("Number of jobs:", nrow(job_grid), "\n")

############################################################
## Thresholds by w
############################################################

threshold_rows <- list()

for (ww in sort(unique(scenario_grid$w_roll))) {
  cat("\nCalibrating thresholds for w =", ww, "with B_CAL =", B_CAL, "\n")
  
  alpha_ad <- adaptive_mcd_alpha(
    w_roll = ww,
    n_spikes = N_SPIKES,
    c_factor = 1.25,
    lower = 0.50,
    upper = 0.95
  )
  
  c_liu <- liu_T2_threshold(p = p, omega = omega, alpha = alpha_level)
  
  c_class <- calibrate_classical_rolling_threshold(
    B = B_CAL,
    w_roll = ww,
    alpha = alpha_level,
    seed = 10000 + ww
  )
  
  c_rrh075 <- calibrate_rrh_threshold(
    B = B_CAL,
    w_roll = ww,
    alpha = alpha_level,
    alpha_mcd_value = 0.75,
    seed = 20000 + ww
  )
  
  c_rrh085 <- calibrate_rrh_threshold(
    B = B_CAL,
    w_roll = ww,
    alpha = alpha_level,
    alpha_mcd_value = 0.85,
    seed = 30000 + ww
  )
  
  c_rrh_ad <- calibrate_rrh_threshold(
    B = B_CAL,
    w_roll = ww,
    alpha = alpha_level,
    alpha_mcd_value = alpha_ad,
    seed = 40000 + ww
  )
  
  threshold_rows[[as.character(ww)]] <- data.frame(
    w_roll = ww,
    alpha_adaptive = alpha_ad,
    c_liu_T2 = c_liu,
    c_class_roll = c_class,
    c_rrh_075 = c_rrh075,
    c_rrh_085 = c_rrh085,
    c_rrh_adaptive = c_rrh_ad
  )
}

threshold_table <- bind_rows(threshold_rows)

write_csv(
  threshold_table,
  file.path(out_table_dir, "thresholds_rrhfd_main.csv")
)

print(threshold_table)

cat("\nGrid and thresholds prepared successfully.\n")