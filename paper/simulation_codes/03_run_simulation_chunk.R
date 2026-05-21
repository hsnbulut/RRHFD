############################################################
## run_rrhfd_main_chunk.R
## Runs one row of job_grid_rrhfd_main.csv on TRUBA
############################################################

rm(list = ls())

req_pkgs <- c("MASS", "robustbase", "dplyr", "readr", "parallel")

for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(MASS)
library(robustbase)
library(dplyr)
library(readr)
library(parallel)

############################################################
## Environment variables
############################################################

JOB_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("JOB_ID", "1")))

GRID_FILE <- Sys.getenv("GRID_FILE", "results/tables/job_grid_rrhfd_main.csv")
THRESHOLD_FILE <- Sys.getenv("THRESHOLD_FILE", "results/tables/thresholds_rrhfd_main.csv")
OUTDIR <- Sys.getenv("OUTDIR", "results/raw_chunks/rrhfd_main")

BASE_SEED <- as.integer(Sys.getenv("BASE_SEED", "20260516"))
MC_CORES <- as.integer(Sys.getenv("MC_CORES", "1"))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

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

spike_index_main <- c(480, 485, 490, 495, 498)

############################################################
## RSDTS functions
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

generate_offline_residuals <- function(
    omega = 10,
    offline_contamination = c("none", "random"),
    gamma_offline = 1.5,
    eps_random = 0.05,
    seed = NULL
) {
  offline_contamination <- match.arg(offline_contamination)
  
  if (!is.null(seed)) set.seed(seed)
  
  Rarr <- array(NA_real_, dim = c(k_final, p, omega))
  
  for (s in 1:omega) {
    dat <- generate_rsdt_series(
      fault_type = "step",
      step_delta = 0.06,
      online_contamination = "none",
      gamma_online = 0,
      seed = ifelse(is.null(seed), NULL, seed + s)
    )
    
    res <- dat$residual
    
    ## Remove artificial step fault effect to approximate fault-free offline residuals
    for (k in k_fault:k_final) {
      res[k, ] <- res[k, ] - F_vec * fault_signal(k, "step", step_delta = 0.06)
    }
    
    if (offline_contamination == "random") {
      pre_pool <- 50:(k_fault - 1)
      n_cont <- max(1, floor(eps_random * length(pre_pool)))
      idx <- sample(pre_pool, size = n_cont, replace = FALSE)
      
      for (kk in idx) {
        res[kk, ] <- res[kk, ] + gamma_offline * c(1, 0)
      }
    }
    
    Rarr[, , s] <- res
  }
  
  Rarr
}

############################################################
## Estimators and metrics
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

logdet_safe <- function(S) {
  as.numeric(
    determinant(
      S + 1e-8 * diag(ncol(S)),
      logarithm = TRUE
    )$modulus[1]
  )
}

run_detectors_one <- function(
    residual,
    offline_Rarr,
    thresholds,
    w_roll
) {
  n <- nrow(residual)
  
  methods <- c(
    "Liu_F_offline",
    "Classical_Rolling",
    "RRH_FD_075",
    "RRH_FD_085",
    "RRH_FD_adaptive"
  )
  
  stat <- matrix(NA_real_, nrow = n, ncol = length(methods))
  colnames(stat) <- methods
  
  alarm <- matrix(0L, nrow = n, ncol = length(methods))
  colnames(alarm) <- methods
  
  logBIR_075 <- rep(NA_real_, n)
  logBIR_085 <- rep(NA_real_, n)
  logBIR_adaptive <- rep(NA_real_, n)
  
  alpha_ad <- thresholds$alpha_adaptive
  
  for (k in (w_roll + 1):n) {
    r_k <- residual[k, ]
    
    ## Liu-type offline Hotelling
    Rk_off <- t(offline_Rarr[k, , ])
    S_off <- (t(Rk_off) %*% Rk_off) / nrow(Rk_off)
    
    T2_liu <- safe_inv_quad(S_off, r_k)
    stat[k, "Liu_F_offline"] <- T2_liu
    alarm[k, "Liu_F_offline"] <- as.integer(T2_liu > thresholds$c_liu_T2)
    
    ## Classical rolling Hotelling
    W <- residual[(k - w_roll):(k - 1), , drop = FALSE]
    mu_C <- colMeans(W)
    S_C <- safe_cov(W)
    
    T2_C <- safe_inv_quad(S_C, r_k - mu_C)
    stat[k, "Classical_Rolling"] <- T2_C
    alarm[k, "Classical_Rolling"] <- as.integer(T2_C > thresholds$c_class_roll)
    
    ## RRH 0.75
    est_075 <- rmcd_est(W, alpha_mcd = 0.75)
    T2_075 <- safe_inv_quad(est_075$cov, r_k - est_075$center)
    stat[k, "RRH_FD_075"] <- T2_075
    alarm[k, "RRH_FD_075"] <- as.integer(T2_075 > thresholds$c_rrh_075)
    
    ## RRH 0.85
    est_085 <- rmcd_est(W, alpha_mcd = 0.85)
    T2_085 <- safe_inv_quad(est_085$cov, r_k - est_085$center)
    stat[k, "RRH_FD_085"] <- T2_085
    alarm[k, "RRH_FD_085"] <- as.integer(T2_085 > thresholds$c_rrh_085)
    
    ## RRH adaptive
    est_ad <- rmcd_est(W, alpha_mcd = alpha_ad)
    T2_ad <- safe_inv_quad(est_ad$cov, r_k - est_ad$center)
    stat[k, "RRH_FD_adaptive"] <- T2_ad
    alarm[k, "RRH_FD_adaptive"] <- as.integer(T2_ad > thresholds$c_rrh_adaptive)
    
    logBIR_075[k] <- logdet_safe(S_C) - logdet_safe(est_075$cov)
    logBIR_085[k] <- logdet_safe(S_C) - logdet_safe(est_085$cov)
    logBIR_adaptive[k] <- logdet_safe(S_C) - logdet_safe(est_ad$cov)
  }
  
  list(
    stat = stat,
    alarm = alarm,
    logBIR_075 = logBIR_075,
    logBIR_085 = logBIR_085,
    logBIR_adaptive = logBIR_adaptive
  )
}

compute_metrics <- function(
    det,
    w_roll,
    k_fault = 500,
    k_final = 2000
) {
  alarm <- det$alarm
  methods <- colnames(alarm)
  
  pre_idx <- (w_roll + 1):(k_fault - 1)
  post_idx <- k_fault:k_final
  bir_idx <- (k_fault - 40):(k_fault + 80)
  
  out <- list()
  
  for (m in methods) {
    A <- alarm[, m]
    
    nf <- sum(A[pre_idx] == 1, na.rm = TRUE)
    nm <- sum(A[post_idx] == 0, na.rm = TRUE)
    
    pf <- nf / length(pre_idx)
    pm <- nm / length(post_idx)
    
    det_times <- post_idx[which(A[post_idx] == 1)]
    
    if (length(det_times) == 0) {
      td <- NA_real_
      detected <- 0L
    } else {
      td <- min(det_times) - k_fault
      detected <- 1L
    }
    
    det_H25 <- as.integer(any(A[k_fault:min(k_final, k_fault + 25)] == 1, na.rm = TRUE))
    det_H50 <- as.integer(any(A[k_fault:min(k_final, k_fault + 50)] == 1, na.rm = TRUE))
    det_H100 <- as.integer(any(A[k_fault:min(k_final, k_fault + 100)] == 1, na.rm = TRUE))
    
    if (m == "RRH_FD_075") {
      bir_mean <- mean(det$logBIR_075[bir_idx], na.rm = TRUE)
      bir_max <- max(det$logBIR_075[bir_idx], na.rm = TRUE)
    } else if (m == "RRH_FD_085") {
      bir_mean <- mean(det$logBIR_085[bir_idx], na.rm = TRUE)
      bir_max <- max(det$logBIR_085[bir_idx], na.rm = TRUE)
    } else {
      bir_mean <- mean(det$logBIR_adaptive[bir_idx], na.rm = TRUE)
      bir_max <- max(det$logBIR_adaptive[bir_idx], na.rm = TRUE)
    }
    
    out[[m]] <- data.frame(
      method = m,
      pf = pf,
      pm = pm,
      td = td,
      detected = detected,
      det_H25 = det_H25,
      det_H50 = det_H50,
      det_H100 = det_H100,
      mean_logBIR = bir_mean,
      max_logBIR = bir_max
    )
  }
  
  do.call(rbind, out)
}

run_one_rep <- function(rep_id, job_row, thresholds, offline_Rarr) {
  seed_i <- BASE_SEED + job_row$scenario_id * 100000 + rep_id
  
  dat <- generate_rsdt_series(
    fault_type = job_row$fault_type,
    step_delta = job_row$step_delta,
    online_contamination = job_row$online_contamination,
    gamma_online = job_row$gamma_online,
    spike_index = spike_index_main,
    seed = seed_i
  )
  
  det <- run_detectors_one(
    residual = dat$residual,
    offline_Rarr = offline_Rarr,
    thresholds = thresholds,
    w_roll = job_row$w_roll
  )
  
  met <- compute_metrics(
    det = det,
    w_roll = job_row$w_roll,
    k_fault = k_fault,
    k_final = k_final
  )
  
  met$rep_id <- rep_id
  met
}

############################################################
## Read job row and thresholds
############################################################

job_grid <- read_csv(GRID_FILE, show_col_types = FALSE)
threshold_table <- read_csv(THRESHOLD_FILE, show_col_types = FALSE)

if (JOB_ID < 1 || JOB_ID > nrow(job_grid)) {
  stop("Invalid JOB_ID / SLURM_ARRAY_TASK_ID.")
}

job_row <- job_grid[JOB_ID, ]

thresholds_df <- threshold_table %>%
  filter(w_roll == job_row$w_roll)

if (nrow(thresholds_df) != 1) {
  stop("Threshold row not found.")
}

thresholds <- as.list(thresholds_df[1, ])

cat("Running JOB_ID =", JOB_ID, "\n")
print(job_row)
print(thresholds)

############################################################
## Generate offline data once for this scenario/block
############################################################

offline_seed <- BASE_SEED + job_row$scenario_id * 999

offline_Rarr <- generate_offline_residuals(
  omega = omega,
  offline_contamination = job_row$offline_contamination,
  gamma_offline = job_row$gamma_online,
  eps_random = job_row$eps_offline,
  seed = offline_seed
)

rep_ids <- job_row$rep_start:job_row$rep_end

cat("Replicates:", min(rep_ids), "to", max(rep_ids), "\n")
cat("MC_CORES:", MC_CORES, "\n")

if (MC_CORES > 1) {
  res_list <- parallel::mclapply(
    rep_ids,
    run_one_rep,
    job_row = job_row,
    thresholds = thresholds,
    offline_Rarr = offline_Rarr,
    mc.cores = MC_CORES
  )
} else {
  res_list <- lapply(
    rep_ids,
    run_one_rep,
    job_row = job_row,
    thresholds = thresholds,
    offline_Rarr = offline_Rarr
  )
}

res <- bind_rows(res_list)

for (nm in names(job_row)) {
  res[[nm]] <- job_row[[nm]]
}

res$alpha_adaptive <- thresholds$alpha_adaptive
res$c_liu_T2 <- thresholds$c_liu_T2
res$c_class_roll <- thresholds$c_class_roll
res$c_rrh_075 <- thresholds$c_rrh_075
res$c_rrh_085 <- thresholds$c_rrh_085
res$c_rrh_adaptive <- thresholds$c_rrh_adaptive

outfile <- file.path(
  OUTDIR,
  sprintf(
    "rrhfd_main_s%03d_b%02d_job%04d.rds",
    job_row$scenario_id,
    job_row$block_id,
    JOB_ID
  )
)

saveRDS(res, outfile)

cat("Saved:", outfile, "\n")