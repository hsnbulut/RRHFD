############################################################
## combine_rrhfd_main_results.R
############################################################

rm(list = ls())

req_pkgs <- c("dplyr", "readr", "tidyr")

for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(dplyr)
library(readr)
library(tidyr)

INDIR <- Sys.getenv("INDIR", "results/raw_chunks/rrhfd_main")
OUTDIR <- Sys.getenv("OUTDIR_TABLES", "results/tables")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(
  INDIR,
  pattern = "\\.rds$",
  full.names = TRUE
)

cat("Number of RDS files:", length(files), "\n")

if (length(files) == 0) {
  stop("No RDS files found.")
}

all_res <- bind_rows(lapply(files, readRDS))

write_csv(
  all_res,
  file.path(OUTDIR, "rrhfd_main_all_results.csv")
)

summary_res <- all_res %>%
  group_by(
    scenario_id,
    study,
    fault_type,
    step_delta,
    w_roll,
    gamma_online,
    online_contamination,
    offline_contamination,
    eps_offline,
    alpha_adaptive,
    method
  ) %>%
  summarise(
    n_reps = n_distinct(rep_id),
    pf = mean(pf, na.rm = TRUE),
    pm = mean(pm, na.rm = TRUE),
    td_mean = mean(td, na.rm = TRUE),
    td_median = median(td, na.rm = TRUE),
    detection_rate = mean(detected, na.rm = TRUE),
    det_H25_rate = mean(det_H25, na.rm = TRUE),
    det_H50_rate = mean(det_H50, na.rm = TRUE),
    det_H100_rate = mean(det_H100, na.rm = TRUE),
    mean_logBIR = mean(mean_logBIR, na.rm = TRUE),
    max_logBIR = mean(max_logBIR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario_id, method)

write_csv(
  summary_res,
  file.path(OUTDIR, "rrhfd_main_summary_results.csv")
)

compact <- summary_res %>%
  select(
    scenario_id,
    study,
    fault_type,
    step_delta,
    w_roll,
    gamma_online,
    online_contamination,
    offline_contamination,
    eps_offline,
    alpha_adaptive,
    method,
    pf,
    pm,
    td_mean,
    det_H25_rate,
    det_H50_rate,
    det_H100_rate,
    mean_logBIR,
    max_logBIR
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = c(
      pf,
      pm,
      td_mean,
      det_H25_rate,
      det_H50_rate,
      det_H100_rate,
      mean_logBIR,
      max_logBIR
    )
  ) %>%
  mutate(
    delta_td_adaptive_minus_classical =
      td_mean_RRH_FD_adaptive - td_mean_Classical_Rolling,
    delta_pm_adaptive_minus_classical =
      pm_RRH_FD_adaptive - pm_Classical_Rolling,
    delta_H25_adaptive_minus_classical =
      det_H25_rate_RRH_FD_adaptive - det_H25_rate_Classical_Rolling,
    delta_H50_adaptive_minus_classical =
      det_H50_rate_RRH_FD_adaptive - det_H50_rate_Classical_Rolling,
    delta_H100_adaptive_minus_classical =
      det_H100_rate_RRH_FD_adaptive - det_H100_rate_Classical_Rolling,
    
    delta_td_adaptive_minus_liu =
      td_mean_RRH_FD_adaptive - td_mean_Liu_F_offline,
    delta_pm_adaptive_minus_liu =
      pm_RRH_FD_adaptive - pm_Liu_F_offline,
    delta_H50_adaptive_minus_liu =
      det_H50_rate_RRH_FD_adaptive - det_H50_rate_Liu_F_offline
  ) %>%
  arrange(scenario_id)

write_csv(
  compact,
  file.path(OUTDIR, "rrhfd_main_compact_comparison.csv")
)

method_overall <- summary_res %>%
  group_by(study, method) %>%
  summarise(
    pf = mean(pf, na.rm = TRUE),
    pm = mean(pm, na.rm = TRUE),
    td_mean = mean(td_mean, na.rm = TRUE),
    det_H25_rate = mean(det_H25_rate, na.rm = TRUE),
    det_H50_rate = mean(det_H50_rate, na.rm = TRUE),
    det_H100_rate = mean(det_H100_rate, na.rm = TRUE),
    mean_logBIR = mean(mean_logBIR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(study, method)

write_csv(
  method_overall,
  file.path(OUTDIR, "rrhfd_main_method_overall.csv")
)

cat("Combined results saved in:", OUTDIR, "\n")