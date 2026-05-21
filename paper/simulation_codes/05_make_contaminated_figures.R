############################################################
# Figures for contaminated RRH-FD simulation results
# Output format: PNG
# Run from project root: RobHotellingFault/
############################################################

rm(list = ls())

required_pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "scales")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them first using install.packages()."
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

############################################################
# Paths
############################################################

tables_dir <- "results/tables/contaminated"
fig_dir    <- "results/figures/contaminated"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

overall_file <- file.path(tables_dir, "rrhfd_main_method_overall.csv")
summary_file <- file.path(tables_dir, "rrhfd_main_summary_results.csv")
compact_file <- file.path(tables_dir, "rrhfd_main_compact_comparison.csv")
all_file     <- file.path(tables_dir, "rrhfd_main_all_results.csv")

if (!file.exists(overall_file)) stop("File not found: ", overall_file)
if (!file.exists(summary_file)) stop("File not found: ", summary_file)
if (!file.exists(compact_file)) stop("File not found: ", compact_file)

overall <- read_csv(overall_file, show_col_types = FALSE)
summary <- read_csv(summary_file, show_col_types = FALSE)
compact <- read_csv(compact_file, show_col_types = FALSE)

############################################################
# Helper functions
############################################################

pick_col <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop("None of these columns were found: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

clean_method_labels <- function(x) {
  dplyr::recode(
    x,
    "Classical_Rolling" = "Classical rolling",
    "Liu_F_offline" = "Liu-type offline",
    "RRH_FD_075" = "RRH-FD 0.75",
    "RRH_FD_085" = "RRH-FD 0.85",
    "RRH_FD_adaptive" = "RRH-FD adaptive",
    .default = x
  )
}

method_order <- c(
  "Classical rolling",
  "Liu-type offline",
  "RRH-FD 0.75",
  "RRH-FD 0.85",
  "RRH-FD adaptive"
)

save_png <- function(plot, filename, width = 8, height = 5.2, dpi = 300) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

############################################################
# Standardize key columns
############################################################

td_col_overall <- pick_col(overall, c("td_mean", "td", "mean_td"))
pf_col_overall <- pick_col(overall, c("pf"))
pm_col_overall <- pick_col(overall, c("pm"))

dr25_col_overall <- pick_col(overall, c("det_H25_rate", "DR25", "DR_25", "det_H25"))
dr50_col_overall <- pick_col(overall, c("det_H50_rate", "DR50", "DR_50", "det_H50"))
dr100_col_overall <- pick_col(overall, c("det_H100_rate", "DR100", "DR_100", "det_H100"))

overall2 <- overall %>%
  mutate(
    Method = clean_method_labels(method),
    Method = factor(Method, levels = method_order),
    td_value = .data[[td_col_overall]],
    pf_value = .data[[pf_col_overall]],
    pm_value = .data[[pm_col_overall]],
    DR25 = .data[[dr25_col_overall]],
    DR50 = .data[[dr50_col_overall]],
    DR100 = .data[[dr100_col_overall]]
  ) %>%
  arrange(Method)

############################################################
# Figure 1: Overall average detection delay by method
############################################################

fig1 <- ggplot(overall2, aes(x = Method, y = td_value, fill = Method)) +
  geom_col(width = 0.72, alpha = 0.90, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.1f", td_value)),
    vjust = -0.35,
    size = 3.5
  ) +
  labs(
    x = NULL,
    y = expression("Average detection delay " * t[d])
  ) +
  expand_limits(y = max(overall2$td_value, na.rm = TRUE) * 1.12) +
  theme_paper()

save_png(fig1, "fig_contaminated_delay_by_method.png", width = 8.2, height = 5.2)

############################################################
# Figure 2: Early detection rates by method
############################################################

dr_long <- overall2 %>%
  select(Method, DR25, DR50, DR100) %>%
  pivot_longer(
    cols = starts_with("DR"),
    names_to = "Horizon",
    values_to = "DetectionRate"
  ) %>%
  mutate(
    Horizon = factor(
      Horizon,
      levels = c("DR25", "DR50", "DR100"),
      labels = c("DR[25]", "DR[50]", "DR[100]")
    )
  )

fig2 <- ggplot(dr_long, aes(x = Horizon, y = DetectionRate, group = Method, color = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_discrete(labels = scales::parse_format()) +
  labs(
    x = "Post-fault detection horizon",
    y = "Detection rate"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_png(fig2, "fig_contaminated_detection_rates.png", width = 8.2, height = 5.2)
############################################################
# Figure 3: Detection delay by fault type
############################################################

td_col_summary <- pick_col(summary, c("td_mean", "td", "mean_td"))
dr50_col_summary <- pick_col(summary, c("det_H50_rate", "DR50", "DR_50", "det_H50"))

summary_fault <- summary %>%
  mutate(
    Method = clean_method_labels(method),
    FaultType = str_to_title(fault_type),
    td_value = .data[[td_col_summary]],
    DR50 = .data[[dr50_col_summary]]
  ) %>%
  filter(Method %in% c("Classical rolling", "RRH-FD adaptive", "RRH-FD 0.85")) %>%
  group_by(FaultType, Method) %>%
  summarise(
    td_mean = mean(td_value, na.rm = TRUE),
    DR50_mean = mean(DR50, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("Classical rolling", "RRH-FD adaptive", "RRH-FD 0.85")),
    FaultType = factor(FaultType, levels = c("Step", "Oscillation", "Gradual", "Intermittent"))
  )

fig3 <- ggplot(summary_fault, aes(x = FaultType, y = td_mean, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, alpha = 0.90) +
  labs(
    x = "Fault type",
    y = expression("Average detection delay " * t[d])
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_png(fig3, "fig_contaminated_fault_type_delay.png", width = 8.4, height = 5.2)

############################################################
# Figure 4: Boundary inflation vs detection-delay improvement
############################################################

logbir_col <- pick_col(
  compact,
  c(
    "mean_logBIR_Classical_Rolling",
    "mean_logBIR_RRH_FD_adaptive",
    "mean_logBIR",
    "logBIR"
  )
)

delta_td_col <- pick_col(
  compact,
  c(
    "delta_td_adaptive_minus_classical",
    "delta_td_RRH_adaptive_minus_Classical",
    "delta_td"
  )
)

compact2 <- compact %>%
  mutate(
    logBIR = .data[[logbir_col]],
    delay_improvement = -.data[[delta_td_col]],
    w_roll = factor(w_roll),
    gamma_online = factor(gamma_online)
  )

fig4 <- ggplot(
  compact2,
  aes(
    x = logBIR,
    y = delay_improvement,
    color = gamma_online,
    shape = w_roll
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  labs(
    x = expression(log(BIR)),
    y = expression("Delay improvement over classical rolling " * (-Delta * t[d])),
    color = expression(gamma[online]),
    shape = "Window length"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_png(fig4, "fig_contaminated_bir_improvement.png", width = 8.2, height = 5.4)

############################################################
# Figure 5: Grouped boundary inflation summary by w and gamma
############################################################

bir_group <- compact2 %>%
  group_by(w_roll, gamma_online) %>%
  summarise(
    logBIR_mean = mean(logBIR, na.rm = TRUE),
    delay_improvement_mean = mean(delay_improvement, na.rm = TRUE),
    DR25_improvement = mean(delta_H25_adaptive_minus_classical, na.rm = TRUE),
    DR50_improvement = mean(delta_H50_adaptive_minus_classical, na.rm = TRUE),
    .groups = "drop"
  )

fig5 <- ggplot(
  bir_group,
  aes(x = gamma_online, y = logBIR_mean, fill = w_roll)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, alpha = 0.90) +
  geom_text(
    aes(label = sprintf("%.2f", logBIR_mean)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.4
  ) +
  labs(
    x = expression(gamma[online]),
    y = expression("Average " * log(BIR)),
    fill = "Window length"
  ) +
  expand_limits(y = max(bir_group$logBIR_mean, na.rm = TRUE) * 1.15) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_png(fig5, "fig_contaminated_boundary_inflation_groups.png", width = 7.6, height = 5.2)

############################################################
# Optional Figure 6: Distribution of detection delay from all replications
# This requires rrhfd_main_all_results.csv
############################################################

if (file.exists(all_file)) {
  all_results <- read_csv(all_file, show_col_types = FALSE)
  
  td_col_all <- pick_col(all_results, c("td", "td_mean"))
  detected_col_all <- pick_col(all_results, c("detected", "detection_rate"))
  
  all_plot <- all_results %>%
    mutate(
      Method = clean_method_labels(method),
      td_value = .data[[td_col_all]],
      detected_value = .data[[detected_col_all]]
    ) %>%
    filter(
      Method %in% method_order,
      detected_value == 1,
      is.finite(td_value)
    ) %>%
    mutate(Method = factor(Method, levels = method_order))
  
  fig6 <- ggplot(all_plot, aes(x = Method, y = td_value, fill = Method)) +
    geom_boxplot(outlier.alpha = 0.15, width = 0.65, show.legend = FALSE) +
    labs(
      x = NULL,
      y = expression("Detection delay " * t[d])
    ) +
    theme_paper()
  
  save_png(fig6, "fig_contaminated_delay_distribution_boxplot.png", width = 8.5, height = 5.4)
}

############################################################
# Save figure-ready summaries
############################################################

write_csv(overall2, file.path(fig_dir, "plot_data_overall_method.csv"))
write_csv(dr_long, file.path(fig_dir, "plot_data_detection_rates.csv"))
write_csv(summary_fault, file.path(fig_dir, "plot_data_fault_type_delay.csv"))
write_csv(compact2, file.path(fig_dir, "plot_data_bir_improvement.csv"))
write_csv(bir_group, file.path(fig_dir, "plot_data_boundary_inflation_groups.csv"))

cat("\nFigures saved to:\n")
cat(normalizePath(fig_dir), "\n\n")
cat("Created PNG files:\n")
print(list.files(fig_dir, pattern = "\\.png$", full.names = FALSE))