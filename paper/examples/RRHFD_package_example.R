############################################################
# Simple RRHFD package application for the manuscript
############################################################

rm(list = ls())

############################################################
# 1. Load package
############################################################

# If the package has not been installed yet, use one of the following:
# install.packages("devtools")
# devtools::install_github("hsnbulut/RRHFD")
#
# or, if using a local zip/source package:
# install.packages("RRHFD.zip", repos = NULL, type = "source")

library(RRHFD)

############################################################
# 2. Settings
############################################################

set.seed(2026)

n          <- 300
p          <- 2
k0         <- 150
w_roll     <- 50
alpha_lvl  <- 0.05
n_spikes   <- 6
B_cal      <- 500

out_dir <- "results/software_example"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Generate clean calibration residuals
############################################################

null_residuals <- simulate_residuals(
  n = 1000,
  p = p,
  sigma = diag(p),
  seed = 1001
)

############################################################
# 4. Calibrate RRHFD threshold
############################################################

cal <- rrh_calibrate_threshold(
  null_residuals = null_residuals,
  w_roll = w_roll,
  B = B_cal,
  alpha_level = alpha_lvl,
  n_spikes = n_spikes,
  seed = 2026,
  use_robustbase = TRUE
)

############################################################
# 5. Generate online residual stream
############################################################

# Transient pre-fault casewise contamination
spike_index <- 120:125
spike_shift <- c(3.5, 2.5)

# Persistent post-fault mean shift
fault_shift <- c(1.0, 0.6)

residual_stream <- simulate_residuals(
  n = n,
  p = p,
  sigma = diag(p),
  fault_time = k0,
  fault_mean = fault_shift,
  contamination = list(
    index = spike_index,
    shift = spike_shift
  ),
  seed = 2027
)

############################################################
# 6. Apply RRHFD detector
############################################################

det <- rrh_detect(
  residuals = residual_stream,
  w_roll = w_roll,
  threshold = cal$threshold,
  n_spikes = n_spikes,
  use_robustbase = TRUE
)

############################################################
# 7. Compute detection metrics
############################################################

met <- rrh_metrics(
  detections = det,
  fault_time = k0,
  horizon = 50
)

first_alarm_time <- if (met$detected == 1) {
  k0 + met$delay
} else {
  NA_integer_
}

summary_table <- data.frame(
  window_length = w_roll,
  alpha_level = alpha_lvl,
  alpha_mcd = cal$alpha_mcd,
  threshold = cal$threshold,
  fault_time = k0,
  first_alarm_time = first_alarm_time,
  detection_delay = met$delay,
  false_alarm_rate = met$pf,
  missed_detection_rate = met$pm,
  detected_within_50 = met$detected_in_horizon,
  mean_logBIR = mean(det$logBIR, na.rm = TRUE)
)

print(summary_table)

write.csv(
  summary_table,
  file.path(out_dir, "rrhfd_software_example_summary.csv"),
  row.names = FALSE
)

write.csv(
  det,
  file.path(out_dir, "rrhfd_software_example_detection_path.csv"),
  row.names = FALSE
)

############################################################
# 8. Plot detection path
############################################################

png(
  filename = file.path(out_dir, "rrhfd_software_example_detection_path.png"),
  width = 1800,
  height = 1100,
  res = 220
)

plot(
  det$time,
  det$stat_rrh,
  type = "l",
  lwd = 1.4,
  xlab = "Time index",
  ylab = expression("RRHFD statistic " * T[R]^2),
  main = ""
)

abline(
  h = cal$threshold,
  lty = 2,
  lwd = 1.3,
  col = "red"
)

abline(
  v = k0,
  lty = 3,
  lwd = 1.3,
  col = "black"
)

rect(
  xleft = min(spike_index),
  ybottom = par("usr")[3],
  xright = max(spike_index),
  ytop = par("usr")[4],
  col = rgb(1, 0.65, 0, 0.15),
  border = NA
)

lines(
  det$time,
  det$stat_rrh,
  lwd = 1.4
)

if (!is.na(first_alarm_time)) {
  points(
    first_alarm_time,
    det$stat_rrh[det$time == first_alarm_time],
    pch = 21,
    bg = "yellow",
    col = "black",
    cex = 1.4
  )
}

legend(
  "topright",
  legend = c(
    "RRHFD statistic",
    "Threshold",
    "Fault onset",
    "Transient pre-fault spikes",
    "First post-fault alarm"
  ),
  lty = c(1, 2, 3, NA, NA),
  pch = c(NA, NA, NA, 15, 21),
  pt.bg = c(NA, NA, NA, rgb(1, 0.65, 0, 0.15), "yellow"),
  col = c("black", "red", "black", rgb(1, 0.65, 0, 0.35), "black"),
  bty = "n",
  cex = 0.85
)

dev.off()

############################################################
# 9. Console message
############################################################

cat("\nFiles created:\n")
cat(file.path(out_dir, "rrhfd_software_example_summary.csv"), "\n")
cat(file.path(out_dir, "rrhfd_software_example_detection_path.csv"), "\n")
cat(file.path(out_dir, "rrhfd_software_example_detection_path.png"), "\n")