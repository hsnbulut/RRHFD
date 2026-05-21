# RRHFD

`RRHFD` is a small R package for the adaptive robust rolling Hotelling T-squared
fault-detection method used in the RRH-FD stochastic monitoring project.

The package is intentionally compact for the first version. It provides:

- adaptive MCD alpha selection,
- robust rolling center and covariance estimation,
- Hotelling T-squared statistics,
- empirical threshold calibration under a clean/null residual stream,
- online alarm generation,
- simple false-alarm, missed-detection, delay, and horizon metrics.

## Minimal Example

```r
library(RRHFD)

set.seed(1)
clean <- simulate_residuals(n = 600, p = 2)

thr <- rrh_calibrate_threshold(
  null_residuals = clean,
  w_roll = 50,
  B = 300,
  alpha_level = 0.05,
  seed = 99
)

stream <- simulate_residuals(
  n = 700,
  p = 2,
  fault_time = 350,
  fault_mean = c(0.5, 0)
)

det <- rrh_detect(
  residuals = stream,
  w_roll = 50,
  threshold = thr$threshold,
  n_spikes = 5
)

rrh_metrics(det, fault_time = 350)
```

If `robustbase` is installed, `rrh_cov_mcd()` uses `robustbase::covMcd()`.
Otherwise it uses a lightweight trimmed robust covariance fallback so the package
can run without extra dependencies.
