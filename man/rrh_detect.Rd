\name{rrh_detect}
\alias{rrh_detect}
\title{Run Online RRH-FD Detection}
\usage{
rrh_detect(residuals, w_roll = 50, threshold, n_spikes = 5,
  alpha_mcd = NULL, ridge = 1e-08, use_robustbase = NULL)
}
\arguments{
\item{residuals}{Numeric residual matrix, one row per time point.}
\item{w_roll}{Rolling window length.}
\item{threshold}{Alarm threshold.}
\item{n_spikes}{Expected number of spikes for adaptive alpha.}
\item{alpha_mcd}{Optional fixed MCD alpha. If \code{NULL}, adaptive alpha is used.}
\item{ridge}{Ridge added to covariance matrices.}
\item{use_robustbase}{Passed to \code{\link{rrh_cov_mcd}}.}
}
\value{A data frame with \code{time}, \code{stat_rrh}, \code{alarm},
\code{alpha_mcd}, \code{logBIR}, and estimator \code{method}.}
\description{Compute rolling robust Hotelling statistics and alarms for a residual matrix.}
\examples{
clean <- simulate_residuals(200, 2, seed = 1)
cal <- rrh_calibrate_threshold(clean, B = 50, w_roll = 30)
det <- rrh_detect(clean, w_roll = 30, threshold = cal$threshold)
head(det)
}
