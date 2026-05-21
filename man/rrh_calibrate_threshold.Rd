\name{rrh_calibrate_threshold}
\alias{rrh_calibrate_threshold}
\title{Empirical RRH-FD Threshold Calibration}
\usage{
rrh_calibrate_threshold(null_residuals = NULL, simulate_fn = NULL,
  B = 1000, w_roll = 50, alpha_level = 0.05, n_spikes = 5,
  alpha_mcd = NULL, ridge = 1e-08, seed = NULL,
  use_robustbase = NULL)
}
\arguments{
\item{null_residuals}{Clean residual matrix. If supplied, random time points are
sampled and their preceding windows are used.}
\item{simulate_fn}{Optional function taking one integer \code{n} and returning an
\code{n} by \code{p} residual matrix.}
\item{B}{Number of calibration replicates.}
\item{w_roll}{Rolling window length.}
\item{alpha_level}{Target false alarm level.}
\item{n_spikes}{Expected number of spikes for adaptive alpha.}
\item{alpha_mcd}{Optional fixed MCD alpha. If \code{NULL}, adaptive alpha is used.}
\item{ridge}{Ridge added to covariance matrices.}
\item{seed}{Optional random seed.}
\item{use_robustbase}{Passed to \code{\link{rrh_cov_mcd}}.}
}
\value{A list with threshold, null statistics, and settings.}
\description{Calibrate an RRH-FD threshold under a clean/null residual stream.}
\examples{
clean <- simulate_residuals(200, 2, seed = 1)
rrh_calibrate_threshold(clean, B = 50, w_roll = 30)
}
