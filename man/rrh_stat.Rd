\name{rrh_stat}
\alias{rrh_stat}
\title{One-Step RRH-FD Statistic}
\usage{
rrh_stat(rk, window, alpha_mcd = 0.75, ridge = 1e-08,
  use_robustbase = NULL)
}
\arguments{
\item{rk}{Current residual vector.}
\item{window}{Historical residual window, one row per time point.}
\item{alpha_mcd}{MCD alpha.}
\item{ridge}{Ridge added to covariance matrices.}
\item{use_robustbase}{Passed to \code{\link{rrh_cov_mcd}}.}
}
\value{A list with statistic, center, covariance, alpha, and estimator method.}
\description{Compute the robust rolling Hotelling statistic for one residual vector.}
