\name{simulate_residuals}
\alias{simulate_residuals}
\title{Simulate Residuals}
\usage{
simulate_residuals(n, p = 2, sigma = diag(p), fault_time = NULL,
  fault_mean = rep(0, p), contamination = NULL, seed = NULL)
}
\arguments{
\item{n}{Number of observations.}
\item{p}{Residual dimension.}
\item{sigma}{Covariance matrix. Defaults to identity.}
\item{fault_time}{Optional fault start time.}
\item{fault_mean}{Mean shift after \code{fault_time}.}
\item{contamination}{Optional list with \code{index} and \code{shift} entries.}
\item{seed}{Optional random seed.}
}
\value{An \code{n} by \code{p} residual matrix.}
\description{Generate a simple multivariate residual stream for examples and calibration.}
