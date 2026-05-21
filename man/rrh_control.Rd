\name{rrh_control}
\alias{rrh_control}
\title{RRH-FD Control Settings}
\usage{
rrh_control(w_roll = 50, alpha_level = 0.05, n_spikes = 5,
  c_factor = 1.25, alpha_mcd = NULL, alpha_lower = 0.5,
  alpha_upper = 0.95, ridge = 1e-08, use_robustbase = NULL)
}
\arguments{
\item{w_roll}{Rolling window length.}
\item{alpha_level}{Target false alarm level.}
\item{n_spikes}{Expected number of pre-fault spikes inside the window.}
\item{c_factor}{Inflation factor for adaptive alpha.}
\item{alpha_mcd}{Optional fixed MCD alpha. If \code{NULL}, adaptive alpha is used.}
\item{alpha_lower}{Lower bound for adaptive alpha.}
\item{alpha_upper}{Upper bound for adaptive alpha.}
\item{ridge}{Ridge added to covariance matrices.}
\item{use_robustbase}{Use \code{robustbase::covMcd()} when available.}
}
\value{A list of control settings.}
\description{Build a small list of detector settings.}
