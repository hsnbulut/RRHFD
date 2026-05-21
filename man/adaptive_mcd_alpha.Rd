\name{adaptive_mcd_alpha}
\alias{adaptive_mcd_alpha}
\title{Adaptive MCD Alpha}
\usage{
adaptive_mcd_alpha(w_roll, n_spikes = 5, c_factor = 1.25,
  lower = 0.5, upper = 0.95)
}
\arguments{
\item{w_roll}{Rolling window length.}
\item{n_spikes}{Expected number of pre-fault spikes inside the window.}
\item{c_factor}{Inflation factor applied to the estimated contamination rate.}
\item{lower}{Lower bound for alpha.}
\item{upper}{Upper bound for alpha.}
}
\value{A numeric alpha value between \code{lower} and \code{upper}.}
\description{Compute the adaptive MCD alpha used by the RRH-FD detector.}
\examples{
adaptive_mcd_alpha(50, n_spikes = 5)
}
