\name{rrh_cov_mcd}
\alias{rrh_cov_mcd}
\title{Robust MCD Covariance Estimate}
\usage{
rrh_cov_mcd(x, alpha = 0.75, ridge = 1e-08, use_robustbase = NULL)
}
\arguments{
\item{x}{Numeric matrix whose rows are residual observations.}
\item{alpha}{MCD alpha value.}
\item{ridge}{Ridge added to covariance matrices.}
\item{use_robustbase}{If \code{TRUE}, require \code{robustbase::covMcd()}. If
\code{FALSE}, always use the trimmed fallback. If \code{NULL}, use
\code{robustbase} when installed.}
}
\value{A list with \code{center}, \code{cov}, \code{alpha}, \code{h},
\code{subset}, and \code{method}.}
\description{Estimate a robust center and covariance for a rolling residual window.}
\examples{
x <- simulate_residuals(80, 2, seed = 1)
rrh_cov_mcd(x, alpha = 0.85)
}
