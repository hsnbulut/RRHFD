\name{rrh_metrics}
\alias{rrh_metrics}
\title{Fault Detection Metrics}
\usage{
rrh_metrics(detections, fault_time, horizon = 50)
}
\arguments{
\item{detections}{Data frame returned by \code{\link{rrh_detect}}.}
\item{fault_time}{Fault start time.}
\item{horizon}{Optional horizon for early detection rate.}
}
\value{A one-row data frame with false alarm rate, missed detection rate,
detection delay, detection flag, and optional horizon detection flag.}
\description{Compute simple pre/post fault metrics from an RRH-FD detection result.}
