% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_utils.r
\name{lambda}
\alias{lambda}
\title{Lambda
 
Seasonal model for the intensity function, \eqn{\lambda(t) = a + b + b \sin(c\,2 \pi/365 + x 2 \,\pi / 365)}.}
\usage{
lambda(x, params, names = F)
}
\arguments{
\item{x}{\code{numeric} vector of time points.}

\item{params}{\code{numeric} vector of parameters of length = 4.}

\item{names}{A vector. Can be used to label}
}
\value{
a \code{numeric} vector of the same length as \code{x}.
}
\description{
Lambda
 
Seasonal model for the intensity function, \eqn{\lambda(t) = a + b + b \sin(c\,2 \pi/365 + x 2 \,\pi / 365)}.
}
\examples{
lambda(c(1,2,3,4,5,6), c(1,1,1,1) )
}
