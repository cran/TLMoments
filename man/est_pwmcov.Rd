% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_pwmcov.R
\name{est_pwmcov}
\alias{est_pwmcov}
\alias{est_pwmcov.numeric}
\alias{est_pwmcov.matrix}
\title{Estimate the covariance matrix of PWM estimations}
\usage{
est_pwmcov(x, order = 0:3, distr = NULL, distr.trim = c(0, 0))

\method{est_pwmcov}{numeric}(x, order = 0:3, distr = NULL, distr.trim = c(0, 0))

\method{est_pwmcov}{matrix}(x, order = 0:3, distr = NULL, distr.trim = c(0, 0))
}
\arguments{
\item{x}{numeric vector or matrix of data.}

\item{order}{numeric vector giving the orders that are returned.}

\item{distr}{character of length 1 which indicates a distribution if a
parametric assumption should be used.}

\item{distr.trim}{integer vector of length 2 indicating the trimming used
to calculate parameters if a parametric assumption is used (i.e. \code{distr} is set).}
}
\value{
numeric matrix
}
\description{
Internal function. Use \link{est_cov}. Description not done yet.
}
\examples{
### Numeric vectors
x <- rgev(500, shape = .2)
est_pwmcov(x)
est_pwmcov(x, distr = "gev")

### Numeric matrices
x <- matrix(rgev(600, shape = .2), nc = 3)
est_pwmcov(x, order = 0:2)
est_pwmcov(x, order = 0:2, distr = "gev")

}
