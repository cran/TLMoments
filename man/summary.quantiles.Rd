% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.quantiles.R
\name{summary.quantiles}
\alias{summary.quantiles}
\title{Summary quantiles}
\usage{
\method{summary}{quantiles}(object, ci.level = 0.9, ...)
}
\arguments{
\item{object}{object of quantiles.}

\item{ci.level}{numeric vector of length 1 giving the confidence level (default is 0.9).}

\item{...}{additional arguments submitted to \code{est_cov}.}
}
\value{
A \code{summary.quantiles}-object, a list with dimensions \itemize{
 \item \code{q}
 \item \code{ci.level}
 \item \code{ci}
 \item \code{cov}
}
It is printed with \code{print.summary.quantiles}.
}
\description{
Calculating and printing of summary statistics to a given quantiles-object.
}
\examples{
x <- cbind(rgev(100, shape = .2), rgev(100, shape = .2))

q <- quantiles(parameters(TLMoments(x[, 1]), "gev"), c(.9, .95, .99))
summary(q)
summary(q, select = c(.9, .99))

q <- quantiles(parameters(TLMoments(x[, 1], rightrim = 1), "gev"), .95)
summary(q)

q <- quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
summary(q)
summary(q, select = .95)

q <- quantiles(as.parameters(loc = 10, scale = 5, shape = .3, distr = "gev"), c(.9, .99))
summary(q)
summary(q, rightrim = 1, set.n = 250)

}
\seealso{
\code{\link{quantiles}}, \code{\link{est_cov}}
}
