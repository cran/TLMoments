% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PWMs.R
\name{PWMs}
\alias{PWMs}
\alias{PWMs.numeric}
\alias{PWMs.matrix}
\alias{PWMs.list}
\alias{PWMs.data.frame}
\alias{PWMs.TLMoments}
\title{Probability weighted moments}
\usage{
PWMs(x, ...)

\method{PWMs}{numeric}(x, max.order = 4L, na.rm = FALSE, ...)

\method{PWMs}{matrix}(x, max.order = 4L, na.rm = FALSE, ...)

\method{PWMs}{list}(x, max.order = 4L, na.rm = FALSE, ...)

\method{PWMs}{data.frame}(x, formula, max.order = 4L, na.rm = FALSE, ...)

\method{PWMs}{TLMoments}(x, ...)
}
\arguments{
\item{x}{numeric vector or matrix, list, or data.frame of data OR an object of TLMoments.}

\item{...}{additional arguments.}

\item{max.order}{integer, maximal order of PWMs.}

\item{na.rm}{logical, indicates if NAs should be removed.}

\item{formula}{if x is of type data.frame a formula has to be submitted.}
}
\value{
numeric vector, matrix, list, or data.frame consisting of the PWMs and
with class \code{PWMs}.
The object contains the following attributes: \itemize{
 \item \code{order}: a integer vector with corresponding PWM orders
 \item \code{source}: a list with background information (used function, data, n, formula;
 mainly for internal purposes)
}
The attributes are hidden in the print-function for a clearer presentation.
}
\description{
Calculates probability weighted moments up to a specific order. Note that PWMs start with
order 0. Acceptable input types are numeric vectors, matrices, lists, and data.frames.
}
\examples{
# Generating data sets:
xmat <- matrix(rnorm(100), nc = 4)
xvec <- xmat[, 3]
xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
xdat <- data.frame(
 station = rep(letters[1:2], each = 50),
 season = rep(c("S", "W"), 50),
 hq = as.vector(xmat)
)

# Calculating PWMs from data:
PWMs(xvec)
PWMs(xmat)
PWMs(xlist)
PWMs(xdat, formula = hq ~ station)
PWMs(xdat, formula = hq ~ season)
PWMs(xdat, formula = hq ~ .)
PWMs(xdat, formula = . ~ station + season)

# Calculating PWMs from L-moments:
PWMs(TLMoments(xvec))
PWMs(TLMoments(xmat))
PWMs(TLMoments(xlist))
PWMs(TLMoments(xdat, hq ~ station))
PWMs(TLMoments(xdat, hq ~ season))
PWMs(TLMoments(xdat, hq ~ .))
PWMs(TLMoments(xdat, . ~ station + season))

# In data.frame-mode invalid names are preceded by "."
xdat <- data.frame(
 beta0 = rep(letters[1:2], each = 50),
 beta1 = as.vector(xmat)
)
PWMs(xdat, formula = beta1 ~ beta0)

}
\references{
Greenwood, J. A., Landwehr, J. M., Matalas, N. C., & Wallis, J. R. (1979). Probability weighted moments: definition and relation to parameters of several distributions expressable in inverse form. Water Resources Research, 15(5), 1049-1054.
}
