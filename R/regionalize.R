#' @title
#' regionalize - Calculating regionalized TL-moments
#' @description
#' regionalize takes the result of TLMoments.matrix and calculates a weighted mean
#' of TL-moments and TL-moment ratios.
#' @param x object returned by TLMoments.matrix. TLMoments.list and
#' TLMoments.data.frame are not supported by now.
#' @param w numeric vector giving the weights. The function ensures
#' that it adds up to 1.
#' @param reg.lambdas logical, if TRUE (default) regionalization is
#' based upon TL-moments. If false it's based on TL-moment-ratios
#' @return a list of two dimensions: \code{lambdas}/\code{ratios} are a matrix
#' consisting of the regionalized TL-moments/TL-moment-ratios. The list has
#' the class \code{TLMoments}. The object contains the following attributes: \itemize{
#'  \item \code{leftrim}: a numeric giving the used leftrim-argument
#'  \item \code{rightrim}: a numeric giving the used rightrim-argument
#'  \item \code{order}: a integer vector with corresponding TL-moment orders
#'  \item \code{computation.method} a character giving the used computation method
#'  \item \code{source}: a list with background information (used function, data,
#'  n, formula; mainly for internal purposes)
#' }
#' @examples
#' tlm <- TLMoments(
#'   matrix(evd::rgev(200, loc = 10, scale = 5, shape = .3), nc = 5),
#'   leftrim = 0, rightrim = 1)
#' regionalize(tlm)
#' parameters(regionalize(tlm), "gev")
#' quantiles(parameters(regionalize(tlm), "gev"), c(.99, .999))
#' quantiles(parameters(regionalize(tlm, reg.lambdas = FALSE), "gev"), c(.99, .999))
#'
#' # With magrittr
#' library(magrittr)
#' matrix(evd::rgev(200, shape = .3), nc = 5) %>%
#'  TLMoments(rightrim = 1) %>%
#'  regionalize %>%
#'  parameters("gev") %>%
#'  quantiles(c(.99, .999))
#' @export
regionalize <- function(x, w = rep(1, ncol(x$lambdas)), reg.lambdas = TRUE) {
  if (!("TLMoments" %in% class(x)))
    stop("x must be object of class TLMoments. ")
  if (!is.matrix(x$lambdas))
    stop("only matrix types are permitted by now. ")

  # Ensure that the sum of weights is 1.
  if (sum(w) != 1) w <- w / sum(w)

  # Two versions:
  # 1) Regionalize lambdas and calculate taus (default)
  # 2) Regionalize taus (and l1) und calculate lambdas
  if (reg.lambdas) {
    l <- apply(x$lambdas, 1, weighted.mean, w = w)
    out <- as.TLMoments(l, leftrim = attr(x, "leftrim"), rightrim = attr(x, "rightrim"))
  } else {
    r <- apply(x$ratios, 1, weighted.mean, w = w)
    l1 <- weighted.mean(x$lambdas[1, ], w)
    out <- as.TLMoments(c(l1, r[-1]), leftrim = attr(x, "leftrim"), rightrim = attr(x, "rightrim"), ratios = TRUE)
  }

  # Adjust attributes
  attr(out, "source") <- attr(x, "source")
  attr(out, "source")$func <- c(attr(out, "source")$func, "regionalize")
  out
}
