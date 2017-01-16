#' @title
#' Calculation of the covariance matrix of estimations of PWMs, TLMoments, parameters, or quantiles
#' @description
#' description not done yet
#' @param x object of \code{PWMs}, \code{TLMoments}, \code{parameters},
#'        or \code{quantiles} constructed using the same-named functions.
#' @param order numeric vector giving the orders of PWMs or TLMoments, if \code{est_cov} is calculated from
#' \code{PWMs} or \code{TLMoments}.
#' @param leftrim,rightrim lower and upper trimming parameter, if \code{est_cov} is calculated from \code{TLMoments}.
#' @param distr character giving the distribution, if \code{est_cov} is calculated with parametric assumption.
#' @param p numeric vector giving the quantile estimates from which the covariance should be calculated, if \code{est_cov}
#' is calculated from \code{quantiles}.
#' @param ... additional arguments given to the sub-functions. Normally not necessary.
#' @return a numeric matrix (if \code{x} is of class \code{PWMs}, \code{parameters}, or
#'          \code{quantiles}) or a list of two matrices (\code{lambdas} and \code{ratios}, if
#'          \code{x} is of class \code{TLMoments}).
#'
#' @examples
#' ### 1: PWMs:
#'
#' x <- evd::rgev(100, shape = .1)
#'
#' est_cov(PWMs(x))
#' est_cov(PWMs(x), order = 1:2)
#'
#' est_cov(PWMs(x), distr = "gev") # parametric with gev assumption
#'
#' x <- cbind(evd::rgev(100, shape = .1), evd::rgev(100, shape = .3))
#' est_cov(PWMs(x))
#' est_cov(PWMs(x), order = 0:1)
#' #cov(t(replicate(100000,
#' # as.vector(PWMs(cbind(evd::rgev(100, shape = .1), evd::rgev(100, shape = .3)), 1)))
#' #))
#'
#' ### 2. TLMoments:
#'
#' x <- evd::rgev(100, shape = .1)
#'
#' est_cov(TLMoments(x))
#' est_cov(TLMoments(x, rightrim = 1))
#' est_cov(TLMoments(x), order = 1:2)
#'
#' est_cov(TLMoments(x), distr = "gev") # parametric with gev assumption
#'
#' x <- cbind(evd::rgev(100, shape = .1), evd::rgev(100, shape = .3))
#' est_cov(TLMoments(x))
#'
#' para <- as.parameters(loc = 10, scale = 5, shape = .2, distr = "gev")
#' est_cov(TLMoments(para, rightrim = 0), distr = "gev", set.n = 100)
#' est_cov(TLMoments(para, rightrim = 1), distr = "gev", set.n = 100)
#'
#' ### 3. Parameters:
#'
#' x <- evd::rgev(100, shape = .1)
#'
#' est_cov(parameters(TLMoments(x), "gev")) # default parametric
#' est_cov(parameters(TLMoments(x, rightrim = 1), "gev"))
#'
#' est_cov(parameters(TLMoments(x), "gev"), np.cov = TRUE) # explicit non-parametric
#' est_cov(parameters(TLMoments(x, rightrim = 1), "gev"), np.cov = TRUE)
#'
#' x <- cbind(evd::rgev(100, shape = .1), evd::rgev(100, shape = .3))
#' est_cov(parameters(TLMoments(x), "gev"))
#'
#' para <- as.parameters(loc = 10, scale = 5, shape = .2, distr = "gev")
#' est_cov(para, set.n = 100)
#' est_cov(para, rightrim = 1, set.n = 100)
#'
#' # var(replicate(1000, parameters(TLMoments(evd::rgev(100, 10, 5, .2)), "gev")[1]))
#' # var(replicate(1000, parameters(TLMoments(evd::rgev(100, 10, 5, .2), rightrim = 1), "gev")[1]))
#'
#' ### 4. Quantiles:
#'
#' library(evd)
#' x <- evd::rgev(100, shape = .2)
#'
#' q <- quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
#' est_cov(q) # default parametric
#' est_cov(q, np.cov = TRUE) # explicit non-parametric
#'
#' q <- quantiles(parameters(TLMoments(c(rep(NA, 1000), x), na.rm = TRUE), "gev"), c(.9, .95, .99))
#' est_cov(q)
#'
#' q <- quantiles(parameters(TLMoments(x, rightrim = 1), "gev"), c(.9, .95, .99))
#' est_cov(q) # default parametric
#' est_cov(q, np.cov = TRUE) # explicit non-parametric
#'
#' x <- evd::rgev(1000, shape = .1)
#' q <- quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
#' est_cov(q) # default parametric
#' est_cov(q, np.cov = TRUE) # explicit non-parametric
#'
#' q <- quantiles(as.parameters(loc = 10, scale = 5, shape = .3, distr = "gev"), c(.9, .99))
#' est_cov(q)
#' est_cov(q, leftrim = 0, rightrim = 1)
#' est_cov(q, leftrim = 0, rightrim = 1, set.n = 100)
#'
#' @rdname est_cov
#' @export
est_cov <- function(x, ...) {
  # TODO: General error catches
  if (!any(c("PWMs", "TLMoments", "parameters", "quantiles") %in% class(x)))
    stop("x must be of class PWMs, TLMoments, parameters, or quantiles. ")

  UseMethod("est_cov")
}

#' @rdname est_cov
#' @method est_cov PWMs
#' @export
est_cov.PWMs <- function(x,
                         order = attr(x, "order"),
                         ...) {
  # TODO: General error catches
  if (!any(c("numeric", "matrix") %in% class(x)))
    stop("To date, only numeric and matrix types of PWMs are supported. ")

  est_pwmcov(attr(x, "source")$data,
             order = order,
             ...)
}

#' @rdname est_cov
#' @method est_cov TLMoments
#' @export
est_cov.TLMoments <- function(x,
                              leftrim = attr(x, "leftrim"),
                              rightrim = attr(x, "rightrim"),
                              order = attr(x, "order"),
                              ...) {
  # TODO: General error catches
  if (!any(c("numeric", "matrix") %in% class(x$lambdas)))
    stop("To date, only numeric and matrix types of TLMoments are supported. ")

  if (attr(x, "source")$func[1] %in% c("as.PWMs", "as.TLMoments", "as.parameters")) { # theoretical values
    est_tlmcov(x, ...)
  } else { # empirical values
    est_tlmcov(attr(x, "source")$data,
               leftrim = leftrim,
               rightrim = rightrim,
               order = order,
               ...)
  }
}

#' @rdname est_cov
#' @method est_cov parameters
#' @export
est_cov.parameters <- function(x,
                               distr = attr(x, "distribution"),
                               leftrim = attr(x, "source")$trimmings[1],
                               rightrim = attr(x, "source")$trimmings[2],
                               ...) {
  # TODO: General error catches
  if (!any(c("numeric", "matrix") %in% class(x)))
    stop("To date, only numeric and matrix types of PWMs are supported. ")

  if (attr(x, "source")$func[1] %in% c("as.PWMs", "as.TLMoments", "as.parameters")) {  # theoretical values
    # if (!are.integer.like(attr(x, "source")$trimmings)) {
    #   attr(x, "source")$trimmings <- c(leftrim, rightrim)
    #   #warning("Calculation based on TL(0,0)-moments. Overwrite with \"leftrim\" and \"rightrim\". ")
    # }
    est_paramcov(x, ...)
  } else { # empirical values
    est_paramcov(attr(x, "source")$data,
                 distr = distr,
                 leftrim = leftrim,
                 rightrim = rightrim,
                 ...)
  }
}

#' @rdname est_cov
#' @method est_cov quantiles
#' @export
est_cov.quantiles <- function(x,
                              distr = attr(x, "distribution"),
                              p = attr(x, "p"),
                              leftrim = attr(x, "source")$trimmings[1],
                              rightrim = attr(x, "source")$trimmings[2],
                              ...) {
  # TODO: General error catches
  if (!any(c("numeric", "matrix") %in% class(x)))
    stop("To date, only numeric and matrix types of PWMs are supported. ")

  if (attr(x, "source")$func[1] %in% c("as.PWMs", "as.TLMoments", "as.parameters")) {  # theoretical values
    # if (!are.integer.like(attr(x, "source")$trimmings)) {
    #   attr(x, "source")$trimmings <- c(leftrim, rightrim)
    #   #warning("Calculation based on TL(0,0)-moments. Overwrite with \"leftrim\" and \"rightrim\". ")
    # }
    est_quancov(x, ...)
  } else { # empirical values
    est_quancov(attr(x, "source")$data,
                distr = distr,
                p = p,
                leftrim = leftrim,
                rightrim = rightrim,
                ...)
  }
}
