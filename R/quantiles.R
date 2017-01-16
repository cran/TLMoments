#' @title
#' Calculating quantiles from distribution parameters
#' @description
#' Calculates quantiles from distribution parameters received by parameters or
#' from a named vector.
#' @param x object returned by parameters or a named vector. In the latter
#' you have to specify the \code{distr}-argument.
#' @param p numeric vector giving the quantiles to calculate.
#' @param distr character object defining the distribution. Supported types are
#' "gev", "gum" and "gpd". You do not need to set this, if \code{x} is from parameters.
#' @return a numeric vector, matrix, list, or data.frame of the quantiles and with class
#' \code{quantiles}. The object contains the following attributes: \itemize{
#'  \item \code{distribution}: a character indicating the used distribution
#'  \item \code{p}: a vector with the calculated quantiles
#'  \item \code{source}: a list with background information (used function, data, n,
#'  formula, trimmings; mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#' @examples
#' # Generating data sets:
#' xmat <- matrix(rnorm(100), nc = 4)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:2], each = 50),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' # Calculating quantiles from parameters-object
#' tlm <- TLMoments(xvec, leftrim = 0, rightrim = 1)
#' quantiles(parameters(tlm, "gev"), c(.9, .99))
#' tlm <- TLMoments(xmat, leftrim = 1, rightrim = 1)
#' quantiles(parameters(tlm, "gumbel"), c(.9, .95, .999))
#' tlm <- TLMoments(xlist)
#' quantiles(parameters(tlm, "gpd"), .999)
#' tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2)
#' quantiles(parameters(tlm, "gev"), seq(.1, .9, .1))
#' tlm <- TLMoments(xdat, hq ~ station + season, leftrim = 0, rightrim = 2)
#' quantiles(parameters(tlm, "gev"), seq(.1, .9, .1))
#'
#' # Distribution can be overwritten (but parameters have to fit)
#' tlm <- TLMoments(xvec, leftrim = 0, rightrim = 1)
#' params <- parameters(tlm, "gev")
#' quantiles(params, c(.9, .99))
#' quantiles(params[1:2], c(.9, .99), distr = "evd::gumbel")
#' evd::qgumbel(c(.9, .99), loc = params[1], scale = params[2])
#'
#'
#' # Using magrittr
#' library(magrittr)
#' evd::rgev(50, shape = .3) %>%
#'   TLMoments(leftrim = 0, rightrim = 1) %>%
#'   parameters("gev") %>%
#'   quantiles(c(.99, .999))
#'
#' # Calculating quantiles to given parameters for arbitrary functions
#' quantiles(c(mean = 10, sd = 3), c(.95, .99), "norm")
#' qnorm(c(.95, .99), mean = 10, sd = 3)
#' quantiles(c(loc = 10, scale = 5, shape = .3), c(.95, .99), "evd::gev")
#' evd::qgev(c(.95, .99), loc = 10, scale = 5, shape = .3)
#' quantiles(c(loc = 10, scale = 5, shape = .3), c(.95, .99), "evd::gpd")
#' evd::qgpd(c(.95, .99), loc = 10, scale = 5, shape = .3)
#'
#' # These give errors:
#' #quantiles(c(loc = 10, scale = 5, shape = .3), c(.95, .99), "notexistingdistribution")
#' #quantiles(c(loc = 10, scale = 5, shpe = .3), c(.95, .99), "evd::gev") # wrong arguments
#' @export
quantiles <- function(x, p, distr = attr(x, "distribution")) {
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")

  UseMethod("quantiles")
}

#' @method quantiles matrix
#' @export
quantiles.matrix <- function(x, p,
                             distr = attr(x, "distribution")) {
  out <- apply(x, 2, quantiles.numeric, p = p, distr = distr)

  class(out) <- c("quantiles", "matrix")
  attr(out, "distribution") <- distr
  attr(out, "p") <- p
  attr(out, "source") <- attr(x, "source")
  attr(out, "source")$func <- c(attr(x, "source")$func, "quantiles")
  attr(out, "source")$parameters <- removeAttributes(x)
  out
}

#' @method quantiles list
#' @export
quantiles.list <- function(x, p,
                           distr = attr(x, "distribution")) {
  out <- lapply(x, quantiles.numeric, p = p, distr = distr)

  # Delete attributes...
  for (i in 1:length(out)) {
    attr(out[[i]], "source") <- NULL
    attr(out[[i]], "distribution") <- NULL
    attr(out[[i]], "p") <- NULL
  }
  # ...and add global attributes
  class(out) <- c("quantiles", "list")
  attr(out, "distribution") <- distr
  attr(out, "p") <- p
  attr(out, "source") <- attr(x, "source")
  attr(out, "source")$func <- c(attr(x, "source")$func, "quantiles")
  attr(out, "source")$parameters <- removeAttributes(x)
  out
}

#' @method quantiles data.frame
#' @export
quantiles.data.frame <- function(x, p,
                                 distr = attr(x, "distribution")) {

  formula <- attr(x, "source")$formula
  nam <- getFormulaSides(formula, names(x))

  out <- apply(x[!(names(x) %in% nam$rhs)], 1, quantiles, p = p, distr = distr)

  if (length(dim(out)) == 2) {
    out <- cbind(x[nam$rhs], as.data.frame(t(out)))
  } else {
    out <- cbind(x[nam$rhs], as.data.frame(out))
    names(out)[-seq_along(nam$rhs)] <- as.character(p)
  }

  class(out) <- c("quantiles", "data.frame")
  attr(out, "distribution") <- distr
  attr(out, "p") <- p
  attr(out, "source") <- attr(x, "source")
  attr(out, "source")$func <- c(attr(x, "source")$func, "quantiles")
  attr(out, "source")$parameters <- removeAttributes(x)
  out
}

#' @method quantiles numeric
#' @export
quantiles.numeric <- function(x, p,
                              distr = attr(x, "distribution")) {
  if (is.null(distr)) stop("Argument distr defining the distribution must be submitted.")

  if (!("parameters" %in% class(x))) {
    x <- as.parameters(x, distr = distr)
  }

  q <- do.call(getQ, c(x = distr, as.list(x)))
  out <- setNames(q(p), p)
  class(out) <- c("quantiles", "numeric")
  attr(out, "distribution") <- distr
  attr(out, "p") <- p
  attr(out, "source") <- attr(x, "source")
  attr(out, "source")$func <- c(attr(x, "source")$func, "quantiles")
  attr(out, "source")$parameters <- removeAttributes(x)
  out
}


#' @export
print.quantiles <- function(x, ...) {
  if ("data.frame" %in% class(x)) {
    print.data.frame(x)
    return(invisible(x))
  }

  tmp <- x
  attributes(tmp) <- NULL
  dim(tmp) <- dim(x)

  names(tmp) <- names(x)
  dimnames(tmp) <- dimnames(x)
  print(tmp)
  invisible(x)
}
