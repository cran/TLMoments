#' @title
#' Trimmed L-moments
#' @description
#' Calculates empirical or theoretical Trimmed L-moments and -ratios up to a specific order.
#' If empirical moments should be calculated, acceptable input types are numeric vectors,
#' matrices, lists, data.frames. TLMoments is type-preservative, so the input type is also
#' the output type. If theoretical moments should be calculated, the input type has to be
#' of class parameters or PWMs, so an object returned by parameters, as.parameters or
#' PWMs, as.PWMs.
#' @param x numeric data in form of vector, matrix, list, or data.frame OR an object
#' of class parameters or PWMs.
#' @param formula if \code{x} is data.frame. See examples.
#' @param leftrim integer indicating lower trimming parameter, has to be greater than 0.
#' @param rightrim integer indicating upper trimming parameter, has to be greater than 0.
#' @param max.order integer, maximum order of Trimmed L-moments/ratios, has to be
#' greater than 1.
#' @param na.rm logical, indicates if NAs should be removed. Only if empirical moments
#' are calculated.
#' @param computation.method character, indicating if the computation is performed via
#' PWMs, direct, recursive, or recurrence (see References Hosking & Balakrishnan, 2015).
#' Possible values are \code{auto} (default, automatically choose appropriate method), \code{pwm},
#' \code{direct}, \code{recursive}, or \code{recurrence}. Only if empirical moments are calculated.
#' @param ... additional arguments
#' @return a list of two dimensions: \code{lambdas}/\code{ratios} are a numeric vector, matrix,
#' list, or data.frame consisting of the TL-moments/TL-moment-ratios. The list has the class
#' \code{TLMoments}.
#' The object contains the following attributes: \itemize{
#'  \item \code{leftrim}: a numeric giving the used leftrim-argument
#'  \item \code{rightrim}: a numeric giving the used rightrim-argument
#'  \item \code{order}: a integer vector with corresponding TL-moment orders
#'  \item \code{computation.method} a character giving the used computation method
#'  \item \code{source}: a list with background information (used function, data, n, formula;
#'  mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#' @references Elamir, E. A., & Seheult, A. H. (2003). Trimmed L-moments. Computational Statistics & Data Analysis, 43(3), 299-314.
#' @references Hosking, J. R. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 105-124.
#' @references Hosking, J. R. M. (2007). Some theory and practical uses of trimmed L-moments. Journal of Statistical Planning and Inference, 137(9), 3024-3039.
#' @references Hosking, J. R. M., & Balakrishnan, N. (2015). A uniqueness result for L-estimators, with applications to L-moments. Statistical Methodology, 24, 69-80.
#' @seealso \code{\link{PWMs}}, \code{\link{parameters}}, \code{\link{quantiles}}, \code{\link{summary.TLMoments}}, \code{\link{as.TLMoments}}
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
#' # Calculating TL-moments from data:
#' TLMoments(xvec, leftrim = 0, rightrim = 1)
#' TLMoments(xmat, leftrim = 1, rightrim = 1)
#' TLMoments(xlist)
#' TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2)
#' TLMoments(xdat, hq ~ season, leftrim = 0, rightrim = 2)
#' TLMoments(xdat, hq ~ ., leftrim = 0, rightrim = 2)
#'
#' # Calculating TL-moments from PWMs:
#' TLMoments(PWMs(xvec))
#' TLMoments(PWMs(xmat), rightrim = 1)
#' TLMoments(PWMs(xlist), leftrim = 1, rightrim = 1)
#' TLMoments(PWMs(xdat, hq ~ station), leftrim = 0, rightrim = 2)
#' TLMoments(PWMs(xdat, hq ~ station + season), leftrim = 0, rightrim = 2)
#'
#' # Calculating TL-moments from parameters:
#' (tlm <- TLMoments(xmat, leftrim = 0, rightrim = 1))
#' TLMoments(parameters(tlm, "gev"))
#'
#' p <- as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' TLMoments(p, rightrim = 1)
#' p <- as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' TLMoments(p, max.order = 6)
#' p <- as.parameters(list(
#'  list(loc = 3, scale = 2, shape = .4),
#'  list(loc = 3, scale = 2, shape = .2)
#' ), distr = "gev")
#' TLMoments(p)
#' p <- as.parameters(data.frame(
#'  station = letters[1:2],
#'  loc = c(2, 3),
#'  scale = c(2, 2),
#'  shape = c(.4, .2)
#' ), .~station, distr = "gev")
#' TLMoments(p)
#' @export
TLMoments <- function(x, ...) UseMethod("TLMoments")

#' @describeIn TLMoments TLMoments for numeric vector of data
#' @method TLMoments numeric
#' @export
TLMoments.numeric <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L,
                              na.rm = FALSE, computation.method = "auto",
                              ...) {
  if (!are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  if (!are.integer.like(leftrim, rightrim) | any(c(leftrim, rightrim) < 0))
    stop("leftrim and rightrim must be positive integers. ")
  if (!is.logical(na.rm))
    stop("na.rm must be TRUE or FALSE. ")
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- TLMoment(x, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)

  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)
  attr(out, "computation.method") <- computation.method
  attr(out, "source") <- list(
    func = "TLMoments",
    data = x,
    n = sum(!is.na(x)),
    formula = NA
  )
  out
}

#' @describeIn TLMoments TLMoments for numeric matrix of data
#' @method TLMoments matrix
#' @export
TLMoments.matrix <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L,
                             na.rm = FALSE, computation.method = "auto",
                             ...) {
  if (!are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  # if (!are.numeric(leftrim, rightrim))
  #   stop("leftrim and rightrim must be numeric. ")
  if (!are.integer.like(leftrim, rightrim) | any(c(leftrim, rightrim) < 0))
    stop("leftrim and rightrim must be positive integers. ")
  if (!is.logical(na.rm))
    stop("na.rm must be TRUE or FALSE. ")
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- apply(x, 2, TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  if (max.order == 1) {
    dim(ls) <- c(1, ncol(x))
  }

  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)
  attr(out, "computation.method") <- computation.method
  attr(out, "source") <- list(
    func = "TLMoments",
    data = x,
    n = apply(x, 2, function(y) sum(!is.na(y))),
    formula = NA
  )
  out
}

#' @describeIn TLMoments TLMoments for numeric list of data
#' @method TLMoments list
#' @export
TLMoments.list <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L, na.rm = FALSE,
                           computation.method = "auto",
                           ...) {
  if (!are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  if (!are.numeric(leftrim, rightrim))
    stop("leftrim and rightrim must be numeric. ")
  if (!is.logical(na.rm))
    stop("na.rm must be TRUE or FALSE. ")
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- lapply(x, TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)

  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)
  attr(out, "computation.method") <- computation.method
  attr(out, "source") <- list(
    func = "TLMoments",
    data = x,
    n = vapply(x, function(y) sum(!is.na(y)), numeric(1)),
    formula = NA
  )
  out
}

#' @describeIn TLMoments TLMoments for numeric data.frame of data
#' @method TLMoments data.frame
#' @export
TLMoments.data.frame <- function(x, formula, leftrim = 0L, max.order = 4L, rightrim = 0L,
                                 na.rm = FALSE, computation.method = "auto",
                                 ...) {
  if (!are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  if (!are.numeric(leftrim, rightrim))
    stop("leftrim and rightrim must be numeric. ")
  if (!is.logical(na.rm))
    stop("na.rm must be TRUE or FALSE. ")
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  nam <- getFormulaSides(formula, names(x))
  agg <- aggregate(nam$new.formula, data = x, FUN = TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  ls <- cbind(agg[-length(agg)], as.data.frame(agg[[length(agg)]]))

  out <- as.TLMoments(ls, as.formula(paste0(". ~ ", paste0(nam$rhs, collapse = "+"))), leftrim = leftrim, rightrim = rightrim)

  attr(out, "computation.method") <- computation.method
  attr(out, "source") <- list(
    func = "TLMoments",
    data = x,
    n = aggregate(nam$new.formula, x, function(y) sum(!is.na(y)))$y,
    formula = nam$new.formula
  )
  out
}



#' @describeIn TLMoments TLMoments for PWMs-object
#' @method TLMoments PWMs
#' @export
TLMoments.PWMs <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  if (!are.integer.like(leftrim, rightrim))
    stop("leftrim and rightrim must be integers")
  if (any(diff(attr(x, "order")) != 1))
    stop("PWM order must not have gaps")

  UseMethod("TLMoments.PWMs")
}

#' @method TLMoments.PWMs numeric
#' @export
TLMoments.PWMs.numeric <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- .Call('TLMoments_PWM_to_TLMoments', PACKAGE = 'TLMoments', x, leftrim, rightrim)
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "computation.method") <- "pwms"
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.PWMs")
  attr(out, "source")$pwms <- removeAttributes(x)
  out
}

#' @method TLMoments.PWMs matrix
#' @export
TLMoments.PWMs.matrix <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- apply(x, 2, function(xx) {
    .Call('TLMoments_PWM_to_TLMoments', PACKAGE = 'TLMoments', xx, leftrim, rightrim)
  })
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "computation.method") <- "pwms"
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.PWMs")
  attr(out, "source")$pwms <- removeAttributes(x)
  out
}

#' @method TLMoments.PWMs list
#' @export
TLMoments.PWMs.list <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- lapply(x, function(xx) {
    .Call('TLMoments_PWM_to_TLMoments', PACKAGE = 'TLMoments', xx, leftrim, rightrim)
  })
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "computation.method") <- "pwms"
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.PWMs")
  attr(out, "source")$pwms <- removeAttributes(x)
  out
}

#' @method TLMoments.PWMs data.frame
#' @export
TLMoments.PWMs.data.frame <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  pwms <- x[, grep("beta[0-9]*", names(x)), drop = FALSE]
  fac <- x[, !grepl("beta[0-9]*", names(x)), drop = FALSE]
  ls <- apply(pwms, 1, function(xx) {
    .Call('TLMoments_PWM_to_TLMoments', PACKAGE = 'TLMoments', xx, leftrim, rightrim)
  })
  ls <- as.data.frame(t(ls))
  names(ls) <- paste0("L", 1:ncol(ls))

  rhs <- dimnames(attr(terms(attr(x, "source")$formula), "factors"))[[2]]
  out <- as.TLMoments(cbind(fac, ls), as.formula(paste0(". ~ ", paste0(rhs, collapse = "+"))), leftrim = leftrim, rightrim = rightrim)

  attr(out, "computation.method") <- "pwms"
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.PWMs")
  attr(out, "source")$pwms <- removeAttributes(x)
  out
}



#' @describeIn TLMoments TLMoments for parameters-object
#' @method TLMoments parameters
#' @export
TLMoments.parameters <- function(x,
                                 leftrim = attr(x, "source")$trimmings[1],
                                 rightrim = attr(x, "source")$trimmings[2],
                                 max.order = 4L,
                                 ...) {

  if (is.na(leftrim)) leftrim <- 0L
  if (is.na(rightrim)) rightrim <- 0L
  if (is.null(max.order) | max.order == 0) max.order <- 4L

  if (!are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  if (!are.numeric(leftrim, rightrim))
    stop("leftrim and rightrim must be numeric. ")

  UseMethod("TLMoments.parameters")
}

#' @method TLMoments.parameters numeric
#' @export
TLMoments.parameters.numeric <- function(x,
                                         leftrim = attr(x, "source")$trimmings[1],
                                         rightrim = attr(x, "source")$trimmings[2],
                                         max.order = attr(x, "source")$max.order,
                                         ...) {

  if (is.na(leftrim)) leftrim <- 0L
  if (is.na(rightrim)) rightrim <- 0L
  if (is.null(max.order) || max.order == 0) max.order <- 4L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      attr(x, "source")$max.order != max.order) { # calculate new

    ls <- calcTLMom(max.order, leftrim, rightrim, qfunc = getQ(x))
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.parameters")
  #attr(out, "source")$trimmings <- NULL
  out
}

#' @method TLMoments.parameters matrix
#' @export
TLMoments.parameters.matrix <- function(x,
                                        leftrim = attr(x, "source")$trimmings[1],
                                        rightrim = attr(x, "source")$trimmings[2],
                                        max.order = attr(x, "source")$max.order,
                                        ...) {

  if (is.na(leftrim)) leftrim <- 0L
  if (is.na(rightrim)) rightrim <- 0L
  if (is.null(max.order) || max.order == 0) max.order <- 4L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      attr(x, "source")$max.order != max.order) { # calculate new

    ls <- apply(x, 2, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.parameters")
  #attr(out, "source")$trimmings <- NULL
  out
}

#' @method TLMoments.parameters list
#' @export
TLMoments.parameters.list <- function(x,
                                      leftrim = attr(x, "source")$trimmings[1],
                                      rightrim = attr(x, "source")$trimmings[2],
                                      max.order = attr(x, "source")$max.order,
                                      ...) {

  if (is.na(leftrim)) leftrim <- 0L
  if (is.na(rightrim)) rightrim <- 0L
  if (is.null(max.order) || max.order == 0) max.order <- 4L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      attr(x, "source")$max.order != max.order) { # calculate new

    ls <- lapply(x, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }
  out <- as.TLMoments(ls, leftrim = leftrim, rightrim = rightrim)

  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.parameters")
  #attr(out, "source")$trimmings <- NULL
  out
}

#' @method TLMoments.parameters data.frame
#' @export
TLMoments.parameters.data.frame <- function(x,
                                            leftrim = attr(x, "source")$trimmings[1],
                                            rightrim = attr(x, "source")$trimmings[2],
                                            max.order = attr(x, "source")$max.order,
                                            ...) {

  if (is.na(leftrim)) leftrim <- 0L
  if (is.na(rightrim)) rightrim <- 0L
  if (is.null(max.order) || max.order == 0) max.order <- 4L

  nam <- getFormulaSides(attr(x, "source")$formula)
  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      attr(x, "source")$max.order != max.order) { # calculate new

    ls <- apply(x[!(names(x) %in% nam$rhs)], 1, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
    ls <- as.data.frame(t(ls))
    names(ls) <- paste0("L", 1:max.order)
    ls <- cbind(x[nam$rhs], ls)

  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }

  out <- as.TLMoments(ls, as.formula(paste0(".~", paste0(nam$rhs, collapse = "+"))), leftrim = leftrim, rightrim = rightrim)

  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "TLMoments.parameters")
  #attr(out, "source")$trimmings <- NULL
  out
}


#' @export
print.TLMoments <- function(x, ...) {
  # if ("data.frame" %in% class(x$lambdas)) {
  #   print.data.frame(x)
  #   return(invisible(x))
  # }

  tmp <- x
  attributes(tmp) <- NULL

  dim(tmp) <- dim(x)
  names(tmp) <- names(x)
  dimnames(tmp) <- dimnames(x)

  print(tmp)
  invisible(x)
}

