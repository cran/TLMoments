#' @title
#' Probability weighted moments
#' @description
#' Calculates probability weighted moments up to a specific order. Note that PWMs start with
#' order 0. Acceptable input types are numeric vectors, matrices, lists, and data.frames.
#' @param x numeric vector or matrix, list, or data.frame of data OR an object of TLMoments
#' @param formula if x is of type data.frame a formula has to be submitted
#' @param max.order integer, maximal order of PWMs
#' @param na.rm logical, indicates if NAs should be removed
#' @param ... additional arguments
#' @return a numeric vector, matrix, list, or data.frame consisting of the PWMs and
#' with class \code{PWMs}.
#' The object contains the following attributes: \itemize{
#'  \item \code{order}: a integer vector with corresponding PWM orders
#'  \item \code{source}: a list with background information (used function, data, n, formula;
#'  mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#' @references Greenwood, J. A., Landwehr, J. M., Matalas, N. C., & Wallis, J. R. (1979). Probability weighted moments: definition and relation to parameters of several distributions expressable in inverse form. Water Resources Research, 15(5), 1049-1054.
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
#' # Calculating PWMs from data:
#' PWMs(xvec)
#' PWMs(xmat)
#' PWMs(xlist)
#' PWMs(xdat, formula = hq ~ station)
#' PWMs(xdat, formula = hq ~ season)
#' PWMs(xdat, formula = hq ~ .)
#' PWMs(xdat, formula = . ~ station + season)
#'
#' # Calculating PWMs from L-moments:
#' PWMs(TLMoments(xvec))
#' PWMs(TLMoments(xmat))
#' PWMs(TLMoments(xlist))
#' PWMs(TLMoments(xdat, hq ~ station))
#' PWMs(TLMoments(xdat, hq ~ season))
#' PWMs(TLMoments(xdat, hq ~ .))
#' PWMs(TLMoments(xdat, . ~ station + season))
#'
#' @rdname PWMs
#' @export
PWMs <- function(x, ...) UseMethod("PWMs")

#' @rdname PWMs
#' @method PWMs numeric
#' @export
PWMs.numeric <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- setNames(
    PWM(x, order = 0L:max.order, na.rm = na.rm),
    paste0("beta", 0:max.order)
  )

  attr(out, "order") <- 0L:max.order
  attr(out, "source") <- list(func = "PWMs",
                              data = x,
                              n = length(x),
                              formula = NA)
  class(out) <- c("PWMs", "numeric")
  out
}

#' @rdname PWMs
#' @method PWMs matrix
#' @export
PWMs.matrix <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- apply(x, 2, PWM, order = 0L:max.order, na.rm = na.rm)

  attr(out, "order") <- 0L:max.order
  attr(out, "source") <- list(func = "PWMs",
                              data = x,
                              n = apply(x, 2, function(y) sum(!is.na(y))),
                              formula = NA)
  class(out) <- c("PWMs", "matrix")
  out
}

#' @rdname PWMs
#' @method PWMs list
#' @export
PWMs.list <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- lapply(x, PWM, order = 0L:max.order, na.rm = na.rm)

  attr(out, "order") <- 0L:max.order
  attr(out, "source") <- list(func = "PWMs",
                              data = x,
                              n = vapply(x, length, numeric(1)),
                              formula = NA)
  class(out) <- c("PWMs", "list")
  out
}

#' @rdname PWMs
#' @method PWMs data.frame
#' @export
PWMs.data.frame <- function(x, formula, max.order = 4L, na.rm = FALSE, ...) {

  nam <- getFormulaSides(formula, names(x))
  r <- aggregate(nam$new.formula, data = x, FUN = PWM, order = 0L:max.order, na.rm = na.rm)
  out <- cbind(r[-length(r)], as.data.frame(r[[length(r)]]))

  attr(out, "order") <- 0L:max.order
  attr(out, "source") <- list(func = "PWMs",
                              data = x,
                              n = aggregate(nam$new.formula, x, length)$y,
                              formula = nam$new.formula)
  class(out) <- c("PWMs", "data.frame")
  out
}


#' @rdname PWMs
#' @method PWMs TLMoments
#' @export
PWMs.TLMoments <- function(x, ...) {
  if (attr(x, "leftrim") != 0 | attr(x, "rightrim") != 0) stop("Transformation to PWMs only runs for L-moments. ")
  if (any(diff(attr(x, "order"))) != 1) stop("Transformation to PWMs only runs for sequent L-moments. ")

   UseMethod("PWMs.TLMoments", x$lambdas)
}

#' @method PWMs.TLMoments numeric
#' @export
PWMs.TLMoments.numeric <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  out <- as.numeric(solve(.Call('TLMoments_Z_C', PACKAGE = 'TLMoments', max.order, 0, 0)) %*% x$lambdas)

  names(out) <- paste0("beta", 0:(max.order-1))
  attr(out, "order") <- 0L:(max.order-1)
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "PWMs")
  class(out) <- c("PWMs", "numeric")
  out
}

#' @method PWMs.TLMoments matrix
#' @export
PWMs.TLMoments.matrix <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  out <- solve(.Call('TLMoments_Z_C', PACKAGE = 'TLMoments', max.order, 0, 0)) %*% x$lambdas
  rownames(out) <- paste0("beta", 0:(max.order-1))
  attr(out, "order") <- 0L:(max.order-1)
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "PWMs")
  class(out) <- c("PWMs", "matrix")
  out
}

#' @method PWMs.TLMoments list
#' @export
PWMs.TLMoments.list <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  A <- solve(.Call('TLMoments_Z_C', PACKAGE = 'TLMoments', max.order, 0, 0))
  out <- lapply(x$lambdas, function(x) setNames(as.numeric(A %*% x), paste0("beta", 0:(max.order-1))))
  attr(out, "order") <- 0L:(max.order-1)
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "PWMs")
  class(out) <- c("PWMs", "list")
  out
}

#' @method PWMs.TLMoments data.frame
#' @export
PWMs.TLMoments.data.frame <- function(x, ...) {
  max.order <- max(attr(x, "order"))

  ls <- x$lambdas[, grep("L[0-9]*", names(x$lambdas))]
  fac <- x$lambdas[, !grepl("L[0-9]*", names(x$lambdas)), drop = FALSE]

  out <- as.data.frame(
    t(solve(.Call('TLMoments_Z_C', PACKAGE = 'TLMoments', max.order, 0, 0)) %*% t(as.matrix(ls)))
  )
  names(out) <- paste0("beta", 0:(max.order-1))
  out <- cbind(fac, out)
  attr(out, "order") <- 0L:(max.order-1)
  attr(out, "source") <- attributes(x)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "PWMs")
  class(out) <- c("PWMs", "data.frame")
  out
}

#' @export
print.PWMs <- function(x, ...) {
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
