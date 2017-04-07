#' @title
#' Convert vector or matrix to a TLMoments-object
#' @description
#' Convert a vector or matrix of TL-moments or TL-moment ratios or a PWMs-object to
#' a TLMoments-object in order to be used with TLMoments-functions.
#' The first position of a vector or the first row of a matrix is
#' always used as the L1-moment. The \code{ratios} argument determines if the
#' following positions or rows are used as TL-moments oder TL-moments ratios.
#' The trimming has to be given using the \code{leftrim} and \code{rightrim} arguments.
#' @param x vector or matrix of TL-moments (or TL-moment ratios if ratios is TRUE) or a PWMs-object.
#' The first position or row is always used as the L1-moment, if vector or matrix are given.
#' @param formula if \code{x} is data.frame. See examples.
#' @param ratios boolean, if TRUE the non-first positions or rows of x give
#' L-moment ratios, if FALSE (default) they give L-moments.
#' If ratios are used and the first position or row is NA, L1 is assumed to be 1!
#' @param leftrim,rightrim integer, order of trimmed L-moments
#' @param ... additional arguments
#' @return an object with class TLMoments, ready to commit to other TLMoments-functions
#' @seealso \code{\link{TLMoments}}
#' @examples
#' ### Vector or matrix as input
#' xmat <- cbind(c(1, .2, .05), c(1, .2, .04), c(1.3, .3, .1))
#' xvec <- xmat[, 1]
#' xlist <- lapply(1:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:3], each = 1),
#'  season = c("S", "W", "S"),
#'  L1 = c(1, 1, 1.3),
#'  L2 = c(.2, .2, .3),
#'  L3 = c(.05, .04, .1)
#' )
#'
#' as.TLMoments(xvec, rightrim = 1)
#' as.TLMoments(xmat, rightrim = 1)
#' as.TLMoments(xlist, rightrim = 1)
#' as.TLMoments(xdat, cbind(L1, L2, L3) ~ station)
#' as.TLMoments(xdat, .~station+season)
#' as.TLMoments(xdat, cbind(L1, L2, L3) ~ .)
#'
#' parameters(as.TLMoments(xvec, rightrim = 0), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(1, .2, .25)), "gev")$para
#'
#' xmat <- cbind(c(NA, .2, -.05), c(NA, .2, .2))
#' xvec <- xmat[, 1]
#'
#' as.TLMoments(xvec, ratios = TRUE)
#' as.TLMoments(xmat, ratios = TRUE)
#' parameters(as.TLMoments(xvec, ratios = TRUE), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(1, .2, -.05)), "gev")$para
#'
#' xmat <- cbind(c(10, .2, -.05), c(10, .2, .2))
#' xvec <- xmat[, 1]
#'
#' as.TLMoments(xvec, ratios = TRUE)
#' as.TLMoments(xmat, ratios = TRUE)
#' parameters(as.TLMoments(xvec, ratios = TRUE), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(10, .2, -.05)), "gev")$para
#'
#' @rdname as.TLMoments
#' @export
as.TLMoments <- function(x, ..., leftrim, rightrim, ratios) UseMethod("as.TLMoments")

#' @describeIn as.TLMoments as.TLMoments for numeric data vectors
#' @method as.TLMoments numeric
#' @export
as.TLMoments.numeric <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {
  dim(x) <- c(length(x), 1)
  out <- as.TLMoments(x, ratios = ratios, leftrim = leftrim, rightrim = rightrim)
  out$lambdas <- drop(out$lambdas)
  out$ratios <- drop(out$ratios)
  out
}

#' @describeIn as.TLMoments as.TLMoments for numeric data matrices
#' @method as.TLMoments matrix
#' @export
as.TLMoments.matrix <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {
  if (!ratios) {
    ls <- x
    taus <- apply(ls, 2, function(y) {
      if (length(y) > 2L) {
        c(NA, y[2]/y[1], y[3:length(y)]/y[2])
      } else if (length(y) == 2L) {
        c(NA, y[2]/y[1])
      } else {
        NA
      }
    })
    if (is.null(dim(taus))) dim(taus) <- c(1, length(taus))
  } else {
    taus <- x
    l1 <- ifelse(is.na(taus[1]), 1, taus[1])
    taus[1] <- NA

    ls <- apply(taus, 2, function(y) {
      if (length(y) > 2L) {
        c(l1, y[2]*l1, y[3:length(y)]*(y[2]*l1))
      } else if (length(y) == 2L) {
        c(l1, y[2]*l1)
      } else {
        l1
      }
    })
  }

  rownames(ls) <- paste0("L", 1L:nrow(ls))
  rownames(taus) <- paste0("T", 1L:nrow(taus))

  out <- list(
    lambdas = ls,
    ratios = taus
  )
  attr(out, "leftrim") <- leftrim
  attr(out, "rightrim") <- rightrim
  attr(out, "order") <- 1L:nrow(ls)
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.TLMoments", #c("as.TLMoments", attr(x, "source")$func),
    data = x, #attr(x, "source")$data,
    n = NA, #attr(x, "source")$n,
    formula = NA #attr(x, "source")$formula
  )

  class(out) <- c("TLMoments", "list")
  out
}

#' @describeIn as.TLMoments as.TLMoments for numeric data lists
#' @method as.TLMoments list
#' @export
as.TLMoments.list <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {
  ls <- lapply(x, as.TLMoments, ratios = ratios, leftrim = leftrim, rightrim = rightrim)

  out <- list(
    lambdas = lapply(ls, getElement, "lambdas"),
    ratios = lapply(ls, getElement, "ratios")
  )
  attr(out, "leftrim") <- leftrim
  attr(out, "rightrim") <- rightrim
  attr(out, "order") <- attr(ls[[1]], "order")
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.TLMoments", #c("as.TLMoments.list", attr(x, "source")$func),
    data = x, #attr(x, "source")$data,
    n = NA, #attr(x, "source")$n,
    formula = NA #attr(x, "source")$formula
  )

  class(out) <- c("TLMoments", "list")
  out
}

#' @describeIn as.TLMoments as.TLMoments for numeric data.frames
#' @method as.TLMoments data.frame
#' @export
as.TLMoments.data.frame <- function(x, formula, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {
  if (!ratios) {
    i <- 2
    while (paste0("L", i) %in% names(x)) {
      if (i == 2) {
        x$.__new_var <- x$L2 / x$L1
      } else {
        x$.__new_var <- x[[paste0("L", i)]] / x$L2
      }
      names(x)[length(names(x))] <- paste0("T", i)
      i <- i+1
    }
  } else {
    print("todo")
  }

  nam <- getFormulaSides(formula, names(x))

  out <- list(
    lambdas = x[, c(which(nam$rhs %in% names(x)), grep("L[0-9]*", names(x)))],
    ratios = x[, c(which(nam$rhs %in% names(x)), grep("T[0-9]*", names(x)))]
  )
  attr(out, "leftrim") <- leftrim
  attr(out, "rightrim") <- rightrim
  attr(out, "order") <- as.integer(gsub(".([0-9]*)", "\\1", nam$lhs))
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.TLMoments", #c("as.TLMoments.data.frame", attr(x, "source")$func),
    data = x, #attr(x, "source")$data,
    n = NA, #attr(x, "source")$n,
    formula = formula
  )

  class(out) <- c("TLMoments", "list")
  out
}
