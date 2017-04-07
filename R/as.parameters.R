#' @title
#' Converting parameter-vectors to a parameters-objects
#' @description
#' description not done yet.
#' @param ... parameters of distribution
#' @param x numeric vector, matrix, list, or data.frame of parameters.
#' @param formula if \code{x} is data.frame a formula has to be given.
#' @param distr character giving the distribution. Function of name
#' q\"distr\" has to be available.
#' @return Named vector of parameters of class parameters.
#' @seealso \code{\link{parameters}}
#' @examples
#' as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' as.parameters(c(loc = 3, scale = 2, shape = .4), distr = "gev")
#'
#' as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' as.parameters(rbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' #as.parameters(matrix(1:9, nr = 3), distr = "gev") # doesn't work!
#'
#' as.parameters(list(list(mean = 2, sd = 1), list(mean = 0, sd = 1)), distr = "norm")
#' as.parameters(list(c(mean = 2, sd = 1), c(mean = 0, sd = 1)), distr = "norm")
#'
#' xdat <- data.frame(station = c(1, 2), mean = c(2, 0), sd = c(1, 1))
#' as.parameters(xdat, cbind(mean, sd) ~ station, distr = "norm")
#' as.parameters(xdat, . ~ station, distr = "norm")
#' as.parameters(xdat, cbind(mean, sd) ~ ., distr = "norm")
#'
#' quantiles(as.parameters(xdat, cbind(mean, sd) ~ ., distr = "norm"), c(.99))
#'
#' # quantile estimation
#' p <- as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' quantiles(p, c(.9, .95))
#' p <- as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' quantiles(p, c(.9, .95))
#' p <- as.parameters(list(list(mean = 2, sd = 1), list(mean = 0, sd = 1)), distr = "norm")
#' quantiles(p, c(.95, .975))
#'
#' # With magrittr
#' library(magrittr)
#' as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev") %>% quantiles(c(.9, .99))
#' @rdname as.parameters
#' @export
as.parameters <- function(..., distr = NULL) {
  if (is.null(distr)) stop("distr of distribution has to be submitted")

  # q <- try(get(paste0("q", distr), mode = "function"), silent = TRUE)
  # if (!is.function(q)) stop(paste0("Found no q-function for ", distr))
  if (grepl("::", x = distr)) { # Falls pkg::func
    f <- sub("^([a-zA-Z0-9]*)::([a-zA-Z0-9]*)$", "\\1::q\\2", x = distr)
    q <- eval(parse(text = paste0("match.fun(", f,")")))
  } else { # falls nur func
    q <- eval(parse(text = paste0("match.fun(q", distr, ")")))
  }
  if (!is.function(q)) stop(paste0("Found no q-function for ", distr))

  UseMethod("as.parameters")
}

#' @describeIn as.parameters as.parameters for numeric data vectors
#' @method as.parameters numeric
#' @export
as.parameters.numeric <- function(..., distr) {
  out <- unlist(list(...))

  attr(out, "distribution") <- distr
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.parameters",
    trimmings = c(NA, NA),
    data = out,
    n = NA,
    formula = NA
  )
  class(out) <- c("parameters", "numeric")
  out
}

#' @describeIn as.parameters as.parameters for numeric data matrices
#' @method as.parameters matrix
#' @export
as.parameters.matrix <- function(x, distr, ...) {
  if (is.null(dimnames(x))) stop("matrix must have rownames or colnames indicating the parameter")

  i <- sapply(dimnames(x), is.null)
  if (sum(i) == 1 & which(i) == 1) {
    x <- t(x)
  }

  attr(x, "distribution") <- distr
  attr(x, "computation.method") <- "input"
  attr(x, "source") <- list(
    func = "as.parameters",
    trimmings = c(NA, NA),
    data = x,
    n = NA,
    formula = NA
  )
  class(x) <- c("parameters", "matrix")
  x
}

#' @describeIn as.parameters as.parameters for numeric data lists
#' @method as.parameters list
#' @export
as.parameters.list <- function(x, distr, ...) {
  out <- lapply(x, function(y) {
    do.call(as.parameters, args = list(unlist(y), distr = distr))
  })

  # Delete attributes...
  for (i in 1:length(out)) {
    attr(out[[i]], "source") <- NULL
    attr(out[[i]], "distribution") <- NULL
  }
  # ...and add global attributes
  attr(out, "distribution") <- distr
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.parameters",
    trimmings = c(NA, NA),
    data = x,
    n = NA,
    formula = NA
  )
  class(out) <- c("parameters", "list")
  out
}

#' @describeIn as.parameters as.parameters for numeric data.frames
#' @method as.parameters data.frame
#' @export
as.parameters.data.frame <- function(x, formula, distr, ...) {

  nam <- getFormulaSides(formula, names(x))
  out <- cbind(x[nam$rhs], x[nam$lhs])

  attr(out, "distribution") <- distr
  attr(out, "computation.method") <- "input"
  attr(out, "source") <- list(
    func = "as.parameters",
    trimmings = c(NA, NA),
    data = x,
    n = NA,
    formula = nam$new.formula
  )
  class(out) <- c("parameters", "data.frame")
  out
}
