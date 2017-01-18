# automatically select computation.method:
#   L-moments: recursive
#   TL-moments: recurrence
#   NI-TL-moments: direct
select_computation <- function(leftrim, rightrim) {
  stopifnot(is.numeric(leftrim), is.numeric(rightrim))

  if (isTRUE(all.equal(leftrim, 0L)) & isTRUE(all.equal(rightrim, 0L))) {
    return("recursive")
  } else if (are.integer.like(leftrim, rightrim)) {
    return("recurrence")
  } else {
    return("direct")
  }
}

# calculate pseudo observations
pseudo <- function(x, k) {
  .Call('TLMoments_pseudo_C', PACKAGE = 'TLMoments', sort(na.omit(x)), k)[rank(x, na.last = "keep"), ]
}

# Functions for general error catching
are.integer.like <- function(...) {
  all(vapply(list(...), function(a) isTRUE(all.equal(a, as.integer(a))), logical(1)))
}
are.numeric <- function(...) {
  all(vapply(list(...), function(a) is.numeric(a), logical(1)))
}

# infix for intervals
`%-+%` <- function(a, b) { cbind(a - b, a + b) }


calcTLMom <- function(maxr, s, t, qfunc) {
  if (s < 0 | t < 0 | maxr < 1) stop("s and t have to be >= 0, maxr has to be >= 1")
  if (!is.function(qfunc)) stop("Wrong qfunc")

  vapply(1L:maxr, function(r) {
    sum(vapply(0:(r-1), function(j) {
      i <- integrate(function(u) u^(s+j) * (1-u)^(t+r-j-1) * qfunc(u), lower = 0, upper = 1)
      if (i$message != "OK") stop("Error occurred while integrating. ")
      (-1)^(r-j-1) * factorial(r-1)*factorial(r+s+t) /r /factorial(j) /factorial(r+t-j-1) /factorial(r-j-1) /factorial(s+j) * i$value
    }, numeric(1)))
  }, numeric(1))
}
#
# calcTLMom(4, 0, 2, function(u) evd::qgumbel(u, loc = 0, scale = 1))
# lmomco::theoTLmoms(lmomco::vec2par(c(0, 1), type = "gum"), nmom = 4, leftrim = 0, rightrim = 2)$lambdas
#
# calcTLMom(4, 2, 1, function(u) evd::qgev(u, loc = 0, scale = 1, shape = .2))
# lmomco::theoTLmoms(lmomco::vec2par(c(0, 1, -.2), type = "gev"), nmom = 4, leftrim = 2, rightrim = 1)$lambdas
#
# microbenchmark::microbenchmark(
#   calcTLMom(4, 2, 1, function(u) evd::qgev(u, loc = 0, scale = 1, shape = .2)),
#   lmomco::theoTLmoms(lmomco::vec2par(c(0, 1, -.2), type = "gev"), nmom = 4, leftrim = 2, rightrim = 1)$lambdas
# )


getQ <- function(x, ...) {
  if (!("parameters" %in% class(x)) & !("character" %in% class(x)))
    stop("x must be of class parameters or character vector")

  UseMethod("getQ")
}

getQ.character <- function(x, ...) {
  distr <- x
  args <- list(...)

  if (grepl("::", x = distr)) { # Falls pkg::func
    f <- sub("^([a-zA-Z0-9]*)::([a-zA-Z0-9]*)$", "\\1::q\\2", x = distr)
    q <- eval(parse(text = paste0("match.fun(", f,")")))
  } else { # falls nur func
    q <- eval(parse(text = paste0("match.fun(q", distr, ")")))
  }
  if (!is.function(q)) stop(paste0("Found no q-function for ", distr))

  if (any(!(names(args) %in% names(formals(q))))) stop("Wrong arguments given.")
  formals(q)[names(args)] <- args
  q
}

getQ.parameters <- function(x) {
  if (!("numeric" %in% class(x))) stop("ATM only for parameters, numeric!")

  distr <- attr(x, "distribution")
  args <- as.list(x)

  do.call(getQ.character, c(x = distr, args))
}
#
# getQ(as.parameters(loc = 9, scale = 5, shape = .3, type = "evd::gev"))
# getQ.character("evd::gev", loc = 10, scale = 4, shape = .2)


getFormulaSides <- function(formula, names = NULL) {

  lhs <- all.vars(update(formula, .~0))
  all <- all.vars(formula)
  rhs <- all[!(all %in% lhs)]

  if (!is.null(names)) {
    if (length(lhs) == 1 && lhs == ".") {
      lhs <- names[!(names %in% rhs)]
    }
    if (length(rhs) == 1 && rhs == ".") {
      rhs <- names[!(names %in% lhs)]
    }
    if ("." %in% all) all <- names

    #if (!all(lhs %in% names)) stop("Formula error")
    #if (!all(rhs %in% names)) stop("Formula error")
  }

  list(lhs = lhs, rhs = rhs, all = all,
       new.formula = as.formula(paste0("cbind(", paste0(lhs, collapse = ","), ") ~ ", paste0(rhs, collapse = "+"))))
}
# getFormulaSides(z ~ x + y)
# getFormulaSides(cbind(z1, z2) ~ x + y)
# getFormulaSides(. ~ x + y)
# getFormulaSides(cbind(z1, z2) ~ .)
# getFormulaSides(. ~ x + y, names = c("z1", "z2", "x", "y"))
# getFormulaSides(cbind(z1, z2) ~ ., names = c("z1", "z2", "x", "y"))


blockdiag <- function(x, j, back = NULL) {
  d <- dim(x)
  if (!is.null(back) & (dim(back)[1] != d[1]*j || dim(back)[2] != d[2]*j)) {
    warning("Wrong dimensions of background matrix. Set to Zero-Matrix. ")
    back <- NULL
  }
  if (is.null(back)) {
    X <- matrix(0, nrow = d[1] * j, ncol = d[2] * j)
  } else {
    X <- back
  }
  for (i in 0:(j-1)) {
    X[(i*d[1]+1):((i+1)*d[1]), (i*d[2]+1):((i+1)*d[2])] <- x
  }
  X
}
# blockdiag(matrix(1:4, 2), 3)
# blockdiag(matrix(1:8, 2), 3)
# blockdiag(matrix(1:4, 2), 3, matrix(NA, nr = 6, nc = 6))
# blockdiag(matrix(1:4, 2), 3, matrix(NA, nr = 6, nc = 7))
blockdiag_list <- function(x, back = NULL) {
  dims <- lapply(x, dim)
  dim <- c(sum(sapply(dims, getElement, 1)), sum(sapply(dims, getElement, 2)))

  if (!is.null(back) & (dim(back)[1] != dim[1] || dim(back)[2] != dim[2])) {
    warning("Wrong dimensions of background matrix. Set to Zero-Matrix. ")
    back <- NULL
  }
  if (is.null(back)) {
    X <- matrix(0, nrow = dim[1], ncol = dim[2])
  } else {
    X <- back
  }
  pos <- list(
    c(0, cumsum(sapply(dims, getElement, 1))),
    c(0, cumsum(sapply(dims, getElement, 2)))
  )
  for (i in 1:length(x)) {
    X[(pos[[1]][i]+1):pos[[1]][i+1], (pos[[2]][i]+1):pos[[2]][i+1]] <- x[[i]]
  }
  X
}
# blockdiag_list(list(matrix(1:4, 2), matrix(1:9, 3)))
# blockdiag_list(list(matrix(1:4, 2), matrix(1:8, 2)), back = matrix(NA, nr = 4, nc = 6))


removeAttributes <- function(x) {
  attr(x, "source") <- NULL
  attr(x, "class") <- NULL
  attr(x, "order") <- NULL
  attr(x, "distribution") <- NULL
  attr(x, "leftrim") <- NULL
  attr(x, "rightrim") <- NULL
  attr(x, "computation.method") <- NULL
  x
}