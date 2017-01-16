#' @title
#' Estimate the covariance matrix of parameter estimations
#' @description
#' description not done yet
#' @param x numeric vector or matrix containing data OR an object of parameters.
#' @param distr character indicating the distribution from which the parameters are calculated. If x is parameters-object, distr does not need to be specified.
#' @param leftrim,rightrim lower and upper trimming parameter used for parameter calculation, have to be non-negative integers.
#' @param np.cov boolean, if TRUE no parametric assumptions are used to calculate the covariance matrix (default FALSE).
#' @param set.n hypothetical data length n if theoretical values are given.
#' @param ... additional arguments, ignored.
#' @return numeric matrix
#' @examples
#' ### Numeric vectors
#' x <- evd::rgev(500, shape = .2)
#'
#' parameters(TLMoments(x), "gev")
#' est_paramcov(x, "gev", 0, 0)
#' #cov(t(replicate(10000, parameters(TLMoments(evd::rgev(500, shape = .2)), "gev"))))
#'
#' parameters(TLMoments(x, rightrim = 1), "gev")
#' est_paramcov(x, "gev", 0, 1)
#' #cov(t(replicate(10000, parameters(TLMoments(evd::rgev(500, shape = .2), rightrim = 1), "gev"))))
#'
#' parameters(TLMoments(x, rightrim = 2), "gev")
#' est_paramcov(x, "gev", 0, 2)
#' #cov(t(replicate(10000, parameters(TLMoments(evd::rgev(500, shape = .2), rightrim = 2), "gev"))))
#'
#' ### Numeric matrices
#' x <- matrix(evd::rgev(600, shape = .2), nc = 3)
#'
#' parameters(TLMoments(x), "gev")
#' est_paramcov(x, "gev", 0, 0)
#' #cov(t(replicate(5000,
#' #  as.vector(parameters(TLMoments(matrix(evd::rgev(600, shape = .2), nc = 3)), "gev")))
#' #))
#'
#' ### parameters-object
#' x <- as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' est_paramcov(x)
#' est_paramcov(x, leftrim = 0, rightrim = 0)
#' est_paramcov(x, leftrim = 0, rightrim = 0, set.n = 100)
#' # distr-argument can be neglected
#'
#' @rdname est_paramcov
#' @export
est_paramcov <- function(x, distr, leftrim = 0L, rightrim = 0L, ...) {
  if ("parameters" %in% class(x)) {
    distr <- attr(x, "distribution")
    if (!any(is.na(attr(x, "source")$trimming))) {
      leftrim <- attr(x, "source")$trimming[1]
      rightrim <- attr(x, "source")$trimming[2]
    }
  }

  # TODO: General error catches
  if (!are.integer.like(leftrim, rightrim))
    stop("leftrim and rightrim have to be integer types!")

  if (!(distr %in% c("evd::gev", "gev", "gumbel", "gpd")))
    stop("distr has to be \"gev\", \"gumbel\", or \"gpd\"")

  if (!(distr %in% c("evd::gev", "gev")))
    stop("only GEV for now")

  if (leftrim != 0 || !(rightrim %in% c(0, 1, 2)))
    stop("only (0,0), (0,1), or (0,2) for now")

  UseMethod("est_paramcov")
}

#' @rdname est_paramcov
#' @method est_paramcov numeric
#' @export
est_paramcov.numeric <- function(x, distr, leftrim = 0L, rightrim = 0L, np.cov = FALSE, ...) {

  if (np.cov) {
    tlmcov <- est_tlmcov(x, leftrim = leftrim, rightrim = rightrim, order = 1:3, ratio.cov = FALSE)
  } else {
    tlmcov <- est_tlmcov(x, leftrim = leftrim, rightrim = rightrim, order = 1:3, distr = distr, ratio.cov = FALSE)
  }
  tlm <- TLMoments(x, leftrim = leftrim, rightrim = rightrim, max.order = 3, na.rm = TRUE)
  A <- CovTLtoPara(distr, l2 = tlm$lambdas[2], l3 = tlm$lambdas[3], leftrim = leftrim, rightrim = rightrim)
  r <- A %*% tlmcov %*% t(A)
  rownames(r) <- colnames(r) <- c("loc", "scale", "shape")

  #attr(r, "source") <- list(call = match.call(), distr = distr, leftrim = leftrim, rightrim = rightrim, np.cov = np.cov)
  r
}

#' @rdname est_paramcov
#' @method est_paramcov matrix
#' @export
est_paramcov.matrix <- function(x, distr, leftrim = 0L, rightrim = 0L, np.cov = FALSE, ...) {

  if (np.cov) {
    tlmcov <- est_tlmcov(x, leftrim = leftrim, rightrim = rightrim, order = 1:3, ratio.cov = FALSE)
  } else {
    tlmcov <- est_tlmcov(x, leftrim = leftrim, rightrim = rightrim, order = 1:3, distr = distr, ratio.cov = FALSE)
  }
  tlm <- TLMoments(x, leftrim = leftrim, rightrim = rightrim, max.order = 3, na.rm = TRUE)
  A <- lapply(1:ncol(x), function(i) CovTLtoPara(distr, l2 = tlm$lambdas[2, i], l3 = tlm$lambdas[3, i], leftrim = leftrim, rightrim = rightrim))
  A <- blockdiag_list(A)

  r <- A %*% tlmcov %*% t(A)
  rownames(r) <- colnames(r) <- paste0(rep(c("loc", "scale", "shape"), ncol(x)), "_", rep(1:ncol(x), each = 3))

  #attr(r, "source") <- list(call = match.call(), distr = distr, leftrim = leftrim, rightrim = rightrim, np.cov = np.cov)
  r
}

#' @rdname est_paramcov
#' @method est_paramcov parameters
#' @export
est_paramcov.parameters <- function(x,
                                    distr = attr(x, "distribution"),
                                    leftrim = attr(x, "source")$trimmings[1],
                                    rightrim = attr(x, "source")$trimmings[2],
                                    set.n = NA,
                                    ...) {
  if (!("numeric" %in% class(x)))
    stop("x must be a numeric parameters-object. ")
  if (!(attr(x, "source")$func[1] %in% c("as.PWMs", "as.TLMoments", "as.parameters")))
    stop("est_paramcov.parameters only for theoretical values. ")

  if (are.integer.like(leftrim) && is.na(rightrim)) {
    rightrim <- 0
  }
  if (are.integer.like(rightrim) && is.na(leftrim)) {
    leftrim <- 0
  }
  if (any(is.na(c(leftrim, rightrim)))) {
    warning("No trimmings found, set to leftrim = 0, rightrim = 0. ")
    leftrim <- rightrim <- 0
  }

  if (is.na(set.n) | !is.numeric(set.n)) {
    warning("Missing or invalid set.n argument. Giving results for n = 100. ")
    n <- 100
  } else { n <- set.n }

  tlm <- TLMoments(x, leftrim = leftrim, rightrim = rightrim, max.order = 3)
  tlmcov <- est_tlmcov(tlm, leftrim = leftrim, rightrim = rightrim, distr = distr, set.n = n, ratio.cov = FALSE)

  A <- CovTLtoPara(distr, tlm$lambdas[2], tlm$lambdas[3], leftrim, rightrim)
  r <- A %*% tlmcov %*% t(A)
  rownames(r) <- colnames(r) <- c("loc", "scale", "shape")

  #attr(r, "source") <- list(call = match.call(), distr = distr, leftrim = leftrim, rightrim = rightrim)
  r
}
