#' @title
#' Converting TL-moments to distribution parameters
#' @description
#' Converts TL-moments to distribution parameters. By now, conversions for gev, gumbel, gpd, and ln3
#' are available. Important trimming options are calculated through known formulas (see references for
#' some of them), other options are calculated through a numerical optimization. In this case there's a
#' warning informing about the experimental nature of this feature.
#' @param tlm object returned by TLMoments
#' @param distr character object defining the distribution. Supported types are
#' "gev", "gum", "gpd", and "ln3"
#' @param ... additional arguments
#' @return Numeric vector, matrix, list, or data.frame of parameter estimates with
#' class \code{parameters}. The object contains the following attributes: \itemize{
#'  \item \code{distribution}: a character indicating the used distribution
#'  \item \code{source}: a list with background information (used function, data, n, formula, trimmings; mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#' @references Elamir, E. A. H. (2010). Optimal choices for trimming in trimmed L-moment method. Applied Mathematical Sciences, 4(58), 2881-2890.
#' @references Fischer, S., Fried, R., & Schumann, A. (2015). Examination for robustness of parametric estimators for flood statistics in the context of extraordinary extreme events. Hydrology and Earth System Sciences Discussions, 12, 8553-8576.
#' @references Hosking, J. R. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 105-124.
#' @references Hosking, J. R. M. (2007). Some theory and practical uses of trimmed L-moments. Journal of Statistical Planning and Inference, 137(9), 3024-3039.
#' @seealso \code{\link{PWMs}}, \code{\link{TLMoments}}, \code{\link{quantiles}}, \code{\link{summary.parameters}}, \code{\link{as.parameters}}. Built-in distributions: \code{\link{pgev}}, \code{\link{pgum}}, \code{\link{pgpd}}, \code{\link{pln3}}.
#' @examples
#' xmat <- matrix(rgev(100, shape = .2), nc = 4)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:2], each = 50),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' tlm <- TLMoments(xvec, leftrim = 0, rightrim = 0)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(xmat, leftrim = 1, rightrim = 1)
#' parameters(tlm, "gum")
#'
#' tlm <- TLMoments(xlist)
#' parameters(tlm, "gpd")
#'
#' tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(xdat, hq ~ station + season, leftrim = 0, rightrim = 2)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(xdat, hq ~ station + season, leftrim = 0, rightrim = 2)
#' parameters(tlm, "ln3")
#'
#'
#' # Numerical calculation
#'
#' tlm <- TLMoments(rgum(200, loc = 5, scale = 2), leftrim = 1, rightrim = 4)
#' parameters(tlm, "gum")
#'
#' tlm <- TLMoments(rgev(200, loc = 10, scale = 5, shape = .4), leftrim = 2, rightrim = 2)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(rln3(200, loc = 3, scale = 1.5, shape = 2), leftrim = 0, rightrim = 1)
#' parameters(tlm, "ln3")
#'
#' # It's A LOT slower:
#' # system.time(replicate(500,
#' #   parameters(TLMoments(rgum(100, loc = 5, scale = 2), 1, 1), "gum")
#' # ))[3]
#' # system.time(replicate(500,
#' #   parameters(TLMoments(rgum(100, loc = 5, scale = 2), 1, 2), "gum")
#' # ))[3]
#'
#'
#' # Using magrittr
#' library(magrittr)
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 0) %>%
#'  parameters("gpd")
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 0) %>%
#'  parameters("gpd", u = 10)
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 1) %>%
#'  parameters("gpd")
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 2) %>%
#'  parameters("gpd")
#'
#' @export
parameters <- function(tlm, distr, ...) {
  if (!("TLMoments" %in% class(tlm))) stop("tlm has to be of class TLMoments")

  valid.distr <- c("gev", "gum", "gpd", "ln3")
  if (!(distr %in% valid.distr))
    stop(paste0("distr has to be one of ", paste(valid.distr, collapse = ", "), "."))

  # Check for sufficient data:
  if (distr == "gev") {
    if (!all(c(1, 2, 3) %in% attr(tlm, "order"))) stop("Parameter calculation of \"gev\" needs at least the first three TL-moments")
  } else if (distr == "gum") {
    if (!all(c(1, 2) %in% attr(tlm, "order"))) stop("Parameter calculation of \"gum\" needs at least the first two TL-moments")
  } else if (distr == "gpd" & "u" %in% names(list(...))) {
    if (!all(c(1, 2) %in% attr(tlm, "order"))) stop("Parameter calculation of \"gpd\" needs at least the first two TL-moments")
  } else if (distr == "gpd") {
    if (!all(c(1, 2, 3) %in% attr(tlm, "order"))) stop("Parameter calculation of \"gpd\" needs at least the first three TL-moments")
  } else if (distr == "ln3") {
    if (!all(c(1, 2, 3) %in% attr(tlm, "order"))) stop("Parameter calculation of \"ln3\" needs at least the first three TL-moments")
  }

  UseMethod("parameters", tlm$lambdas)
}

#' @method parameters matrix
#' @export
parameters.matrix <- function(tlm, distr, ...) {
  out <- sapply(1L:ncol(tlm$lambdas), function(i) {
    tmp <- list(lambdas = tlm$lambdas[, i], ratios = tlm$ratios[, i])
    attributes(tmp) <- attributes(tlm)
    parameters.numeric(tmp, distr, ...)
  })

  attr(out, "distribution") <- distr
  attr(out, "source") <- attributes(tlm)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "parameters")
  attr(out, "source")$trimmings <- c(attr(tlm, "leftrim"), attr(tlm, "rightrim"))
  attr(out, "source")$lambdas <- tlm$lambdas
  attr(out, "source")$max.order <- max(attr(tlm, "order"))
  class(out) <- c("parameters", "matrix")
  out
}

#' @method parameters list
#' @export
parameters.list <- function(tlm, distr, ...) {
  out <- lapply(1L:length(tlm$lambdas), function(i) {
    tmp <- list(lambdas = tlm$lambdas[[i]], ratios = tlm$ratios[[i]])
    attributes(tmp) <- attributes(tlm)
    parameters.numeric(tmp, distr, ...)
  })

  attr(out, "distribution") <- distr
  attr(out, "source") <- attributes(tlm)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "parameters")
  attr(out, "source")$trimmings <- c(attr(tlm, "leftrim"), attr(tlm, "rightrim"))
  attr(out, "source")$lambdas <- tlm$lambdas
  attr(out, "source")$max.order <- max(attr(tlm, "order"))
  class(out) <- c("parameters", "list")
  out
}

#' @method parameters data.frame
#' @export
parameters.data.frame <- function(tlm, distr, ...) {

  out <- sapply(1L:nrow(tlm$lambdas), function(i) {
    tmp <- list(lambdas = unlist(tlm$lambdas[i, grep("L[0-9]*", names(tlm$lambdas))]),
                ratios = unlist(tlm$ratios[i, grep("T[0-9]*", names(tlm$ratios))]))
    attributes(tmp)[c("class", "order", "leftrim", "rightrim")] <- attributes(tlm)[c("class", "order", "leftrim", "rightrim")]
    parameters(tmp, distr, ...)
  })
  rhs <- dimnames(attr(terms(attr(tlm, "source")$formula), "factors"))[[2]]
  out <- cbind(tlm$lambdas[rhs], as.data.frame(t(out)))

  attr(out, "distribution") <- distr
  attr(out, "source") <- attributes(tlm)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "parameters")
  attr(out, "source")$trimmings <- c(attr(tlm, "leftrim"), attr(tlm, "rightrim"))
  attr(out, "source")$lambdas <- tlm$lambdas
  attr(out, "source")$max.order <- max(attr(tlm, "order"))
  class(out) <- c("parameters", "data.frame")
  out
}

#' @method parameters numeric
#' @export
parameters.numeric <- function(tlm, distr, ...) {

  leftrim <- attr(tlm, "leftrim")
  rightrim <- attr(tlm, "rightrim")

  out <- switch(
    distr,
    gev = gev.est(tlm$lambdas["L1"], tlm$lambdas["L2"], tlm$ratios["T3"], leftrim, rightrim),
    gum = gum.est(tlm$lambdas["L1"], tlm$lambdas["L2"], leftrim, rightrim),
    gpd = gpd.est(tlm$lambdas["L1"], tlm$lambdas["L2"], tlm$ratios["T3"], leftrim, rightrim, ...),
    ln3 = ln3.est(tlm$lambdas["L1"], tlm$lambdas["L2"], tlm$ratios["T3"], leftrim, rightrim),
    stop("Not implemented. ")
  )
  out <- unlist(out)

  attr(out, "distribution") <- distr
  attr(out, "source") <- attributes(tlm)$source
  attr(out, "source")$func <- c(attr(out, "source")$func, "parameters")
  attr(out, "source")$trimmings <- c(attr(tlm, "leftrim"), attr(tlm, "rightrim"))
  attr(out, "source")$lambdas <- tlm$lambdas
  attr(out, "source")$max.order <- max(attr(tlm, "order"))
  class(out) <- c("parameters", "numeric")
  out
}

#' @export
print.parameters <- function(x, ...) {
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


#### GUM ####
gum.est <- function(l1, l2, s, t) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gum.tl00.est(l1, l2),
         `0-1` = gum.tl01.est(l1, l2),
         `1-0` = gum.tl10.est(l1, l2),
         `1-1` = gum.tl11.est(l1, l2),
         gum.numerical.est(l1, l2, s, t))
}
# TL(0,0)=L
gum.tl00.est <- function(l1, l2) {
  scale <- l2 / log(2)
  loc <- l1 - 0.577 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
# TL(0,1), TL(1,0), TL(1,1) from Elamir, 2010: Optimal Choices for Trimming in Trimmed L-moment Method
gum.tl01.est <- function(l1, l2) {
  scale <- l2 / 0.431
  loc <- l1 + 0.116 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
gum.tl10.est <- function(l1, l2) {
  scale <- l2 / 0.607
  loc <- l1 - 1.269 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
gum.tl11.est <- function(l1, l2) {
  scale <- l2 / 0.353
  loc <- l1 - .459 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
gum.numerical.est <- function(l1, l2, s, t) {
  #warning("Calculating numerical solution. ")

  # scale
  f1 <- function(scale) calcTLMom(2, s, t, qgum, loc = 0, scale = scale)[2] - l2
  i <- 0; while (sign(f1(10^-i)) == sign(f1(10^i))) i <- i+1
  u <- uniroot(f1, c(10^-i, 10^i))
  scale <- u$root
  # loc
  f2 <- function(loc) calcTLMom(2, s, t, qgum, loc = loc, scale = scale)[1] - l1
  i <- 0; while (sign(f2(-10^i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "computation.method") <- "numerical"
  return(out)
}


#### GEV ####
gev.est <- function(l1, l2, t3, s, t) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gev.tl00.est(l1, l2, t3),
         `0-1` = gev.tl01.est(l1, l2, t3),
         `0-2` = gev.tl02.est(l1, l2, t3),
         `1-1` = gev.tl11.est(l1, l2, t3),
         `1-0` = gev.tl10.est(l1, l2, t3),
         gev.numerical.est(l1, l2, t3, s, t))
}
# TL(0,0)=L der GEV
gev.tl00.est <- function(l1, l2, t3) {
  z <- 2/(3+t3) - log(2)/log(3)
  k <- 7.8590 * z + 2.9554 * z^2

  # Hosking (1990)
  gk <- gamma(1+k)
  al <- l2*k/((1-2^-k)*gk)
  xi <- l1 + al*(gk-1)/k

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
# TL(0,1) der GEV
gev.tl01.est <- function(l1, l2, t3) {
  z <- 1/(2+t3) * 10/9 - (2*log(2)-log(3))/(3*log(3)-2*log(4))
  k <- 8.567394*z - 0.675969*z^2

  gk <- gamma(k)
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2/3 * 1/gk * 1/(V3 - 2*V2 + 1)
  xi <- l1 - al/k - al*gk*(V2 - 2)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
# TL(0,2) der GEV
gev.tl02.est <- function(l1, l2, t3) {

  z <- (t3+5/3) * 6/5 - (3*log(5)-8*log(4)+6*log(3))/(log(4)-3*log(3)+3*log(2))
  k <- -2.468959*z + 1.130074*z^2 - 0.635912*z^3

  gk <- gamma(k)
  V4 <- (1/4)^k
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2 * 1/gk * 1/(-4*V4 + 12*V3 - 12*V2 + 4)
  xi <- l1 - al/k - al*gk*(-V3 + 3*V2 - 3)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
# TL(1,0) der GEV
gev.tl10.est <- function(l1, l2, t3) {

  z <- (t3-40/3) * 9/20 - (log(3)-log(4))/(log(2)-log(3))
  k <- -9.066941*z - 3.374925*z^2 - 0.303208*z^3

  gk <- gamma(k)
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2 * 1/gk * 1/(3*V2 - 3*V3)
  xi <- l1 - al/k + al*gk*V2

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
# TL(1,1) der GEV
gev.tl11.est <- function(l1, l2, t3) {

  z <- t3 * 9/20 + (log(3)-2*log(4)+log(5))/(log(2)-2*log(3)+log(4))
  k <- 25.3171*z - 91.5507*z^2 + 110.0626*z^3 - 46.5518*z^4

  gk <- gamma(k)
  V4 <- (1/4)^k
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 1/gk * 1/(3*V2 - 6*V3 + 3*V4)
  xi <- l1 - al/k - al*gk*(-3*V2 + 2*V3)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}
gev.numerical.est <- function(l1, l2, t3, s, t) {
  #warning("Calculating numerical solution. ")

  # shape
  f1 <- function(shape) {
    a <- calcTLMom(3, s, t, qgev, loc = 0, scale = 1, shape = shape)
    a[3]/a[2] - t3
  }
  i <- 0; while (sign(f1(-10^i)) == sign(f1(10^i))) i <- i+1
  u <- uniroot(f1, c(-10^i, 10^i))
  shape <- u$root
  # scale
  f2 <- function(scale) calcTLMom(2, s, t, qgev, loc = 0, scale = scale, shape = shape)[2] - l2
  i <- 0; while (sign(f2(10^-i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(10^-i, 10^i))
  scale <- u$root
  # loc
  f3 <- function(loc) calcTLMom(1, s, t, qgev, loc = loc, scale = scale, shape = shape)[1] - l1
  i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
  u <- uniroot(f3, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "numerical"
  return(out)
}

#### LN3 ####
ln3.est <- function(l1, l2, t3, s, t) {
  #warning("Calculating numerical solution. ")

  # 1. step, T3 <=> shape
  f1 <- function(shape) {
    a <- calcTLMom(3, s, t, qln3, loc = 0, scale = 1, shape = shape)
    a[3]/a[2] - t3
  }
  i <- 0; while (sign(f1(10^-i)) == sign(f1(10^i))) i <- i+1
  u <- uniroot(f1, c(10^-i, 10^i))
  shape <- u$root
  # 2. step, L2 <=> scale
  f2 <- function(scale) calcTLMom(2, s, t, qln3, loc = 0, scale = scale, shape = shape)[2] - l2
  i <- 0; while (sign(f2(10^-i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(10^-i, 10^i))
  scale <- u$root
  # 3. step L1 <=> loc
  f3 <- function(loc) calcTLMom(1, s, t, qln3, loc = loc, scale = scale, shape = shape)[1] - l1
  i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
  u <- uniroot(f3, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "numerical"
  return(out)
}


#### GPD ####
gpd.est <- function(l1, l2, t3, s, t, ...) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gpd.tl00.est(l1, l2, t3, ...),
         `0-1` = gpd.tl01.est(l1, l2, t3, ...),
         `0-2` = gpd.tl02.est(l1, l2, t3, ...),
         gpd.numerical.est(l1, l2, t3, s, t, ...))
}
# GPD
gpd.tl00.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    shape <- -(l1-u)/l2 + 2
    scale <- (l1-u) * (-shape+1)
    loc <- u
  } else {
    shape <- (3*t3 - 1)/(t3 + 1)
    scale <- l2 * (shape-1) * (shape-2)
    loc <- l1 + scale/(shape-1)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}

gpd.tl01.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    # from Hosking, 2007.
    shape <- -3/2 * (l1-u)/l2 + 3
    scale <- (l1-u) * (-shape+2)
    loc <- u
  } else {
    # from Fischer, 2015.
    shape <- (36*t3 - 8)/(9*t3 + 8)
    scale <- 2/3 * l2 * (shape-2) * (shape-3)
    loc <- l1 + scale/(shape-2)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}

gpd.tl02.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    # from Hosking, 2007.
    shape <- -2 * (l1-u)/l2 + 4
    scale <- (l1-u) * (-shape+3)
    loc <- u
  } else {
    # from Fischer, 2015.
    shape <- (30*t3 - 5)/(6*t3 + 5)
    scale <- 1/2 * l2 * (shape-3) * (shape-4)
    loc <- l1 + scale/(shape-3)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "formula"
  return(out)
}

gpd.numerical.est <- function(l1, l2, t3, s, t, u = NULL) {
  #warning("Calculating numerical solution. ")

  if (!is.null(u) & is.numeric(u)) {

    # shape
    f1u <- function(shape, u) {
      a <- calcTLMom(2, s, t, qgpd, loc = u, scale = 1, shape = shape)
      a[2] - l2
    }
    i <- 0; while (sign(f1u(-10^i, u)) == sign(f1u(10^i, u))) i <- i+1
    ur <- uniroot(f1u, c(-10^i, 10^i), u = u)
    shape <- ur$root
    # scale
    f2u <- function(scale, u) calcTLMom(1, s, t, qgpd, loc = u, scale = scale, shape = shape)[1] - l1
    i <- 0; while (sign(f2u(10^-i, u)) == sign(f2u(10^i, u))) i <- i+1
    ur <- uniroot(f2u, c(10^-i, 10^i), u = u)
    scale <- ur$root
    # loc
    loc <- u

  } else {

    # shape
    f1 <- function(shape) {
      a <- calcTLMom(3, s, t, qgpd, loc = 0, scale = 1, shape = shape)
      a[3]/a[2] - t3
    }
    i <- 0; while (sign(f1(-10^i)) == sign(f1(10^i))) i <- i+1
    ur <- uniroot(f1, c(-10^i, 10^i))
    shape <- ur$root
    # scale
    f2 <- function(scale) calcTLMom(2, s, t, qgpd, loc = 0, scale = scale, shape = shape)[2] - l2
    i <- 0; while (sign(f2(10^-i)) == sign(f2(10^i))) i <- i+1
    ur <- uniroot(f2, c(10^-i, 10^i))
    scale <- ur$root
    # loc
    f3 <- function(loc) calcTLMom(1, s, t, qgpd, loc = loc, scale = scale, shape = shape)[1] - l1
    i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
    ur <- uniroot(f3, c(-10^i, 10^i))
    loc <- ur$root

  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "computation.method") <- "numerical"
  return(out)
}
