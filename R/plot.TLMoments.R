#' @title
#' L-Moment-ratio diagram
#' @description
#' Generates a ggplot2 object containing a scatterplot of TL skewness and TL kurtosis
#' as well as the theoretical curves and points of several distributions (today: GEV, GPD,
#' LN3, GUM, EXP, NORM).
#'
#' @param x object of TLMoments.
#' @param distr character indicating the plotted theoretical distributions, see details.
#' @param ... additional arguments, not used at the moment.
#' @details
#' \code{distr}: this can either be a vector containing the abbreviations of the
#' theoretical distributions (gev, gpd, ln3, gum, exp, or norm) or one of the
#' shortcuts \"all\" (default), \"only-lines\", or \"only-points\" that indicate
#' all distributions, all distributions displayed as lines (i.e. gev, gpd, ln3),
#' or all distributions displayed as points (ie. gum, exp, norm), respectively.
#'
#' Values of theoretical distributions are pre-calculated for several trimmings.
#' If other trimmings are selected this results in a (small) delay for calculation.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' xmat <- matrix(rgev(1000, shape = .3), nc = 10)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:10], each = 100),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' library(ggplot2)
#' plot(TLMoments(xvec))
#' plot(TLMoments(xlist), distr = c("gev", "gum"))
#' plot(TLMoments(xmat), distr = "only-points")
#' plot(TLMoments(xmat), distr = "only-lines") + scale_colour_viridis_d()
#' plot(TLMoments(xmat, 0, 1))
#' plot(TLMoments(xmat, 0, 1)) + coord_cartesian(xlim = c(-.05, .4), ylim = c(.05, .2))
#' plot(TLMoments(xdat, hq ~ station, 1, 0))
#' }
#' @method plot TLMoments
#' @export
plot.TLMoments <- function(x, distr, ...) {
  if (!("TLMoments" %in% class(x))) stop("First argument has to be of class TLMoments. ")
  if (!all(c(3, 4) %in% attr(x, "order"))) stop("Object must contain T2 and T3. ")
  UseMethod("plot.TLMoments", x$lambdas)
}

#' @method plot.TLMoments numeric
#' @export
plot.TLMoments.numeric <- function(x, distr = "all", ...) {
  lmrdiagram(x$ratios[3], x$ratios[4],
             trim = c(attr(x, "leftrim"), attr(x, "rightrim")),
             distr = distr)
}

#' @method plot.TLMoments matrix
#' @export
plot.TLMoments.matrix <- function(x, distr = "all", ...) {
  lmrdiagram(x$ratios[3, ], x$ratios[4, ],
             trim = c(attr(x, "leftrim"), attr(x, "rightrim")),
             distr = distr)
}

#' @method plot.TLMoments list
#' @export
plot.TLMoments.list <- function(x, distr = "all", ...) {
  lmrdiagram(vapply(x$ratios, getElement, "T3", FUN.VALUE = numeric(1)),
             vapply(x$ratios, getElement, "T4", FUN.VALUE = numeric(1)),
             trim = c(attr(x, "leftrim"), attr(x, "rightrim")),
             distr = distr)
}

#' @method plot.TLMoments data.frame
#' @export
plot.TLMoments.data.frame <- function(x, distr = "all", ...) {
  lmrdiagram(x$ratios$T3, x$ratios$T4,
             trim = c(attr(x, "leftrim"), attr(x, "rightrim")),
             distr = distr)
}

lmrdiagram <- function(t3, t4, trim = c(0, 0), distr = c("all")) {

  if (length(distr) == 1 && distr == "all")
    distr <- c("gum", "exp", "norm", "gev", "gpd", "ln3")
  if (length(distr) == 1 && distr == "only-lines")
    distr <- c("gev", "gpd", "ln3")
  if (length(distr) == 1 && distr == "only-points")
    distr <- c("gum", "exp", "norm")

  tlmr <- tlmomentratios[tlmomentratios$leftrim == trim[1] & tlmomentratios$rightrim == trim[2], ]
  if (nrow(tlmr) == 0) {
    tlmr <- getTLMomsByDistr(distr, trim)
  }
  tlmr_points <- tlmr[tlmr$distr %in% intersect(c("gum", "exp", "norm"), distr), ]
  tlmr_lines <- tlmr[tlmr$distr %in% intersect(c("gev", "gpd", "ln3"), distr), ]

  # ggplot-mode
  lab_pref <- ifelse(all(trim == c(0, 0)), "L", paste0("TL(", paste(trim, collapse = ","), ")"))

  ggplot2::ggplot(
    data.frame(T3 = t3, T4 = t4),
    ggplot2::aes_(~T3, ~T4)
  ) +
    ggplot2::labs(x = paste(lab_pref, "skewness"), y = paste(lab_pref, "kurtosis")) +
    ggplot2::coord_cartesian(xlim = range(t3), ylim = range(t4)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes_(~T3, ~T4, colour = ~distr, linetype = ~distr), data = tlmr_lines) +
    ggplot2::geom_point(ggplot2::aes_(~T3, ~T4, shape = ~distr), data = tlmr_points, size = 4) #+
    #ggplot2::geom_point(ggplot2::aes(T3, T4, colour = distr), data = ticks) +
    #ggplot2::annotate("line", x = a0$T3+c(-.002, .002), y = a0$T4-d*c(-.002, .002), colour = "#E41A1C")

  # ticks <- subset(tlmr_lines, shape %in.equal% seq(-.9, .9, .1))
  #
  #
  # a0 <- subset(tlmr_lines, distr == "gev" & shape %in.equal% .3)
  # a <- subset(tlmr_lines, distr == "gev" & shape %in.equal% c(.29, .31))
  # d <- (a$T3[2]-a$T3[1])/(a$T4[2]-a$T4[1])

}
