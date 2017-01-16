#' @title
#' Summary TLMoments
#' @description
#' Calculating and printing of summary statistics to a given TLMoments-object.
#' @param object object of TLMoments
#' @param ci.level numeric vector of length 1 giving the confidence level (default is 0.9).
#' @param ... additional arguments submitted to \code{est_cov}.
#' @return A \code{summary.TLMoments}-object, a list with dimensions \itemize{
#'  \item \code{tlm}
#'  \item \code{ci.level}
#'  \item \code{lambda.ci}
#'  \item \code{lambda.cov}
#'  \item \code{ratio.ci}
#'  \item \code{ratio.cov}
#' }
#' It is printed with \code{print.summary.TLMoments}.
#' @examples
#' tlm <- TLMoments(evd::rgev(100, shape = .2))
#' summary(tlm)
#'
#' tlm <- TLMoments(evd::rgev(100, shape = .2), rightrim = 1)
#' summary(tlm)
#'
#' tlm <- TLMoments(evd::rgev(100, shape = .2), max.order = 2, rightrim = 1)
#' summary(tlm)
#'
#' tlm <- TLMoments(matrix(evd::rgev(100, shape = .2), nc = 2))
#' summary(tlm)
#'
#' tlm <- TLMoments(matrix(evd::rgev(100, shape = .2), nc = 2), max.order = 3)
#' summary(tlm, ci = .95, distr = "gev")
#'
#' tlm <- as.TLMoments(c(15, 5, 1.3))
#' summary(tlm, distr = "gev", n = 100)
#'
#' @method summary TLMoments
#' @export
summary.TLMoments <- function(object, ci.level = .9, ...) {
  if (length(ci.level) != 1 | !is.numeric(ci.level)) stop("ci must be a numeric vector of length 1. ")
  if (!("TLMoments" %in% class(object))) stop("First argument has to be of class TLMoments. ")

  UseMethod("summary.TLMoments", object$lambdas)
}

#' @method summary.TLMoments numeric
#' @export
summary.TLMoments.numeric <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # covs
  cov <- est_cov(object, order = attr(object, "order"), ...)

  # lambda ci
  lambda_ci <- object$lambdas %-+% (u * sqrt(diag(cov$lambdas)))
  lambda_ci <- cbind(lambda_ci[, 1], object$lambdas, lambda_ci[, 2])
  colnames(lambda_ci) <- c("LCL", "lambda_hat", "UCL")

  out <- list(
    tlm = object,
    ci.level = ci.level,
    lambda.ci = lambda_ci,
    lambda.cov = cov$lambdas
  )

  # tau ci
  if (length(attr(object, "order")) >= 2) {
    ratio_ci <- object$ratios[-1] %-+% (u * sqrt(diag(as.matrix(cov$ratios))))
    out$ratio.ci <- cbind(ratio_ci[, 1], object$ratios[-1], ratio_ci[, 2])
    colnames(out$ratio.ci) <- c("LCL", "tau_hat", "UCL")
    out$ratios.cov <- cov$ratios
  }

  class(out) <- c("summary.TLMoments")
  out
}

#' @method summary.TLMoments matrix
#' @export
summary.TLMoments.matrix <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # covs
  cov <- est_cov(object, order = attr(object, "order"), ...)

  # lambda ci
  lambda_ci <- as.vector(object$lambdas) %-+% (u * sqrt(diag(cov$lambdas)))
  lambda_ci <- cbind(lambda_ci[, 1], as.vector(object$lambdas), lambda_ci[, 2])
  colnames(lambda_ci) <- c("LCL", "lambda_hat", "UCL")

  out <- list(
    tlm = object,
    ci.level = ci.level,
    lambda.ci = lambda_ci,
    lambda.cov = cov$lambdas
  )

  # tau ci
  if (length(attr(object, "order")) >= 2) {
    ratio_ci <- as.numeric(object$ratios[-1, , drop = FALSE]) %-+% (u * sqrt(diag(cov$ratios)))
    out$ratio.ci <- cbind(ratio_ci[, 1], as.numeric(object$ratios[-1, ]), ratio_ci[, 2])
    colnames(out$ratio.ci) <- c("LCL", "tau_hat", "UCL")
    out$ratios.cov <- cov$ratios
  }

  class(out) <- c("summary.TLMoments")
  out
}

#' @export
print.summary.TLMoments <- function(x, ...) {
  if (attr(x$tlm, "source")$func[1] %in% c("as.PWMs", "as.TLMoments", "as.parameters")) {
    # Theoretical data
  } else {
    # Real data
    ns <- attr(x$tlm, "source")$n
    cat(length(ns), " data row(s) with n = ", paste0(ns, collapse = ", "), ".\n", sep = "")
    cat("TL(", attr(x$tlm, "leftrim"), ",", attr(x$tlm, "rightrim"), ") calculated. \n", sep = "")
  }
  cat("\n")
  cat("Approximate ", x$ci.level, "% confidence interval of TL moments: \n", sep = "")
  print(x$lambda.ci)
  #cat("\n")
  if (!is.null(x$ratio.ci)) {
    cat("Approximate ", x$ci.level, "% confidence interval of TL moment ratios: \n", sep = "")
    print(x$ratio.ci)
    #cat("\n")
  }
  #cat("Covariance matrix of TL-Moments estimates: \n", sep = "")
  #print(x$cov)
}
