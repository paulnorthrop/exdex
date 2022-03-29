# ============================== confint.dgaps ============================== #

#' Confidence intervals for the extremal index \eqn{\theta}
#'
#' \code{confint} method for objects of class \code{c("dgaps", "exdex")}.
#' Computes confidence intervals for \eqn{\theta} based on an object returned
#' from \code{\link{dgaps}}.  Two types of interval may be returned:
#' (a) intervals based on approximate large-sample normality of the estimator
#' of \eqn{\theta}, which are symmetric about the point estimate,
#' and (b) likelihood-based intervals.
#'
#' @param object An object of class \code{c("dgaps", "exdex")}, returned by
#'   \code{\link{dgaps}}.
#' @param parm Specifies which parameter is to be given a confidence interval.
#'   Here there is only one option: the extremal index \eqn{\theta}.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param interval_type A character scalar: \code{"norm"} for intervals of
#'   type (a), \code{"lik"} for intervals of type (b).
#' @param conf_scale A character scalar.  If \code{interval_type = "norm"} then
#'   \code{conf_scale} determines the scale on which we use approximate
#'   large-sample normality of the estimator to estimate confidence intervals.
#'
#'   If \code{conf_scale = "theta"}
#'   then confidence intervals are estimated for \eqn{\theta} directly.
#'   If \code{conf_scale = "log"} then confidence intervals are first
#'   estimated for \eqn{\log\theta}{log\theta} and then transformed back
#'   to the \eqn{\theta}-scale.
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any confidence limits that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  Otherwise,
#'   limits that are greater than 1 may be obtained.
#'   If \code{constrain = TRUE} then any lower confidence limits that are
#'   less than 0 are set to 0.
#' @param ... Further arguments. None are used currently.
#' @details Two type of interval are calculated: (a) an interval based on the
#'   approximate large sample normality of the estimator of \eqn{\theta}
#'   (if \code{conf_scale = "theta"}) or of \eqn{\log\theta}{log\theta}
#'   (if \code{conf_scale = "log"}) and (b) a likelihood-based interval,
#'   based on the approximate large sample chi-squared, with 1 degree of
#'   freedom, distribution of the log-likelihood ratio statistic.
#' @return A list of class c("confint_dgaps", "exdex") containing the
#'   following components.
#'   \item{cis}{A matrix with columns giving the lower and upper confidence
#'   limits. These are labelled as (1 - level)/2 and 1 - (1 - level)/2 in
#'   \% (by default 2.5\% and 97.5\%).
#'   The row names indicate the type of interval:
#'   \code{norm} for intervals based on large sample normality and \code{lik}
#'   for likelihood-based intervals.}
#'   \item{call}{The call to \code{spm}.}
#'   \item{object}{The input object \code{object}.}
#'   \item{level}{The input \code{level}.}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197–213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' theta <- dgaps(newlyn, u)
#' cis <- confint(theta)
#' cis
#' plot(cis)
#' @export
confint.dgaps <- function (object, parm = "theta", level = 0.95,
                           interval_type = c("both", "norm", "lik"),
                           conf_scale = c("theta", "log"), constrain = TRUE,
                           ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  Call <- match.call(expand.dots = TRUE)
  parm <- match.arg(parm)
  if (level <= 0 || level >= 1) {
    stop("''level'' must be in (0, 1)")
  }
  interval_type <- match.arg(interval_type)
  conf_scale <- match.arg(conf_scale)
  theta <- coef(object)
  se <- sqrt(vcov(object))
  if (interval_type == "norm" || interval_type == "both") {
    # Symmetric confidence intervals, based on large sample normal theory
    # The intervals are (initially) centred on the unconstrained estimate of
    # theta, which may be greater than 1
    z_val <- stats::qnorm(1 - (1 - level) / 2)
    if (conf_scale == "theta") {
      lower <- theta - z_val * se
      upper <- theta + z_val * se
    } else {
      lower <- exp(log(theta) - z_val * se / theta)
      upper <- exp(log(theta) + z_val * se / theta)
    }
  } else {
    lower <- upper <- NULL
  }
  if (interval_type == "lik" || interval_type == "both") {
    # Likelihood-based confidence intervals.
    temp <- dgaps_conf_int(theta_mle = theta, ss = object$ss,
                           conf = 100 * level)
    lower <- c(lower, temp[1])
    upper <- c(upper, temp[2])
  }
  # Constrain the intervals to (0, 1] if required
  if (constrain) {
    lower <- pmin(lower, 1)
    lower <- pmax(lower, 0)
    upper <- pmin(upper, 1)
  }
  temp <- cbind(lower, upper)
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1), "%")
  colnames(temp) <- pct
  row_names <- interval_type
  if (row_names == "both") {
    row_names <- c("norm", "lik")
  }
  rownames(temp) <- row_names
  temp <- list(cis = temp, call = Call, object = object, level = level)
  class(temp) <- c("confint_dgaps", "exdex")
  return(temp)
}

# ========================= Method for confint_dgaps ======================== #

# ---------------------------- plot.confint_dgaps --------------------------- #

#' Plot diagnostics for a confint_dgaps object
#'
#' \code{plot} method for an objects of class
#' \code{c("confint_dgaps", "exdex")}.
#'
#' @param x an object of class \code{c("confint_dgaps", "exdex")}, a result of
#'   a call to \code{\link{confint.dgaps}}.
#' @param y Not used.
#' @param ... Further arguments to be passed to
#'   \code{\link[chandwich]{plot.confint}}.
#' @return Nothing is returned.
#' @seealso \code{\link{confint.dgaps}}: \code{confint} method for
#'   class \code{c("dgaps", "exdex")}.
#' @section Examples:
#' See the examples in \code{\link{confint.dgaps}}.
#' @export
plot.confint_dgaps <- function(x, y = NULL, ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  if (!("lik" %in% rownames(x$cis))) {
    stop("Plot method not available when interval_type = ''norm''")
  }
  prof_ci <- x$cis[rownames(x$cis) == "lik"]
  theta <- seq(0.99 * prof_ci[1], 1.01 * prof_ci[2], length = 100)
  kloglik <- function(theta) {
    do.call(dgaps_loglik, c(list(theta = theta), x$object$ss))
  }
  prof_lik <- vapply(theta, kloglik, 0.0)
  my_plot <- function(xx, y, ..., type = "l", xlab = expression(theta),
                      ylab = "log-likelihood",
                      main = paste0(100 * x$level, "% confidence interval")) {
    graphics::plot(xx, y, ..., type = type, xlab = xlab, ylab = ylab,
                   main = main)
  }
  my_plot(theta, prof_lik, ...)
  cutoff <- x$object$max_loglik - stats::qchisq(x$level, 1) / 2
  graphics::abline(h = cutoff)
  graphics::axis(1, at = prof_ci,  labels = round(prof_ci, 3), tick = FALSE,
                 mgp = c(3, 0.15, 0))
  return(invisible())
}

# --------------------------- print.confint_dgaps --------------------------- #

#' Print method for a confint_dgaps object
#'
#' \code{print} method for class \code{c("confint_dgaps", "exdex")}.
#'
#' @param x an object of class \code{c("confint_dgaps", "exdex")}, a result of
#'   a call to \code{\link{confint.dgaps}}.
#' @param ... Additional optional arguments to be passed to
#'   \code{\link{print.default}}
#' @details Prints the matrix of confidence intervals for \eqn{\theta}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{dgaps}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.dgaps}}: \code{confint} method for
#'   class \code{"dgaps"}.
#' @export
print.confint_dgaps <- function(x, ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  print(x$cis, ...)
  return(invisible(x))
}
