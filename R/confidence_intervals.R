#' Confidence intervals for the extremal index \eqn{\theta}
#'
#' \code{confint} method for objects of class \code{"exdex"}.
#' Computes confidence intervals for \eqn{theta} based on an object returned
#' from \code{\link{spm}}.
#'
#' @param object An object of class \code{"exdex"}, returned by
#'   \code{\link{spm}}.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any confidence limits that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  Otherwise,
#'   limits that are greater than 1 may be obtained.
#'   If \code{constrain = TRUE} then any lower confidence limits that are
#'   less than 0 are set to 0.
#' @param conf_scale A character scalar.   If \code{conf_scale = "theta"}
#'   then confidence intervals are estimated for \eqn{\theta} directly.
#'   If \code{conf_scale = "log_theta"} then confidence intervals are first
#'   estimated for \eqn{log\theta} and then transformed back to the
#'   \eqn{\theta}-scale.
#' @param ... Further arguments.  None are used at present.
#' @details Add details.
#' @return A matrix with columns giving lower and upper confidence limits for
#'   each parameter. These will be labelled as (1 - level)/2 and
#'   1 - (1 - level)/2 in \% (by default 2.5\% and 97.5\%).
#'   The row names are the names of the model parameters,
#'   if these are available.
#' @examples
#' res <- spm(newlyn, 20)
#' confint(res)
#' @export
confint.exdex <- function (object, level = 0.95, constrain = TRUE,
                           conf_scale = c("theta", "log_theta"), ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  conf_scale <- match.arg(conf_scale)
  # Symmetric confidence intervals, based on large sample normal theory
  # The intervals are (initially) centred on the unconstrained estimate of
  # theta, which may be greater than 1
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  if (conf_scale == "theta") {
    lower <- res$unconstrained_theta - z_val * res$se
    upper <- res$unconstrained_theta + z_val * res$se
  } else {
    lower <- exp(log(res$unconstrained_theta) - z_val * res$se / res$theta)
    upper <- exp(log(res$unconstrained_theta) + z_val * res$se / res$theta)
  }
  # Constrain to (0, 1] if required
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
  return(temp)
}
