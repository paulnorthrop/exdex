#' Semiparametric maxima estimator of the extremal index
#'
#' Calculates the semiparametric maxima estimator of the extremal index
#' \eqn{\theta} based on sliding or disjoint block maxima based on Northrop (2015).
#'
#' @param data A numeric vector of raw data.
#' @param b A numeric scalar.  The block size.
#' @param sliding A logical scalar indicating whether use sliding blocks
#' (\code{TRUE}) or disjoint blocks (\code{FALSE}).
#' @param constrain A logical scalar indicating whether or not to constrain the
#' mle to lie in the interval (0, 1].
#' @param conf  A numeric scalar.  If \code{conf} is supplied then the
#'   the standard error and \code{conf}\% likelihood-based confidence interval
#'   for \eqn{\theta}, and the standard error of the estimator of \eqn{\theta},
#'   are estimated using block bootstrapping, implemented by
#'   \code{\link[boot]{tsboot}}.  See Northrop (2015) and references therein
#'   for more information.
#' @param R An integer scalar.  The number of bootstrap resamples used to
#'   estimate the standard error and confidence intervals.
#'   See \code{\link[boot]{tsboot}}.
#' @details The extremal index \eqn{\theta} is estimated using the semiparametric
#' maxima estimator of Northrop (2015).  If \code{sliding = TRUE} then the
#' function uses sliding block maxima, that is, the largest value observed in
#' \emph{all} blocks of \code{b} observations, whereas if \code{sliding = FALSE}
#' then disjoint block maxima, that is, the largest values in non-overlapping
#' blocks of \code{b} observations, are used.  If \code{constrain = TRUE} then
#' if the raw estimate of the extremal index is greater than one then a value of
#' 1 is returned. Otherwise (\code{constrain = FALSE}) the raw estimate is
#' returned, even if it is greater than 1.
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta_mle} : } {The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item {\code{theta_se} : } {The estimated standard error of the MLE.}
#'     \item {\code{theta_ci} : } {(If \code{conf} is supplied) a numeric
#'       vector of length two giving lower and upper confidence limits for
#'       \eqn{\theta}.}
#'   }
#'   If \code{conf} is not supplied then only the MLE \code{theta_mle}
#'   is returned.
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @references Northrop, P. J. (2015) \emph{An efficient semiparametric maxima
#' estimator of the extremal index} Extremes, \strong{18}(4), 585-603.
#' \url{http://dx.doi.org/10.1007/s10687-015-0221-5}
#' @examples
#' spm_mle(newlyn, 20)
#' spm_mle(newlyn, 20, sliding = FALSE)
#'
#' spm_mle(newlyn, 20, sliding = FALSE, conf = 95)
#' # When sliding = TRUE bootstrapping is slow
#' \dontrun{
#' spm_mle(newlyn, 20, conf = 95)
#' }
#' @export
spm_mle <- function(data, b, sliding = TRUE, constrain = TRUE, conf = NULL,
                    R = 1000){
  #
  # Function whose returned value depends on for_boot
  #   FALSE: estimate of theta (constrained to (0, 1] if constrain = TRUE))
  #   TRUE: estimate of log(theta) (not constrained to (0, 1])
  spm_estimates <- function(data, for_boot = FALSE) {
    # Calculate the required block maxima
    if (sliding) {
      temp <- sliding_maxima(data, b)
    } else {
      temp <- disjoint_maxima(data, b)
    }
    # Extract x ~ F (only xs contributing to y are included) and y ~ G
    x <- temp$x
    y <- temp$y
    # Empirical c.d.f. of raw (`daily') values
    Fhat <- stats::ecdf(x)
    # Evaluate Fx at y
    Fhaty <- Fhat(y)
    # `Bias-adjust' the empirical c.d.f. of based on the Xs: by subtracting b in
    # numerator and denominator we remove Xs that are in the same block as Y
    # We use of m-b in the denominator rather than the m-b+1 in Northrop (2015)
    m <- length(x)
    Fhaty <- (m * Fhaty - b) / (m - b)
    # Calculate the estimate of theta
    log_v <- b * log(Fhaty)
    theta_mle <- -1 / mean(log_v)
    # For the bootstrap return only the log of the (unconstrained) estimate
    if (for_boot) {
      return(log(theta_mle))
    }
    # Constrain to (0, 1] if required
    if (constrain) {
      theta_mle <- pmin(theta_mle, 1)
    }
    return(theta_mle)
  }
  # Find the point estimate of theta and the raw data that contribute to it
  theta_mle <- spm_estimates(data = data, for_boot = FALSE)
  # If conf is NULL don't do any boostrapping.  Only return the MLE.
  if (is.null(conf)) {
    return(list(theta_mle = theta_mle))
  }
  # Find the optimal value of the mean block length to use in boot::tsboot().
  l_opt <- np::b.star(data)[1]
  # Do the block bootstrap
  theta_boot <- boot::tsboot(tseries = data, statistic = spm_estimates, R = R,
                             l = l_opt, sim = "geom", for_boot = TRUE)
  theta_se <- stats::sd(exp(theta_boot$t))
  # Calculate bootstrap CIs for log(theta)
  ci_log_theta <- boot::boot.ci(theta_boot, conf = conf / 100, type = "basic",
                                index = 1)
  boot_low <- ci_log_theta$basic[,4]
  boot_up <- ci_log_theta$basic[,5]
  conf_int <- c(pmin(exp(boot_low), 1), pmin(exp(boot_up), 1))
  return(list(theta_mle = theta_mle, theta_se = theta_se, theta_ci = conf_int))
}
