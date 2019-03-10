#' Confidence intervals for the extremal index \eqn{\theta}
#'
#' \code{confint} method for objects of class \code{"exdex"}.
#' Computes confidence intervals for \eqn{\theta} based on an object returned
#' from \code{\link{spm}}.  Two types of interval are returned:
#' (a) intervals based on approximate large-sample normality of the estimators
#' of \eqn{\theta}, which are symmetric about the respective point estimates,
#' and (b) likelihood-based intervals based on an adjustment of a naive
#' (pseudo-) loglikelihood, using the \code{\link[chandwich]{adjust_loglik}}
#' function in the \code{\link[chandwich]{chandwich}} package.
#'
#' @param object An object of class \code{"exdex"}, returned by
#'   \code{\link{spm}}.
#' @param parm A character scalar specifying whether to estimate
#'   confidence intervals based on sliding maxima or disjoint maxima.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any confidence limits that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  Otherwise,
#'   limits that are greater than 1 may be obtained.
#'   If \code{constrain = TRUE} then any lower confidence limits that are
#'   less than 0 are set to 0.
#' @param conf_scale A character scalar.  Determines the scale on which
#'   we use approximate large-sample normality of the estimators to
#'   estimate confidence intervals.
#'
#'   If \code{conf_scale = "theta"}
#'   then confidence intervals are estimated for \eqn{\theta} directly.
#'   If \code{conf_scale = "log_theta"} then confidence intervals are first
#'   estimated for \eqn{\log\theta}{log\theta} and then transformed back
#'   to the \eqn{\theta}-scale.
#'
#'   Any bias-adjustment requested in the original call to \code{\link{spm}},
#'   using it's \code{bias_adjust} argument, is automatically applied here.
#' @param bias_adjust A logical scalar.  If \code{bias_adjust = TRUE} then,
#'   if appropriate, bias-adjustment is also applied to the loglikelihood
#'   before it is adjusted using \code{\link[chandwich]{adjust_loglik}}.
#'   This is performed only if, in the call to
#'   \code{\link{spm}}, \code{bias_adjust = "BB3"} or
#'   \code{"BB1"} was specified, that is, we have
#'   \code{object$bias_adjust = "BB3"}
#'   or \code{"BB1"}.  In these cases the relevant
#'   component of \code{object$bias_val} is used to scale \eqn{\theta} so
#'   that the location of the maximum of the loglikelihood lies at the
#'   bias-adjusted estimate of \eqn{\theta}.
#'
#'   If \code{bias_adjust = FALSE} or \code{object$bias_adjust = "none"}
#'   or \code{"N"} then no bias-adjustment of the
#'   intervals is performed.  In the latter case this is because the
#'   bias-adjustment is applied in the creation of the data in
#'   \code{object$N2015_data} and \code{object$BB2018_data}, on which the
#'   naive likelihood is based.
#' @param type A character scalar.  The argument \code{type} to be passed to
#'   \code{\link[chandwich]{conf_intervals}} in order to estimate the
#'   likelihood-based intervals.  See \strong{Details}.
#'   Using \code{type = "none"} is \emph{not} advised because then the
#'   intervals are based on naive estimated standard errors.  In particular,
#'   if (the default) \code{sliding = TRUE} was used in the call to
#'   \code{\link{spm}} then the likelihood-based confidence intervals provide
#'   \emph{vast} underestimates of uncertainty.
#' @param plot A logical scalar.  If \code{plot = TRUE} then a plot of the
#'   adjust loglikelihood is produced, using
#'   \code{\link[chandwich]{plot.confint}}.
#' @param ndec An integer scalar.  The legend (if included on the plot)
#'   contains the confidence limits rounded to \code{ndec} decimal places.
#' @param ... Further arguments to be passed to
#'   \code{\link[chandwich]{plot.confint}}, used only if
#'   \code{plot = TRUE}.
#' @details The likelihood-based intervals are estimated using the
#'   \code{\link[chandwich]{adjust_loglik}} function in the
#'   \code{\link[chandwich]{chandwich}} package, followed by a call to
#'   \code{\link[chandwich]{conf_intervals}}.
#'
#'   If \code{object$se} contains \code{NA}s, because the block size \code{b}
#'   was too small or too large in the call to \code{\link{spm}} then
#'   confidence intervals cannot be estimated and an error will be thrown.
#'   See the \strong{Details} section of the \code{\link{spm}} documentation
#'   for more information.
#' @return A matrix with columns giving the lower and upper confidence limits.
#'   These will be labelled as (1 - level)/2 and 1 - (1 - level)/2 in \%
#'   (by default 2.5\% and 97.5\%).
#'   The row names are a concatentation of the variant of the estimator
#'   (N2015 for Northrop (2015), BB2018 for Berghaus and Bucher (2018))
#'   and the type of interval
#'   (sym for symmetric and lik for likelihood-based).
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @examples
#' res <- spm(newlyn, 20)
#' # I can't include these examples until new chandwich is on CRAN
#' #confint(res)
#' #confint(res, plot = TRUE)
#' @export
confint.exdex <- function (object, parm = c("sliding", "disjoint"),
                           level = 0.95, constrain = TRUE,
                           conf_scale = c("theta", "log_theta"),
                           bias_adjust = TRUE,
                           type = c("vertical", "cholesky", "spectral",
                                    "none"),
                           plot = FALSE, ndec = 2, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  if (is.na(object$se_sl[1])) {
    temp <- matrix(NA, nrow = 4, ncol = 2)
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    colnames(temp) <- pct
    rownames(temp) <- c("N2015sym", "BB2018sym", "N2015lik", "BB0218lik")
    return(temp)
  }
  parm <- match.arg(parm)
  if (parm == "sliding" && type == "none") {
    warning("The likelihood-based CIs are vast underestimates of uncertainty!")
  }
  conf_scale <- match.arg(conf_scale)
  type <- match.arg(type)
  # Set the components that we need, based on argument maxima
  if (parm == "sliding") {
    uncon <- object$uncon_theta_sl
    se <- object$se_sl
    theta <- object$theta_sl
    yz_data <- object$data_sl
    bias_val <- object$bias_sl
  } else {
    uncon <- object$uncon_theta_dj
    se <- object$se_dj
    theta <- object$theta_dj
    yz_data <- object$data_dj
    bias_val <- object$bias_dj
  }
  # Symmetric confidence intervals, based on large sample normal theory
  # The intervals are (initially) centred on the unconstrained estimate of
  # theta, which may be greater than 1
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  if (conf_scale == "theta") {
    lower <- uncon - z_val * se
    upper <- uncon + z_val * se
  } else {
    lower <- exp(log(uncon) - z_val * se / theta)
    upper <- exp(log(uncon) + z_val * se / theta)
  }
  names(lower) <- paste0(names(lower), "sym")
  names(upper) <- paste0(names(upper), "sym")
  #
  # Likelihood-based confidence intervals.  We use the chandwich package
  # to adjust the independence loglikelihood so that it curvature at its
  # maximum corresponds to the estimated standard errors that result from
  # the theory in Berghaus and Bucher (2018). Note that:
  #
  # 1. we must use raw MLEs, not the bias-adjusted versions, because these
  #    give the location of the maximum of the independence loglikelihood
  #
  # 2. we give the option, via the argument bias_adjust, to bias-adjust
  #    the intervals based on the bias_adjust argument supplied in the
  #    original call to spm(): if bias_adjust = "BB3" or "BB1" in that call
  #    then the bias-adjustment is subtracted from the limits of the interval;
  #    if bias_adjust = "N" in that call then no adjustment in performed
  #    because this bias-adjustment is applied in the creation of the data
  #    in object$N2015_data and object$BB2018_data.
  #
  exponential_loglik <- function(theta, data) {
    if (theta <= 0) return(-Inf)
    return(log(theta) - theta * data)
  }
  # Bias-adjust, if requested and if appropriate
  # We can achieve this by scaling the data by the ratio of the naive MLE
  # to the bias-adjusted MLE
  mleN <- 1 / mean(yz_data[, "N2015"])
  mleBB <- 1 / mean(yz_data[, "BB2018"])
  if (bias_adjust && object$bias_adjust %in% c("BB3", "BB1")) {
    scaleN <- mleN / (mleN - bias_val["N2015"])
    scaleBB <- mleBB / (mleBB - bias_val["BB2018"])
  } else {
    scaleN <- 1
    scaleBB <- 1
  }
  # Northrop (2015)
  n <- length(yz_data[, "N2015"])
  H <- as.matrix(-n / mleN ^ 2)
  V <- as.matrix(H ^ 2 * se["N2015"] ^ 2)
  # Note the multiplication of the data by scaleN and the division of the
  # naive estimate by scaleN, to obtain the bias_adjusted estimate
  adjN <- chandwich::adjust_loglik(loglik = exponential_loglik,
                                   data = yz_data[, "N2015"] * scaleN,
                                   p = 1, par_names = "theta",
                                   mle = mleN / scaleN, H = H, V = V)
  # Berghaus and Bucher (2018)
  n <- length(yz_data[, "BB2018"])
  H <- as.matrix(-n / mleBB ^ 2)
  V <- as.matrix(H ^ 2 * se["BB2018"] ^ 2)
  # Note the multiplication of the data by scaleBB and the division of the
  # naive estimate by scaleBB, to obtain the bias_adjusted estimate
  adjBB <- chandwich::adjust_loglik(loglik = exponential_loglik,
                                    data = yz_data[, "BB2018"] * scaleBB,
                                    p = 1, par_names = "theta",
                                    mle = mleBB / scaleBB, H = H, V = V)
  # Avoid chandwich::conf_intervals()'s profiling messages
  tempN <- suppressMessages(chandwich::conf_intervals(adjN, conf = 100
                                                      * level, type = type,
                                                      lower = 0))
  tempBB <- suppressMessages(chandwich::conf_intervals(adjBB, conf = 100 *
                                                         level, type = type,
                                                       lower = 0))
  # Add the likelihood-based intervals to the symmetric ones
  lower <- c(lower, tempN$prof_CI[1], tempBB$prof_CI[1])
  upper <- c(upper, tempN$prof_CI[2], tempBB$prof_CI[2])
  names(lower)[3:4] <- c("N2015lik", "BB2018lik")
  names(upper)[3:4] <- c("N2015lik", "BB2018lik")
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
  # Produce a plot of the adjusted loglikelihood, if requested
  # Owing to the different scales of the loglikelihoods for N2015 and BB2018 we
  # shoof them to have a maximum of 0, so that we can display them on one plot
  if (plot) {
    shoofN <- max(tempN$prof_loglik_vals)
    shoofBB <- max(tempBB$prof_loglik_vals)
    tempN$prof_loglik_vals <- tempN$prof_loglik_vals - shoofN
    tempBB$prof_loglik_vals <- tempBB$prof_loglik_vals - shoofBB
    tempN$max_loglik <- tempN$max_loglik - shoofN
    tempBB$max_loglik <- tempBB$max_loglik - shoofBB
    # Round confidence limits for inclusion in the legend
    fmt <- paste0("%.", ndec, "f")
    roundN <- sprintf(fmt, round(temp["N2015lik", ], ndec))
    roundBB <- sprintf(fmt, round(temp["BB2018lik", ], ndec))
    my_leg <- NULL
    my_leg[1] <- paste("N2015:    (", roundN[1], ",", roundN[2], ")")
    my_leg[2] <- paste("BB2018: (", roundBB[1], ",", roundBB[2], ")")
    # A clunky way to avoid conflict between my choice of legend and
    # (legend) title and those that the user might supply via ...
    user_args <- list(...)
    if (is.null(user_args$legend_pos)) {
      if (is.null(user_args$legend) && is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, legend = my_leg,
             title = "estimator", legend_pos = "bottom", ...)
      } else if (is.null(user_args$legend) && !is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, legend = my_leg, legend_pos = "bottom",
             ...)
      } else if (!is.null(user_args$legend) && is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, title = "estimator", legend_pos = "bottom",
             ...)
      } else {
        plot(x = tempN, y = tempBB, legend_pos = "bottom", ...)
      }
    } else {
      # Make user_args$legend_pos NULL because otherwise
      # is.null(user_args$legend) returns FALSE
      user_args$legend_pos <- NULL
      if (is.null(user_args$legend) && is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, legend = my_leg,
             title = "estimator", ...)
      } else if (is.null(user_args$legend) && !is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, legend = my_leg, ...)
      } else if (!is.null(user_args$legend) && is.null(user_args$title)) {
        plot(x = tempN, y = tempBB, title = "estimator", ...)
      } else {
        plot(x = tempN, y = tempBB, ...)
      }
    }
#    abline(v = temp["N2015lik", ])
#    abline(v = temp["BB2018lik", ], lty = 2)
  }
  return(temp)
}
