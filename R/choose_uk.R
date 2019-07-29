# ================================= choose_uk =================================

#' Threshold \eqn{u} and runs parameter \eqn{K} diagnostic for the \eqn{K}-gaps
#' estimator
#'
#' Creates data for a plot to aid the choice of the threshold and
#' run parameter \eqn{K} for the \eqn{K}-gaps estimator (see
#' \code{\link{kgaps}}).  \code{\link{plot.choose_uk}} creates the plot.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u,k Numeric vectors.  \code{u} is a vector of
#'   extreme value thresholds applied to data.  \code{k} is a vector of values
#'   of the run parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps}} for more details.
#' @details For each combination of threshold in \code{u} and \eqn{K}
#'   in \code{k} the functions \code{\link{kgaps}} and \code{\link{kgaps_imt}}
#'   are called in order to estimate \eqn{\theta} and to perform the
#'   information matrix test of Suveges and Davison (2010).
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292"}
#' @return An object (a list) of class \code{c("choose_uk", "exdex")}
#'   containing
#'   \item{imt }{an object of class \code{c("kgaps_imt", "exdex")} returned
#'     from \code{\link{kgaps_imt}}.}
#'   \item{theta }{a \code{length(u)} by \code{length(k)} matrix.
#'     Element (i,j) of \code{theta} contains an object (a list) of class
#'     \code{c("kgaps", "exdex")}, a result of a call
#'     \code{kgaps(data, u[j], k[i])} to \code{\link{kgaps}}.}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{kgaps_imt}} for the information matrix test under the
#'   \eqn{K}-gaps model
#' @seealso \code{\link{plot.choose_uk}} to produce the diagnostic plot.
#' @examples
#' # One run parameter K, many thresholds u
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_uk(newlyn, u = u, k = 1)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#'
#' # One threshold u, many run parameters K
#' u <- quantile(newlyn, probs = 0.9)
#' imt_theta <- choose_uk(newlyn, u = u, k = 1:5)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#' @export
choose_uk <- function(data, u, k = 1) {
  n_u <- length(u)
  n_k <- length(k)
  theta <- matrix(rep(list(), n_u * n_k), n_k, n_u)
  # Function to set the correct element of the matrix of lists theta
  # i indexes k, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  for (i in 1:n_k) {
    for (j in 1:n_u) {
      theta[[comp(i, j)]] <- kgaps(data, u[j], k[i])
    }
  }
  imt <- kgaps_imt(data, u, k)
  res <- list(imt = imt, theta = theta)
  class(res) <- c("choose_uk", "exdex")
  return(res)
}

# ============================= plot.choose_uk ===============================

#' Plot block length diagnostic for the semiparametric maxima estimator
#'
#' \code{plot} method for objects inheriting from class \code{"choose_uk"},
#' returned from \code{\link{choose_uk}}
#'
#' @param x an object of class \code{c("choose_uk", "exdex")}, a result of a
#'   call to \code{\link{choose_uk}}.
#' @param y A character scalar indicating what should be plotted.
#'  This is only relevant if, in the call to \code{\link{choose_uk}} that
#'  produced \code{x} either only one value threshold was supplied via
#'  \code{u} or only one value of \eqn{K} was supplied via \code{k}.
#'  In that event, information matrix test statistics are plotted if
#'  \code{y = "imts"} and estimates of, and confidence intervals for,
#'  \eqn{\theta} is plotted if \code{y = "theta"}.
#'  Otherwise, information matrix test statistics are plotted.
#' @param level A numeric scalar in (0, 1).  The confidence level used in calculating confidence intervals
#'   for \eqn{\theta}.  Only relevant if \code{y = "theta"}.
#' @param interval_type A character scalar.  The type of confidence interval
#'   to be plotted, if \code{y = "theta"}.  See \code{\link{confint.kgaps}}.
#' @param conf_scale A character scalar.  If \code{interval_type = "norm"} then
#'   \code{conf_scale} determines the scale on which we use approximate
#'   large-sample normality of the estimator to estimate confidence intervals.
#'   See \code{\link{confint.kgaps}}.
#' @param alpha A numeric scalar in (0, 1). The significance level to be used
#'   for the information matrix test.
#' @param constrain A logical scalar.  The argument \code{constrain} to
#'  \code{\link{confint.kgaps}}.
#' @param for_abline Only relevant when \code{y = "imt"} and at one of
#'   \code{u} or \code{k} is scalar. A list of graphical parameters to be
#'   passed to \code{\link{abline}} to indicate the critical value of the
#'   information matrix test implied by \code{alpha}.
#' @param digits An integer. Used for formatting the value of the threshold
#'   with \code{\link[base:Round]{signif}} before adding its value to a plot.
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}}.
#' @details Add details
#' @return Nothing is returned.
#' @seealso \code{\link{choose_uk}}.
#' @section Examples:
#' See the examples in \code{\link{choose_uk}}.
#' @export
plot.choose_uk <- function(x, y = c("imts", "theta"), level = 0.95,
                           interval_type = c("norm", "lik"),
                           conf_scale = c("theta", "log"), alpha = 0.05,
                           constrain = TRUE,
                           for_abline = list(lty = 2, lwd = 1, col = 1),
                           digits = 3, ...) {
  y <- match.arg(y)
  interval_type <- match.arg(interval_type)
  conf_scale <- match.arg(conf_scale)
  # Extract the values of k and u
  k <- x$imt$k
  u <- x$imt$u
  n_k <- length(k)
  n_u <- length(u)
  if (n_k == 1 && n_u == 1) {
    stop("Object contains only 1 threshold and one value of K")
  }
  # Function to set the correct element of the matrix of lists theta
  # i indexes k, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  def_par <- graphics::par(no.readonly = TRUE)
  # My plotting functions: to give defaults but allow the user to override
  my_matplot <- function(x, y, ..., type = "l", lty = 1, col = 1,
                         xlab = my_xlab, ylab = my_ylab, ylim = my_ylim) {
    graphics::matplot(x = x, y = y, ..., type = type, lty = lty, col = col,
                      xlab = xlab, ylab = ylab, ylim = ylim)
  }
  my_title <- function(..., main = my_main) {
    graphics::title(..., main = main)
  }
  # Critical value for the IMT (in case we need it)
  crit <- stats::qchisq(alpha, df = 1, lower.tail = FALSE)
  # One of k or u is scalar
  cond1 <- n_k == 1 && n_u > 1
  cond2 <- n_k > 1 && n_u == 1
  if (cond1 || cond2)  {
     max_uk <- max(n_k, n_u)
     ymat <- matrix(NA, ncol= 3, nrow = max_uk)
     if (cond1) {
       xvec <- u
     } else {
       xvec <- k
     }
     loop_vec <- 1:max_uk
     if (y == "theta") {
       for (ij in loop_vec) {
         if (cond1) {
           kgaps_object <- x$theta[[comp(1, ij)]]
         } else {
           kgaps_object <- x$theta[[comp(ij, 1)]]
         }
         temp <- confint(kgaps_object, level = level,
                         interval_type = interval_type,
                         conf_scale = conf_scale, constrain = constrain)
         ymat[ij, 1] <- kgaps_object$theta
         ymat[ij, 2:3] <- temp
       }
       my_ylab <- "theta"
       my_xlab <- ifelse(cond1, "threshold u", "run parameter K")
       my_ylim <- c(0, 1)
       my_matplot(xvec, ymat, ...)
       my_main <- ifelse(cond1, paste0("run parameter K = ", k),
                         paste0("threshold u = ", signif(u, digits = digits)))
       my_title(...)
     } else {
       my_ylab <- "IMT"
       my_xlab <- ifelse(cond1, "threshold u", "run parameter K")
       my_ylim <- c(0, max(x$imt$imt))
       if (cond1) {
         xvec <- x$imt$u
         ymat <- x$imt$imt
       } else {
         xvec <- x$imt$k
         ymat <- t(x$imt$imt)
       }
       my_matplot(xvec, ymat)
       my_main <- ifelse(cond1, paste0("run parameter K = ", k),
                         paste0("threshold u = ", signif(u, digits = digits)))
       my_title(...)
       for_abline <- c(for_abline, h = crit)
       do.call(graphics::abline, for_abline)
     }
  }
  # Both k and u are vectors
  graphics::par(def_par)
  return(invisible())
}
