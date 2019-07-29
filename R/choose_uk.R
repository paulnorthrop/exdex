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
#'   in \code{k} the functions \code{\link{kgaps}} and \code{link{kgaps_imt}}
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
#' @examples
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_uk(newlyn, u = u, k = 1:5)
#' plot(imt_theta)
#' @export
choose_uk <- function(data, u, k = 1) {
  n_u <- length(u)
  n_k <- length(k)
  theta <- matrix(rep(list(), n_u * n_k), n_k, n_u)
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
#' @param interval_type A character scalar.  The type of confidence interval
#'   to be plotted, if \code{y = "theta"}.  See \code{\link{confint.spm}}.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#'   Only relevant if \code{y = "theta"}.
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}} and/or \code{\link[graphics]{axis}}.
#' @details Add details
#' @return Nothing is returned.
#' @seealso \code{\link{choose_uk}}.
#' @section Examples:
#' See the examples in \code{\link{choose_uk}}.
#' @export
plot.choose_uk <- function(x, y = c("imts", "theta"), interval_type, level,
                           ...) {
  y <- match.arg(y)
  # Extract the values of k and u
  k <- x$imt$k
  u <- x$imt$u
  n_k <- length(k)
  n_u <- length(u)
  if (n_k == 1 && n_u == 1) {
    stop("Object contains only 1 threshold and one value of K")
  }
  def_par <- graphics::par(no.readonly = TRUE)
  # 1. k is a scalar and u is a vector
  if (n_k == 1 && n_u > 1) {
     if (y == "theta") {
       graphics::matplot()
     } else {
       graphics::matplot(x$imt$u, x$imt$imt, type = "l")
     }
  }
  # 2. k is a scalar and u is a vector
  # 3. Both k and u are vectors
  graphics::par(def_par)
  return()
}
