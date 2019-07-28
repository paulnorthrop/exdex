# ================================= choose_uk =================================

#' Threshold and runs parameter diagnostic for the K-gaps estimator
#'
#' Creates data for a plot to aid the choice of the threshold and
#' run parameter \eqn{K} for the K-gaps estimator (see \code{\link{kgaps}}).
#' \code{\link{plot.choose_uk}} creates the plot.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh,k Numeric vectors.  \code{thresh} is a vector of
#'   extreme value thresholds applied to data.  \code{k} is a vector of values
#'   of the run parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps}} for more details.
#'   The information matrix test is performed a over grid of all
#'   combinations of threshold and \eqn{K} in the vectors \code{thresh}
#'   and \code{k}.
#' @details For each combination of threshold in \code{thresh} and \eqn{K}
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
#'   \item{theta }{a \code{length(thresh)} by \code{length(k)} matrix.
#'     Element (i,j) of \code{theta} contains an object (a list) of class
#'     \code{c("kgaps", "exdex")}, a result of a call
#'     \code{kgaps(data, thresh[j], k[i])} to \code{\link{kgaps}}.}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @seealso \code{\link{kgaps_imt}} for the information matrix test under the
#'   K-gaps model
#' @examples
#' thresh <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' res <- choose_uk(newlyn, thresh = thresh, k = 1:5)
#' @export
choose_uk <- function(data, thresh, k = 1) {
  n_thresh <- length(thresh)
  n_k <- length(k)
  theta <- matrix(rep(list(), n_thresh * n_k), n_k, n_thresh)
  comp <- function(i, j) {
    return((i - 1) * n_thresh + j)
  }
  for (i in 1:n_k) {
    for (j in 1:n_thresh) {
      theta[[comp(i, j)]] <- kgaps(data, thresh[j], k[i])
    }
  }
  imt <- kgaps_imt(data, thresh, k)
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
#' @param y Not used.
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}} and/or \code{\link[graphics]{axis}}.
#' @details Add details
#' @return Nothing is returned.
#' @seealso \code{\link{choose_uk}}.
#' @section Examples:
#' See the examples in \code{\link{choose_uk}}.
#' @export
plot.choose_uk <- function(x, y, ...) {
  return()
}
