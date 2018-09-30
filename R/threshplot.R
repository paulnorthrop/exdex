# ================================= threshplot =================================
#
#' Plot of Theta against Quantile Level of threshold
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on the K-gaps model for threshold inter-exceedances times of
#' Suveges and Davison (2010).
#'
#'
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param tmin A numeric scalar.  Minimum quantile level of threshold applied
#'   to data.
#' @param tmax A numeric scalar.  Maximum quantile level of threshold applied
#'   to data.
#' @param conf A numeric scalar in (0,100). Confidence level required for the
#'   confidence intervals.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param ... Arguments to be passed to methods, such as
#' @details The maximum likelihood estimate of the extremal index \eqn{\theta}
#'   under the K-gaps model of Suveges and Davison (2010) is calculated.
#'   The form of the log-likelihood is given in the \strong{Details} section of
#'    \code{\link{kgaps_stats}}. The plot produced presents the movement of the
#'    extremal index \eqn{\theta} along the quantile level of threshold
#'    indicated by the user. Confidence intervals are also presented on the
#'    graph. The confidence level can be chosen by the user, using the argument
#'    \code{conf}.
#' @references Coles, S.G. An Introduction to Statistical Modeling of Extreme
#' Values (2001)
#' \url{http://dx.doi.org/10.1007/978-1-4471-3675-0}
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta_mle} : } {The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item {\code{theta_se} : } {The estimated standard error of the MLE.}
#'     \item {\code{theta_ci} : } {(If \code{conf} is supplied) a numeric
#'       vector of length two giving lower and upper confidence limits for
#'       \eqn{\theta}.}
#'     \item {\code{ss} : } {The list of summary statistics returned from
#'       \code{\link{kgaps_stats}}.}
#'   }
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation for the
#'   K-gaps model.
#' @examples
#' threshplot(newlyn, 40, 80)
#' @export
threshplot <- function (data, tmin, tmax, conf = 95, k = 1, ... ){

  quant <- seq ( tmin , tmax , by = 1)

  x <- (tmax - tmin) + 1

  p <- quant/100

  thresh <- quantile( data, probs = p )

  mplot<- matrix ( 0 ,nrow = x, ncol = 2, byrow = TRUE, dimnames = list( NULL ,
                   c("Theta", "Quantile")))
  ci <- matrix ( 0 ,nrow = x, ncol = 2, byrow = TRUE, dimnames = list( NULL ,
                c("Lower", "Upper")))

  for( i in 1:x){
    theta <- kgaps_mle( data, thresh[i], conf = conf, k = k)
    mplot[ i, 2 ] <- p[i]
    mplot[ i, 1 ] <- theta$theta_mle
    ci[ i, 1 ] <- theta$theta_ci[1]
    ci[ i, 2 ] <- theta$theta_ci[2]
  }

  y <- cbind(mplot[, "Theta"], ci)
  user_args <- list(...)
  if (is.null(user_args$lty)) {
    user_args$lty <- c(1, 2, 2)
  }
  if (is.null(user_args$col)) {
    user_args$col <- 1
  }
  if (is.null(user_args$xlab)) {
    user_args$xlab <- "Quantile level of threshold"
  }
  if (is.null(user_args$ylab)) {
    user_args$ylab <- "Theta"
  }
  if (is.null(user_args$type)) {
    user_args$type = "l"
  }
  if (is.null(user_args$pch)) {
    user_args$pch = 16
  }
  user_args$x <- mplot[, "Quantile"]
  user_args$y <- y
  do.call(graphics::matplot, user_args)
  return_matrix <- cbind(mplot[, 2:1], ci)
  return(invisible(return_matrix))

}
