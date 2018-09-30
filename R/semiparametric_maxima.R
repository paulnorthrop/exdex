#' Semiparametric maxima estimator of the extremal index
#'
#' Add description
#'
#' @param data A numeric vector of raw data.
#' @param b A numeric scalar.  The block size.
#' @param sliding A logical scalar.
#' @param F_adjust A logical scalar.
#' @param constrain A logical scalar.
#' @param N_or_BB A character scalar.
#' @details Explain re treatment of missings
#' @return Add details.
#' @seealso \code{\link{disjoint_maxima}} for the calculation of the maxima
#'   over disjoint blocks.
#' @seealso \code{\link{sliding_maxima}} for the calculation of the maxima
#'   over sliding blocks.
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index \emph{Extremes}, \strong{18}(4),
#'   585-603. \url{http://dx.doi.org/10.1007/s10687-015-0221-5}
#' @examples
#' spmax_mle(newlyn, 20)
#' spmax_mle(newlyn, 20, sliding = FALSE)
#' @export
spmax_mle <- function(data, b, sliding = TRUE, F_adjust = TRUE,
                      constrain = TRUE, N_or_BB = "N"){
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
  if (F_adjust) {
    m <- length(x)
    Fhaty <- (m * Fhaty - b) / (m - b)
  }
  # Calculate the estimate of theta
  if (N_or_BB == "N") {
    log_v <- b * log(Fhaty)
    theta_mle <- -1 / mean(log_v)
  } else {
    theta_mle <- 1 / (b * mean(1 - Fhaty))
  }
  # Constrain to (0, 1] if required
  if (constrain) {
    theta_mle <- pmin(theta_mle, 1)
  }
  return(theta_mle)
}
