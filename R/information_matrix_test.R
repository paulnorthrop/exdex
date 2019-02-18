# ================================= kgaps_imt =================================
#
#' Information matrix test under the K-gaps model
#'
#' Performs the information matrix test of Suveges and Davison (2010) to
#' diagnose misspecification of the K-gaps model
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh,k Numeric vectors.  \code{thresh} is a vector of
#'   extreme value thresholds applied to data.  \code{k} is a vector of values
#'   of the run parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps_mle}} for more details.
#'   The information matrix test is performed a over grid of all
#'   combinations of threshold and \eqn{K} in the vectors \code{thresh}
#'   and \code{k}.
#' @details Add details
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{http://dx.doi.org/10.1214/09-AOAS292"}
#' @return A list containing
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @examples
#' thresh <- quantile(newlyn, probs = 0.90)
#' pjn <- kgaps_imt(newlyn, thresh)
#' @export
kgaps_imt <- function(data, thresh, k = 1) {
  # Function to return only the MLE of theta
  mle_only <- function(k, data, thresh) {
    return(kgaps_mle(data, thresh, k, inc_cens = FALSE)$theta_mle)
  }
  theta <- T_mat <- p_mat <- NULL
  n_u <- length(thresh)
  n_k <- length(k)
  # Beginning of loop over all thresholds ----------
  for (iu in 1:n_u) {
    u <- thresh[iu]
    # Calculate the MLE of theta for each value of k
    thetahat <- vapply(k, mle_only, 0, data = data, thresh = u)
    # sample size of x
    nx <- length(data)
    # positions of exceedances of u
    exc_u <- (1:nx)[data > u]
    # number of exceedances
    n_u <- length(exc_u)
    # proportion of values that exceed u
    q_u <- n_u / nx
#    q_u <- (n_u - 1) / nx # mev
    # inter-exceedance times
    T_u <- diff(exc_u)
    #
    # Create a length(T) by length(k) matrix with column i containing
    # the values of S_k for k=k[i] ...
    #
    # One column for each value of k
    S_k <- sapply(k, function(k) pmax(T_u - k, 0))
    c_mat <- q_u * S_k
    theta_mat <- matrix(thetahat, ncol = n_k, nrow = nrow(S_k), byrow = TRUE)
    ld <- ifelse(c_mat == 0, -1 / (1 - theta_mat), 2 / theta_mat) - c_mat
    neg_ldd <- ifelse(c_mat == 0, 1 / (1 - theta_mat) ^ 2, 2 / theta_mat ^ 2)
    In <- colMeans(neg_ldd)
    Jn <- colMeans(ld ^ 2)
    Dn <- Jn - In
    dc <- ld ^ 2 - neg_ldd
    dcd <- 4 * c_mat / theta_mat ^ 2 + ifelse(c_mat == 0, 0, -4 / theta_mat ^ 3)
    Dnd <- colMeans(dcd)
    # Multiply the columns of ld by the corresponding elements of Dnd / In
    temp <- ld * rep(Dnd / In, rep(nrow(ld), ncol(ld)))
    Vn <- colMeans((dc - temp) ^ 2)
    test_stats <- (n_u - 1) * Dn ^ 2 / Vn
#    test_stats <- n_u * Dn ^ 2 / Vn # mev
    pvals <- stats::pchisq(test_stats, df = 1, lower.tail = FALSE)
    theta <- rbind(theta, thetahat)
    T_mat <- rbind(T_mat, test_stats)
    p_mat <- rbind(p_mat, pvals)
  }
  # End of loop over thresholds ----------
  colnames(T_mat) <- colnames(p_mat) <- colnames(theta) <- k
  u_ps <- as.numeric(substr(names(thresh), 1,
                            nchar(names(thresh), type = "c") - 1))
  rownames(T_mat) <- rownames(p_mat) <- rownames(theta) <- u_ps
  return(list(IMT = T_mat, p = p_mat, theta = theta))
}


