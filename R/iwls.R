# =================================== iwls ====================================
#
#' Iterated weighted least squares estimation of the extremal index
#'
#' Estimates the extremal index \eqn{\theta} using the iterated weighted least
#' squares method of Suveges (2007)
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param conf  A numeric scalar.  If \code{conf} is supplied then a
#'   \code{conf}\% confidence interval for \eqn{\theta} is estimated using
#'   bootstrapping.
#' @param maxit A numeric scalar.  The maximum number of iterations.
#' @details The iterated weighted least squares algorithm on page 46 of
#'   Suveges (2007) is used to estimate the value of the extremal index.
#'   This approach uses the time \emph{gaps} between successive exceedances
#'   in the data \code{data} of the threshold \code{thresh}.  The \eqn{i}th
#'   gap is defined as \eqn{T_i - 1}, where \eqn{T_i} is the difference in
#'   the occurrence time of exceedance \eqn{i} and exceedance \eqn{i + 1}.
#'   Therefore, threshold exceedances at adjacent time points produce a gap
#'   of zero.
#'
#'   The model underlying this approach is an exponential-point mas mixture
#'   for \emph{scaled gaps}, that is, gaps multiplied by the proportion of
#'   values in  \code{data} that exceed \code{thresh}.  Under this model
#'   scaled gaps are zero (`within-cluster' interexceedance times) with
#'   probability \eqn{1 - \theta} and otherwise (`between-cluster'
#'   interexceedance times) follow an exponential distribution with mean
#'   \eqn{1 / \theta}.
#'   The estimation method is based on fitting the `broken stick' model of
#'   Ferro (2003) to an exponential quantile-quantile plot of all of the
#'   scaled gaps.  Specifically, the broken stick is a horizontal line
#'   and a line with gradient \eqn{1 / \theta} which intersect at
#'   \eqn{(-\log\theta, 0)}{(-log \theta, 0)}.  The algorithm on page 46 of
#'   Suveges (2007) uses a weighted least squares minimization applied to
#'   the exponential
#'   part of this model to seek a compromise between the role of \eqn{\theta}
#'   as the proportion of interexceedance times that are between-cluster
#'   and the reciprocal of the mean of an exponential distribution for these
#'   interexceedance times.  The weights (see Ferro (2003)) are based on the
#'   variances of order statistics of a standard exponential sample: larger
#'   order statistics have larger sampling variabilities and therefore
#'   receive smaller weight than smaller order statistics.
#'
#'   Note that in step (1) of the algorithm on page 46 of Suveges there is a
#'   typo: \eqn{N_c + 1} should be \eqn{N}, where \eqn{N} is the number of
#'   threshold exceedances.  Also, the gaps are scaled as detailed above,
#'   not by their mean.
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta} : } {The estimate of \eqn{\theta}.}
#'     \item {\code{se} : } {(If \code{conf} is supplied) the estimated
#'       standard error of the estimate.}
#'     \item {\code{theta_ci} : } {(If \code{conf} is supplied) a numeric
#'       vector of length two giving lower and upper confidence limits for
#'       \eqn{\theta}.}
#'     \item {\code{conv} : } {A convergence indicator: 0 indicates successful
#'       convergence; 1 indicates that \code{maxit} has been reached.}
#'     \item {\code{niter} : } {The number of iterations performed.}
#'   }
#' @references Suveges, M. (2007) Likelihood estimation of the extremal
#'   index. \emph{Extremes}, \strong{10}, 41-55.
#'   \url{https://doi.org/10.1007/s10687-007-0034-2}
#' @references Ferro, C.A.T. (2003) Statistical methods for clusters of
#'   extreme values. Ph.D. thesis, Lancaster University.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @examples
#' thresh <- quantile(newlyn, probs = 0.90)
#' iwls(newlyn, thresh)
#' @export
iwls <- function(data, thresh, conf = NULL, maxit = 100) {
  if (!is.numeric(thresh) || length(thresh) != 1) {
    stop("thresh must be a numeric scalar")
  }
  if (thresh >= max(data)) {
    stop("thresh must be less than max(data)")
  }
  # Calculate the quantities required to call iwls_fun()
  #
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > thresh]
  N <- length(exc_u)
  # Inter-exceedances times, (largest first) 1-gaps, number of non-zero 1-gaps
  T_u <- diff(exc_u)
  S_1 <- pmax(T_u - 1, 0)
  # Initial value of n_wls (the number of non-zero 1-gaps)
  n_wls <- length(S_1 > 0)
  # Sort the 1-gaps (largest to smallest) and scale by the sample proportion
  # of values that exceed u
  # [Bottom of page 45 of Suveges (2007), but not mentioned in the algorithm]
  S_1_sort <- sort(S_1, decreasing = TRUE)
  qhat <- N / nx
  S_1_sort <- S_1_sort * qhat
  # Standard exponential quantiles (based on ALL the N-1 1-gaps)
  exp_qs <- -log(1:(N - 1) / N)
  # Weights for the LS fit.  Larger value have large sampling variability and
  # therefore have smaller weights
  ws <-  rev(1 / cumsum(1 / (N:1) ^ 2))
  old_n_wls <- n_wls
  diff_n_wls <- 1
  niter <- 1L
  while (diff_n_wls != 0 & niter < maxit) {
    temp <- iwls_fun(n_wls, N, S_1_sort, exp_qs, ws, nx)
    n_wls <- temp$n_wls
    diff_n_wls <- n_wls - old_n_wls
    old_n_wls <- n_wls
    niter <- niter + 1L
  }
  conv <- ifelse(diff_n_wls > 0, 1, 0)
  n_wls <- temp$n_wls
  theta <- temp$theta
  theta_se <- NULL
  conf_int <- NULL
  return(list(theta = theta, se = theta_se, theta_ci = conf_int,
              conv = conv, niter = niter))
}

# ================================== iwls_fun =================================

iwls_fun <- function(n_wls, N, S_1_sort, exp_qs, ws, nx) {
  #
  # This function implements the algorithm on page 46 of Suveges (2007).
  # [In step (1) there is a typo in the paper: in x_i the N_C+1 should be N.]
  #
  # Args:
  # n_wls    : A numeric scalar.  The number of the largest 1-gaps to include
  #            in the weighted least squares estimation.
  # N        : A numeric scalar.  The number of threshold excesses.
  # S_1_sort : A numeric N-vector.  Sorted (largest to smallest) scaled 1-gaps.
  #            The scaling multiplies the raw 1-gaps by the sample proportion
  #            of values that exceed u.
  # exp_qs   : A numeric N-vector.  Standard exponential quantiles (order
  #            statistics) for a sample of size N.
  # ws       : A numeric N-vector.  Weights for the least squares fit.
  # nx       : A numeric scalar.  The number of raw observations.
  #
  # Returns: A list with components
  #    theta : A numeric scalar.  The new estimate of theta.
  #    n_wls : A numeric scalar.  The new value of n_wls.
  #
  # Extract the values corresponding to the largest n_wls 1-gaps
  # Extract the largest n_wls scaled 1-gaps (ordered largest to smallest)
  chi_i <- S_1_sort[1:n_wls]
  # Standard exponential quantiles, based on N 1-gaps (largest to smallest)
  x_i <- exp_qs[1:n_wls]
  # Extract the weights for the values in chi_i
  ws <- ws[1:n_wls]
  # Weighted least squares for (chi_i, x_i)
  temp <- stats::lm(chi_i ~ x_i, weights = ws)
  ab <- temp$coefficients
  # Estimate theta
  theta <- min(exp(ab[1] / ab[2]), 1)
  # Update n_wls
  n_wls <- floor(theta * (N - 1))
  return(list(theta = theta, n_wls = n_wls))
}

