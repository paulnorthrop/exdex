# =================================== dgaps ===================================
#
#' Maximum likelihood estimation using left-censored inter-exceedances times
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on a model for threshold inter-exceedances times of
#' Holesovsky and Fusek (2020).  We refer to this as the \eqn{D}-gaps model,
#' because it uses a tuning parameter \eqn{D}, whereas the related \eqn{K}-gaps
#' model of Suveges and Davison (2010) has a tuning parameter \eqn{K}.
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then \code{\link{split_by_NAs}} is
#'   used to divide the data further into sequences of non-missing values,
#'   stored in different columns in a matrix.  Again, the log-likelihood
#'   is constructed as a sum of contributions from different columns.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param D A numeric scalar.  The censoring parameter \eqn{D}. Threshold
#'   inter-exceedances times that are not larger than \code{D} units are
#'   left-censored, occurring with probability
#'   \eqn{\log(1 - \theta e^{-\theta d})}{log(1 - exp(-\theta d))},
#'   where \eqn{d = q D} and \eqn{q} is the probability with which the
#'   threshold \eqn{u} is exceeded.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations. It is known that these times are greater
#'   than or equal to the time observed.
#' @details If \code{inc_cens = FALSE} then the maximum likelihood estimate of
#'   the extremal index \eqn{\theta} under the \eqn{D}-gaps model of
#'   Holesovsky and Fusek (2020) is calculated.
#'   If \code{inc_cens = TRUE} then information from right-censored
#'   inter-exceedance times is also included in the likelihood to be maximized.
#'   For an explanation of the idea see Attalides (2015).  The form of the
#'   log-likelihood is given in the \strong{Details} section of
#'   \code{\link{dgaps_stat}}.
#'
#'   It is possible that the estimate of \eqn{\theta} is equal to 1, and also
#'   possible that it is equal to 0. \code{\link{dgaps_stat}} explains the
#'   respective properties of the data that cause these events to occur.
#' @return An object (a list) of class \code{c("dgaps", "exdex")} containing
#'     \item{\code{theta} }{The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item{\code{se} }{The estimated standard error of the MLE, calculated
#'       using an algebraic expression for the observed information.}
#'     \item{\code{se_check} }{The estimated standard error of the MLE,
#'       estimated using \code{\link[stats:optim]{optimHess}}.}
#'     \item{\code{ss} }{The list of summary statistics returned from
#'       \code{\link{kgaps_stat}}.}
#'     \item{\code{d, u, inc_cens} }{The input values of \code{d},
#'       \code{u} and \code{inc_cens}.}
#'     \item{\code{call }}{The call to \code{dgaps}.}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{kgaps_stat}} for the calculation of sufficient
#'   statistics for the \eqn{K}-gaps model.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197–213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @examples
#' ### S&P 500 index
#'
#' u <- quantile(sp500, probs = 0.60)
#' theta <- dgaps(sp500, u = u, D = 1)
#' theta
#' summary(theta)
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = 0.60)
#' theta <- dgaps(newlyn, u = u, D = 2)
#' theta
#' summary(theta)
#' @export
dgaps <- function(data, u, D = 1, inc_cens = TRUE) {
  Call <- match.call(expand.dots = TRUE)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data, na.rm = TRUE)) {
    stop("u must be less than max(data)")
  }
  if (!is.numeric(D) || length(D) != 1) {
    stop("d must be a numeric scalar")
  }
  # If there are missing values then use split_by_NAs to extract sequences
  # of non-missing values
  if (anyNA(data) && is.null(attr(data, "split_by_NAs_done"))) {
    data <- split_by_NAs(data)
  }
  # Calculate sufficient statistics for each column in data and then sum
  stats_list <- apply(as.matrix(data), 2, dgaps_stat, u = u, D = D,
                      inc_cens = inc_cens)
  ss <- Reduce(f = function(...) Map("+", ...), stats_list)
  # If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
  # and the likelihood is maximized at theta = 1.
  N0 <- ss$N0
  # If N1 = 0 then we are in the degenerate case where there is one cluster
  # (all K-gaps are zero) and the likelihood is maximized at theta = 0.
  N1 <- ss$N1
  theta_se_check <- NA
  if (N1 == 0) {
    theta_mle <- 0L
  } else if (N0 == 0) {
    theta_mle <- 1L
  } else {
    # Use the K-gaps estimate as an initial estimate
    theta_init <- kgaps(data = data, u = u, k = D, inc_cens = inc_cens)$theta
    dgaps_negloglik <- function(theta) {
      return(-do.call(dgaps_loglik, c(list(theta = theta), ss)))
    }
    temp <- stats::optim(theta_init, dgaps_negloglik, method = "Brent",
                         lower = 0, upper = 1)
    theta_mle <- temp$par
    theta_se_check <- 1 / sqrt(optimHess(theta_mle, dgaps_negloglik))
  }
  # Estimate standard error
  obs_info <- 0
  if (N0 > 0) {
    obs_info <- obs_info + N0 * gdd_theta(theta_mle, q_u = ss$q_u, D = ss$D)
  }
  if (N1 > 0) {
    obs_info <- obs_info + 2 * N1 / theta_mle ^ 2
  }
  theta_se <- sqrt(1 / obs_info)
  res <- list(theta = theta_mle, se = theta_se, se_check = theta_se_check,
              ss = ss, D = D, u = u, inc_cens = inc_cens, call = Call)
  class(res) <- c("dgaps", "exdex")
  return(res)
}

# ================================ dgaps_stat =================================

#' Sufficient statistics for the left-censored inter-exceedances time model
#'
#' Calculates sufficient statistics for the the left-censored inter-exceedances
#' time \eqn{D}-gaps model for the extremal index \eqn{\theta}.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param D A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @details The sample \eqn{K}-gaps are
#'   \eqn{S_0, S_1, ..., S_{N-1}, S_N}{S_0, S_1, ..., S_(N-1), S_N},
#'   where \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)} are uncensored and
#'   \eqn{S_0} and \eqn{S_N} are censored.  Under the assumption that the
#'   \eqn{K}-gaps are independent, the log-likelihood of the \eqn{K}-gaps
#'   model is given by
#'   \deqn{l(\theta; S_0, \ldots, S_N) = N_0 \log(1 - \theta) +
#'     2 N_1 \log \theta - \theta q (S_0 + \cdots + S_N),}{%
#'     l(\theta; S_0, ..., S_N) = N_0 log(1 - \theta) + 2 N_1 log \theta -
#'     \theta q (S_0 + ... + S_N),}
#'    where \eqn{q} is the threshold exceedance probability,
#'    \eqn{N_0} is the number of sample \eqn{K}-gaps that are equal to zero and
#'    (apart from an adjustment for the contributions of \eqn{S_0} and
#'    \eqn{S_N}) \eqn{N_1} is the number of positive sample \eqn{K}-gaps.
#'    Specifically, \eqn{N_1} is equal to the number of
#'    \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)}
#'    that are positive plus \eqn{(I_0 + I_N) / 2}, where \eqn{I_0 = 1} if
#'    \eqn{S_0} is greater than zero and similarly for \eqn{I_N}.
#'    The differing treatment of uncensored and censored \eqn{K}-gaps reflects
#'    differing contributions to the likelihood.
#'    For full details see Suveges and Davison (2010) and Attalides (2015).
#'
#'    If \eqn{N_1 = 0} then we are in the degenerate case where there is one
#'    cluster (all \eqn{K}-gaps are zero) and the likelihood is maximized at
#'    \eqn{\theta = 0}.
#'
#'    If \eqn{N_0 = 0} then all exceedances occur singly (all \eqn{K}-gaps are
#'    positive) and the likelihood is maximized at \eqn{\theta = 1}.
#' @return A list containing the sufficient statistics, with components
#'     \item{\code{N0} }{the number of zero \eqn{K}-gaps.}
#'     \item{\code{N1} }{contribution from non-zero \eqn{K}-gaps (see
#'       \strong{Details}).}
#'     \item{\code{sum_qs} }{the sum of the (scaled) \eqn{K}-gaps, i.e.
#'       \eqn{q (S_0 + \cdots + S_N)}{q (S_0 + ... + S_N)}, where \eqn{q}
#'       is estimated by the proportion of threshold exceedances.}
#'     \item{\code{n_kgaps} }{the number of \eqn{K}-gaps, including 2
#'       censored \eqn{K}-gaps if \code{inc_cens = TRUE}.}
#'     \item{\code{q_u} }{the sample proportion of values that exceed the
#'       threshold.}
#'     \item{\code{d} }{the input value of \code{d}.}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197–213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' dgaps_stat(newlyn, u = u, D = 1)
#' @export
dgaps_stat <- function(data, u, D = 1, inc_cens = TRUE) {
  data <- stats::na.omit(data)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (!is.numeric(D) || length(D) != 1) {
    stop("k must be a numeric scalar")
  }
  # If all the data are smaller than the threshold then return null results
  if (u >= max(data, na.rm = TRUE)) {
    return(list(N0 = 0, N1 = 0, sum_qtd = 0, n_dgaps = 0, q_u = 0, D = D))
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  left_censored <- T_u <= D
  # N0, N1, sum of scaled inter-exceedance times that are greater than d,
  # that is, not left-censored
  N1 <- sum(!left_censored)
  N0 <- N_u - 1 - N1
  T_gt_D <- T_u[!left_censored]
  sum_qtd <- sum(q_u * T_gt_D)
  # Store the number of K-gaps, for use by nobs.kgaps()
  n_dgaps <- N0 + N1
  # Include censored inter-exceedance times?
  if (inc_cens) {
    # censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    # T_u_cens <= d adds no information, because we have no idea to which part
    # of the log-likelihood they would contribute
    left_censored_cens <- T_u_cens <= D
    # N0, N1, sum of scaled inter-exceedance times that are greater than d,
    # that is, not left-censored
    N1_cens <- sum(T_u_cens[!left_censored_cens])
    n_gaps <- n_dgaps + N1_cens
    T_gt_D_cens <- T_u_cens[!left_censored_cens]
    sum_qtd_cens <- sum(q_u * T_gt_D_cens)
    # Add contributions.
    # Note: we divide N1_cens by two because a right-censored inter-exceedance
    # times that is not left-censored at d (i.e. is greater than d) contributes
    # theta exp(-theta q_u T_u) to the D-gaps likelihood, but an uncensored
    # observation contributes theta^2 exp(-theta q_u T_u).
    N1 <- N1 + N1_cens / 2
    sum_qtd <- sum_qtd + sum_qtd_cens
  }
  return(list(N0 = N0, N1 = N1, sum_qtd = sum_qtd, n_dgaps = n_dgaps,
              q_u = q_u, D = D))
}
