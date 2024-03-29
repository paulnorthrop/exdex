# =================================== kgaps ===================================
#
#' Maximum likelihood estimation for the \eqn{K}-gaps model
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on the \eqn{K}-gaps model for threshold inter-exceedances times of
#' Suveges and Davison (2010).
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
#' @param k A non-negative numeric scalar. Run parameter \eqn{K}, as defined in
#'   Suveges and Davison (2010).  Threshold inter-exceedances times that are
#'   not larger than \code{k} units are assigned to the same cluster, resulting
#'   in a \eqn{K}-gap equal to zero. Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.  In practice, \eqn{k} should
#'   be no smaller than 1, because when \eqn{k} is less than 1 the estimate
#'   of \eqn{\theta} is always equal to 1.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations.  It is known that these times are greater
#'   than or equal to the time observed. See Attalides (2015) for details.
#'   If \code{data} has multiple columns then there will be right-censored
#'   first and last inter-exceedance times for each column.
#' @details If \code{inc_cens = FALSE} then the maximum likelihood estimate of
#'   the extremal index \eqn{\theta} under the \eqn{K}-gaps model of
#'   Suveges and Davison (2010) is calculated.
#'
#'   If \code{inc_cens = TRUE} then information from right-censored
#'   first and last inter-exceedance times is also included in the likelihood
#'   to be maximized, following Attalides (2015).  The form of the
#'   log-likelihood is given in the \strong{Details} section of
#'   \code{\link{kgaps_stat}}.
#'
#'   It is possible that the estimate of \eqn{\theta} is equal to 1, and also
#'   possible that it is equal to 0. \code{\link{kgaps_stat}} explains the
#'   respective properties of the data that cause these events to occur.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @return An object (a list) of class \code{c("kgaps", "exdex")} containing
#'     \item{\code{theta} }{The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item{\code{se} }{The estimated standard error of the MLE, calculated
#'       using an algebraic expression for the observed information.
#'       If \code{k = 0} then \code{se} is returned as \code{0}.}
#'     \item{\code{se_exp} }{The estimated standard error of the MLE,
#'       calculated using an algebraic expression for the expected information.
#'       If the estimate of \eqn{\theta} is 0 or 1 then \code{se_exp} is
#'       \code{NA}.}
#'     \item{\code{ss} }{The list of summary statistics returned from
#'       \code{\link{kgaps_stat}}.}
#'     \item{\code{k, u, inc_cens} }{The input values of \code{k},
#'       \code{u} and \code{inc_cens}.}
#'     \item{\code{max_loglik }}{The value of the log-likelihood at the MLE.}
#'     \item{\code{call }}{The call to \code{kgaps}.}
#' @seealso \code{\link{kgaps_confint}} to estimate confidence intervals
#'   for \eqn{\theta}.
#' @seealso \code{\link{kgaps_methods}} for S3 methods for \code{"kgaps"}
#'   objects.
#' @seealso \code{\link{kgaps_imt}} for the information matrix test, which
#'   may be used to inform the choice of the pair (\code{u, k}).
#' @seealso \code{\link{choose_uk}} for a diagnostic plot based on
#'   \code{\link{kgaps_imt}}.
#' @seealso \code{\link{kgaps_stat}} for the calculation of sufficient
#'   statistics for the \eqn{K}-gaps model.
#' @seealso \code{\link[revdbayes]{kgaps_post}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for Bayesian inference
#'   about \eqn{\theta} using the \eqn{K}-gaps model.
#' @examples
#' ### S&P 500 index
#'
#' u <- quantile(sp500, probs = 0.60)
#' theta <- kgaps(sp500, u)
#' theta
#' summary(theta)
#' coef(theta)
#' nobs(theta)
#' vcov(theta)
#' logLik(theta)
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = 0.60)
#' theta <- kgaps(newlyn, u, k = 2)
#' theta
#' summary(theta)
#'
#' ### Cheeseboro wind gusts
#'
#' theta <- kgaps(cheeseboro, 45, k = 3)
#' theta
#' summary(theta)
#' @export
kgaps <- function(data, u, k = 1, inc_cens = TRUE) {
  Call <- match.call(expand.dots = TRUE)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data, na.rm = TRUE)) {
    stop("'u' must be less than 'max(data, na.rm = TRUE)'")
  }
  if (!is.numeric(k) || k < 0 || length(k) != 1) {
    stop("k must be a non-negative scalar")
  }
  # If there are missing values then use split_by_NAs to extract sequences
  # of non-missing values
  if (anyNA(data) && is.null(attr(data, "split_by_NAs_done"))) {
    data <- split_by_NAs(data)
  }
  # Estimate the marginal exceedance probability q_u
  q_u <- mean(data > u, na.rm = TRUE)
  # Calculate sufficient statistics for each column in data and then sum
  stats_list <- apply(as.matrix(data), 2, kgaps_stat, u = u, q_u = q_u, k = k,
                      inc_cens = inc_cens)
  ss <- Reduce(f = function(...) Map("+", ...), stats_list)
  # If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
  # and the likelihood is maximized at theta = 1.
  N0 <- ss$N0
  # If N1 = 0 then we are in the degenerate case where there is one cluster
  # (all K-gaps are zero) and the likelihood is maximized at theta = 0.
  N1 <- ss$N1
  if (N1 == 0) {
    theta_mle <- 0L
  } else if (N0 == 0) {
    theta_mle <- 1L
  } else {
    sum_qs <- ss$sum_qs
    theta_mle <- kgaps_quad_solve(N0, N1, sum_qs)
  }
  # Estimate standard error
  # For completeness add an estimate based on the expected information
  # If N1 = 0 then the estimate of theta is 0 and we return NA for se_exp
  # If N0 = 0 then the estimate of theta is 1 and the expected information is
  # not defined and we return NA for se_exp
  if (N1 > 0 && N0 > 0) {
    exp_info <- kgaps_exp_info(theta = theta_mle, ss = ss, inc_cens = inc_cens)
  } else {
    exp_info <- NA
  }
  se_exp <- 1 / sqrt(exp_info)
  # Based on the observed information
  obs_info <- 0
  if (N0 > 0) {
    obs_info <- obs_info + N0 / (1 - theta_mle) ^ 2
  }
  if (N1 > 0) {
    obs_info <- obs_info + 2 * N1 / theta_mle ^ 2
  }
  theta_se <- sqrt(1 / obs_info)
  # If K = 0 then the estimate of theta is 1 by default
  # We return a SE equal to 0 so that the estimation of SEs for return level
  # estimates in the lite package work in this case
  if (k == 0) {
    theta_se <- 0
  }
  max_loglik <- do.call(kgaps_loglik, c(list(theta = theta_mle), ss))
  res <- list(theta = theta_mle, se = theta_se, se_exp = se_exp, ss = ss,
              k = k, u = u, inc_cens = inc_cens, max_loglik = max_loglik,
              call = Call)
  class(res) <- c("kgaps", "exdex")
  return(res)
}

# ================================ kgaps_stat =================================

#' Sufficient statistics for the \eqn{K}-gaps model
#'
#' Calculates sufficient statistics for the \eqn{K}-gaps model for the extremal
#' index \eqn{\theta}. Called by \code{\link{kgaps}}.
#'
#' @param data A numeric vector of raw data.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param q_u A numeric scalar.  An estimate of the probability with which
#'   the threshold \code{u} is exceeded.  If \code{q_u} is missing then it is
#'   calculated using \code{mean(data > u, na.rm = TRUE)}.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times relating to the
#'   first and last observation.  It is known that these times are greater
#'   than or equal to the time observed. See Attalides (2015) for details.
#' @details The sample \eqn{K}-gaps are
#'   \eqn{S_0, S_1, ..., S_{N-1}, S_N}{S_0, S_1, ..., S_(N-1), S_N},
#'   where \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)} are uncensored and
#'   \eqn{S_0} and \eqn{S_N} are right-censored.  Under the assumption that the
#'   \eqn{K}-gaps are independent, the log-likelihood of the \eqn{K}-gaps
#'   model is given by
#'   \deqn{l(\theta; S_0, \ldots, S_N) = N_0 \log(1 - \theta) +
#'     2 N_1 \log \theta - \theta q (S_0 + \cdots + S_N),}{%
#'     l(\theta; S_0, ..., S_N) = N_0 log(1 - \theta) + 2 N_1 log \theta -
#'     \theta q (S_0 + ... + S_N),}
#'    where
#'     \itemize{
#'       \item \eqn{q} is the threshold exceedance probability, estimated by
#'         the proportion of threshold exceedances,
#'       \item \eqn{N_0} is the number of uncensored sample \eqn{K}-gaps that
#'         are equal to zero,
#'       \item (apart from an adjustment for the contributions of \eqn{S_0}
#'         and \eqn{S_N}) \eqn{N_1} is the number of positive sample
#'         \eqn{K}-gaps,
#'       \item specifically, if \code{inc_cens = TRUE} then \eqn{N_1} is equal
#'         to the number of \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)}
#'         that are positive plus \eqn{(I_0 + I_N) / 2}, where \eqn{I_0 = 1} if
#'         \eqn{S_0} is greater than zero and \eqn{I_0 = 0} otherwise, and
#'         similarly for \eqn{I_N}.
#'     }
#'    The differing treatment of uncensored and right-censored \eqn{K}-gaps
#'    reflects differing contributions to the likelihood.  Right-censored
#'    \eqn{K}-gaps that are equal to zero add no information to the likelihood.
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
#'     \item{\code{sum_qs} }{the sum of the (scaled) \eqn{K}-gaps, that is,
#'       \eqn{q (S_0 + \cdots + S_N)}{q (S_0 + ... + S_N)}, where \eqn{q}
#'       is estimated by the proportion of threshold exceedances.}
#'     \item{\code{n_kgaps} }{the number of \eqn{K}-gaps that contribute to the
#'       log-likelihood.}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' kgaps_stat(newlyn, u)
#' @export
kgaps_stat <- function(data, u, q_u, k = 1, inc_cens = TRUE) {
  if (missing(q_u)) {
    q_u <- mean(data > u, na.rm = TRUE)
  }
  data <- stats::na.omit(data)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (!is.numeric(k) || k < 0 || length(k) != 1) {
    stop("k must be a non-negative scalar")
  }
  # If all the data are smaller than the threshold then return null results
  if (u >= max(data, na.rm = TRUE)) {
    return(list(N0 = 0, N1 = 0, sum_qs = 0, n_kgaps = 0))
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k, 0)
  # N0, N1, sum of scaled K-gaps
  N1 <- sum(S_k > 0)
  N0 <- N_u - 1 - N1
  sum_qs <- sum(q_u * S_k)
  # Store the number of K-gaps, for use by nobs.kgaps()
  n_kgaps <- N0 + N1
  # Include right-censored inter-exceedance times?
  if (inc_cens) {
    # Right-censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    S_k_cens <- pmax(T_u_cens - k, 0)
    # N0, N1, sum of scaled K-gaps
    # S_k_cens = 0 adds no information, because P(S >= 0) = 1
    N1_cens <- sum(S_k_cens > 0)
    n_kgaps <- n_kgaps + N1_cens
    # Remove the right-censored K-gaps that are equal to zero
    # (This is is not necessary here, but we do it for clarity)
    S_k_cens <- S_k_cens[S_k_cens > 0]
    sum_s_cens <- sum(q_u * S_k_cens)
    # Add contributions.
    # Note: we divide N1_cens by two because a right-censored non-zero K-gap
    # S_c contributes theta exp(-theta q_u S_c) to the K-gaps likelihood,
    # whereas a non-censored non-zero K-gap contributes
    # theta^2 exp(-theta q_u S_c).
    # See equation (4.3) of Attalides (2015)
    N1 <- N1 + N1_cens / 2
    sum_qs <- sum_qs + sum_s_cens
  }
  return(list(N0 = N0, N1 = N1, sum_qs = sum_qs, n_kgaps = n_kgaps))
}
