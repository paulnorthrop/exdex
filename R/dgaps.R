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
#'   If \code{data} has multiple columns then there will be right-censored
#'   first and last inter-exceedance times for each column.
#' @details If \code{inc_cens = FALSE} then the maximum likelihood estimate of
#'   the extremal index \eqn{\theta} under the \eqn{D}-gaps model of
#'   Holesovsky and Fusek (2020) is calculated.  Under this model
#'   inter-exceedance times that are less than or equal to \eqn{D} are
#'   left-censored, as a strategy to mitigate model mis-specification resulting
#'   from the fact that inter-exceedance times that are equal to 0 are expected
#'   under the model but only positive inter-exceedance times can be observed
#'   in practice.
#'
#'   If \code{inc_cens = TRUE} then information from the right-censored
#'   first and last inter-exceedance times are also included in the likelihood
#'   to be maximized.
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
#'       using an algebraic expression for the observed information.  If the
#'       estimate of \eqn{\theta} is 0 then \code{se} is \code{NA}.}
#'     \item{\code{se_exp} }{The estimated standard error of the MLE,
#'       calculated using an algebraic expression for the expected information.
#'       If the estimate of \eqn{\theta} is 0 then \code{se_exp} is \code{NA}.
#'       This is provided because cases may be encountered where the observed
#'       information is not positive.}
#'     \item{\code{ss} }{The list of summary statistics returned from
#'       \code{\link{dgaps_stat}}.}
#'     \item{\code{D, u, inc_cens} }{The input values of \code{D},
#'       \code{u} and \code{inc_cens}.}
#'     \item{\code{max_loglik }}{The value of the log-likelihood at the MLE.}
#'     \item{\code{call }}{The call to \code{dgaps}.}
#' @seealso \code{\link{dgaps_confint}} to estimate confidence intervals
#'   for \eqn{\theta}.
#' @seealso \code{\link{dgaps_methods}} for S3 methods for \code{"dgaps"}
#'   objects.
#' @seealso \code{\link{dgaps_imt}} for the information matrix test, which
#'   may be used to inform the choice of the pair (\code{u, D}).
#' @seealso \code{\link{choose_ud}} for a diagnostic plot based on
#'   \code{\link{dgaps_imt}}.
#' @seealso \code{\link{dgaps_stat}} for the calculation of sufficient
#'   statistics for the \eqn{D}-gaps model.
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
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
#' coef(theta)
#' nobs(theta)
#' vcov(theta)
#' logLik(theta)
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = 0.60)
#' theta <- dgaps(newlyn, u = u, D = 2)
#' theta
#' summary(theta)
#'
#' ### Uccle July temperatures
#'
#' # Using vector input, which merges data from different years
#' u <- quantile(uccle720$temp, probs = 0.9, na.rm = TRUE)
#' theta <- dgaps(uccle720$temp, u = u, D = 2)
#' theta
#'
#' # Using matrix input to separate data from different years
#' u <- quantile(uccle720m, probs = 0.9, na.rm = TRUE)
#' theta <- dgaps(uccle720m, u = u, D = 2)
#' theta
#' @export
dgaps <- function(data, u, D = 1, inc_cens = TRUE) {
  Call <- match.call(expand.dots = TRUE)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data, na.rm = TRUE)) {
    stop("'u' must be less than 'max(data, na.rm = TRUE)'")
  }
  if (!is.numeric(D) || D < 0 || length(D) != 1) {
    stop("D must be a non-negative scalar")
  }
  # If there are missing values then use split_by_NAs to extract sequences
  # of non-missing values
  if (anyNA(data) && is.null(attr(data, "split_by_NAs_done"))) {
    data <- split_by_NAs(data)
  }
  # Estimate the marginal exceedance probability q_u
  q_u <- mean(data > u, na.rm = TRUE)
  # Calculate sufficient statistics for each column in data and then sum
  stats_list <- apply(as.matrix(data), 2, dgaps_stat, u = u, q_u = q_u, D = D,
                      inc_cens = inc_cens)
  ss <- Reduce(f = function(...) Map("+", ...), stats_list)
  # Add q_u and D to the list of arguments to be passed to functions that
  # calculated quantities based on the log-likelihood
  ss$q_u <- q_u
  ss$D <- D
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
    # Use the K-gaps estimate as an initial estimate
    theta_init <- kgaps(data = data, u = u, k = D, inc_cens = inc_cens)$theta
    dgaps_negloglik <- function(theta) {
      return(-do.call(dgaps_loglik, c(list(theta = theta), ss)))
    }
    temp <- stats::optim(theta_init, dgaps_negloglik, method = "Brent",
                         lower = 0, upper = 1)
    theta_mle <- temp$par
  }
  # Estimate standard error.  In some cases there will be problems with the
  # observed information.  Therefore, also calculate estimated standard errors
  # based on the expected information, using a modified (for inc_cens= TRUE)
  # version of equation (11) on page 202 of Holesovsky and Fusek (2020).
  # If N1 = 0 then the estimate of theta is 0 and we return NA for se_exp
  if (N1 > 0) {
    exp_info <- dgaps_exp_info(theta = theta_mle, ss = ss, inc_cens = inc_cens)
  } else {
    exp_info <- NA
  }
  #
  # If N1 = 0 then the estimate of theta is 0. The contribution to obs_info
  # from the N0 > 0 case is not constrained to be positive unless D = 0 (when
  # the calculation is the same as K-gaps).  Therefore, we set the SE to NA
  # if N1 = 0 unless D = 0. Note: at least one of N0 and N1 must be positive.
  obs_info <- 0
  if (N0 > 0) {
    if (N1 > 0 || D == 0) {
      obs_info <- obs_info - N0 * gdd_theta(theta_mle, q_u = q_u, D = D)
    } else {
      obs_info <- NA
    }
  }
  if (N1 > 0) {
    obs_info <- obs_info + 2 * N1 / theta_mle ^ 2
  }
  # The observed information is not guaranteed to be positive
  # If it is not positive then return NA for the estimated SE
  # Similarly for the expected information
  if (!is.na(obs_info) && obs_info <= 0) {
    theta_se <- NA
    se_exp <- NA
  } else {
    theta_se <- sqrt(1 / obs_info)
    se_exp <- 1 / sqrt(exp_info)
  }
  max_loglik <- do.call(dgaps_loglik, c(list(theta = theta_mle), ss))
  res <- list(theta = theta_mle, se = theta_se, se_exp = se_exp, ss = ss,
              D = D, u = u, inc_cens = inc_cens, max_loglik = max_loglik,
              call = Call)
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
#' @param q_u A numeric scalar.  An estimate of the probability with which
#'   the threshold \code{u} is exceeded.  If \code{q_u} is missing then it is
#'   calculated using \code{mean(data > u, na.rm = TRUE)}.
#' @param D A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times relating to the
#'   first and last observation.  It is known that these times are greater
#'   than or equal to the time observed. See Attalides (2015) for details.
#' @details The sample inter-exceedance times are
#'   \eqn{T_0, T_1, ..., T_{N-1}, T_N}{T_0, T_1, ..., T_(N-1), T_N},
#'   where \eqn{T_1, ..., T_{N-1}}{T_1, ..., T_(N-1)} are uncensored and
#'   \eqn{T_0} and \eqn{T_N} are right-censored.  Under the assumption that the
#'   inter-exceedance times are independent, the log-likelihood of the
#'   \eqn{D}-gaps model is given by
#'   \deqn{l(\theta; T_0, \ldots, T_N) = N_0 \log(1 - \theta e^{-\theta d}) +
#'     2 N_1 \log \theta - \theta q (I_0 T_0 + \cdots + I_N T_N),}{%
#'     l(\theta; T_0, ..., T_N) = N_0 log(1 - exp(-\theta d)) +
#'     2 N_1 log \theta - \theta q (I_0 T_0 + ... + I_N T_N),}
#'    where
#'     \itemize{
#'       \item{\eqn{q} is the threshold exceedance probability, estimated by
#'         the proportion of threshold exceedances,}
#'       \item{\eqn{d = q D},}
#'       \item{\eqn{I_j = 1} if \eqn{T_j > D} and \eqn{I_j = 0} otherwise,}
#'       \item{\eqn{N_0} is the number of sample inter-exceedance times that
#'         are left-censored, that is, are less than or equal to \eqn{D},}
#'       \item{(apart from an adjustment for the contributions of \eqn{T_0} and
#'         \eqn{T_N}) \eqn{N_1} is the number of inter-exceedance times that
#'         are uncensored, that is, are greater than \eqn{D},}
#'       \item{specifically, if \code{inc_cens = TRUE} then \eqn{N_1} is equal
#'         to the number of \eqn{T_1, ..., T_{N-1}}{T_1, ..., T_(N-1)} that are
#'         uncensored plus \eqn{(I_0 + I_N) / 2}.}
#'     }
#'    The differing treatment of uncensored and censored \eqn{K}-gaps reflects
#'    differing contributions to the likelihood. Right-censored
#'    inter-exceedance times whose observed values are less than or equal to
#'    \eqn{D} add no information to the likelihood because we do not know to
#'    which part of the likelihood they should contribute.
#'
#'    If \eqn{N_1 = 0} then we are in the degenerate case where there is one
#'    cluster (all inter-exceedance times are left-censored) and the likelihood
#'    is maximized at \eqn{\theta = 0}.
#'
#'    If \eqn{N_0 = 0} then all exceedances occur singly (no inter-exceedance
#'    times are left-censored) and the likelihood is maximized at
#'    \eqn{\theta = 1}.
#' @return A list containing the sufficient statistics, with components
#'     \item{\code{N0} }{the number of left-censored inter-exceedance times.}
#'     \item{\code{N1} }{contribution from inter-exceedance times that are not
#'       left-censored (see \strong{Details}).}
#'     \item{\code{sum_qtd} }{the sum of the (scaled) inter-exceedance times
#'       that are not left-censored, that is,
#'       \eqn{q (I_0 T_0 + \cdots + I_N T_N)}{q (I_0 T_0 + ... + I_N T_N)},
#'       where \eqn{q} is estimated by the proportion of threshold
#'       exceedances.}
#'     \item{\code{n_dgaps} }{the number of inter-exceedances that contribute
#'       to the log-likelihood.}
#'     \item{\code{q_u} }{the sample proportion of values that exceed the
#'       threshold.}
#'     \item{\code{D} }{the input value of \code{D}.}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{dgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{D}-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' dgaps_stat(newlyn, u = u, D = 1)
#' @export
dgaps_stat <- function(data, u, q_u, D = 1, inc_cens = TRUE) {
  if (missing(q_u)) {
    q_u <- mean(data > u, na.rm = TRUE)
  }
  data <- stats::na.omit(data)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (!is.numeric(D) || D < 0 || length(D) != 1) {
    stop("D must be a non-negative scalar")
  }
  # If all the data are smaller than the threshold then return null results
  if (u >= max(data, na.rm = TRUE)) {
    return(list(N0 = 0, N1 = 0, sum_qtd = 0, n_dgaps = 0))
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  # Inter-exceedances times and left-censoring indicator
  T_u <- diff(exc_u)
  left_censored <- T_u <= D
  # N0, N1, sum of scaled inter-exceedance times that are greater than d,
  # that is, not left-censored
  N1 <- sum(!left_censored)
  N0 <- N_u - 1 - N1
  T_gt_D <- T_u[!left_censored]
  sum_qtd <- sum(q_u * T_gt_D)
  # Store the number of D-gaps, for use by nobs.dgaps()
  n_dgaps <- N0 + N1
  # Include censored inter-exceedance times?
  if (inc_cens) {
    # censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    # T_u_cens <= d adds no information, because we have no idea to which part
    # of the log-likelihood they would contribute
    left_censored_cens <- T_u_cens <= D
    # N0, N1, sum of scaled inter-exceedance times that are greater than D,
    # that is, not left-censored
    N1_cens <- sum(!left_censored_cens)
    n_dgaps <- n_dgaps + N1_cens
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
  return(list(N0 = N0, N1 = N1, sum_qtd = sum_qtd, n_dgaps = n_dgaps))
}
