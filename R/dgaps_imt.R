# ================================= dgaps_imt =============================== #

#' Information matrix test under the \eqn{D}-gaps model
#'
#' Performs an information matrix test (IMT) to diagnose misspecification of
#' the \eqn{D}-gaps model of Holesovsky and Fusek (2020).
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then \code{\link{split_by_NAs}} is
#'   used to divide the data into sequences of non-missing values.
#' @param u,D Numeric vectors.  \code{u} is a vector of extreme value
#'   thresholds applied to data.  \code{D} is a vector of values of the
#'   left-censoring parameter \eqn{D}, as defined in Holesovsky and Fusek
#'   (2020). See \code{\link{dgaps}}.
#'
#'   Any values in \code{u} that are greater than all the observations in
#'   \code{data} will be removed without a warning being given.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations. See \code{\link{dgaps}}.
#' @details The general approach follows Suveges and Davison (2010).
#'   The \eqn{D}-gaps IMT is performed a over grid of all
#'   combinations of threshold and \eqn{D} in the vectors \code{u}
#'   and \code{D}.  If the estimate of \eqn{\theta} is 0 then the
#'   IMT statistic, and its associated \eqn{p}-value is \code{NA}.
#' @return An object (a list) of class \code{c("dgaps_imt", "exdex")}
#'   containing
#'   \item{imt }{A \code{length(u)} by \code{length(D)} numeric matrix.
#'     Column i contains, for \eqn{D} = \code{D[i]}, the values of the
#'     information matrix test statistic for the set of thresholds in
#'     \code{u}.  The column names are the values in \code{D}.
#'     The row names are the approximate empirical percentage quantile levels
#'     of the thresholds in \code{u}.}
#'   \item{p }{A \code{length(u)} by \code{length(D)} numeric matrix
#'     containing the corresponding \eqn{p}-values for the test.}
#'   \item{theta }{A \code{length(u)} by \code{length(D)} numeric matrix
#'     containing the corresponding estimates of \eqn{\theta}.}
#'   \item{u,D }{The input \code{u} and \code{D}.}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{dgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{D}-gaps model.
#' @examples
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt <- dgaps_imt(newlyn, u = u, D = 1:5)
#'
#' ### S&P 500 index
#'
#' u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
#' imt <- dgaps_imt(sp500, u = u, D = 1:5)
#'
#' ### Cheeseboro wind gusts (a matrix containing some NAs)
#'
#' probs <- c(seq(0.5, 0.98, by = 0.025), 0.99)
#' u <- quantile(cheeseboro, probs = probs, na.rm = TRUE)
#' imt <- dgaps_imt(cheeseboro, u = u, D = 1:5)
#'
#' ### Uccle July temperatures
#'
#' probs <- c(seq(0.7, 0.98, by = 0.025), 0.99)
#' u <- quantile(uccle720m, probs = probs, na.rm = TRUE)
#' imt <- dgaps_imt(uccle720m, u = u, D = 1:5)
#' @export
dgaps_imt <- function(data, u, D = 1, inc_cens = TRUE) {
  if (any(D < 0)) {
    stop("D must be non-negative")
  }
  # Remove any thresholds that are greater than all the observations
  u_ok <- vapply(u, function(u) any(data > u), TRUE)
  u <- u[u_ok]
  # If there are missing values then use split_by_NAs to extract sequences
  # of non-missing values
  if (anyNA(data) && is.null(attr(data, "split_by_NAs_done"))) {
    data <- split_by_NAs(data)
  } else if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  # Function to perform the IMT for individual values of u and k
  imt_by_uD <- function(u, D) {
    # If the threshold is not high enough then return NAs
    if (u >= max(data, na.rm = TRUE)) {
      return(c(Tn = NA, pvalue = NA))
    }
    # Estimate theta.  If the estimated SE is missing (the observed information
    # was not positive) then return NAs for the test statistic and p-value
    temp <- dgaps(data, u, D, inc_cens = inc_cens)
    theta <- temp$theta
    if (is.na(temp$se)) {
      return(c(imt = NA, p = NA, theta = theta))
    }
    # Contributions to the test statistic from each observation, returning a list
    # with a list of (ldj, Ij, Jj, dj, Ddj, n_dgaps) for each column in data
    imt_stats_list <- apply(data, 2, dgaps_imt_stat, theta = theta, u = u,
                            D = D, inc_cens = inc_cens)
    # Concatenate the results from different columns
    sc <- Reduce(f = function(...) Map(c, ...), imt_stats_list)
    # Calculate the components of the test statistic
    ndgaps <- sum(sc$n_dgaps)
    In <- sum(sc$Ij) / ndgaps
    Jn <- sum(sc$Jj) / ndgaps
    Dn <- Jn - In
    Dnd <- sum(sc$Ddj) / ndgaps
    Vn <- sum((sc$dj - Dnd * sc$ldj / In) ^ 2) / ndgaps
    Tn <- ndgaps * Dn ^ 2 / Vn
    pvalue <- stats::pchisq(Tn, df = 1, lower.tail = FALSE)
    return(c(imt = Tn, p = pvalue, theta = theta))
  }
  uD <- expand.grid(u = u, D = D)
  res <- mapply(imt_by_uD, u = uD$u, D = uD$D)
  imt <- matrix(res[1, ], nrow = length(u), ncol = length(D), byrow = FALSE)
  p <- matrix(res[2, ], nrow = length(u), ncol = length(D), byrow = FALSE)
  theta <- matrix(res[3, ], nrow = length(u), ncol = length(D), byrow = FALSE)
  colnames(imt) <- colnames(p) <- colnames(theta) <- D
  if (is.null(names(u))) {
    u_ps <- round(100 * sapply(u, function(x) mean(data < x)))
  } else {
    u_ps <- as.numeric(substr(names(u), 1, nchar(names(u), type = "c") - 1))
  }
  rownames(imt) <- rownames(p) <- rownames(theta) <- u_ps
  res <- list(imt = imt, p = p, theta = theta, u = u, D = D)
  class(res) <- c("dgaps_imt", "exdex")
  return(res)
}

# ============================= dgaps_imt_stat ============================== #

#' Statistics for the \eqn{D}-gaps information matrix test
#'
#' Calculates the components required to calculate the value of the information
#' matrix test under the \eqn{D}-gaps model, using vector data input.
#' Called by \code{\link{dgaps_imt}}.
#'
#' @param data A numeric vector of raw data.  Missing values are allowed, but
#'   they should not appear between non-missing values, that is, they only be
#'   located at the start and end of the vector.  Missing values are omitted
#'   using \code{\link[stats]{na.omit}}.
#' @param theta A numeric scalar. An estimate of the extremal index
#'  \eqn{\theta}, produced by \code{\link{dgaps}}.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param D A numeric scalar.  The censoring parameter \eqn{D}. Threshold
#'   inter-exceedances times that are not larger than \code{D} units are
#'   left-censored, occurring with probability
#'   \eqn{\log(1 - \theta e^{-\theta d})}{log(1 - exp(-\theta d))},
#'   where \eqn{d = q D} and \eqn{q} is the probability with which the
#'   threshold \eqn{u} is exceeded.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations. See \code{\link{dgaps}}.
#' @return A list
#' relating the quantities given on pages 18-19 of
#'   Suveges and Davison (2010).  All but the last component are vectors giving
#'   the contribution to the quantity from each \eqn{D}-gap, evaluated at the
#'   input value \code{theta} of \eqn{\theta}.
#'   \item{\code{ldj} }{the derivative of the log-likelihood with respect to
#'     \eqn{\theta} (the score)}
#'   \item{\code{Ij} }{the observed information}
#'   \item{\code{Jj} }{the square of the score}
#'   \item{\code{dj} }{\code{Jj} - \code{Ij}}
#'   \item{\code{Ddj} }{the derivative of \code{Jj} - \code{Ij} with respect
#'     to \eqn{\theta}}
#'   \item{\code{n_dgaps} }{the number of \eqn{D}-gaps that contribute to the
#'     log-likelihood.}
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @export
dgaps_imt_stat <- function(data, theta, u, D = 1, inc_cens = TRUE) {
  data <- stats::na.omit(data)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (!is.numeric(D) || D < 0 || length(D) != 1) {
    stop("D must be a non-negative scalar")
  }
  #
  # Calculate the statistics in log-likelihood, as in dgaps_stat()
  #
  # If all the data are smaller than the threshold then return null results
  if (u >= max(data, na.rm = TRUE)) {
    return(list(ldj = 0, Ij = 0, Jj = 0, dj = 0, Ddj = 0, n_dgaps = 0))
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and left-censoring indicator
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
  # Multipliers for terms in ld, ...
  mldj <- rep_len(2, n_dgaps)
  mIj <- rep_len(2, n_dgaps)
  mDdj1 <- rep_len(4, n_dgaps)
  mDdj2 <- rep_len(4, n_dgaps)
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
    # Supplement the inter-exceedance times with the right-censored
    # inter-exceedance times that are greater than D
    T_u <- c(T_u, T_gt_D_cens)
    # Supplement the multipliers
    mldj <- c(mldj, rep_len(1, N1_cens))
    mIj <- c(mIj, rep_len(1, N1_cens))
    mDdj1 <- c(mDdj1, rep_len(2, N1_cens))
    mDdj2 <- c(mDdj2, rep_len(0, N1_cens))
  }
  # Calculate the statistics in the IMT
  # An estimate theta = 1 occurs if none of the inter-exceedance times are
  # left-censored, that is, if they are all greater than D so that N0 = 0
  # in this case we never attempt to divide by 0.
  # An estimate theta = 0 occurs only if all the inter-exceedance times are
  # left-censored (T_u <= d):
  #  in this case we divide by zero in calculating Ddj.
  # If this happens then we convert the NaN to NA.
  #
  # Note: all the right-censored inter-exceedance times that are included in
  # T_u satisfy T_u > D, so the T_u == 0 terms have not contribution from the
  # right-censored observations
  gd <- gd_theta(theta, q_u, D)
  gdd <- gdd_theta(theta, q_u, D)
  ldj <- ifelse(T_u > D, mldj / theta - q_u * T_u, gd)
  Ij <- ifelse(T_u > D, mIj / theta ^ 2, -gdd)
  # If N1 = 0 and D > 0 then the estimate of theta is 0 and the observed
  # information is not well-behaved : not even constrained to be positive.
  # Therefore, we set Ij to NA in this case so that the IMT will also be NA.
  if (N1 == 0 && D > 0) {
    Ij <- NA
  }
  Jj <- ldj ^ 2
  dj <- Jj - Ij
  DdjTgtD <- 2 * gd_theta(theta, q_u, D) * gdd_theta(theta, q_u, D) +
    gddd_theta(theta, q_u, D)
  Ddj <- ifelse(T_u > D, mDdj1 * q_u * T_u / theta ^ 2 - mDdj2 / theta ^ 3,
                DdjTgtD)
  res <- list(ldj = ldj, Ij = Ij, Jj = Jj, dj = dj, Ddj = Ddj,
              n_dgaps = n_dgaps)
  return(res)
}
