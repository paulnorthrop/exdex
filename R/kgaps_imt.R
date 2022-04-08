# ================================= kgaps_imt =============================== #

#' Information matrix test under the \eqn{K}-gaps model
#'
#' Performs the information matrix test (IMT) of Suveges and Davison (2010) to
#' diagnose misspecification of the \eqn{K}-gaps model.
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then \code{\link{split_by_NAs}} is
#'   used to divide the data into sequences of non-missing values.
#' @param u,k Numeric vectors.  \code{u} is a vector of extreme value
#'   thresholds applied to data.  \code{k} is a vector of values of the run
#'   parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps}} for more details.
#'
#'   Any values in \code{u} that are greater than all the observations in
#'   \code{data} will be removed without a warning being given.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times, relating to the
#'   first and last observations.  See Attalides (2015) for details.
#' @details The \eqn{K}-gaps IMT is performed a over grid of all
#'   combinations of threshold and \eqn{K} in the vectors \code{u}
#'   and \code{k}.  If the estimate of \eqn{\theta} is 0 then the
#'   IMT statistic, and its associated \eqn{p}-value is \code{NA}.
#'
#'   For details of the IMT see Suveges and Davison
#'   (2010).  There are some typing errors on pages 18-19 that have been
#'   corrected in producing the code: the penultimate term inside \code{{...}}
#'   in the middle equation on page 18 should be \eqn{(c_j(K))^2}, as should
#'   the penultimate term in the first equation on page 19; the \code{{...}}
#'   bracket should be squared in the 4th equation on page 19; the factor
#'   \eqn{n} should be \eqn{N-1} in the final equation on page 19.
#' @return An object (a list) of class \code{c("kgaps_imt", "exdex")}
#'   containing
#'   \item{imt }{A \code{length(u)} by \code{length(k)} numeric matrix.
#'     Column i contains, for \eqn{K} = \code{k[i]}, the values of the
#'     information matrix test statistic for the set of thresholds in
#'     \code{u}.  The column names are the values in \code{k}.
#'     The row names are the approximate empirical percentage quantile levels
#'     of the thresholds in \code{u}.}
#'   \item{p }{A \code{length(u)} by \code{length(k)} numeric matrix
#'     containing the corresponding \eqn{p}-values for the test.}
#'   \item{theta }{A \code{length(u)} by \code{length(k)} numeric matrix
#'     containing the corresponding estimates of \eqn{\theta}.}
#'   \item{u,k }{The input \code{u} and \code{k}.}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{choose_uk}} for graphical diagnostic to aid the choice
#'   of the threshold \eqn{u} and the run parameter \eqn{K}.
#' @examples
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt <- kgaps_imt(newlyn, u = u, k = 1:5)
#'
#' ### S&P 500 index
#'
#' u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
#' imt <- kgaps_imt(sp500, u = u, k = 1:5)
#'
#' ### Cheeseboro wind gusts (a matrix containing some NAs)
#'
#' probs <- c(seq(0.5, 0.98, by = 0.025), 0.99)
#' u <- quantile(cheeseboro, probs = probs, na.rm = TRUE)
#' imt <- kgaps_imt(cheeseboro, u = u, k = 1:5)
#' @export
kgaps_imt <- function(data, u, k = 1, inc_cens = TRUE) {
  if (any(k < 0)) {
    stop("k must be non-negative")
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
  imt_by_uk <- function(u, k) {
    # If the threshold is not high enough then return NAs
    if (u >= max(data, na.rm = TRUE)) {
      return(c(Tn = NA, pvalue = NA))
    }
    # Estimate theta
    theta <- kgaps(data, u, k, inc_cens = inc_cens)$theta
    # Contributions to the test statistic from each observation, returning a list
    # with a list of (ldj, Ij, Jj, dj, Ddj, n_kgaps) for each column in data
    imt_stats_list <- apply(data, 2, kgaps_imt_stat, theta = theta, u = u,
                            k = k, inc_cens = inc_cens)
    # Concatenate the results from different columns
    sc <- Reduce(f = function(...) Map(c, ...), imt_stats_list)
    # Calculate the components of the test statistic
    nkgaps <- sum(sc$n_kgaps)
    In <- sum(sc$Ij) / nkgaps
    Jn <- sum(sc$Jj) / nkgaps
    Dn <- Jn - In
    Dnd <- sum(sc$Ddj) / nkgaps
    Vn <- sum((sc$dj - Dnd * sc$ldj / In) ^ 2) / nkgaps
    Tn <- nkgaps * Dn ^ 2 / Vn
    pvalue <- stats::pchisq(Tn, df = 1, lower.tail = FALSE)
    return(c(imt = Tn, p = pvalue, theta = theta))
  }
  uk <- expand.grid(u = u, k = k)
  res <- mapply(imt_by_uk, u = uk$u, k = uk$k)
  imt <- matrix(res[1, ], nrow = length(u), ncol = length(k), byrow = FALSE)
  p <- matrix(res[2, ], nrow = length(u), ncol = length(k), byrow = FALSE)
  theta <- matrix(res[3, ], nrow = length(u), ncol = length(k), byrow = FALSE)
  colnames(imt) <- colnames(p) <- colnames(theta) <- k
  if (is.null(names(u))) {
    u_ps <- round(100 * sapply(u, function(x) mean(data < x)))
  } else {
    u_ps <- as.numeric(substr(names(u), 1, nchar(names(u), type = "c") - 1))
  }
  rownames(imt) <- rownames(p) <- rownames(theta) <- u_ps
  res <- list(imt = imt, p = p, theta = theta, u = u, k = k)
  class(res) <- c("kgaps_imt", "exdex")
  return(res)
}

# ============================= kgaps_imt_stat ============================== #

#' Statistics for the information matrix test
#'
#' Calculates the components required to calculate the value of the information
#' matrix test under the \eqn{K}-gaps model, using vector data input.
#' Called by \code{\link{kgaps_imt}}.
#'
#' @param data A numeric vector of raw data.  Missing values are allowed, but
#'   they should not appear between non-missing values, that is, they only be
#'   located at the start and end of the vector.  Missing values are omitted
#'   using \code{\link[stats]{na.omit}}.
#' @param theta A numeric scalar. An estimate of the extremal index
#'  \eqn{\theta}, produced by \code{\link{kgaps}}.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @return A list relating the quantities given on pages 18-19 of
#'   Suveges and Davison (2010).  All but the last component are vectors giving
#'   the contribution to the quantity from each \eqn{K}-gap, evaluated at the
#'   input value \code{theta} of \eqn{\theta}.
#'   \item{\code{ldj} }{the derivative of the log-likelihood with respect to
#'     \eqn{\theta} (the score)}
#'   \item{\code{Ij} }{the observed information}
#'   \item{\code{Jj} }{the square of the score}
#'   \item{\code{dj} }{\code{Jj} - \code{Ij}}
#'   \item{\code{Ddj} }{the derivative of \code{Jj} - \code{Ij} with respect
#'     to \eqn{\theta}}
#'   \item{\code{n_kgaps} }{the number of \eqn{K}-gaps that contribute to the
#'     log-likelihood.}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{https://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @export
kgaps_imt_stat <- function(data, theta, u, k = 1, inc_cens = TRUE) {
  data <- stats::na.omit(data)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (!is.numeric(k) || k < 0 || length(k) != 1) {
    stop("k must be a non-negative scalar")
  }
  #
  # Calculate the statistics in log-likelihood, as in kgaps_stat()
  #
  # If all the data are smaller than the threshold then return null results
  if (u >= max(data, na.rm = TRUE)) {
    return(list(ldj = 0, Ij = 0, Jj = 0, dj = 0, Ddj = 0, n_kgaps = 0))
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k, 0)
  # N0, N1, sum of scaled K-gaps
  N1 <- sum(S_k > 0)
  N0 <- N_u - 1 - N1
  sum_qs <- sum(q_u * S_k)
  # Store the number of K-gaps
  n_kgaps <- N0 + N1
  # Values of c^(K) = q_u * S^(K)
  qS <- q_u * S_k
  # Multipliers for terms in ld, ...
  mldj <- rep_len(2, n_kgaps)
  mIj <- rep_len(2, n_kgaps)
  mDdj1 <- rep_len(4, n_kgaps)
  mDdj2 <- rep_len(4, n_kgaps)
  # Include censored inter-exceedance times?
  if (inc_cens) {
    # censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    S_k_cens <- pmax(T_u_cens - k, 0)
    # N0, N1, sum of scaled K-gaps
    # S_k_cens = 0 adds no information, because P(S >= 0) = 1
    N1_cens <- sum(S_k_cens > 0)
    n_kgaps <- n_kgaps + N1_cens
    # Remove the censored K-gaps that are equal to zero
    S_k_cens <- S_k_cens[S_k_cens > 0]
    sum_s_cens <- sum(q_u * S_k_cens)
    # Add contributions.
    # Note: we divide N1_cens by two because a censored non-zero K-gap S_c
    # contributes theta exp(-theta q_u S_c) to the K-gaps likelihood,
    # whereas a non-censored non-zero K-gap contributes
    # theta^2 exp(-theta q_u S_c).
    # See equation (4.3) of Attalides (2015)
    N1 <- N1 + N1_cens / 2
    sum_qs <- sum_qs + sum_s_cens
    qS_cens <- q_u * S_k_cens
    # Supplement the c^(K) terms
    qS <- c(qS, qS_cens)
    # Supplement the multipliers
    mldj <- c(mldj, rep_len(1, N1_cens))
    mIj <- c(mIj, rep_len(1, N1_cens))
    mDdj1 <- c(mDdj1, rep_len(2, N1_cens))
    mDdj2 <- c(mDdj2, rep_len(0, N1_cens))
  }
  # Calculate the statistics in the IMT
  # An estimate theta = 1 occurs if all the K-gaps are positive (qS > 0):
  #  in this case we never attempt to divide by 0.
  # An estimate theta = 0 occurs only if all the K-gaps are zero (qS = 0):
  #  in this case we divide by zero in calculating Ddj (mDdj1 * qS / theta ^ 2)
  # If this happens then we convert the NaN to NA.
  # Note: all the right-censored K-gaps have qS > 0, so the qS == 0 terms
  #  have no contribution from the right-censored observations
  ldj <- ifelse(qS == 0, -1 / (1 - theta), mldj / theta) - qS
  Ij  <- ifelse(qS == 0, 1 / (1 - theta) ^ 2, mIj / theta ^ 2)
  Jj <- ldj ^ 2
  dj <- Jj - Ij
  Ddj <- mDdj1 * qS / theta ^ 2 - ifelse(qS == 0, 0, mDdj2 / theta ^ 3)
  Ddj[is.nan(Ddj)] <- NA
  res <- list(ldj = ldj, Ij = Ij, Jj = Jj, dj = dj, Ddj = Ddj,
              n_kgaps = n_kgaps)
  return(res)
}
