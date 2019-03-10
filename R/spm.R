# ============================== spm() ====================================== #

#' Semiparametric maxima estimator of the extremal index
#'
#' Estimates the extremal index \eqn{\theta} using a semiparametric block
#' maxima estimator of Northrop (2015) and a variant of this estimator
#' studied by Berghaus and Bucher (2018).  Estimates of uncertainty are
#' calculated using the asymptotic theory developed by Berghaus and
#' Bucher (2018).
#'
#' @param data A numeric vector of raw data.
#' @param b A numeric scalar.  The block size.
#' @param bias_adjust A character scalar.  Is bias-adjustment of the
#'   raw estimate of \eqn{\theta} performed using the bias-reduced
#'   estimator (\code{bias_adjust = "BB3"}), derived in Section 5 of
#'   Berghaus and Bucher (2018); or a simpler version
#'   (\code{bias_adjust = "BB1"}), in which the raw estimate is multiplied by
#'   \eqn{(k-1) / k}, where \eqn{k} is the number of blocks; or the
#'   bias-adjustment of the empirical distribution function used to calculate
#'   the estimate, as detailed in Section 2 of Northrop (2015).  When disjoint
#'   maxima are used \code{bias_adjust = "BB1"} and \code{bias_adjust = "N"}
#'   give identical estimates of the Berghaus and Bucher (2018) variant,
#'   as explained at the end of Section 5 of Berghaus and Bucher (2018).
#'   If \code{bias_adjust = "none"} then no bias-adjustment is performed.
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any estimates that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  This is carried out
#'   \emph{after} any bias-adjustment.  Otherwise,
#'   estimates that are greater than 1 may be obtained.
#' @param varN A logical scalar.  If \code{varN = TRUE} then the estimation
#'   of the sampling variance of the Northrop (2015) estimator is tailored
#'   to that estimator.  Otherwise, the sampling variance derived in
#'   Berghaus and Bucher (2018) is used.
#'   See \strong{Details} for further information.
#' @param which_dj A character scalar.  Determines which set of disjoint
#'   maxima are used to calculate an estimate of \eqn{\theta}: \code{"first"},
#'   only the set whose first block starts on the first observation in
#'   \code{x}; \code{"last"}, only the set whose last block end on the last
#'   observation in \code{x}.
#' @details The extremal index \eqn{\theta} is estimated using the
#'   semiparametric maxima estimator of Northrop (2015) and variant
#'   of this studied by Berghaus and Bucher (2018).  In each case a sample
#'   of 'data' is derived from the input data \code{data}, based on the
#'   empirical distribution function of these data evaluated at the
#'   maximum values of  of blocks of \code{b} contiguous values in \code{data}.
#'
#'   The estimators are based on an assumption that these 'data' are sampled
#'   approximately from an exponential distribution with mean \eqn{1/\theta}.
#'   For details see page 2309 of Berghaus and Bucher (2018), where the
#'   'data' for the Northrop (2015) estimator are denoted \eqn{Y} and
#'   those for the Berghaus and Bucher (2018) are denoted \eqn{Z}.
#'   For convenience, we will refer to these values as the
#'   \eqn{Y}-data and the \eqn{Z}-data.
#'
#'   If \code{sliding = TRUE} then the function uses sliding block maxima,
#'   that is, the largest value observed in \emph{all}
#'   (\code{length(data) - b + 1}) blocks of \code{b} observations.
#'   If \code{sliding = FALSE} then disjoint block maxima, that is, the
#'   largest values in (\code{floor(n / b)}) disjoint blocks of \code{b}
#'   observations, are used.
#'
#'   Estimation of the sampling variances of the estimators is based
#'   on Proposition 4.1 of Berghaus and Bucher (2018).
#'   For the Northrop (2015) variant the user has the choice either to
#'   use the sampling variance based on the Berghaus and Bucher (2018)
#'   estimator, i.e. the \eqn{Z}-data (\code{varN = FALSE}) or an analogous
#'   version tailored to the Northrop (2015) estimator that uses the
#'   \eqn{Y}-data (\code{varN = TRUE}).
#'
#'   A condition imposed in Proposition 4.1 of Berghaus and Bucher (2018)
#'   means that \code{b} must be no smaller than \eqn{k^{1/2}} and no larger
#'   than \eqn{k^2}, where \eqn{k} is \code{floor(length(data) / b)},
#'   i.e. \eqn{k} is the number of complete blocks.  If this is not the case
#'   then a warning will be given and
#'     \itemize{
#'       \item{estimated standard errors will be missing from the returned
#'             object,}
#'       \item{if \code{bias_adjust == "BB3"} then bias-adjustment
#'             based on \code{bias_adjust == "BB1"} will be performed instead,
#'             because the former relies on the estimated variances of the
#'             estimators.}
#'     }
#' @return A list of class \code{c("exdex", "spm")} containing the
#'   components listed below.  The components that are vectors are
#'   labelled to indicate the estimator to which the constituent values
#'   relate: N2015 for Northrop (2015) and BB2018 for
#'   Berghaus and Bucher (2018).
#'   \item{theta_sl, theta_dj}{ Vectors containing the estimates of
#'     \eqn{\theta} resulting from sliding maxima and disjoint maxima
#'     respectively.}
#'   \item{se_sl, se_dj}{The estimated standard errors associated
#'     with the estimates in \code{theta_sl} and \code{theta_dj}.}
#'   \item{bias_sl, bias_dj}{The respective values of the
#'       bias-adjustment applied to the raw estimates.  This is only
#'       relevant if \code{bias_adjust} is "BB3" or "BB1".  Otherwise,
#'       \code{bias_sl} and \code{bias_dj} are \code{c(NA, NA)}.}
#'   \item{uncon_theta_sl, uncon_theta_dj}{The estimates of \eqn{\theta}
#'     are the constraint that they lie in (0, 1] has been applied.}
#'   \item{data_sl, data_dj}{Vectors of the values of the \eqn{Y}-data and
#'     the \eqn{Z}-data.}
#'   \item{sigma2dj, sigma2dj_for_sl}{Estimates of the variance
#'    \eqn{\sigma_{{\rm dj}}^2}{\sigma^2_dj}
#'    defined on pages 2314-2315 of Berghaus and Bucher (2018).
#'    The form of the estimates is given on page 2319.
#'    \code{sigma2dj} is used in estimating the standard error \code{se_dj},
#'    \code{sigma2dj_for_sl} in estimating the standard error \code{se_sl}.}
#'   \item{sigma2sl}{Estimates of the variance
#'     \eqn{\sigma_{{\rm sl}}^2}{\sigma^2_sl}.
#'     defined on pages 2314-2315 of Berghaus and Bucher (2018).
#'     The form of the estimates is given on page 2319.
#'     \code{sigma2dj_for_sl} is used to estimate
#'     \eqn{\sigma_{{\rm dj}}^2}{\sigma^2_dj} for this purpose.}
#'   \item{b}{The input value of \code{b}.}
#'   \item{bias_adjust}{The input value of \code{bias_adjust}.}
#'   \item{call}{The call to \code{spm}.}
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @seealso \code{\link{confint.exdex}} to estimate confidence intervals
#'   for \eqn{theta}.
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @examples
#' temp <- spm(-as.vector(sp500[2:6550]), 250)
#'
#' temp <- spm(newlyn, 20)
#' @export
spm <- function(data, b, bias_adjust = c("BB3", "BB1", "N", "none"),
                constrain = TRUE, varN = TRUE,
                which_dj = c("last", "first")) {
  Call <- match.call(expand.dots = TRUE)
  #
  # Check inputs
  #
  if (missing(data) || length(data) == 0L || mode(data) != "numeric")
    stop("'data' must be a non-empty numeric vector")
  if (any(!is.finite(data))) {
    stop("'data' contains missing or infinite values")
  }
  if (is.matrix(data)) stop("'data' must be a vector")
  data <- as.vector(data)
  if (!is.numeric(b) || length(b) != 1) {
    stop("'b' must be a numeric scalar (specifically, a positive integer)")
  }
  if (b < 1) {
    stop("'b' must be no smaller than 1")
  }
  bias_adjust <- match.arg(bias_adjust)
  if (!is.logical(constrain) || length(constrain) != 1) {
    stop("'constrain' must be a logical scalar")
  }
  if (!is.logical(varN) || length(varN) != 1) {
    stop("'varN' must be a logical scalar")
  }
  #
  # Check that the value of b satisfies the inequality in Proposition 4.1
  # of Berghaus and Bucher (2018).  If not then we don't calculate estimates
  # of uncertainty and we cannot use the BB3 bias-adjustment.
  #
  k_n <- floor(length(data) / b)
  if (b < sqrt(k_n)) {
    b_ok <- FALSE
    warn1 <- "b is too small"
  } else if (b > k_n ^ 2) {
    b_ok <- FALSE
    warn1 <- "b is too large"
  } else {
    b_ok <- TRUE
  }
  if (!b_ok) {
    warn2 <- "Estimates of uncertainty will be missing"
    if (bias_adjust == "BB3") {
      bias_adjust <- "BB1"
      warn3 <- "'bias_adjust' has been changed to ''BB1''"
      warning("\n", warn1, "\n", warn2, "\n", warn3)
    } else {
      warning("\n", warn1, "\n", warn2)
    }
  }
  which_dj <- match.arg(which_dj)
  #
  # Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
  # We require the disjoint maxima to do this.  If sliding = TRUE then
  # pass these to spm_sigmahat_dj using the dj_maxima argument
  # Only do this is b_ok = TRUE.
  # Otherwise, just calculate point estimates of theta
  # At this point these estimates have not been bias-adjusted, unless
  # bias_adjust = "N".
  if (b_ok) {
    # Find all sets of maxima of disjoint blocks of length b
    all_max <- all_maxima(data, b)
    res <- ests_sigmahat_dj(all_max, b, which_dj, bias_adjust)
  } else {
    all_max <- all_maxima(data, b, which_dj)
    # Disjoint maxima
    Fhaty <- ecdf2(all_max$xd, all_max$yd)
    k_n <- length(all_max$yd)
    m <- length(all_max$xd)
    const <- -log(m - b + k_n)
    if (bias_adjust == "N") {
      Fhaty <- (m * Fhaty - b) / (m - b)
    }
    res <- list()
    res$theta_dj <- c(-1 / mean(b * log0const(Fhaty, const)),
                      1 / (b * mean(1 - Fhaty)))
    names(res$theta_dj) <- c("N2015", "BB2018")
  }
  # Sliding maxima
  Fhaty <- ecdf2(all_max$xs, all_max$ys)
  # Avoid over-writing the `disjoint' sample size k_n: it is needed later
  k_n_sl <- length(all_max$ys)
  m <- length(all_max$xs)
  const <- -log(m - b + k_n_sl)
  if (bias_adjust == "N") {
    Fhaty <- (m * Fhaty - b) / (m - b)
  }
  res$theta_sl <- c(-1 / mean(b * log0const(Fhaty, const)),
                    1 / (b * mean(1 - Fhaty)))
  names(res$theta_sl) <- c("N2015", "BB2018")
  #
  # Add the values of the Y-data and the Z-data to the output
  res$data_sl <- cbind(N2015 = -b * log(Fhaty), BB2018 = b * (1 - Fhaty))
  #
  # Estimate the sampling variances of the estimators
  #
  if (b_ok) {
    res$sigma2sl <- res$sigma2dj_for_sl - (3 - 4 * log(2)) / res$theta_sl ^ 2
    indexN <- ifelse(varN, 2, 1)
    if (varN) {
      index <- 1:2
    } else {
      index <- c(2, 2)
    }
    res$se_dj <- res$theta_dj ^ 2 * sqrt(res$sigma2dj[index] / k_n)
    res$se_sl <- res$theta_sl ^ 2 * sqrt(res$sigma2sl[index] / k_n)
  } else {
    res$se_dj <- c(N2015 = NA, BB2018 = NA)
    res$se_sl <- c(N2015 = NA, BB2018 = NA)
  }
  #
  # Perform BB2018 bias-adjustment if required
  #
  if (bias_adjust == "BB3") {
    res$bias_dj <- res$theta_dj / k_n + res$theta_dj ^ 3 * res$sigma2dj / k_n
    res$bias_sl <- res$theta_sl / k_n + res$theta_sl ^ 3 * res$sigma2sl / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    res$theta_sl <- res$theta_sl - res$bias_sl
  } else if (bias_adjust == "BB1") {
    res$bias_dj <- res$theta_dj / k_n
    res$bias_sl <- res$theta_sl / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    res$theta_sl <- res$theta_sl - res$bias_sl
  } else {
    res$bias_dj <- res$bias_sl <- c(N2015 = NA, BB2018 = NA)
  }
  #
  # Save the unconstrained estimates, so that they can be returned
  res$uncon_theta_dj <- res$theta_dj
  res$uncon_theta_sl <- res$theta_sl
  #
  # Constrain to (0, 1] if required
  if (constrain) {
    res$theta_dj <- pmin(res$theta_dj, 1)
    res$theta_sl <- pmin(res$theta_sl, 1)
  }
  #
  res$bias_adjust <- bias_adjust
  res$b <- b
  res$call <- Call
  class(res) <- c("exdex", "spm")
  return(res)
}

ests_sigmahat_dj <- function(all_max, b, which_dj, bias_adjust){
  # Which of the raw values in x are <= each of the values in y?
  # For each of the block maxima in y calculate the numbers of the raw
  # values in each block that are <= the block maximum
  # k_n is the number of blocks
  k_n <- nrow(all_max$yd)
  # m is the number of raw observations
  m <- nrow(all_max$xd)
  # Value to replace log(0), in the unlikely event that this happens
  const <- -log(m - b + k_n)
  block <- rep(1:k_n, each = b)
  sum_fun <- function(x, y) {
    return(vapply(y, function(y) sum(x <= y), 0))
  }
  # This returns a list with k_n elements.  The ith element of the list
  # (a of length vector k_n) contains the numbers of values in the ith block
  # that are <= each of the block maxima in y
  #
  # Function to calculate Fhaty and UsumN for each set of disjoint block maxima
  UsumN_fn <- function(i) {
    y <- all_max$yd[, i]
    x <- all_max$xd[, i]
    nums_list <- tapply(x, block, sum_fun, y = y)
    # Make this into a matrix
    # Column j contains the numbers of values in the ith block that are <= each
    # of the block maxima in y
    # Row i contains the numbers of values in each block that are <=
    # block maximum i in y
    nums_mat <- do.call(cbind, nums_list)
    # Therefore, the row sums contain the total number of values that are <=
    # each block maximum in y
    # The corresponding proportion is Fhaty in spm(): ecdf of x evaluated at y,
    # in the disjoint (sliding = FALSE) case
    Fhaty <- rowSums(nums_mat) / m
    # For each block, we want an equivalent vector obtained when we delete that
    # block
    fun <- function(i, x) {
      rowSums(x[, -i])
    }
    # Column j contains the numbers of values outside of block j that are <= each
    # of the block maxima in y
    # Row i contains the number of values that are outside block 1, ..., k_n
    # and <= block maximum i in y
    # The proportions are Fhat_{-j}(M_{ni}), i, j = 1, ..., k_n
    FhatjMni <- vapply(1:k_n, fun, numeric(k_n), x = nums_mat) / (m - b)
    # Column j enables us to calculate Yhatni(j) and/or Zhatni(j)
    UsumN <- -b * colMeans(log0const(FhatjMni, const))
    Usum <- b * (1 - colMeans(FhatjMni))
    return(list(Nhat = Fhaty, UsumN = UsumN, Usum = Usum))
  }
  fun_value <- list(numeric(k_n), numeric(k_n), numeric(k_n))
  # Use all sets of disjoint maxima to estimate sigmahat2_dj for sliding maxima
  which_vals <- 1:ncol(all_max$yd)
  temp <- vapply(which_vals, UsumN_fn, fun_value)
  Nhat <- do.call(cbind, temp[1, ])
  # BB2018
  Zhat <- b * (1 - Nhat)
  That <- colMeans(Zhat)
  Usum <- do.call(cbind, temp[3, ])
  Usum <-  t(k_n * That - (k_n - 1) * t(Usum))
  Bhat <- t(t(Zhat + Usum) - 2 * That)
  # N2015
  ZhatN <- -b * log(Nhat)
  ThatN <- colMeans(ZhatN)
  UsumN <- do.call(cbind, temp[2, ])
  UsumN <-  t(k_n * ThatN - (k_n - 1) * t(UsumN))
  # Bhat was mean-centred, but BhatN isn't (quite)
  BhatN <- t(t(ZhatN + UsumN) - 2 * ThatN)
  BhatN <- t(t(BhatN) - colMeans(BhatN))
  # Estimate sigma2_dj.
  # First calculate an estimate for each set of disjoint block maxima
  sigmahat2_dj <- apply(Bhat, 2, function(x) sum(x ^ 2) / length(x))
  sigmahat2_djN <- apply(BhatN, 2, function(x) sum(x ^ 2) / length(x))
  # Then, for sliding maxima, take the mean value over all sets of maxima
  sigmahat2_dj_for_sl <- sum(sigmahat2_dj) / length(sigmahat2_dj)
  sigmahat2_dj_for_slN <- sum(sigmahat2_djN) / length(sigmahat2_djN)
  sigma2dj_for_sl <- c(sigmahat2_dj_for_slN, sigmahat2_dj_for_sl)
  # For disjoint maxima pick either the first or last value, based on which_dj
  which_dj <- switch(which_dj, first = 1, last = length(sigmahat2_dj))
  sigma2dj <- c(sigmahat2_djN[which_dj], sigmahat2_dj[which_dj])
  #
  # Point estimates: disjoint maxima. Component i of ThatN (N) and That (BB)
  # contains (the reciprocal of) point estimates of theta based on set of
  # disjoint maxima i.  Use the mean of these estimates as an overall estimate.
  # Perform the Northrop (2015) `bias-adjustment' of Fhaty (Nhat here) if
  # requested and recalculate That and ThatN
  #
  if (bias_adjust == "N") {
    Nhat <- (m * Nhat - b) / (m - b)
    That <- colMeans(b * (1 - Nhat))
    ThatN <- colMeans(-b * log0const(Nhat, const))
  }
  theta_dj <- 1 / c(ThatN[which_dj], That[which_dj])
  names(theta_dj) <- names(sigma2dj) <- names(sigma2dj_for_sl) <-
    c("N2015", "BB2018")
  return(list(sigma2dj = sigma2dj, sigma2dj_for_sl = sigma2dj_for_sl,
              theta_dj = theta_dj,
              data_dj = cbind(N2015 = -b * log(Nhat[, which_dj]),
                              BB2018 = b * (1 - Nhat[, which_dj]))))
}
