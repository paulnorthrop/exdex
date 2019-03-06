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
#' @param sliding A logical scalar.  If \code{sliding = TRUE} then the
#'   estimates are based on sliding block maxima.  Otherwise, the estimates
#'   are based on disjoint block maxima.
#' @param bias_adjust A character scalar.  Is bias-adjustment of the
#'   raw estimate of \eqn{\theta} performed using the bias-reduced
#'   estimator (\code{bias_adjust = "BB3"}), derived in Section 5 of
#'   Berghaus and Bucher (2018); or a simpler version of this
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
#'   \eqn{Z}-data (\code{varN = TRUE}).
#'
#'   A condition imposed in Proposition 4.1 of Berghaus and Bucher (2018)
#'   means that \code{b} must be no smaller than \eqn{k^{1/2}} and no larger
#'   than \eqn{k^2}, where \eqn{k} is \code{floor(length(data) / b)},
#'   i.e. \eqn{k} is the number of complete blocks.  If this is not the case
#'   then a warning will be given and
#'     \itemize{
#'       \item{estimated standard errors will be missing from the returned
#'             object,}
#'       \item{if \code{bias_adjust == "BB3"} then this bias-adjustment
#'             based on \code{bias_adjust == "BB1"} will be performed instead,
#'             because the former relies on the estimated variances of the
#'             estimators.}
#'     }
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta} : } {A vector containing the estimates \eqn{\theta}
#'       resulting from the two variants of the semiparametric estimator,
#'       labelled N2015 for Northrop (2015) and BB2018 for
#'       Berghaus and Bucher (2018).}
#'     \item {\code{se} : } {The respective estimated standard errors.}
#'     \item {\code{uncontrained_theta} : } {The estimates of \eqn{\theta}
#'       without the constraint that they lie in (0, 1])}
#'     \item {\code{N2015_data} : } {The values of the \eqn{Y}-data.}
#'     \item {\code{BB2018_data} : } {The values of the \eqn{Z}-data.}
#'     \item {\code{bias_val} : } {The respective values of the
#'       bias-adjustment applied to the raw estimates.  This is only relevant
#'       if \code{bias_adjust} is "BB3" or "BB1".  Otherwise, \code{bias_val}
#'       is \code{NA}.}
#'     \item {\code{bias_adjust} : } {The input value of \code{bias_adjust}.}
#'     \item {\code{b} : } {The input value of \code{b}.}
#'     \item {\code{sliding} : } {The input value of \code{sliding}.}
#'     \item {\code{call} : } {The call to \code{spm}.}
#'   }
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
#' spm(-as.vector(sp500[2:6550]), 250)
#'
#' spm(newlyn, 20)
#' spm(newlyn, 20, sliding = FALSE)
#' @export
spm <- function(data, b, sliding = TRUE,
                bias_adjust = c("BB3", "BB1", "N", "none"), constrain = TRUE,
                varN = TRUE) {
  Call <- match.call(expand.dots = TRUE)
  # Check inputs
  if (missing(data) || length(data) == 0L || mode(data) != "numeric") {
    stop("'data' must be a non-empty numeric vector")
  }
  if (any(!is.finite(data))) {
    stop("'data' contains missing or infinite values")
  }
  if (is.matrix(data)) {
    stop("'data' must be a vector")
  }
  data <- as.vector(data)
  if (!is.numeric(b) || length(b) != 1) {
    stop("'b' must be a numeric scalar (specifically, a positive integer)")
  }
  if (b < 1) {
    stop("'b' must be no smaller than 1")
  }
  if (!is.logical(sliding) || length(sliding) != 1) {
    stop("'sliding' must be a logical scalar")
  }
  bias_adjust <- match.arg(bias_adjust)
  if (!is.logical(constrain) || length(constrain) != 1) {
    stop("'sliding' must be a logical scalar")
  }
  if (!is.logical(varN) || length(varN) != 1) {
    stop("'varN' must be a logical scalar")
  }
  # Check that the value of b satisfies the inequality in Proposition 4.1
  # of Berghaus and Bucher (2018)
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
    }
    warning("\n", warn1, "\n", warn2, "\n", warn3)
  }
  # A function that returns N2015 and BB2018 estimates of theta
  # (constrained to (0, 1] if constrain = TRUE))
  spm_estimates <- function(data) {
    # Calculate the block maxima
    if (sliding) {
      temp <- sliding_maxima(data, b)
    } else{
      temp <- disjoint_maxima(data, b)
    }
    # Extract x ~ F (only xs contributing to y are included) and y ~ G
    x <- temp$x
    y <- temp$y
    # Empirical c.d.f. of raw (`daily') values
    Fhat <- stats::ecdf(x)
    # Evaluate Fx at y
    Fhaty <- Fhat(y)
    if (bias_adjust == "N") {
      # `Bias-adjust' the empirical c.d.f. of Y based on the Xs: by subtracting
      # b in numerator and denominator we remove Xs that are in the same block
      # as Y. We use m-b in the denominator rather than the m-b+1 in
      # Northrop (2015)
      m <- length(x)
      Fhaty <- (m * Fhaty - b) / (m - b)
      # In the unlikely event that an element of Fhaty is equal to zero,
      # i.e. a block maximum is less than all the data from outside that
      # block, we force Fhaty to be positive
      Fhaty[Fhaty == 0] <- 1 / (m - b + length(y))
    }
    # Calculate the estimate of theta:
    # theta_N: Northrop (2015) and theta_BB: Berghaus and Bucher (2018)
    theta_N <- -1 / mean(b * log(Fhaty))
    theta_BB <- 1 / (b * mean(1 - Fhaty))
    # Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
    # We require the disjoint maxima to do this.  If sliding = TRUE then
    # pass these to spm_sigmahat_dj using the dj_maxima argument
    # Only do this is b_ok = TRUE
    if (b_ok) {
      if (sliding) {
        sigma2hat_dj <- spm_sigmahat_dj(data = data, b = b,
                                        dj_maxima = disjoint_maxima(data, b))
      } else {
        sigma2hat_dj <- spm_sigmahat_dj(data = data, b = b, dj_maxima = temp)
      }
      # If sliding = TRUE then estimate sigma2hat_sl
      # Otherwise use sigma2hat_dj
      indexN <- ifelse(varN, 2, 1)
      if (sliding) {
        sigma2hat_N <- sigma2hat_dj[indexN] - (3 - 4 * log(2)) / theta_N ^ 2
        sigma2hat_BB <- sigma2hat_dj[1] - (3 - 4 * log(2)) / theta_BB ^ 2
      } else {
        sigma2hat_N <- sigma2hat_dj[indexN]
        sigma2hat_BB <- sigma2hat_dj[1]
      }
    }
    # Estimate the sampling variances of the estimators
    theta <- c(theta_N, theta_BB)
    if (b_ok) {
      vars <- theta ^ 4 * c(sigma2hat_N, sigma2hat_BB) / k_n
    } else {
      vars <- c(NA, NA)
    }
    # Perform BB2018 bias-adjustment if required
    bias_N <- bias_BB <- NA
    if (bias_adjust == "BB3") {
      bias_N <- theta_N / k_n + theta_N ^ 3 * sigma2hat_N / k_n
      theta_N <- theta_N * (1 - 1 / k_n) - theta_N ^ 3 * sigma2hat_N / k_n
      bias_BB <- theta_BB / k_n + theta_BB ^ 3 * sigma2hat_BB / k_n
      theta_BB <- theta_BB * (1 - 1 / k_n) - theta_BB ^ 3 * sigma2hat_BB / k_n
      theta <- c(theta_N, theta_BB)
    } else if (bias_adjust == "BB1") {
      bias_N <- theta_N / k_n
      theta_N <- theta_N * (1 - 1 / k_n)
      bias_BB <- theta_BB / k_n
      theta_BB <- theta_BB * (1 - 1 / k_n)
      theta <- c(theta_N, theta_BB)
    }
    # Save the unconstrained estimates, so that they can be returned
    unconstrained_theta <- theta
    # Constrain to (0, 1] if required
    if (constrain) {
      theta <- pmin(theta, 1)
    }
    se <- sqrt(vars)
    return(list(theta = theta, se = se,
                unconstrained_theta = unconstrained_theta,
                N2015_data = -b * log(Fhaty),
                BB2018_data = b * (1 - Fhaty),
                bias_val = c(bias_N, bias_BB)))
  }
  # End of function spm_estimates() ----------
  #
  # Find the point estimate of theta and the raw data that contribute to it
  res <- spm_estimates(data = data)
  estimator_names <- c("N2015", "BB2018")
  names(res$theta) <- estimator_names
  names(res$se) <- estimator_names
  names(res$unconstrained_theta) <- estimator_names
  names(res$bias_val) <- estimator_names
  res$bias_adjust <- bias_adjust
  res$b <- b
  res$sliding <- sliding
  res$call <- Call
  class(res) <- c("exdex", "spm")
  return(res)
}

#' Estimates \eqn{\sigma_{dj}^2}
#'
#' Estimates the value of \eqn{\sigma_{dj}^2}, a variance involved in
#' Sections 3, 4 and 5 of Berghaus and Bucher (2018).
#' These estimates are required in \code{\link{spm}}, which estimates
#' the extremal index using the methods described in Northrop (2015) and
#' Berghaus and Bucher (2018), in order to perform bias-adjustment and
#' estimation of uncertainty.
#'
#' @param data A numeric vector of raw data.
#' @param b A numeric scalar.  The block size.
#' @param dj_maxima A list returned from a call
#'   \code{disjoint_maxima(data, b)} to \code{\link{disjoint_maxima}}.
#'   If this isn't supplied then this call is made inside
#'   \code{spm_sigmahat_dj}.
#' @param check A logical scalar.  Used only to check this function during
#'   development.  If \code{check = TRUE} then a matrix with 5 columns is
#'   returned, each columns containing the result of 5 different ways to
#'   calculate the middle term of Bhatnj in the second line of the equation
#'   displayed midpage on page 2319 of Berghaus and Bucher (2018).
#' @details The computations follow page 2319 of Berghaus and Bucher (2018),
#'   which relates to the estimator given in equation (1.4) on page 2309
#'   of that paper. Also calculated is an analogous estimate based on the,
#'   asymptotically equivalent, estimator proposed in Northrop (2015),
#'   given in equation (1.3) on page 2309 of Berghaus and Bucher (2018).
#' @return A vector of length two.  The first component is the estimate
#'  given by the first part of the penultimate displayed equation on page
#'  2319 of Berghaus and Bucher (2018). The second component is an analogous
#'  estimate based on the estimator proposed in Northrop (2015).
#' @seealso \code{\link{spm}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @examples
#' spm_sigmahat_dj(newlyn, 20)
#' @export
spm_sigmahat_dj <- function(data, b, dj_maxima, check = FALSE){
#  if (missing(dj_maxima)){
#    dj_maxima <- disjoint_maxima(data, b)
#  }
  all_dj_maxima <- all_disjoint_maxima(data, b)
  # The number of blocks and the number of raw observations that contribute
  k_n <- nrow(all_dj_maxima$y)
  m <- nrow(all_dj_maxima$x)
  lenxx <- m - b
  const <- -log(m - b + k_n)
  # block indicator
  block <- rep(1:k_n, each = b)
  # Set up some functions
#  BB2018_fn <- function(x, y) {
#    return(mean(1 - ecdf2(x, y)))
#  }
  BB2018_fn <- function(x, y) {
    return(1 - sum(ecdf2(x, y)) / k_n)
  }
  #  loobBB2018_fn <- function(x, block, xvec, y) {
#    xx <- xvec[block != x]
#    return(mean(1 - ecdf2(xx, y)))
#  }
  loobBB2018_fn <- function(x, block, xvec, y) {
    xx <- xvec[block != x]
    return(1 - sum(ecdf2(xx, y)) / k_n)
  }
  # In the unlikely event that an element of Fhaty is equal to zero,
  # i.e. a block maximum is less than all the data from outside that
  # block, we force Fhaty to be positive
#  loobN2015_fn <- function(x, block, xvec, y) {
#    xx <- xvec[block != x]
#    return(mean(-log0const(ecdf2(xx, y), const)))
#  }
#  loobN2015_fn <- function(x, block, xvec, y) {
#    xx <- xvec[block != x]
#    return(sum(-log0const(ecdf2(xx, y), const)) / k_n)
#  }
  loobN2015_fn <- function(x, block, xvec, y) {
    xx <- xvec[block != x]
#    return(sum(-log0const(ecdf1(xx, y, lenxx), const)) / k_n)
#    return(sum(-log0const(ecdf2(xx, y), const))) # plus divide UsumN by k_n below
    return(sum(-log0const(ecdf2(xx, y), const)) / k_n)
  }
  ests_fn <- function(i) {
    # y: disjoint block maxima, x: (only) the raw values that contribute to y
    x <- all_dj_maxima$x[, i]
    y <- all_dj_maxima$y[, i]
    # The ecdf of the data evaluated at the block maxima
    Nhat <- ecdf2(x, y)
    # BB2018
    Zhat <- b * (1 - Nhat)
    That <- mean(Zhat)
    # N2015
    ZhatN <- -b * log(Nhat)
    ThatN <- mean(ZhatN)
    Usum <- b * tapply(x, block, BB2018_fn, y = y)
    UsumN <- b * vapply(1:k_n, loobN2015_fn, 0, block = block, xvec = x, y = y)
#    print(UsumN)
#    my_temp <- my_fun(newlyn, b)
#    UsumN <- my_temp$UsumN
#    print(UsumN)
#    print(Usum)
#    Usum <- my_temp$Usum
#    print(k_n * That - (k_n - 1) * Usum)
    UsumN <-  k_n * ThatN - (k_n - 1) * UsumN
    #
    Bhat <- Zhat + Usum - 2 * That
    BhatN <- ZhatN + UsumN - 2 * ThatN
    # Bhat is mean-centred, but BhatN isn't (quite)
    BhatN <- BhatN - mean(BhatN)
    sigmahat2_dj <- mean(Bhat ^ 2)
    sigmahat2_djN <- mean(BhatN ^ 2)
    return(c(sigmahat2_dj, sigmahat2_djN))
  }
  if (!check) {
    temp <- vapply(1:ncol(all_dj_maxima$y), ests_fn, c(0, 0))
    ests <- rowMeans(temp)
    return(ests)
  }
  # 4 (effectively 3) other ways to calculate Usum
  # (Usum1 is commented out because it incurs rounding error from
  # Nhat -> Zhat -> Nhat that can be non-negligible)
  x <- all_dj_maxima$x[, 1]
  y <- all_dj_maxima$y[, 1]
#  Fhat <- stats::ecdf(x)
#  Nhat <- Fhat(y)
  Nhat <- ecdf2(x, y)
  Zhat <- b * (1 - Nhat)
  That <- mean(Zhat)
  Usum <- b * tapply(x, block, BB2018_fn, y = y)
#  Uhats <- Fhat(x)
  Uhats <- ecdf2(x, x)
  # Usum1 <- colSums(vapply(Uhats, function(x) x > 1 - Zhat / b, rep(0, k_n)))
  Usum2 <- colSums(vapply(Uhats, function(x) x > Nhat, rep(0, k_n)))
  Usum3 <- colSums(vapply(x, function(x) x > y, rep(0, k_n)))
  # Usum4 is the analogous calculation to UsumN
  Usum4 <- b * vapply(1:k_n, loobBB2018_fn, 0, block = block, xvec = x, y = y)
  Usum4 <-  k_n * That - (k_n - 1) * Usum4
  # Aggregate the first 3
  # Usum1 <- tapply(Usum1, block, sum) / k_n
  Usum2 <- tapply(Usum2, block, sum) / k_n
  Usum3 <- tapply(Usum3, block, sum) / k_n
  return(cbind(Usum, Usum2, Usum3, Usum4))
}
