#' Internal exdex functions
#'
#' Internal exdex functions
#' @details
#' These functions are not intended to be called by the user.
#' @name exdex-internal
#' @keywords internal
NULL

# =============================== To check spm() ============================ #

#' @keywords internal
#' @rdname exdex-internal
spm_check <- function(data, b, sliding = TRUE,
                      bias_adjust = c("BB3", "BB1", "N", "none"),
                      constrain = TRUE, varN = TRUE,
                      which_dj = c("last", "first")) {
  Call <- match.call(expand.dots = TRUE)
  # We don't check inputs here because this function is only used for
  # testing cases where I have ensured that the input are OK
  which_dj <- match.arg(which_dj)
  # Check that the value of b satisfies the inequality in Proposition 4.1
  # of Berghaus and Bucher (2018)
  k_n <- floor(length(data) / b)
  # Assume that the value of b is OK (for variances to be estimated)
  b_ok <- TRUE
  # If we want the last set of disjoint maxima then reverse the data vector
  if (which_dj == "last") {
    pass_data <- rev(data)
  } else {
    pass_data <- data
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
                                        dj_maxima = disjoint_maxima(data, b),
                                        which_dj = which_dj)
      } else {
        sigma2hat_dj <- spm_sigmahat_dj(data = data, b = b, dj_maxima = temp,
                                        which_dj = which_dj)
      }
      # If sliding = TRUE then estimate sigma2hat_sl
      # Otherwise use sigma2hat_dj
      indexN <- ifelse(varN, 2, 1)
      if (sliding) {
        sigma2hat_N <- sigma2hat_dj[indexN] - (3 - 4 * log(2)) / theta_N ^ 2
        sigma2hat_BB <- sigma2hat_dj[1] - (3 - 4 * log(2)) / theta_BB ^ 2
      } else {
        sigma2hat_N <- sigma2hat_dj[indexN + 2]
        sigma2hat_BB <- sigma2hat_dj[1 + 2]
      }
    }
    # Estimate the sampling variances of the estimators
    theta <- c(theta_N, theta_BB)
    if (b_ok) {
      vars <- theta ^ 4 * c(sigma2hat_N, sigma2hat_BB) / k_n
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
  res <- spm_estimates(data = pass_data)
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

#' @keywords internal
#' @rdname exdex-internal
spm_sigmahat_dj <- function(data, b, dj_maxima, check = FALSE, which_dj){
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
  loobN2015_fn <- function(x, block, xvec, y) {
    xx <- xvec[block != x]
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
    # For dj maxima add the first values.  This is the case because, if
    # which_dj = "last" we used rev(data) to reverse the data
    ests <- c(ests, temp[, 1])
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

# ========================== Functions used in spm() ======================== #

# log(x), but return a constant const for an x = 0

#' @keywords internal
#' @rdname exdex-internal
log0const_slow <- function(x, const) {
  ifelse(x == 0, const, log(x))
}

#' @keywords internal
#' @rdname exdex-internal
log0const <- function(x, const) {
  return(log(x + !x) + const * !x)
}

# =================== Empirical c.d.f. of x, evaluated at y ================= #

#' @keywords internal
#' @rdname exdex-internal
ecdf3 <- function(x, y) {
  return(vapply(y, function(y) mean(x <= y), 0))
}

#' @keywords internal
#' @rdname exdex-internal
ecdf2 <- function(x, y) {
  return(vapply(y, function(y) sum(x <= y) / length(x), 0))
}

#' @keywords internal
#' @rdname exdex-internal
ecdf1 <- function(x, y, lenx) {
  return(vapply(y, function(y) sum(x <= y) / lenx, 0))
}
