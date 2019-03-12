# ============================ coef.spm() =================================== #

#' Extract Model Coefficients
#'
#' \code{coef} method for class \code{c("spm", "exdex")}.
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the estimate
#'   of the extremal index \eqn{\theta} based on sliding maxima or on disjoint
#'   maxima.
#' @param estimator A character scalar specifying which variant of the
#'   semiparametic maxima estimator to use: Northrop (2015) or
#'   Berghaus and Bucher (2018).  If \code{estimator = "both"} then
#'   estimates for both variants are returned.
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any estimates that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  Otherwise,
#'   estimates that are greater than 1 may be obtained.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar (or a vector of length 2 if
#'   \code{estimator = "both"}): the required estimate(s) of the extremal index
#'   \eqn{\theta}.
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @export
coef.spm <- function(object, maxima = c("sliding", "disjoint"),
                     estimator = c("N2015", "BB2018", "both"),
                     constrain = FALSE, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  maxima <- match.arg(maxima)
  estimator <- match.arg(estimator)
  if (constrain) {
    ests <- switch(maxima,
                   sliding = object$uncon_theta_sl,
                   disjoint = object$uncon_theta_dj)
  } else{
    ests <- switch(maxima,
                   sliding = object$theta_sl,
                   disjoint = object$theta_dj)
  }
  if (estimator == "both") {
    estimator <- c("N2015", "BB2018")
    ests <- ests[estimator]
  } else {
    ests <- ests[estimator]
    names(ests) <- "theta"
  }
  return(ests)
}

# ============================ vcov.spm() =================================== #

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' \code{vcov} method for class \code{c("spm", "exdex")}.
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the
#'   estimated variance of the estimator of the extremal index \eqn{\theta}
#'   based on sliding maxima or on disjoint maxima.
#' @param estimator A character scalar specifying which variant of the
#'   semiparametic maxima estimator to use: Northrop (2015) or
#'   Berghaus and Bucher (2018).  If \code{estimator = "both"} then the
#'   estimated variances of both variants are returned.
#' @param ... Further arguments.  None are used.
#' @return A 1 by 1 numeric matrix if \code{estimator = "N2015"} or
#'   \code{"BB2018"} and a vector of length 2 if \code{estimator = "both"},
#'   containing the estimated variance(s) of the estimator(s).
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @export
vcov.spm <- function(object, maxima = c("sliding", "disjoint"),
                     estimator = c("N2015", "BB2018", "both"), ...) {
  maxima <- match.arg(maxima)
  estimator <- match.arg(estimator)
  se <- switch(maxima,
               sliding = object$se_sl,
               disjoint = object$se_dj)
  if (estimator == "both") {
    estimator <- c("N2015", "BB2018")
    vcov <- se[estimator] ^ 2
  } else {
    vcov <- as.matrix(se[estimator] ^ 2)
    dimnames(vcov) <- list("theta", "theta")
  }
  return(vcov)
}

# ============================ nobs.spm() =================================== #

#' Extract the Number of Observations from a Fit
#'
#' \code{nobs} method for class \code{c("spm", "exdex")}.
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the
#'   number of observed sliding maxima or disjoint maxima.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the number of observations used in the fit.
#' @export
nobs.spm <- function(object, maxima = c("sliding", "disjoint"), ...) {
  maxima <- match.arg(maxima)
  n <- switch(maxima,
              sliding = nrow(object$data_sl),
              disjoint = nrow(object$data_dj))
  return(n)
}

# ============================ print.spm() ================================== #

#' Print method for an spm object
#'
#' \code{print} method for class \code{c("spm", "exdex")}.
#'
#' @param x an object of class \code{c("spm", "exdex")}, a result of
#'   a call to \code{\link{spm}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Prints the original call to \code{\link{spm}}
#'   and the estimates of the extremal index \eqn{\theta}, based on both
#'   variants of the semiparametric maxima estimator and both sliding
#'   and disjoint blocks.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.spm}}: \code{confint} method for
#'   class \code{"spm"}.
#' @export
print.spm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimates of the extremal index theta:\n")
  coef_sl <- coef(res, maxima = "sliding", estimator = "both")
  coef_dj <- coef(res, maxima = "disjoint", estimator = "both")
  coefs <- rbind(sliding = coef_sl, disjoint = coef_dj)
  print.default(format(coefs, digits = digits), print.gap = 2L,
                quote = FALSE)
  return(invisible(x))
}

# =============================== summary.spm =============================== #

#' Summary method for an spm object
#'
#' \code{summary} method for class \code{"spm"}
#'
#' @param object an object of class "spm", a result of a call to
#'   \code{\link{spm}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... Additional arguments.  None are used in this function.
#' @return Returns a list containing the list element \code{object$call}
#'   and a numeric matrix \code{summary} giving, for both variants of the
#'   semiparametric estimator and both sliding and disjoint blocks,
#'   the (bias-adjusted) Estimate, the estimated standard error (Std. Error),
#'   and the bias adjustment (Bias adj.) applied to obtain the estimate, i.e.
#'   the value subtracted from the raw estimate.  If any of the
#'   (bias-adjusted) estimates are greater than 1 then a column
#'   containing the unconstrained estimates (Uncon. estimate) is added.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.spm}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{spm}}.
#' @export
summary.spm <- function(object, digits = max(3, getOption("digits") - 3L),
                        ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  res <- object["call"]
  theta <- signif(c(object$theta_sl, object$theta_dj), digits = digits)
  se <- signif(c(object$se_sl, object$se_dj), digits = digits)
  bias_adj <- signif(c(object$bias_sl, object$bias_dj), digits = digits)
  res$matrix <- cbind(`Estimate` = theta, `Std. Error` = se,
                      `Bias adj.` = bias_adj)
  uncon <- c(object$uncon_theta_sl, object$uncon_theta_dj)
  if (any(uncon > 1)) {
    res$matrix <- cbind(res$matrix, `Uncon. estimate` =
                          signif(uncon, digits = digits))
  }
  rownames(res$matrix) <- c("N2015, sliding", "BB2018, sliding",
                            "N2015, disjoint", "BB2018, disjoint")
  class(res) <- "summary.spm"
  return(res)
}

# ============================ print.summary.spm ============================ #

#' Print method for objects of class \code{"summary.spm"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.spm"}.
#'
#' @param x An object of class "summary.pm", a result of a call to
#'   \code{\link{summary.spm}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#' @return Prints the numeric matrix \code{x$summary} returned from
#' \code{\link{summary.spm}}.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.spm}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{spm}}.
#' @export
print.summary.spm <- function(x, ...) {
  if (!inherits(x, "summary.spm")) {
    stop("use only with \"summary.spm\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$summary, ...)
  invisible(x)
}
