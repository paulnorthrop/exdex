#' Methods for objects of class \code{"spm"}
#'
#' Methods for objects of class \code{c("spm", "exdex")} returned from
#' \code{\link{spm}}.
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param x
#'   \code{print.spm}. An object of class \code{c("spm", "exdex")}, a
#'   result of a call to \code{\link{spm}}.
#'
#'   \code{print.summary.spm}. An object of class \code{"summary.spm"}, a
#'   result of a call to \code{\link{summary.spm}}.
#' @param digits
#'   \code{print.spm}. The argument \code{digits} to
#'   \code{\link{print.default}}.
#'
#'   \code{summary.spm}. An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... For \code{print.summary.spm}, additional arguments passed to
#'   \code{\link{print.default}}.
#' @return
#'   \code{coef.spm}. A numeric scalar (or a vector of length 3 if
#'   \code{estimator = "all"}): the required estimate(s) of the extremal index
#'   \eqn{\theta}.
#'
#'   \code{vcov.spm}. A \eqn{1 \times 1}{1 x 1} numeric matrix if
#'   \code{estimator = "N2015"} or \code{"BB2018"} and a vector of length 3 if
#'   \code{estimator = "all"}, containing the estimated variance(s) of the
#'   estimator(s).
#'
#'   \code{nobs.spm}. A numeric scalar: the number of observations used in the
#'   fit.
#'
#'   \code{print.spm}. The argument \code{x}, invisibly.
#'
#'   \code{summary.spm}. Returns an object (a list) of class
#'   \code{"summary.spm"} containing the list element \code{object$call} and a
#'   numeric matrix \code{matrix} giving, for all three variants of the
#'   semiparametric estimator and both sliding and disjoint blocks,
#'   the (bias-adjusted) Estimate of the extremal index \eqn{\theta},
#'   the estimated standard error (Std. Error),
#'   and the bias adjustment (Bias adj.) applied to obtain the estimate, i.e.
#'   the value subtracted from the raw estimate.  If any of the
#'   (bias-adjusted) estimates are greater than 1 then a column
#'   containing the unconstrained estimates (Uncon. estimate) is added.
#'
#'   \code{print.summary.spm}. The argument \code{x}, invisibly.
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @section Examples:
#' See the examples in \code{\link{spm}}.
#' @name spm_methods
NULL
## NULL

# ============================ coef.spm() =================================== #

#' Extract Model Coefficients from an \code{"spm"} object
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the estimate
#'   of the extremal index \eqn{\theta} based on sliding maxima or on disjoint
#'   maxima.
#' @param estimator A character vector specifying which of the three variants
#'   of the semiparametric maxima estimator to use: \code{"N2015", "BB2018"}
#'   or \code{"BB2018b"}.  See \code{\link{spm}} for details.
#'   If \code{estimator = "all"} then all three estimates are returned.
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any estimates that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  Otherwise,
#'   estimates that are greater than 1 may be obtained.
#' @rdname spm_methods
#' @export
coef.spm <- function(object, maxima = c("sliding", "disjoint"),
                     estimator = "all", constrain = FALSE, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  if ("all" %in% estimator) {
    estimator <- c("N2015", "BB2018", "BB2018b")
  }
  if (!all(estimator %in% c("N2015", "BB2018", "BB2018b"))) {
    stop("estimator must be a subset of c(''N2015'', ''BB2018'', ''BB2018b'')")
  }
  maxima <- match.arg(maxima)
  if (constrain) {
    ests <- switch(maxima,
                   sliding = object$uncon_theta_sl,
                   disjoint = object$uncon_theta_dj)
  } else{
    ests <- switch(maxima,
                   sliding = object$theta_sl,
                   disjoint = object$theta_dj)
  }
  ests <- ests[estimator]
  return(ests)
}

# ============================ vcov.spm() =================================== #

#' Calculate Variance-Covariance Matrix for an \code{"spm"} object
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the
#'   estimated variance of the estimator of the extremal index \eqn{\theta}
#'   based on sliding maxima or on disjoint maxima.
#' @param estimator A character vector specifying which of the three variants
#'   of the semiparametric maxima estimator to use: \code{"N2015", "BB2018"}
#'   or \code{"BB2018b"}. See \code{\link{spm}} for details. If
#'   \code{estimator = "all"} then the estimated variances of all variants are
#'   returned.
#' @rdname spm_methods
#' @export
vcov.spm <- function(object, maxima = c("sliding", "disjoint"),
                     estimator = "all", ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  if ("all" %in% estimator) {
    estimator <- c("N2015", "BB2018", "BB2018b")
  }
  if (!all(estimator %in% c("N2015", "BB2018", "BB2018b"))) {
    stop("estimator must be a subset of c(''N2015'', ''BB2018'', ''BB2018b'')")
  }
  maxima <- match.arg(maxima)
  se <- switch(maxima,
               sliding = object$se_sl,
               disjoint = object$se_dj)
  vc <- se[estimator] ^ 2
  if (length(vc) == 1) {
    dim(vc) <- c(1, 1)
    dimnames(vc) <- list(estimator, estimator)
  }
  return(vc)
}

# ============================ nobs.spm() =================================== #

#' Extract the Number of Observations from an \code{"spm"} object
#'
#' @param object and object of class \code{c("spm", "exdex")} returned from
#'   \code{\link{spm}}.
#' @param maxima A character scalar specifying whether to return the
#'   number of observed sliding maxima or disjoint maxima.
#' @rdname spm_methods
#' @export
nobs.spm <- function(object, maxima = c("sliding", "disjoint"), ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  maxima <- match.arg(maxima)
  n <- switch(maxima,
              sliding = nrow(object$data_sl),
              disjoint = nrow(object$data_dj))
  return(n)
}

# ============================ print.spm() ================================== #

#' Print method for an \code{"spm"} object
#'
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @details \code{print.spm} prints the original call to \code{\link{spm}}
#'   and the estimates of the extremal index \eqn{\theta}, based on all three
#'   variants of the semiparametric maxima estimator and both sliding
#'   and disjoint blocks.
#' @rdname spm_methods
#' @export
print.spm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimates of the extremal index theta:\n")
  coef_sl <- coef(x, maxima = "sliding")
  coef_dj <- coef(x, maxima = "disjoint")
  coefs <- rbind(sliding = coef_sl, disjoint = coef_dj)
  print.default(format(coefs, digits = digits), print.gap = 2L,
                quote = FALSE)
  # Add warning that se_sl is missing and the consequences
  if (any(is.na(x$se_sl))) {
    which_na <- names(x$se_sl)[is.na(x$se_sl)]
    cat("\nStd. Errors missing for estimator(s):", which_na)
    if (x$bias_adjust == "BB3") {
      cat("\nBias-adjustment changed from BB3 to BB1 for estimator(s):",
          which_na)
    }
  }
  return(invisible(x))
}

# =============================== summary.spm =============================== #

#' Summary method for an \code{"spm"} object
#'
#' @param object an object of class \code{"spm"}, a result of a call to
#'   \code{\link{spm}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @rdname spm_methods
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
                            "BB2018b, sliding",
                            "N2015, disjoint", "BB2018, disjoint",
                            "BB2018b, disjoint")
  # Add warning that se_sl is missing and the consequences
  if (any(is.na(object$se_sl))) {
    which_na <- names(object$se_sl)[is.na(object$se_sl)]
    if (object$bias_adjust == "BB3") {
      res$warning <- "Bias-adjustment changed from BB3 to BB1 for estimator(s):"
      for (i in 1:length(which_na)) {
        res$warning <- c(res$warning, which_na[i])
      }
    }
  }
  class(res) <- "summary.spm"
  return(res)
}

# ============================ print.summary.spm ============================ #

#' Print method for objects of class \code{"summary.spm"}
#'
#' @rdname spm_methods
#' @export
print.summary.spm <- function(x, ...) {
  if (!inherits(x, "summary.spm")) {
    stop("use only with \"summary.spm\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  if (!is.null(x$warning)) {
    cat("\n")
    cat(x$warning)
  }
  invisible(x)
}
