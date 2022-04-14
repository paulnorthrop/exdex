#' Methods for objects of class \code{"kgaps"}
#'
#' Methods for objects of class \code{c("kgaps", "exdex")} returned from
#' \code{\link{kgaps}}.
#' @param object and object of class \code{c("kgaps", "exdex")} returned from
#'   \code{\link{kgaps}}.
#' @param type A character scalar. Should the estimate of the variance be based
#'   on the observed information or the expected information?
#' @param x
#'   \code{print.kgaps}. An object of class \code{c("kgaps", "exdex")}, a
#'   result of a call to \code{\link{kgaps}}.
#'
#'   \code{print.summary.kgaps}. An object of class \code{"summary.kgaps"}, a
#'   result of a call to \code{\link{summary.kgaps}}.
#' @param se_type A character scalar. Should the estimate of the standard error
#'   be based on the observed information or the expected information?
#' @param digits
#'   \code{print.kgaps}. The argument \code{digits} to
#'   \code{\link{print.default}}.
#'
#'   \code{summary.kgaps}. An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... For \code{print.summary.kgaps}, additional arguments passed to
#'   \code{\link{print.default}}.
#' @return
#'   \code{coef.kgaps}. A numeric scalar: the estimate of the extremal index
#'   \eqn{\theta}.
#'
#'   \code{vcov.kgaps}. A \eqn{1 \times 1}{1 x 1} numeric matrix containing the
#'   estimated variance of the estimator.
#'
#'   \code{nobs.kgaps}. A numeric scalar: the number of inter-exceedance times
#'   used in the fit. If \code{x$inc_cens = TRUE} then this includes up to 2
#'   censored observations.
#'
#'   \code{logLik.kgaps}. An object of class \code{"logLik"}: a numeric scalar
#'   with value equal to the maximised log-likelihood.  The returned object
#'   also has attributes \code{nobs}, the numbers of \eqn{K}-gaps that
#'   contribute to the log-likelihood and \code{"df"}, which is equal to the
#'   number of total number of parameters estimated (1).
#'
#'   \code{print.kgaps}. The argument \code{x}, invisibly.
#'
#'   \code{summary.kgaps}. Returns a list containing the list element
#'   \code{object$call} and a numeric matrix \code{summary} giving the estimate
#'   of the extremal index \eqn{\theta} and the estimated standard error
#'   (Std. Error).
#'
#'   \code{print.summary.kgaps}. The argument \code{x}, invisibly.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{confint.kgaps}} for confidence intervals for
#'   \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{kgaps}}.
#' @name kgaps_methods
NULL
## NULL

# =========================== coef.kgaps() ================================== #

#' Extract Model Coefficients from a \code{"kgaps"} object
#'
#' @rdname kgaps_methods
#' @export
coef.kgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  val <- object$theta
  names(val) <- "theta"
  return(val)
}

# =========================== vcov.kgaps() ================================== #

#' Calculate Variance-Covariance Matrix for a \code{"kgaps"} object
#'
#' @rdname kgaps_methods
#' @export
vcov.kgaps <- function(object, type = c("observed", "expected"), ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  type <- match.arg(type)
  if (type == "observed") {
    vc <- object$se ^ 2
  } else {
    vc <- object$se_exp ^ 2
  }
  dim(vc) <- c(1, 1)
  dimnames(vc) <- list("theta", "theta")
  return(vc)
}

# =========================== nobs.kgaps() ================================== #

#' Extract the Number of Observations from a \code{"kgaps"} object
#'
#' @rdname kgaps_methods
#' @export
nobs.kgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$ss$n_kgaps)
}

# ================================ logLik.kgaps ============================== #

#' Extract log-likelihood for objects of class \code{"kgaps"}
#'
#' @rdname kgaps_methods
#' @export
logLik.kgaps <- function(object, ...) {
  if (!inherits(object, "kgaps")) {
    stop("use only with \"kgaps\" objects")
  }
  val <- object$max_loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- 1
  class(val) <- "logLik"
  return(val)
}

# ============================ print.kgaps() ================================== #

#' Print method for a \code{"kgaps"} object
#'
#' @rdname kgaps_methods
#' @export
print.kgaps <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimate of the extremal index theta:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  return(invisible(x))
}

# =============================== summary.kgaps =============================== #

#' Summary method for a \code{"kgaps"} object
#'
#' @rdname kgaps_methods
#' @export
summary.kgaps <- function(object, se_type = c("observed", "expected"),
                          digits = max(3, getOption("digits") - 3L), ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  se_type <- match.arg(se_type)
  res <- object["call"]
  theta <- signif(object$theta, digits = digits)
  if (se_type == "observed") {
    se <- signif(object$se, digits = digits)
  } else {
    se <- signif(object$se_exp, digits = digits)
  }
  res$matrix <- cbind(`Estimate` = theta, `Std. Error` = se)
  rownames(res$matrix) <- c("theta")
  class(res) <- "summary.kgaps"
  return(res)
}

# ============================ print.summary.kgaps ============================ #

#' Print method for objects of class \code{"summary.kgaps"}
#'
#' @rdname kgaps_methods
#' @export
print.summary.kgaps <- function(x, ...) {
  if (!inherits(x, "summary.kgaps")) {
    stop("use only with \"summary.kgaps\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}
