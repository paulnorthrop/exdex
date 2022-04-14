#' Methods for objects of class \code{"dgaps"}
#'
#' Methods for objects of class \code{c("dgaps", "exdex")} returned from
#' \code{\link{dgaps}}.
#' @param object and object of class \code{c("dgaps", "exdex")} returned from
#'   \code{\link{dgaps}}.
#' @param type A character scalar. Should the estimate of the variance be based
#'   on the observed information or the expected information?
#' @param x
#'   \code{print.dgaps}. An object of class \code{c("dgaps", "exdex")}, a
#'   result of a call to \code{\link{dgaps}}.
#'
#'   \code{print.summary.dgaps}. An object of class \code{"summary.dgaps"}, a
#'   result of a call to \code{\link{summary.dgaps}}.
#' @param se_type A character scalar. Should the estimate of the standard error
#'   be based on the observed information or the expected information?
#' @param digits
#'   \code{print.dgaps}. The argument \code{digits} to
#'   \code{\link{print.default}}.
#'
#'   \code{summary.dgaps}. An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... For \code{print.summary.dgaps}, additional arguments passed to
#'   \code{\link{print.default}}.
#' @return
#'   \code{coef.dgaps}. A numeric scalar: the estimate of the extremal index
#'   \eqn{\theta}.
#'
#'   \code{vcov.dgaps}. A \eqn{1 \times 1}{1 x 1} numeric matrix containing the
#'   estimated variance of the estimator.
#'
#'   \code{nobs.dgaps}. A numeric scalar: the number of inter-exceedance times
#'   used in the fit. If \code{x$inc_cens = TRUE} then this includes up to 2
#'   censored observations.
#'
#'   \code{logLik.dgaps}. An object of class \code{"logLik"}: a numeric scalar
#'   with value equal to the maximised log-likelihood.  The returned object
#'   also has attributes \code{nobs}, the numbers of \eqn{K}-gaps that
#'   contribute to the log-likelihood and \code{"df"}, which is equal to the
#'   number of total number of parameters estimated (1).
#'
#'   \code{print.dgaps}. The argument \code{x}, invisibly.
#'
#'   \code{summary.dgaps}. Returns a list containing the list element
#'   \code{object$call} and a numeric matrix \code{summary} giving the estimate
#'   of the extremal index \eqn{\theta} and the estimated standard error
#'   (Std. Error).
#'
#'   \code{print.summary.dgaps}. The argument \code{x}, invisibly.
#' @seealso \code{\link{dgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{confint.dgaps}} for confidence intervals for
#'   \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{dgaps}}.
#' @name dgaps_methods
NULL
## NULL

# =========================== coef.dgaps() ================================== #

#' Extract Model Coefficients from a \code{"dgaps"} object
#'
#' @rdname dgaps_methods
#' @export
coef.dgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  val <- object$theta
  names(val) <- "theta"
  return(val)
}

# =========================== vcov.dgaps() ================================== #

#' Calculate Variance-Covariance Matrix for a \code{"dgaps"} object
#'
#' @rdname dgaps_methods
#' @export
vcov.dgaps <- function(object, type = c("observed", "expected"), ...) {
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

# =========================== nobs.dgaps() ================================== #

#' Extract the Number of Observations from a \code{"dgaps"} object
#'
#' @rdname dgaps_methods
#' @export
nobs.dgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$ss$n_dgaps)
}

# ================================ logLik.dgaps ============================== #

#' Extract log-likelihood for objects of class \code{"dgaps"}
#'
#' @rdname dgaps_methods
#' @export
logLik.dgaps <- function(object, ...) {
  if (!inherits(object, "dgaps")) {
    stop("use only with \"dgaps\" objects")
  }
  val <- object$max_loglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- 1
  class(val) <- "logLik"
  return(val)
}

# ============================ print.dgaps() ================================== #

#' Print method for a \code{"dgaps"} object
#'
#' @rdname dgaps_methods
#' @export
print.dgaps <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
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

# =============================== summary.dgaps =============================== #

#' Summary method for a \code{"dgaps"} object
#'
#' @rdname dgaps_methods
#' @export
summary.dgaps <- function(object, se_type = c("observed", "expected"),
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
  class(res) <- "summary.dgaps"
  return(res)
}

# ============================ print.summary.dgaps ============================ #

#' Print method for objects of class \code{"summary.dgaps"}
#'
#' @rdname dgaps_methods
#' @export
print.summary.dgaps <- function(x, ...) {
  if (!inherits(x, "summary.dgaps")) {
    stop("use only with \"summary.dgaps\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}
