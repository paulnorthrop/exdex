# =========================== coef.dgaps() ================================== #

#' Extract Model Coefficients from a \code{"dgaps"} object
#'
#' \code{coef} method for class \code{c("dgaps", "exdex")}.
#'
#' @param object and object of class \code{c("kaps", "exdex")} returned from
#'   \code{\link{dgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the estimate of the extremal index \eqn{\theta}.
#' @export
coef.dgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$theta)
}

# =========================== vcov.dgaps() ================================== #

#' Calculate Variance-Covariance Matrix for a \code{"dgaps"} object
#'
#' \code{vcov} method for class \code{c("dgaps", "exdex")}.
#'
#' @param object and object of class \code{c("dgaps", "exdex")} returned from
#'   \code{\link{dgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A 1 by 1 numeric matrix containing the estimated variance of the
#'   estimator.
#' @export
vcov.dgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$se ^ 2)
}

# =========================== nobs.dgaps() ================================== #

#' Extract the Number of Observations from a \code{"dgaps"} object
#'
#' \code{nobs} method for class \code{c("dgaps", "exdex")}.
#'
#' @param object and object of class \code{c("dgaps", "exdex")} returned from
#'   \code{\link{dgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the number of inter-exceedance times used in the
#'   fit. If \code{x$inc_cens = TRUE} then this includes up to 2 censored
#'   observations.
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
#' \code{nobs} method for class \code{c("dgaps", "exdex")}.
#'
#' @param object an object of class \code{"dgaps"}, a result of a call to
#'   \code{\link{dgaps}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return An object of class \code{"logLik"}: a numeric scalar with
#' value equal to the maximised log-likelihood.  The returned object also has
#' attributes \code{nobs}, the numbers of \eqn{K}-gaps that contribute to the
#' log-likelihood and \code{"df"}, which is equal to the number of total number
#' of parameters estimated (1).
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
#' \code{print} method for class \code{c("dgaps", "exdex")}.
#'
#' @param x an object of class \code{c("dgaps", "exdex")}, a result of
#'   a call to \code{\link{dgaps}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Prints the original call to \code{\link{dgaps}}
#'   and the estimate of the extremal index \eqn{\theta}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{dgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{confint.dgaps}}: \code{confint} method for
#'   class \code{"dgaps"}.
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
#' \code{summary} method for class \code{"dgaps"}
#'
#' @param object an object of class "dgaps", a result of a call to
#'   \code{\link{dgaps}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... Additional arguments.  None are used in this function.
#' @return Returns a list containing the list element \code{object$call}
#'   and a numeric matrix \code{summary} giving the estimate of the extremal
#'   index \eqn{\theta} and the estimated standard error (Std. Error).
#' @seealso \code{\link{dgaps}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.dgaps}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{dgaps}}.
#' @export
summary.dgaps <- function(object, digits = max(3, getOption("digits") - 3L),
                        ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  res <- object["call"]
  theta <- signif(object$theta, digits = digits)
  se <- signif(object$se, digits = digits)
  res$matrix <- cbind(`Estimate` = theta, `Std. Error` = se)
  rownames(res$matrix) <- c("theta")
  class(res) <- "summary.dgaps"
  return(res)
}

# ============================ print.summary.dgaps ============================ #

#' Print method for objects of class \code{"summary.dgaps"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.dgaps"}.
#'
#' @param x An object of class \code{"summary.dgaps"}, a result of a call to
#'   \code{\link{summary.dgaps}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#' @return Prints the numeric matrix \code{x$summary} returned from
#' \code{\link{summary.dgaps}}.
#' @seealso \code{\link{dgaps}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.dgaps}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{dgaps}}.
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
