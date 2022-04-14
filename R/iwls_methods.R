#' Methods for objects of class \code{"iwls"}
#'
#' Methods for objects of class \code{c("iwls", "exdex")} returned from
#' \code{\link{iwls}}.
#' @param object and object of class \code{c("iwls", "exdex")} returned from
#'   \code{\link{iwls}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in these functions.
#' @return
#'   \code{coef.iwls}. A numeric scalar: the estimate of the extremal index
#'   \eqn{\theta}.
#'
#'   \code{nobs.iwls}. A numeric scalar: the number of inter-exceedance times
#'   used in the fit.
#'
#'   \code{print.iwls}. The argument \code{x}, invisibly.
#' @seealso \code{\link{iwls}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @section Examples:
#' See the examples in \code{\link{iwls}}.
#' @name iwls_methods
NULL
## NULL

# =========================== coef.iwls() ================================== #

#' Extract Model Coefficients from an \code{"iwls"} object
#'
#' @param object and object of class \code{c("iwls", "exdex")} returned from
#'   \code{\link{iwls}}.
#' @param ... Further arguments.  None are used.
#' @rdname iwls_methods
#' @export
coef.iwls <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$theta)
}

# =========================== nobs.iwls() ================================== #

#' Extract the Number of Observations from an \code{"iwls"} object
#'
#' @param object and object of class \code{c("iwls", "exdex")} returned from
#'   \code{\link{iwls}}.
#' @param ... Further arguments.  None are used.
#' @rdname iwls_methods
#' @export
nobs.iwls <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$n_gaps)
}

# ============================ print.iwls() ================================== #

#' Print method for an \code{"iwls"} object
#'
#' @param x an object of class \code{c("iwls", "exdex")}, a result of
#'   a call to \code{\link{iwls}}.
#' @details \code{print.iwls} prints the original call to \code{\link{iwls}}
#'   and the estimate of the extremal index \eqn{\theta}.
#' @rdname iwls_methods
#' @export
print.iwls <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Convergence (0 means success):", x$conv, "\n\n")
  cat("Estimate of the extremal index theta:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  return(invisible(x))
}
