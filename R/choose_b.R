#' Block length diagnostic plot
#'
#' Description
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @examples
#' \dontrun{
#' # Plot like the top left of Northrop (2015)
#' # Remove the last 14 values because 2880 has lots of factors
#' b_vals <- c(2,3,4,5,6,8,9,10,12,15,16,18,20,24,30,32,36,40,45,48,54,60)
#' res <- choose_b(newlyn[1:2880], b_vals)
#' # Some b are too small for the sampling variance of the sliding blocks
#' # estimator to be estimated
#' plot(res)
#' plot(res, estimator = "BB2018")
#' plot(res, maxima = "disjoint")
#'
#' # S&P 500 Composite Berghaus and Bucher (2018)
#' b_vals <- c(10, seq(from = 25, to = 350, by = 25), 357)
#' res500 <- choose_b(sp500, c(10, b_vals, 357))
#' plot(res500, ylim = c(0, 1))
#' plot(res500, estimator = "BB2018", ylim = c(0, 1))
#' }
#' @export
choose_b <- function(data, b, bias_adjust = c("BB3", "BB1", "N", "none"),
                     constrain = TRUE, varN = TRUE, level = 0.95,
                     conf_scale = c("theta", "log"),
                     interval_type = c("sym", "lik"),
                     type = c("vertical", "cholesky", "spectral", "none")) {
  Call <- match.call(expand.dots = TRUE)
  # All other inputs are checked in the calls to spm() and confint.spm()
  interval_type <- match.arg(interval_type)
  which_ests <- c("N2015", "BB2018")
  which_cis <- paste0(which_ests, interval_type)
  # Objects in which to store the estimates and confidence intervals
  n_b <- length(b)
  empty_matrix <- matrix(0, nrow = n_b, ncol = length(which_ests))
  colnames(empty_matrix) <- which_ests
  rownames(empty_matrix) <- b
  theta_sl <- lower_sl <- upper_sl <- theta_dj <- lower_dj <- upper_dj <-
    empty_matrix
  # Loop over block lengths
  for (i in 1:n_b) {
    ests <- spm(data = data, b = b[i], bias_adjust = bias_adjust,
                constrain = constrain, varN = varN)
    theta_sl[i, ] <- ests$theta_sl[which_ests]
    theta_dj[i, ] <- ests$theta_dj[which_ests]
    # Avoid chandwich::conf_intervals()'s profiling messages
    cis <- suppressMessages(confint(ests, maxima = "sliding",
                                    interval_type = interval_type,
                                    constrain = constrain,
                                    conf_scale = conf_scale,
                                    bias_adjust = TRUE, type = type))
    lower_sl[i, ] <- cis$cis[which_cis, 1]
    upper_sl[i, ] <- cis$cis[which_cis, 2]
    cis <- suppressMessages(confint(ests, maxima = "disjoint",
                                    interval_type = interval_type,
                                    constrain = constrain,
                                    conf_scale = conf_scale,
                                    bias_adjust = TRUE, type = type))
    lower_dj[i, ] <- cis$cis[which_cis, 1]
    upper_dj[i, ] <- cis$cis[which_cis, 2]
  }
  res <- list(theta_sl = theta_sl, lower_sl = lower_sl, upper_sl = upper_sl,
              theta_dj = theta_dj, lower_dj = lower_dj, upper_dj = upper_dj,
              b = b)
  res$call <- Call
  class(res) <- c("choose_b", "exdex")
  return(res)
}


# =========================== plot.choose_b ===========================

#' Plot diagnostics for an exdex object
#'
#' \code{plot} method for objects inheriting from class \code{"choose_b"},
#' returned from \code{\link{choose_b}}
#'
#' @param x an object of class \code{c("choose_b", "exdex")}, a result of a
#'   call to \code{\link{choose_b}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}} and/or \code{\link[graphics]{axis}}.
#' @param estimator Choice of estimator: \code{"N2015"} for Northrop (2015),
#'   \code{"BB2018"} for Berghaus and Bucher (2018).
#'   See \code{\link{spm}} for details.
#' @param maxima Should the estimator be based on sliding or disjoint maxima?
#' @details Produces a simple diagnostic plot to aid the choice of block
#'   length \eqn{b} based on the object returned from \code{\link{choose_b}}.
#'   Estimates of \eqn{b} and approximate \code{conf}\% confidence intervals
#'   are plotted against the value of \eqn{b} used to produce each estimate.
#'   The type of confidence interval is determined by the arguments
#'   \code{interval_type}, \code{conf_scale} and \code{type} provided in the
#'   call to \code{\link{choose_b}}.
#' @return In addition to producing the plot a list of the arguments used
#'   by \code{\link[graphics]{matplot}}, \code{\link[graphics]{axis}} is
#'   returned (invisibly).
#' @seealso \code{\link{choose_b}}.
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @examples
#' b_vals <- seq(from = 25, to = 350, by = 25)
#' res <- choose_b(sp500, b_vals)
#' plot(res)
#' @export
plot.choose_b <- function(x, y, ..., estimator = c("N2015", "BB2018"),
                          maxima = c("sliding", "disjoint")) {
  if (!inherits(x, "choose_b")) {
    stop("use only with \"choose_b\" objects")
  }
  maxima <- match.arg(maxima)
  estimator <- match.arg(estimator)
  if (maxima == "sliding") {
    y_data <- cbind(x$lower_sl[, estimator], x$theta_sl[, estimator],
                    x$upper_sl[, estimator])
  } else {
    y_data <- cbind(x$lower_dj[, estimator], x$theta_dj[, estimator],
                    x$upper_dj[, estimator])
  }
  if (maxima == "sliding") {
    if (estimator == "N2015") {
      y_lab <- expression(hat(theta)[sl]^{N})
    } else {
      y_lab <- expression(hat(theta)[sl]^{BB})
    }
  } else {
    if (estimator == "N2015") {
      y_lab <- expression(hat(theta)[dj]^{N})
    } else {
      y_lab <- expression(hat(theta)[dj]^{BB})
    }
  }
  x_data <- x$b
  x_lab <- "block size, b"
  xy_args <- list(x = x_data, y = y_data)
  # Look for user-supplied arguments to matplot.
  user_args <- list(...)
  m_cond <- names(user_args) %in% methods::formalArgs(graphics::matplot)
  a_cond <- names(user_args) %in% methods::formalArgs(graphics::axis)
  s_cond <- names(user_args) %in% methods::formalArgs(graphics::segments)
  matplot_args <- user_args[!a_cond | m_cond]
  axis_args <- user_args[!m_cond | a_cond]
  segments_args <- user_args[s_cond]
  axis_args$col <- 1
  if (is.null(matplot_args$xlab)) {
    matplot_args$xlab <- x_lab
  }
  if (is.null(matplot_args$ylab)) {
    matplot_args$ylab <- y_lab
  }
  if (is.null(matplot_args$type)) {
    matplot_args$type <- "l"
  }
  if (is.null(matplot_args$col)) {
    matplot_args$col <- 1
  }
  if (is.null(matplot_args$lty)) {
    matplot_args$lty <- c(1, 1, 1)
  }
  if (is.null(matplot_args$mgp)) {
    matplot_args$mgp <- c(1.75, 0.6, 0)
  }
  if (is.null(axis_args$mgp)) {
    axis_args$mgp <- c(1.75, 0.6, 0)
  }
  all_args <- c(xy_args, matplot_args)
  do.call(graphics::matplot, c(all_args, axes = FALSE))
  axis_args$side <- 2
  do.call(graphics::axis, axis_args)
  axis_args$side <- 1
  axis_args$at <- pretty(x_data)
  do.call(graphics::axis, axis_args)
  if (!is.null(axis_args$lwd)) {
    graphics::box(lwd = axis_args$lwd)
  } else {
    graphics::box()
  }
  return(invisible(list(matplot_args = matplot_args, axis_args = axis_args)))
}
