#' Block length diagnostic plot
#'
#' Description
#' @examples
#' #res <- choose_b(newlyn, 20:60)
#' #plot(res)
choose_b <- function(data, b, sliding = TRUE,
                     bias_adjust = c("BB3", "BB1", "N", "none"),
                     constrain = TRUE, varN = TRUE,
                     parm = c("both", "N2015", "BB2018"), level = 0.95,
                     conf_scale = c("theta", "log"),
                     interval_type = c("lik", "sym"),
                     type = c("vertical", "cholesky", "spectral", "none")) {
  Call <- match.call(expand.dots = TRUE)
  # All other inputs are checked in the calls to spm() and confint.spm()
  interval_type <- match.arg(interval_type)
  parm <- match.arg(parm)
  if (parm == "both") {
    which_ests <- c("N2015", "BB2018")
  } else {
    which_ests <- parm
  }
  which_cis <- paste0(which_ests, interval_type)
  # Objects in which to store the estimates and confidence intervals
  n_b <- length(b)
  empty_matrix <- matrix(0, nrow = n_b, ncol = length(which_ests))
  colnames(empty_matrix) <- which_ests
  rownames(empty_matrix) <- b
  theta <- lower <- upper <- empty_matrix
  # Loop over block lengths
  for (i in 1:n_b) {
    ests <- spm(data = data, b = b[i], bias_adjust = bias_adjust,
                constrain = constrain, varN = varN)
    theta[i, ] <- ests$theta[which_ests]
    # Avoid chandwich::conf_intervals()'s profiling messages
    cis <- suppressMessages(confint(ests, parm = parm, constrain = constrain,
                                    conf_scale = conf_scale,
                                    bias_adjust = TRUE, type = type))
    lower[i, ] <- cis[which_cis, 1]
    upper[i, ] <- cis[which_cis, 2]
  }
  res <- list(theta = theta, lower = lower, upper = upper, b = b)
  res$call <- Call
  class(res) <- c("exdex", "block")
  return(res)
}


# =========================== plot.exdex ===========================

#' Plot diagnostics for an exdex object
#'
#' \code{plot} method for objects of class "exdex" returned from
#' \code{\link{choose_b}}
#'
#' @param x an object of clas "exdex", a result of a call to
#'   \code{\link{choose_b}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}}, \code{\link[graphics]{axis}}
#'   and/or \code{\link[graphics]{segments}}.
#' @details Produces a simple threshold diagnostic plot based on the object
#'   returned from \code{\link{stability}}.
#'   The MLEs of the GP shape parameter $\eqn{\xi}$ and
#'   approximate \code{conf}\% confidence intervals
#'   for \eqn{\xi} are plotted against the threshold used to fit the GP model.
#'   This plot is used to choose a threshold above which the underlying GP
#'   shape parameter may be approximately constant. See Chapter 4 of
#'   Coles (2001).  See also the vignette "Introducing threshr".
#'   as described in .
#'   See also the vignette "Introducing threshr".
#' @return In addition to producing the plot a list of the arguments used
#'   by \code{\link[graphics]{matplot}}, \code{\link[graphics]{axis}} is
#'   returned (invisibly).
#' @seealso \code{\link{choose_b}}.
plot.exdex <- function(x, y, ..., vertical = TRUE) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  y_data <- cbind(x$lower, x$theta, x$upper)
  y_lab <- expression(theta)
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
    matplot_args$type <- c("l", "b", "l")
  }
  if (is.null(matplot_args$pch)) {
    matplot_args$pch <- c(0, 16, 0)
  }
  if (is.null(matplot_args$col)) {
    matplot_args$col <- 1
  }
  if (is.null(matplot_args$lty)) {
    if (vertical) {
      matplot_args$lty <- c(0, 1, 0)
    } else {
      matplot_args$lty <- c(2, 1, 2)
    }
  }
  all_args <- c(xy_args, matplot_args)
  do.call(graphics::matplot, c(all_args, axes = FALSE))
  if (vertical) {
    segments_args$x0 <- x_data
    segments_args$x1 <- x_data
    segments_args$y0 <- x$lower
    segments_args$y1 <- x$upper
    do.call(graphics::segments, segments_args)
  }
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
