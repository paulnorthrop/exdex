# ================================= choose_ud =================================

#' Threshold \eqn{u} and runs parameter \eqn{D} diagnostic for the \eqn{D}-gaps
#' estimator
#'
#' Creates data for a plot to aid the choice of the threshold and
#' run parameter \eqn{D} for the \eqn{D}-gaps estimator (see
#' \code{\link{dgaps}}).  \code{\link{plot.choose_ud}} creates the plot.
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then \code{\link{split_by_NAs}} is
#'   used to divide the data into sequences of non-missing values.
#' @param u,D Numeric vectors.  \code{u} is a vector of
#'   extreme value thresholds applied to data.  \code{D} is a vector of values
#'   of the run parameter \eqn{D}, as defined in Holesovsky and Fusek (2020).
#'   See \code{\link{dgaps}} for more details.
#'
#'   Any values in \code{u} that are greater than all the observations in
#'   \code{data} will be removed without a warning being given.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations. It is known that these times are greater
#'   than or equal to the time observed.
#'   If \code{data} has multiple columns then there will be right-censored
#'   first and last inter-exceedance times for each column.  See
#'   \strong{Details} in \code{\link{dgaps}}.
#' @details For each combination of threshold in \code{u} and \eqn{D}
#'   in \code{D} the functions \code{\link{dgaps}} and \code{\link{dgaps_imt}}
#'   are called in order to estimate \eqn{\theta} and to perform the
#'   information matrix test of Holesovsky and Fusek (2020).
#' @references Holesovsky, J. and Fusek, M. Estimation of the extremal index
#'   using censored distributions. Extremes 23, 197-213 (2020).
#'   \doi{10.1007/s10687-020-00374-3}
#' @return An object (a list) of class \code{c("choose_ud", "exdex")}
#'   containing
#'   \item{imt }{an object of class \code{c("dgaps_imt", "exdex")} returned
#'     from \code{\link{dgaps_imt}}.}
#'   \item{theta }{a \code{length(u)} by \code{length(D)} matrix.
#'     Element (i,j) of \code{theta} contains an object (a list) of class
#'     \code{c("dgaps", "exdex")}, a result of a call
#'     \code{dgaps(data, u[j], D[i])} to \code{\link{dgaps}}.}
#' @seealso \code{\link{dgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{D}-gaps model.
#' @seealso \code{\link{dgaps_imt}} for the information matrix test under the
#'   \eqn{D}-gaps model
#' @seealso \code{\link{plot.choose_ud}} to produce the diagnostic plot.
#' @examples
#' ### S&P 500 index
#'
#' # Multiple thresholds and left-censoring parameters
#' u <- quantile(sp500, probs = seq(0.2, 0.9, by = 0.1))
#' imt_theta <- choose_ud(sp500, u = u, D = 1:5)
#' plot(imt_theta)
#' plot(imt_theta, uprob = TRUE)
#' plot(imt_theta, y = "theta")
#'
#' # One left-censoring parameter D, many thresholds u
#' u <- quantile(sp500, probs = seq(0.2, 0.9, by = 0.1))
#' imt_theta <- choose_ud(sp500, u = u, D = 1)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#'
#' # One threshold u, many left-censoring parameters D
#' u <- quantile(sp500, probs = 0.9)
#' imt_theta <- choose_ud(sp500, u = u, D = 1:5)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_ud(newlyn, u = u, D = 1:5)
#' plot(imt_theta, uprob = TRUE)
#'
#' ### Cheeseboro wind gusts (a matrix containing some NAs)
#'
#' probs <- c(seq(0.5, 0.95, by = 0.05), 0.99)
#' u <- quantile(cheeseboro, probs = probs, na.rm = TRUE)
#' imt_theta <- choose_ud(cheeseboro, u, D = 1:6)
#' plot(imt_theta, uprob = FALSE, lwd = 2)
#'
#' ### Uccle July temperatures
#'
#' probs <- c(seq(0.7, 0.95, by = 0.05), 0.99)
#' u <- quantile(uccle720m, probs = probs, na.rm = TRUE)
#' imt_theta <- choose_ud(uccle720m, u, D = 1:5)
#' plot(imt_theta, uprob = TRUE, lwd = 2)
#' @export
choose_ud <- function(data, u, D = 1, inc_cens = TRUE) {
  # If there are missing values then use split_by_NAs to extract sequences
  # of non-missing values
  if (anyNA(data) && is.null(attr(data, "split_by_NAs_done"))) {
    data <- split_by_NAs(data)
  } else if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  # Remove any thresholds that are greater than all the observations
  u_ok <- vapply(u, function(u) any(data > u), TRUE)
  u <- u[u_ok]
  n_u <- length(u)
  n_D <- length(D)
  theta <- matrix(rep(list(), n_u * n_D), n_D, n_u)
  # Function to set the correct element of the matrix of lists theta
  # i indexes D, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  for (i in 1:n_D) {
    for (j in 1:n_u) {
      theta[[comp(i, j)]] <- dgaps(data = data, u = u[j], D = D[i],
                                   inc_cens = inc_cens)
    }
  }
  imt <- dgaps_imt(data = data, u = u, D = D, inc_cens = inc_cens)
  res <- list(imt = imt, theta = theta)
  class(res) <- c("choose_ud", "exdex")
  return(res)
}

# ============================= plot.choose_ud ===============================

#' Plot threshold \eqn{u} and runs parameter \eqn{D} diagnostic for the
#' \eqn{D}-gaps estimator
#'
#' \code{plot} method for objects inheriting from class \code{"choose_ud"},
#' returned from \code{\link{choose_ud}}
#'
#' @param x an object of class \code{c("choose_ud", "exdex")}, a result of a
#'   call to \code{\link{choose_ud}}.
#' @param y A character scalar indicating what should be plotted on the
#'   vertical axes of the plot: information matrix test statistics (IMTS)
#'   if \code{y = "imts"} and estimates of \eqn{\theta} if \code{y = "theta"}.
#'   If \code{y = "theta"}, and either \code{x$u} or \code{x$D} have length
#'   one, then 100\code{level}\% confidence intervals are added to the plot.
#' @param level A numeric scalar in (0, 1).  The confidence level used in
#'   calculating confidence intervals for \eqn{\theta}.  Only relevant if
#'   \code{y = "theta"} and either \code{x$u} or \code{x$D} have length one.
#' @param interval_type A character scalar.  The type of confidence interval
#'   to be plotted, if \code{y = "theta"}.  See \code{\link{confint.dgaps}}.
#' @param conf_scale A character scalar.  If \code{interval_type = "norm"} then
#'   \code{conf_scale} determines the scale on which we use approximate
#'   large-sample normality of the estimator to estimate confidence intervals.
#'   See \code{\link{confint.dgaps}}.
#' @param alpha A numeric vector with entries in (0, 1). The size of the test
#'   to be performed.
#' @param constrain A logical scalar.  The argument \code{constrain} to
#'  \code{\link{confint.dgaps}}.
#' @param for_abline Only relevant when \code{y = "imts"} and at one of
#'   \code{u} or \code{D} is scalar. A list of graphical parameters to be
#'   passed to \code{\link{abline}} to indicate the critical value of the
#'   information matrix test (IMT) implied by \code{alpha}.
#' @param digits An integer. Used for formatting the value of the threshold
#'   with \code{\link[base:Round]{signif}} before adding its value to a plot.
#' @param uprob A logical scalar. Should we plot \code{x$u} on the
#'   horizontal axis (\code{uprob = FALSE}) or the approximate sample quantile
#'   to which \code{x$u} corresponds (\code{uprob = FALSE})?
#' @param leg_pos A character scalar.  The position of any legend added to
#'   a plot.  Only relevant when both the arguments \code{u} and \code{D}
#'   in the call to \code{\link{choose_ud}} have length greater than one.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}.
#' @details The type of plot produced depends mainly on \code{y}.
#'
#'   If \code{y = "imts"} then the values of IMTS are plotted against the
#'   thresholds in \code{x$u} (or their corresponding approximate sample
#'   quantile levels if \code{uprob = TRUE}) for each value of \eqn{D}
#'   in \code{x$D}.  Horizontal lines are added to indicate the critical
#'   values of the IMT for the significance levels in \code{alpha}.
#'   We would not reject at the 100\code{alpha}\% level combinations of
#'   threshold and \eqn{D} corresponding to values of the IMTS that fall
#'   below the line.
#'
#'   If \code{y = "theta"} then estimates of \eqn{\theta} are plotted on the
#'   vertical axis.  If both \code{x$u} and \code{x$D$} have length greater
#'   than one then only these estimates are plotted.  If either \code{x$u}
#'   or \code{x$D} have length one then approximate 100\code{level}\%
#'   confidence intervals are added to the plot and the variable,
#'   \code{x$u} or \code{x$D} that has length greater than one is plotted on
#'   the horizontal axis.
#' @return Nothing is returned.
#' @seealso \code{\link{choose_ud}}.
#' @section Examples:
#' See the examples in \code{\link{choose_ud}}.
#' @export
plot.choose_ud <- function(x, y = c("imts", "theta"), level = 0.95,
                           interval_type = c("norm", "lik"),
                           conf_scale = c("theta", "log"), alpha = 0.05,
                           constrain = TRUE,
                           for_abline = list(lty = 2, lwd = 1, col = 1),
                           digits = 3, uprob = FALSE,
                           leg_pos = if (y == "imts") "topright" else "topleft",
                           ...) {
  y <- match.arg(y)
  interval_type <- match.arg(interval_type)
  conf_scale <- match.arg(conf_scale)
  # Extract the values of D and u
  D <- x$imt$D
  u <- x$imt$u
  n_D <- length(D)
  n_u <- length(u)
  # Approximate sample quantiles of threshold u
  u_ps <- rownames(x$imt$imt)
  if (n_D == 1 && n_u == 1) {
    stop("Object contains only 1 threshold and 1 value of D")
  }
  # Function to set the correct element of the matrix of lists theta
  # i indexes D, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  # My plotting functions: to give defaults but allow the user to override
  my_matplot <- function(x, y, ..., type = "l", lty = my_lty, col = my_col,
                         xlab = my_xlab, ylab = my_ylab, ylim = my_ylim,
                         lwd = my_lwd) {
    graphics::matplot(x = x, y = y, ..., type = type, lty = lty, col = col,
                      xlab = xlab, ylab = ylab, ylim = ylim, lwd = lwd)
  }
  my_title <- function(..., main = my_main) {
    graphics::title(..., main = main)
  }
  # Critical value for the IMT (in case we need it)
  crit <- stats::qchisq(alpha, df = 1, lower.tail = FALSE)
  # One of D or u is scalar
  cond1 <- n_D == 1 && n_u > 1
  cond2 <- n_D > 1 && n_u == 1
  my_lwd <- 1
  if (cond1 || cond2)  {
    my_col <- my_lty <- 1
    max_ud <- max(n_D, n_u)
    ymat <- matrix(NA, ncol = 3, nrow = max_ud)
    if (cond1) {
      if (uprob) {
        xvec <- u_ps
      } else {
        xvec <- u
      }
    } else {
      xvec <- D
    }
    loop_vec <- 1:max_ud
    if (y == "theta") {
      for (ij in loop_vec) {
        if (cond1) {
          dgaps_object <- x$theta[[comp(1, ij)]]
        } else {
          dgaps_object <- x$theta[[comp(ij, 1)]]
        }
        temp <- confint(dgaps_object, level = level,
                        interval_type = interval_type,
                        conf_scale = conf_scale, constrain = constrain)
        ymat[ij, 1] <- dgaps_object$theta
        ymat[ij, 2:3] <- temp$cis
      }
      my_ylab <- "theta"
      my_xlab <- ifelse(cond1, "threshold u", "left-censoring parameter D")
      if (uprob & cond1) {
        my_xlab <- "sample quantile level of threshold u"
      }
      my_ylim <- c(0, 1)
      my_matplot(xvec, ymat, ...)
      my_main <- ifelse(cond1, paste0("left-censoring parameter D = ", D),
                        paste0("threshold u = ", signif(u, digits = digits),
                               " (", u_ps, "% quantile)"))
      my_title(...)
    } else {
      my_ylab <- "IMTS"
      my_xlab <- ifelse(cond1, "threshold u", "left-censoring parameter D")
      if (uprob & cond1) {
        my_xlab <- "sample quantile level of threshold u"
      }
      my_ylim <- c(0, max(x$imt$imt, na.rm = TRUE))
      if (cond1) {
        if (uprob) {
          xvec <- u_ps
        } else {
          xvec <- x$imt$u
        }
        ymat <- x$imt$imt
      } else {
        xvec <- x$imt$D
        ymat <- t(x$imt$imt)
      }
      my_matplot(xvec, ymat, ...)
      my_main <- ifelse(cond1, paste0("left-censoring parameter D = ", D),
                        paste0("threshold u = ", signif(u, digits = digits),
                               " (", u_ps, "% quantile)"))
      my_title(...)
      for_abline <- c(for_abline, list(h = crit))
      do.call(graphics::abline, for_abline)
      graphics::mtext(as.character(alpha), 4, at = crit, las = 1, cex = 0.8,
                      adj = 1, padj = 1)
    }
  } else {
    my_lty <- 1
    my_col <- 1:n_D
    if (uprob) {
      xvec <- u_ps
    } else {
      xvec <- x$imt$u
    }
    if (y == "theta") {
      ymat <- x$imt$theta
      my_ylab <- "theta"
      my_ylim <- c(0, 1)
    } else {
      ymat <- x$imt$imt
      my_ylab <- "IMTS"
      my_ylim <- c(0, max(x$imt$imt, na.rm = TRUE))
    }
    if (uprob) {
      my_xlab <- "sample quantile level of threshold u"
    } else {
      my_xlab <- "threshold u"
    }
    my_matplot(xvec, ymat, ...)
    if (y == "imts") {
      for_abline <- c(for_abline, list(h = crit))
      do.call(graphics::abline, for_abline)
      graphics::mtext(as.character(alpha), 4, at = crit, las = 1, cex = 0.8,
                      adj = 1, padj = 1)
    }
    user_args <- list(...)
    if (is.null(user_args$lty)) {
      leg_lty <- my_lty
    } else {
      leg_lty <- user_args$lty
    }
    if (is.null(user_args$col)) {
      leg_col <- my_col
    } else {
      leg_col <- user_args$col
    }
    if (is.null(user_args$lwd)) {
      leg_lwd <- my_lwd
    } else {
      leg_lwd <- user_args$lwd
    }
    graphics::legend(leg_pos, legend = paste0("D = ", D), lty = leg_lty,
                     col = leg_col, lwd = leg_lwd)
  }
  return(invisible())
}
