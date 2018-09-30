pjn <- function (data, tmin, tmax, ...) {
  quant <- seq(tmin, tmax, by = 1)
  x <- (tmax - tmin) + 1
  mplot <- matrix(0, nrow = x, ncol = 2, byrow = TRUE,
                  dimnames = list(NULL, c("Theta", "Quantile")))
  ci <- matrix(0, nrow = x, ncol = 2, byrow = TRUE,
               dimnames = list(NULL, c("Lower", "Upper")))
  for (i in 1:x) {
    p <- quant[i]/100
    thresh <- quantile(data, probs = p)
    theta <- kgaps_mle(data, thresh, conf = 95)
    mplot[i, 2] <- p
    mplot[i, 1] <- theta$theta_mle
    ci[i, 1] <- theta$theta_ci[1]
    ci[i, 2] <- theta$theta_ci[2]
  }
  y <- cbind(mplot[, "Theta"], ci)
  user_args <- list(...)
  if (is.null(user_args$lty)) {
    user_args$lty <- c(1, 2, 2)
  }
  if (is.null(user_args$col)) {
    user_args$col <- 1
  }
  if (is.null(user_args$xlab)) {
    user_args$xlab <- "quantile level of threshold"
  }
  if (is.null(user_args$ylab)) {
    user_args$ylab <- "theta"
  }
  if (is.null(user_args$type)) {
    user_args$type = "l"
  }
  if (is.null(user_args$pch)) {
    user_args$pch = 16
  }
  user_args$x <- mplot[, "Quantile"]
  user_args$y <- y
  do.call(graphics::matplot, user_args)
  return_matrix <- cbind(mplot[, 2:1], ci)
  return(invisible(return_matrix))
}
