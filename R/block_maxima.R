# =============================== sliding_maxima ==============================

#' Sliding block maxima
#'
#' Calculates the maxima of all blocks of \code{b} contiguous values in the
#' vector \code{x}.
#'
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @details The function \code{\link[zoo]{rollapply}} in the
#' \code{zoo} package is used with argument \code{FUN = max} to calculate
#'   sliding (rolling) maxima.  \code{na.rm = TRUE} is passed to
#'   \code{max} so that blocks containing missing values produce a
#'   non-missing result.
#'
#'   \code{sliding_maxima} is used within \code{\link{spm}} to perform
#'   semiparametric estimation of the extremal index based on block maxima.
#' @return A list containing
#'   \itemize{
#'     \item {\code{y} : } {a numeric vector of the sliding block maxima.}
#'     \item {\code{x} : } {the input \code{x}, included for
#'       consistency with the output from \code{\link{disjoint_maxima}}}.
#'   }
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @seealso \code{\link{disjoint_maxima}} for the calculation of the maxima
#'   over disjoint blocks.
#' @references Zeileis, A. and Grothendieck, G. (2005). zoo: S3
#'   Infrastructure for Regular and Irregular Time Series. \emph{Journal of
#'   Statistical Software}, \strong{14}(6), 1-27.
#'   \url{https://doi.org/10.18637/jss.v014.i06}
#' @export
sliding_maxima <- function(x, b = 1){
  y <- as.numeric(zoo::rollapply(data = zoo::zoo(x), width = b, FUN = max,
                                 na.rm = TRUE))
  return(list(y = y, x = x))
}

# =============================== disjoint_maxima =============================

#' Disjoint block maxima
#'
#' Calculates the maxima of disjoint blocks of \code{b} contiguous values in
#' the vector \code{x}.
#'
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @details If \code{length(x)} is not an integer multiple of \code{b} then
#'   only the first \code{b * floor(n / b)} observations in \code{x} are used.
#'   \code{na.rm = TRUE} is passed to \code{max} so that blocks containing
#'   missing values produce a non-missing result.
#'
#'   \code{disjoint_maxima} is used within \code{\link{spm}} to perform
#'   semiparametric estimation of the extremal index based on block maxima.
#' @return A list containing
#'   \itemize{
#'     \item {\code{y} : } {the sliding block maxima.}
#'     \item {\code{x} : } {the subset of the input \code{x} that contributes
#'       to the values in \code{y}.}
#'   }
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @seealso \code{\link{sliding_maxima}} for the calculation of the maxima
#'   over sliding blocks.
#' @export
disjoint_maxima <- function(x, b = 1){
  n <- length(x)
  # number of maxima of blocks of length b
  n_max <- floor(n / b)
  # take only the first 1 to n_max*b observations
  x <- x[1:(n_max * b)]
  # block indicator: 1, ..., 1, ..., n_max, ..., n_max
  ind <- rep(1:n_max, each = b)
  # calculate block maxima
  y <- as.numeric(tapply(x, ind, max, na.rm = TRUE))
  return(list(y = y, x = x))
}
