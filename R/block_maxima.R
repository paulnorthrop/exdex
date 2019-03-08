# =============================== sliding_maxima ==============================

#' Sliding block maxima
#'
#' Calculates the maxima of all blocks of \code{b} contiguous values in the
#' vector \code{x}. In \code{\link{exdex}} this function is used within
#' \code{\link{all_maxima}} to provide the first step of computations in
#' \code{\link{spm}}.
#'
#' @keywords internal
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @details The function \code{\link[zoo]{rollapply}} in the
#' \code{zoo} package is used with argument \code{FUN = max} to calculate
#'   sliding (rolling) maxima.  \code{na.rm = TRUE} is passed to
#'   \code{max} so that blocks containing missing values produce a
#'   non-missing result.
#' @return A list containing
#'   \itemize{
#'     \item {\code{y} : } {a numeric vector of the sliding block maxima.}
#'     \item {\code{x} : } {the input \code{x}.}
#'   }
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @seealso \code{\link{all_maxima}} for the calculation of both sliding and
#'   disjoint maxima.
#' @references Zeileis, A. and Grothendieck, G. (2005). zoo: S3
#'   Infrastructure for Regular and Irregular Time Series. \emph{Journal of
#'   Statistical Software}, \strong{14}(6), 1-27.
#'   \url{https://doi.org/10.18637/jss.v014.i06}
#' @examples
#' \dontrun{
#' x <- 1:11
#' sliding_maxima(x, 3)
#' }
sliding_maxima <- function(x, b = 1){
  y <- as.numeric(zoo::rollapply(data = zoo::zoo(x), width = b, FUN = max,
                                 na.rm = TRUE))
  return(list(y = y, x = x))
}

# =============================== disjoint_maxima =============================

#' Disjoint block maxima
#'
#' Calculates the maxima of disjoint blocks of \code{b} contiguous values in
#' the vector \code{x}.  In \code{\link{exdex}} this function
#' is used only for the purposes of checking \code{\link{all_maxima}}.
#'
#' @keywords internal
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @param which_dj A character scalar.  Determines Which set of disjoint
#'   maxima are calculated: \code{"first"}, only the set whose first block
#'   starts on the first observation in \code{x}; \code{"last"}, only the
#'   set whose last block end on the last observation in \code{x}.
#' @details If \code{length(x)} is not an integer multiple of \code{b} then
#'   only the first \code{b * floor(n / b)} observations in \code{x} are used.
#'   \code{na.rm = TRUE} is passed to \code{max} so that blocks containing
#'   missing values produce a non-missing result.
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
#' \dontrun{
#' x <- 1:11
#' disjoint_maxima(x, 3)
#' disjoint_maxima(x, 3, which_dj = "last")
#' }
disjoint_maxima <- function(x, b = 1, which_dj = c("first", "last")){
  which_dj <- match.arg(which_dj)
  if (which_dj == "last") {
    x <- rev(x)
  }
  n <- length(x)
  # number of maxima of blocks of length b
  n_max <- floor(n / b)
  # take only the first 1 to n_max*b observations
  x <- x[1:(n_max * b)]
  # block indicator: 1, ..., 1, ..., n_max, ..., n_max
  ind <- rep(1:n_max, each = b)
  # calculate block maxima
  y <- as.numeric(tapply(x, ind, max, na.rm = TRUE))
  if (which_dj == "last") {
    x <- rev(x)
    y <- rev(y)
  }
  return(list(y = y, x = x))
}

# ============================= all_disjoint_maxima ===========================

#' All sets of disjoint block maxima
#'
#' Calculates all sets of the maxima of disjoint blocks of \code{b} contiguous
#' values in the vector \code{x}.  In \code{\link{exdex}} this function
#' is used only for the purposes of checking \code{\link{all_maxima}}.
#'
#' @keywords internal
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @param which_dj A character scalar.  Determines Which sets of disjoint
#'   maxima are calculated: \code{"all"}, all sets; \code{"first"}, only the
#'   set whose first block starts on the first observation in \code{x};
#'   \code{"last"}, only the set whose last block end on the last observation
#'   in \code{x}.
#' @details If \code{length(x)} is not an integer multiple of \code{b} then
#'   only the first \code{b * floor(n / b)} observations in \code{x} are used.
#'   \code{na.rm = TRUE} is passed to \code{max} so that blocks containing
#'   missing values produce a non-missing result.
#' @return A list containing
#'   \itemize{
#'     \item {\code{y} : } {the sliding block maxima.}
#'     \item {\code{x} : } {the subset of the input \code{x} that contributes
#'       to the values in \code{y}.}
#'   }
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @seealso \code{\link{all_maxima}} for the calculation of both sliding and
#'   disjoint maxima.
#' @examples
#' \dontrun{
#' x <- 1:11
#' all_disjoint_maxima(x, 3)
#' all_disjoint_maxima(x, 3, which_dj = "first")
#' all_disjoint_maxima(x, 3, which_dj = "last")
#' }
all_disjoint_maxima <- function(x, b = 1,
                                which_dj = c("all", "first", "last")){
  which_dj <- match.arg(which_dj)
  n <- length(x)
  # The number of maxima of blocks of length b
  n_max <- floor(n / b)
  # Set the possible first indices
  first_value <- switch(which_dj,
                        all = 1:(n - n_max * b + 1),
                        first = 1,
                        last = n - n_max * b + 1)
  # block indicator: 1, ..., 1, ..., n_max, ..., n_max
  ind <- rep(1:n_max, each = b)
  get_maxima <- function(first) {
    last <- first + n_max * b - 1
    # take only the first_value to n_max * b + first_value - 1 observations
    xx <- x[first:last]
    # calculate block maxima
    y <- as.numeric(tapply(xx, ind, max, na.rm = TRUE))
    return(c(y, xx))
  }
  temp <- vapply(first_value, FUN = get_maxima, numeric(n_max * (b + 1)))
  y_mat <- temp[1:n_max, , drop = FALSE]
  x_mat <- temp[-(1:n_max), , drop = FALSE]
  return(list(y_mat = y_mat, x_mat = x_mat))
}

# ================================= all_maxima ================================

#' Sliding and disjoint block maxima
#'
#' Calculates the (sliding) maxima of all blocks of \code{b} contiguous values
#' and all sets of the maxima of disjoint blocks of \code{b} contiguous values
#' in the vector \code{x}.  This provides the first step of computations in
#' \code{\link{spm}}.
#'
#' @keywords internal
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @param which_dj A character scalar.  Determines Which sets of disjoint
#'   maxima are calculated: \code{"all"}, all sets; \code{"first"}, only the
#'   set whose first block starts on the first observation in \code{x};
#'   \code{"last"}, only the set whose last block end on the last observation
#'   in \code{x}.
#' @details \strong{Sliding maxima.} The function \code{\link[zoo]{rollapply}} in the
#' \code{zoo} package is used with argument \code{FUN = max} to calculate
#'   sliding (rolling) maxima.
#'
#'   \strong{Disjoint maxima.}  If \code{n = length(x)} is an integer
#'   multiple of \code{b}, or if \code{which_dj = "first"} or
#'   \code{which_dj = "last"} then only one set of \code{n / b} disjoint
#'   block maxima are returned.
#'   Otherwise, \code{n - floor(n / b) * b + 1} sets of \code{floor(n / b)}
#'   disjoint block maxima are returned.  Set \code{i} are the disjoint maxima
#'   of \code{x[i:(i + floor(n / b) * b - 1)]}.  That is, all possible sets
#'   of contiguous disjoint maxima achieving the maxima length of
#'   \code{floor(n / b)} are calculated.
#'
#'   In both instances \code{na.rm = TRUE} is passed to \code{\link{max}} so
#'   that blocks containing missing values produce a non-missing result.
#'
#'   Also returned are the values in \code{x} that contribute to each set
#'   of block maxima.
#' @return A list containing
#'   \itemize{
#'     \item {\code{ys} : } {a numeric vector containing one set of sliding
#'       block maxima.}
#'     \item {\code{xs} : } {a numeric vector containing the values that
#'       contribute to \code{ys}, that is, the whole input vector \code{x}.}
#'     \item {\code{yd} : } {if \code{which_dj = "all"} a \code{floor(n / b)}
#'       by \code{n - floor(n / b) * b + 1} numeric matrix.  Each column
#'       contains a set of disjoint maxima.  Otherwise, a \code{floor(n / b)}
#'       by 1 numeric matrix containing one set of block maxima.}
#'     \item {\code{xd} : } {if \code{which_dj = "all"} a
#'       \code{floor(n / b) * b} by \code{n - floor(n / b) * b + 1} numeric
#'       matrix.  Each column contains the values in \code{x} that contribute
#'       to the corresponding column in \code{yd}.  Otherwise, a
#'       \code{floor(n / b)} by 1 numeric matrix containing one the one set of
#'       the values in \code{x} that contribute to \code{yd}.}
#'   }
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @seealso \code{\link{sliding_maxima}} for the calculation of the maxima
#'   over sliding blocks.
#' @examples
#' \dontrun{
#' x <- 1:11
#' all_maxima(x, 3)
#' all_maxima(x, 3, which_dj = "first")
#' all_maxima(x, 3, which_dj = "last")
#' }
all_maxima <- function(x, b = 1, which_dj = c("all", "first", "last")){
  which_dj <- match.arg(which_dj)
  # First calculate the sliding block maxima.  All the disjoint maxima that
  # we need are contained in s_max, and we need the sliding maxima anyway
  s_max <- sliding_maxima(x = x, b = b)
  # The number of maxima of blocks of length b
  n <- length(x)
  n_max <- floor(n / b)
  # Set the possible first indices
  first_value <- switch(which_dj,
                        all = 1:(n - n_max * b + 1),
                        first = 1,
                        last = n - n_max * b + 1)
  # A function to return block maxima and contributing values starting from
  # the first value first
  get_maxima <- function(first) {
    s_ind <- seq.int(from = first, by = b, length.out = n_max)
    return(c(s_max$y[s_ind], x[first:(first + n_max * b - 1)]))
  }
  temp <- vapply(first_value, FUN = get_maxima, numeric(n_max * (b + 1)))
  yd <- temp[1:n_max, , drop = FALSE]
  xd <- temp[-(1:n_max), , drop = FALSE]
  return(list(ys = s_max$y, xs = s_max$x, yd = yd, xd = xd))
}
