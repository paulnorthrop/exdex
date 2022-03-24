# ================================ split_by_NAs ============================= #

#' Divides data into parts that contain no missing values
#'
#' Splits the values in a numeric matrix column-wise into sequences of
#' non-missing values.
#'
#' @param x A vector or matrix.
#' @details For each column in \code{x}, \code{split_by_NAs} finds runs of
#'   values that contain no missing values and assigns them to a column in the
#'   matrix that is returned.  Different columns are treated separately.
#'   If there are no missing values in a column then that column appears
#'   unmodified in the output matrix.  Please see the \strong{Examples}
#'   for illustrations.
#' @return A matrix containing a column for each run of non-missing values in
#'   \code{x}.  The number of rows is equal to the longest run of non-missing
#'   values in \code{x} and will therefore be at most \code{nrow{x}}.  The
#'   matrix is padded with \code{NA} values at the end of each column, where
#'   necessary.
#'
#'   The returned object has an attribute called \code{split_by_NAs_done}
#'   whose value is \code{TRUE}, so that in programming one can avoid calling
#'   \code{split_by_NAs} more than once.
#' @examples
#' # Create a simple numeric matrix and insert some NAs
#' x <- matrix(1:50, 10, 5)
#' x[c(3, 8), 1] <- NA
#' x[c(1:2, 5, 10), 3] <- NA
#' x[1:3, 4] <- NA
#' x[7:10, 5] <- NA
#' x
#'
#' res <- split_by_NAs(x)
#' res
#'
#' # An example of a character matrix
#' x <- matrix(c(letters, letters[1:18]), 11, 4)
#' x[c(1:2, 5:11), 2] <- NA
#' x[c(2:4, 6:11), 3] <- NA
#' x[1:10, 4] <- NA
#'
#' res <- split_by_NAs(x)
#' res
#' @export
split_by_NAs <- function(x) {
  #  If x is not a matrix then make it a matrix
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  # temp is a list of lists, one entry in main list for each column in x
  # temp[[i]]$lengths: lengths of successive sequences of NAs and non-NAs
  # temp[[i]]$values: a FALSE indicates a non-NA component in temp$lengths
  temp <- apply(!is.na(x), 2, rle)
  # A vector of the columns in which each sequence lives in x
  column <- rep(1:ncol(x), vapply(X = temp, FUN = function(x) sum(x$values),
                                  FUN.VALUE = 0.0))
  if (length(temp) == 1) {
    max_leng <- max(temp[[1]]$lengths)
    n_seq <- sum(temp[[1]]$values)
  } else {
    # Find the length max_leng of the longest sequence of non-missing values
    max_leng <- Reduce(f = function(...) Map(max, ...), temp)$lengths
    # Find the total number of sequences of non-missing values
    n_seq <- Reduce(f = function(...) Map(sum, ...), temp)$values
  }
  # Create vectors, from and to, containing the starting and ending rows for
  # each sequence, so that sequence i is x[from[i]:to[i], column[i]]
  from_fn <- function(x) {
    return((cumsum(x$lengths) - x$lengths + 1)[x$values])
  }
  to_fn <- function(x) {
    return(cumsum(x$lengths)[x$values])
  }
  from <- unlist(sapply(temp, from_fn), use.names = FALSE)
  to <- unlist(sapply(temp, to_fn), use.names = FALSE)
  # Create the new matrix of data, filling columns with NAs if necessary
  newx_fn <- function(i, from, to, column) {
    na_fill <- max_leng - (to[i] - from[i] + 1)
    c(x[from[i]:to[i], column[i]], rep(NA, na_fill))
  }
  newx <- sapply(X = 1:n_seq, FUN = newx_fn, from = from, to = to,
                 column = column)
  dimnames(newx) <- NULL
  attr(newx, "split_by_NAs_done") <- TRUE
  return(newx)
}
