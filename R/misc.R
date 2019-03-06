# ========================== Miscellaneous functions ======================== #

# log(x), but return a constant const for an x = 0

log0const_slow <- function(x, const) {
  ifelse(x == 0, const, log(x))
}

log0const <- function(x, const) {
  return(log(x + !x) + const * !x)
}

# Empirical c.d.f. of x, evaluated at y

ecdf3 <- function(x, y) {
  return(vapply(y, function(y) mean(x <= y), 0))
}

ecdf2 <- function(x, y) {
  return(vapply(y, function(y) sum(x <= y) / length(x), 0))
}

ecdf1 <- function(x, y, lenx) {
  return(vapply(y, function(y) sum(x <= y) / lenx, 0))
}
