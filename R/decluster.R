# Make notation consistent with iwls.R

#' Declustering
#'
#' Describe
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k An integer scalar
#' @examples
#' thresh <- quantile(newlyn, probs = 0.99)
#' x <- find_clusters(newlyn, thresh)
find_clusters <- function(data, thresh, k = 1) {
  if (any(is.na(data))) {
    stop("No missing values are allowed in ''data''")
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > thresh]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k,0)
  # Threshold inter-exceedances times that are not larger than k units
  # (i.e. a K-gap equal to zero) are assigned to the same cluster
  print(exc_u)
  print(which(S_k == 0))
  print(which(S_k > 0))
  # A cluster starts with an exceedance that is preceded by a positive K-gap
  # and ends with the next exceedance (possibly the exceedance itself) that
  # is followed by a positive K-gap
  clus <- numeric(nx)
  #
  exc_ind <- data > thresh
  print(exc_ind)
  if (k == 1) {
    cluster_durations <- rle(exc_ind)
  } else {
    print("HERE")
    print(zoo::rollsum(exc_ind, k = k))
    print(as.numeric(zoo::rollsum(exc_ind, k = k) > 0))
    cluster_durations <- rle(zoo::rollsum(exc_ind, k = k) > 0)
    which_to_reduce <- cluster_durations$values > 0
    cluster_durations$lengths[which_to_reduce] <-
      cluster_durations$lengths[which_to_reduce] - k + 1
  }
  return(cluster_durations)
}

# What do I want?
# 1. For FW2012: a cluster membership indicator of the same length as data.
#    0 for no cluster, i for cluster i
# 2. For bootstrap estimation of CIs for theta:
