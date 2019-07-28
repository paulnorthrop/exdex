context("choose_uk")

# Check that calling choose_uk() with vector arguments thresh and k gives
# the same results as calling kgaps repeatedly with scalar arguments

thresh <- stats::quantile(newlyn, probs = c(0.1, 0.90))
k_vals <- 1:3
cres <- choose_uk(newlyn, thresh = thresh, k = k_vals)
n_thresh <- length(thresh)
comp <- function(i, j) {
  return((i - 1) * n_thresh + j)
}
for (i in 1:length(k_vals)) {
  for (j in 1:length(thresh)) {
    res <- kgaps(newlyn, thresh[j], k_vals[i])
    temp <- cres$theta[[comp(i, j)]]
    # The calls will be different1
    temp$call <- res$call <- NULL
    test_that("choose_k agrees with kgaps", {
      testthat::expect_equal(temp, res)
    })
  }
}
