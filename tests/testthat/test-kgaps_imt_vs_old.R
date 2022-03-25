#context("kgaps_imt_vs_old")

# Check that the new kgaps_imt() gives the same results as kgaps_imt_old()
# We need to use inc_cens = FALSE because kgaps_imt_old() does not
# allow inc_cens = TRUE.  kgaps_imt_old() can't accommodate the cheeseboro
# data (matrix with missings)

# Newlyn

u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
k <- 1:5
# New
res <- kgaps_imt(newlyn, u = u, k = k, inc_cens = FALSE)
# Old
res2 <- kgaps_imt_old(newlyn, u = u, k = k)

test_that("newlyn: new and old kgaps_imt() agree", {
  testthat::expect_equal(res, res2)
})

# S&P 500

u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
k <- 1:5
# New
res <- kgaps_imt(sp500, u = u, k = k, inc_cens = FALSE)
# Old
res2 <- kgaps_imt_old(sp500, u = u, k = k)

test_that("S&P500: new and old kgaps_imt() agree", {
  testthat::expect_equal(res, res2)
})
