#context("choose_ud")

# For inc_cens = FALSE -----

# Check that calling choose_ud() with vector arguments u and D gives
# the same results as calling dgaps repeatedly with scalar arguments

u <- stats::quantile(newlyn, probs = c(0.1, 0.90))
D_vals <- 1:3
cres <- choose_ud(newlyn, u = u, D = D_vals, inc_cens = FALSE)
n_u <- length(u)
comp <- function(i, j) {
  return((i - 1) * n_u + j)
}
for (i in 1:length(D_vals)) {
  for (j in 1:length(u)) {
    res <- dgaps(newlyn, u[j], D_vals[i], inc_cens = FALSE)
    temp <- cres$theta[[comp(i, j)]]
    # The calls will be different1
    temp$call <- res$call <- NULL
    test_that("choose_ud agrees with dgaps", {
      testthat::expect_equal(temp, res)
    })
  }
}

# =============================== plot.choose_ud = =============================

# Check that plot.choose_ud works

# S&P 500 index
u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
imt_theta <- choose_ud(sp500, u = u, D = 1:5, inc_cens = FALSE)

udplot <- plot(imt_theta)
test_that("plot.choose_ud works", {
  testthat::expect_identical(udplot, NULL)
})
udplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab", lwd = 2,
               col = 1:5)
test_that("plot.choose_ud works, user plot args", {
  testthat::expect_identical(udplot, NULL)
})

# One left-censoring parameter D, many thresholds u
u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
imt_theta <- choose_ud(sp500, u = u, D = 1, inc_cens = FALSE)
udplot <- plot(imt_theta)
test_that("plot.choose_ud works", {
  testthat::expect_identical(udplot, NULL)
})
udplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab",
               lwd = 2, col = 1:5)
test_that("plot.choose_ud works, user plot args", {
  testthat::expect_identical(udplot, NULL)
})

# One threshold u, many left-censoring parameters D
u <- quantile(sp500, probs = 0.9)
imt_theta <- choose_ud(sp500, u = u, D = 1:5, inc_cens = FALSE)
test_that("plot.choose_ud works", {
  testthat::expect_identical(udplot, NULL)
})
udplot <- plot(imt_theta, y = "theta", ylim = c(0, 1), xlab = "my xlab",
               lwd = 2, col = 1:5)
test_that("plot.choose_ud works, user plot args", {
  testthat::expect_identical(udplot, NULL)
})

