#context("summary.kgaps")

# Check that summary.kgaps() returns the correct values

# For inc_cens = FALSE -----

### S&P 500 index
u <- quantile(sp500, probs = 0.60)
theta <- kgaps(sp500, u, inc_cens = FALSE)
res1 <- signif(c(theta$theta, theta$se),
               digits = max(3, getOption("digits") - 3L))
res2 <- summary(theta)$matrix

test_that("Fitted object and summary() agree", {
  testthat::expect_equivalent(res1, res2)
})

# Check that kgaps() gives the same output for vector and matrix data input

theta2 <- kgaps(as.matrix(sp500), u, inc_cens = FALSE)

theta$call <- NULL
theta2$call <- NULL

test_that("kgaps: vector data vs matrix data", {
  testthat::expect_equivalent(theta, theta2)
})
