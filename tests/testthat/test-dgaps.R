#context("summary.dgaps")

### Check that summary.dgaps() returns the correct values

# For inc_cens = FALSE -----

### S&P 500 index
u <- quantile(sp500, probs = 0.60)
theta <- dgaps(sp500, u, inc_cens = FALSE)
res1 <- signif(c(theta$theta, theta$se),
               digits = max(3, getOption("digits") - 3L))
res2 <- summary(theta)$matrix

test_that("Fitted object and summary() agree", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check that dgaps() gives the same output for vector and matrix data input

# For inc_cens = FALSE -----

theta2 <- dgaps(as.matrix(sp500), u, inc_cens = FALSE)

theta$call <- NULL
theta2$call <- NULL

test_that("dgaps: vector data vs matrix data", {
  testthat::expect_equal(theta, theta2, ignore_attr = TRUE)
})

### Check that summary.dgaps() returns the correct values

# For inc_cens = TRUE -----

### S&P 500 index
u <- quantile(sp500, probs = 0.60)
theta <- dgaps(sp500, u, inc_cens = TRUE)
res1 <- signif(c(theta$theta, theta$se),
               digits = max(3, getOption("digits") - 3L))
res2 <- summary(theta)$matrix

test_that("Fitted object and summary() agree", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check that dgaps() gives the same output for vector and matrix data input

# For inc_cens = TRUE -----

theta2 <- dgaps(as.matrix(sp500), u, inc_cens = TRUE)

theta$call <- NULL
theta2$call <- NULL

test_that("dgaps: vector data vs matrix data", {
  testthat::expect_equal(theta, theta2, ignore_attr = TRUE)
})

