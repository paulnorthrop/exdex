#context("kgaps_imt_stat_simple")

# A small example dataset

u <- 0
x <- c(-1, -1, 1, -1, 1, 1, -1, -1, -1)

# 1. Check that kgaps_stat() returns the correct output

k <- 0
res1 <- kgaps_stat(data = x, u, k = k, inc_cens = FALSE)
correct1 <- list(N0 = 0, N1 = 2, sum_qs = 1, n_kgaps = 2)
test_that("simple: kgaps_stat() correct for k = 0, inc_cens = FALSE", {
  testthat::expect_equal(res1, correct1)
})
res2 <- kgaps_stat(data = x, u, k = k, inc_cens = TRUE)
correct2 <- list(N0 = 0, N1 = 3, sum_qs = 8/3, n_kgaps = 4)
test_that("simple: kgaps_stat() correct for k = 0, inc_cens = TRUE", {
  testthat::expect_equal(res2, correct2)
})

k <- 1
res1 <- kgaps_stat(data = x, u, k = k, inc_cens = FALSE)
correct1 <- list(N0 = 1, N1 = 1, sum_qs = 1/3, n_kgaps = 2)
test_that("simple: kgaps_stat() correct for k = 1, inc_cens = FALSE", {
  testthat::expect_equal(res1, correct1)
})
res2 <- kgaps_stat(data = x, u, k = k, inc_cens = TRUE)
correct2 <- list(N0 = 1, N1 = 2, sum_qs = 4/3, n_kgaps = 4)
test_that("simple: kgaps_stat() correct for k = 1, inc_cens = TRUE", {
  testthat::expect_equal(res2, correct2)
})

k <- 2
res1 <- kgaps_stat(data = x, u, k = k, inc_cens = FALSE)
correct1 <- list(N0 = 2, N1 = 0, sum_qs = 0, n_kgaps = 2)
test_that("simple: kgaps_stat() correct for k = 2, inc_cens = FALSE", {
  testthat::expect_equal(res1, correct1)
})
res2 <- kgaps_stat(data = x, u, k = k, inc_cens = TRUE)
correct2 <- list(N0 = 2, N1 = 1/2, sum_qs = 1/3, n_kgaps = 3)
test_that("simple: kgaps_stat() correct for k = 2, inc_cens = TRUE", {
  testthat::expect_equal(res2, correct2)
})

k <- 3
res1 <- kgaps_stat(data = x, u, k = k, inc_cens = FALSE)
correct1 <- list(N0 = 2, N1 = 0, sum_qs = 0, n_kgaps = 2)
test_that("simple: kgaps_stat() correct for k = 3, inc_cens = FALSE", {
  testthat::expect_equal(res1, correct1)
})
res2 <- kgaps_stat(data = x, u, k = k, inc_cens = TRUE)
correct2 <- list(N0 = 2, N1 = 0, sum_qs = 0, n_kgaps = 2)
test_that("simple: kgaps_stat() correct for k = 3, inc_cens = TRUE", {
  testthat::expect_equal(res2, correct2)
})

# 2. Check that the lengths of the vector of statistics returned from
#    kgaps_imt_stat() equal n_kgaps

test_fn <- function(k, inc_cens) {
  theta <- kgaps(data = x, u = u, k = k, inc_cens = inc_cens)$theta
  res <- kgaps_imt_stat(data = x, theta = theta, u = u, k = k,
                        inc_cens = inc_cens)
  # Find all the lengths, except n_kgaps
  lengs <- sapply(res, length)[-6]
  return(res$n_kgaps - unique(lengs))
}

test_that("simple: IMT stats lengths, k = 0, inc_cens = FALSE", {
  testthat::expect_equal(test_fn(k = 0, inc_cens = FALSE), 0)
})
test_that("simple: IMT stats lengths, k = 0, inc_cens = TRUE", {
  testthat::expect_equal(test_fn(k = 0, inc_cens = TRUE), 0)
})
test_that("simple: IMT stats lengths, k = 1, inc_cens = FALSE", {
  testthat::expect_equal(test_fn(k = 1, inc_cens = FALSE), 0)
})
test_that("simple: IMT stats lengths, k = 1, inc_cens = TRUE", {
  testthat::expect_equal(test_fn(k = 1, inc_cens = TRUE), 0)
})
test_that("simple: IMT stats lengths, k = 2, inc_cens = FALSE", {
  testthat::expect_equal(test_fn(k = 2, inc_cens = FALSE), 0)
})
test_that("simple: IMT stats lengths, k = 2, inc_cens = TRUE", {
  testthat::expect_equal(test_fn(k = 2, inc_cens = TRUE), 0)
})
test_that("simple: IMT stats lengths, k = 3, inc_cens = FALSE", {
  testthat::expect_equal(test_fn(k = 3, inc_cens = FALSE), 0)
})
test_that("simple: IMT stats lengths, k = 3, inc_cens = TRUE", {
  testthat::expect_equal(test_fn(k = 3, inc_cens = TRUE), 0)
})

# 3. Check that the new kgaps_imt() gives the same results as kgaps_imt_old()
#    We need to use inc_cens = FALSE because kgaps_imt_old() does not
#    allow inc_cens = TRUE.

k <- 1:5
# New
res <- kgaps_imt(data = x, u = u, k = k, inc_cens = FALSE)
# Old
res2 <- kgaps_imt_old(data = x, u = u, k = k)

test_that("simple: new and old kgaps_imt() agree", {
  testthat::expect_equal(res, res2)
})

# 4. Check nobs.kgaps()

nobs_fn <- function(k, inc_cens) {
  return(nobs(kgaps(data = x, u = u, k = k, inc_cens = inc_cens)))
}
test_that("simple: nobs, k = 0, inc_cens = FALSE", {
  testthat::expect_equal(nobs_fn(k = 0, inc_cens = FALSE), 2)
})
test_that("simple: nobs, k = 0, inc_cens = TRUE", {
  testthat::expect_equal(nobs_fn(k = 0, inc_cens = TRUE), 4)
})
test_that("simple: nobs, k = 0, inc_cens = FALSE", {
  testthat::expect_equal(nobs_fn(k = 1, inc_cens = FALSE), 2)
})
test_that("simple: nobs, k = 0, inc_cens = TRUE", {
  testthat::expect_equal(nobs_fn(k = 1, inc_cens = TRUE), 4)
})
test_that("simple: nobs, k = 0, inc_cens = FALSE", {
  testthat::expect_equal(nobs_fn(k = 2, inc_cens = FALSE), 2)
})
test_that("simple: nobs, k = 0, inc_cens = TRUE", {
  testthat::expect_equal(nobs_fn(k = 2, inc_cens = TRUE), 3)
})
test_that("simple: nobs, k = 0, inc_cens = FALSE", {
  testthat::expect_equal(nobs_fn(k = 3, inc_cens = FALSE), 2)
})
test_that("simple: nobs, k = 0, inc_cens = TRUE", {
  testthat::expect_equal(nobs_fn(k = 3, inc_cens = TRUE), 2)
})
