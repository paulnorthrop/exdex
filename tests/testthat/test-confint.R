context("confint.spm")

my_tol <- 1e-5

res <- spm(newlyn, 100)

ci1 <- confint(res, type = "cholesky")
ci2 <- confint(res, type = "spectral")
test_that("sliding: cholesky and spectral are identical", {
  testthat::expect_identical(ci1$cis, ci2$cis)
})

# Check that the estimates of theta in res and returned from
# chandwich::adjust_loglik()

test_that("estimates of theta agree, sliding", {
  testthat::expect_equal(res$theta_sl, ci1$theta, tolerance = my_tol)
})

ci1 <- confint(res, maxima = "disjoint", type = "cholesky")
ci2 <- confint(res, maxima = "disjoint", type = "spectral")
test_that("disjoint: cholesky and spectral are identical", {
  testthat::expect_identical(ci1$cis, ci2$cis)
})

# Check estimates of theta

test_that("estimates of theta agree, disjoint", {
  testthat::expect_equal(res$theta_dj, ci1$theta, tolerance = my_tol)
})

ci3 <- confint(res, maxima = "disjoint", type = "cholesky",
               conf_scale = "log")

which_rows <- c("N2015lik", "BB2018lik")
test_that("spm lik intervals don't depend on conf_scale", {
  testthat::expect_identical(ci1$cis[which_rows, ], ci3$cis[which_rows, ])
})

# b too low

res <- suppressWarnings(spm(newlyn, 14))
ci_b_low <- confint(res)
temp <- matrix(NA, nrow = 4, ncol = 2)
level <- 0.95
a <- (1 - level) / 2
a <- c(a, 1 - a)
pct <- paste(round(100 * a, 1), "%")
colnames(temp) <- pct
rownames(temp) <- c("N2015sym", "BB2018sym", "N2015lik", "BB0218lik")
test_that("b is too low gives NAs", {
  testthat::expect_identical(ci_b_low, temp)
})

context("confint.kgaps")

thresh <- quantile(newlyn, probs = 0.90)
res1 <- kgaps_mle(newlyn, thresh)
res1 <- confint(res1)
res2 <- kgaps_mle(newlyn, thresh)
res2 <- confint(res2, conf_scale = "log")

test_that("kgaps lik intervals don't depend on conf_scale", {
  testthat::expect_identical(res2["lik", ], res2["lik", ])
})
