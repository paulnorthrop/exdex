context("confint.spm")

my_tol <- 1e-5

# =================================== spm ====================================

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

ci1 <- confint(res, maxima = "disjoint", interval_type = "both",
               type = "cholesky")
ci2 <- confint(res, maxima = "disjoint", interval_type = "both",
               type = "spectral")
test_that("disjoint: cholesky and spectral are identical", {
  testthat::expect_identical(ci1$cis, ci2$cis)
})

# Check estimates of theta

test_that("estimates of theta agree, disjoint", {
  testthat::expect_equal(res$theta_dj, ci1$theta, tolerance = my_tol)
})

ci3 <- confint(res, maxima = "disjoint", type = "cholesky",
               interval_type = "both", conf_scale = "log")

which_rows <- c("N2015lik", "BB2018lik")
test_that("spm lik intervals don't depend on conf_scale", {
  testthat::expect_identical(ci1$cis[which_rows, ], ci3$cis[which_rows, ])
})

# ============================= plot.confint.spm =============================

# Check that plot.confint_spm works

cis <- confint(res, interval_type = "both")
ciplot <- plot(cis)
test_that("plot.confint_spm works, sliding", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = "BB2018", main = "BB2018 only")
test_that("plot.confint_spm works, sliding, BB2018 only, add title", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = c("N2015", "BB2018"),
               main = "N2015 and BB2018", legend = c("cool", "neat"))
test_that("plot.confint_spm works, sliding, 2 ests, user legend", {
  testthat::expect_identical(ciplot, NULL)
})

cis <- confint(res, interval_type = "both", maxima = "disjoint")
ciplot <- plot(cis, xlab = "my xlab", lwd = 2, col = "blue")
test_that("plot.confint_spm works, user plot args, disjoint", {
  testthat::expect_identical(ciplot, NULL)
})

# ================================== kgaps ===================================

context("confint.kgaps")

thresh <- quantile(newlyn, probs = 0.90)

res <- kgaps(newlyn, thresh)
res1 <- confint(res)
res2 <- confint(res, conf_scale = "log")
test_that("kgaps lik intervals don't depend on conf_scale", {
  testthat::expect_identical(res1["lik", ], res2["lik", ])
})

# Repeat for inc_cens = TRUE

res <- kgaps(newlyn, thresh, inc_cens = TRUE)
res1 <- confint(res)
res2 <- confint(res, conf_scale = "log")
test_that("kgaps lik intervals don't depend on conf_scale", {
  testthat::expect_identical(res1["lik", ], res2["lik", ])
})

