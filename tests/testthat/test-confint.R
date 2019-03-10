context("confint.exdex")

res <- spm(newlyn, 100)
ci1 <- confint(res, plot = TRUE, type = "cholesky")
ci2 <- confint(res, plot = TRUE, type = "spectral")

test_that("cholesky and spectral are identical", {
  testthat::expect_identical(ci1, ci2)
})

