#context("dgaps_imt")

# Check that calling dgaps_imt() with vector arguments u and k gives
# the same results as calling dgaps_imt() repeatedly with scalar arguments

# For inc_cens = FALSE -----
inc_cens <- FALSE

u <- stats::quantile(newlyn, probs = c(0.85, 0.90, 0.95))
d_vals <- 1:4
all_res <- dgaps_imt(newlyn, u, d_vals, inc_cens = inc_cens)
all_IMT <- all_res$imt
all_p <- all_res$p
all_theta <- all_res$theta

ind_IMT <- ind_p <- ind_theta <- all_IMT
for (i in 1:length(u)) {
  for (j in 1:length(d_vals)) {
    temp <- dgaps_imt(newlyn, u = u[i], D = d_vals[j], inc_cens = inc_cens)
    ind_IMT[i, j] <- temp$imt
    ind_p[i, j] <- temp$p
    ind_theta[i, j] <- temp$theta
  }
}

my_tol <- 1e-5

test_that("IMT values agree", {
  testthat::expect_equal(all_IMT, ind_IMT, tolerance = my_tol)
})
test_that("p-values agree", {
  testthat::expect_equal(all_p, ind_p, tolerance = my_tol)
})
test_that("MLEs of theta values agree", {
  testthat::expect_equal(all_theta, ind_theta, tolerance = my_tol)
})

# For inc_cens = TRUE -----
inc_cens <- TRUE

u <- stats::quantile(newlyn, probs = c(0.85, 0.90, 0.95))
d_vals <- 1:4
all_res <- dgaps_imt(newlyn, u, d_vals, inc_cens = inc_cens)
all_IMT <- all_res$imt
all_p <- all_res$p
all_theta <- all_res$theta

ind_IMT <- ind_p <- ind_theta <- all_IMT
for (i in 1:length(u)) {
  for (j in 1:length(d_vals)) {
    temp <- dgaps_imt(newlyn, u = u[i], D = d_vals[j], inc_cens = inc_cens)
    ind_IMT[i, j] <- temp$imt
    ind_p[i, j] <- temp$p
    ind_theta[i, j] <- temp$theta
  }
}

my_tol <- 1e-5

test_that("IMT values agree", {
  testthat::expect_equal(all_IMT, ind_IMT, tolerance = my_tol)
})
test_that("p-values agree", {
  testthat::expect_equal(all_p, ind_p, tolerance = my_tol)
})
test_that("MLEs of theta values agree", {
  testthat::expect_equal(all_theta, ind_theta, tolerance = my_tol)
})
