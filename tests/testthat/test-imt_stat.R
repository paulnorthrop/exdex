#context("kgaps_imt_stat")

# Check that the lengths of the vector of statistics returned from
# kgaps_imt_stat() equal n_kgaps

# Newlyn

u <- quantile(newlyn, probs = 0.9)
k <- 1

# inc_cens = TRUE
theta <- kgaps(newlyn, u = u, k = k)$theta
res <- kgaps_imt_stat(newlyn, theta = theta, u = u, k = k)
# Find all the lengths, except n_kgaps
lengs <- sapply(res, length)[-6]

test_that("newlyn: IMT stats lengths equal n_kgaps, inc_cens = TRUE", {
  testthat::expect_equal(res$n_kgaps, unique(lengs))
})

# inc_cens = FALSE
theta <- kgaps(newlyn, u = u, k = k, inc_cens = FALSE)$theta
res <- kgaps_imt_stat(newlyn, theta = theta, u = u, k = k, inc_cens = FALSE)
# Find all the lengths, except n_kgaps
lengs <- sapply(res, length)[-6]

test_that("newlyn: IMT stats lengths equal n_kgaps, inc_cens = FALSE", {
  testthat::expect_equal(res$n_kgaps, unique(lengs))
})

# Cheeseboro

# inc_cens = TRUE
u <- quantile(cheeseboro, probs = 0.9, na.rm = TRUE)
k <- 3
theta <- kgaps(cheeseboro, u = u, k = k)$theta
res <- kgaps_imt_stat(cheeseboro, theta = theta, u = u, k = k)
# Find all the lengths, except n_kgaps
lengs <- sapply(res, length)[-6]

test_that("cheeseboro: IMT stats lengths equal n_kgaps, inc_cens = TRUE", {
  testthat::expect_equal(res$n_kgaps, unique(lengs))
})

# inc_cens = FALSE
u <- quantile(cheeseboro, probs = 0.9, na.rm = TRUE)
k <- 3
theta <- kgaps(cheeseboro, u = u, k = k, inc_cens = FALSE)$theta
res <- kgaps_imt_stat(cheeseboro, theta = theta, u = u, k = k,
                      inc_cens = FALSE)
# Find all the lengths, except n_kgaps
lengs <- sapply(res, length)[-6]

test_that("cheeseboro: IMT stats lengths equal n_kgaps, inc_cens = FALSE", {
  testthat::expect_equal(res$n_kgaps, unique(lengs))
})
