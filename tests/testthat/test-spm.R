context("spm vs spm_check")

# Check that spm(), faster but not very transparent, gives the same results
# as spm_check(), slower but more transparent

# 8 cases:
# bias_adjust in c("BB3", "BB1", "N", "none")
bias_adjust_vec <- c("BB3", "BB1", "N", "none")
# which_dj in c("last", "first")
which_dj_vec <- c("last", "first")
# block size: pick a big one so that the tests aren't slow
# The permitted range of b for these data is 15 - 196
# We must respect this here because spm_check() doesn't check b
# It also doesn't check the format of the data
b <- 180
# Tolerance
my_tol <- 1e-5

for (i in 1:4){
  for (j in 1:2) {
    res <- spm(newlyn, b = b,
               bias_adjust = bias_adjust_vec[i],
               which_dj = which_dj_vec[j])
    res_sl <- spm_check(newlyn, b = b, sliding = TRUE,
                        bias_adjust = bias_adjust_vec[i],
                        which_dj = which_dj_vec[j])
    res_dj <- spm_check(newlyn, b = b, sliding = FALSE,
                        bias_adjust = bias_adjust_vec[i],
                        which_dj = which_dj_vec[j])
    my_text <- paste(bias_adjust_vec[i], which_dj_vec[j])
    test_that(paste(my_text, "sliding, theta"), {
      testthat::expect_equal(res$theta_sl, res_sl$theta, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, theta"), {
      testthat::expect_equal(res$theta_dj, res_dj$theta, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, se"), {
      testthat::expect_equal(res$se_sl, res_sl$se, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, se"), {
      testthat::expect_equal(res$se_dj, res_dj$se, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, bias"), {
      testthat::expect_equal(res$bias_sl, res_sl$bias_val, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, bias"), {
      testthat::expect_equal(res$bias_dj, res_dj$bias_val, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_sl"), {
      testthat::expect_equal(summary(res$data_sl),
                             summary(cbind(N2015 = res_sl$N2015_data,
                                           BB2018 = res_sl$BB2018_data)),
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_dj"), {
      testthat::expect_equal(summary(res$data_dj),
                             summary(cbind(N2015 = res_dj$N2015_data,
                                           BB2018 = res_dj$BB2018_data)),
                             tolerance = my_tol)
    })
  }
}

context("spm when b is too low or too high")

# Check that the results are as expected when b is too low or too high.
# In these cases:
#   (a) the estimated SEs are missing
#   (b) if bias_adjust = "BB3" then bias_adjust is changed to "BB1"

# The permitted range of b for these data is 15 - 196

# b too low
b_low <- 14
res1 <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "BB1"))
res3 <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "BB3"))

test_that(paste("b low, bias_dj"), {
  testthat::expect_equal(res1$bias_dj, res3$bias_dj, tolerance = my_tol)
})
test_that(paste("b low, bias_sl"), {
  testthat::expect_equal(res1$bias_sl, res3$bias_sl, tolerance = my_tol)
})
test_that(paste("b low, se_dj"), {
  testthat::expect_identical(res1$se_dj, c(N2015 = NA, BB2018 = NA))
})
test_that(paste("b low, se_sl"), {
  testthat::expect_identical(res1$se_sl, c(N2015 = NA, BB2018 = NA))
})

# b too high
b_low <- 200
res1 <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "BB1"))
res3 <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "BB3"))

test_that(paste("b high, bias_dj"), {
  testthat::expect_equal(res1$bias_dj, res3$bias_dj, tolerance = my_tol)
})
test_that(paste("b high, bias_sl"), {
  testthat::expect_equal(res1$bias_sl, res3$bias_sl, tolerance = my_tol)
})
test_that(paste("b high, se_dj"), {
  testthat::expect_identical(res1$se_dj, c(N2015 = NA, BB2018 = NA))
})
test_that(paste("b high, se_sl"), {
  testthat::expect_identical(res1$se_sl, c(N2015 = NA, BB2018 = NA))
})

context("spm: equivalence of BB2018 when bias_adjust = ''BB1'' and ''N''")

b <- 100
resBB1 <- spm(newlyn, b = b, bias_adjust = "BB1")
resN <- spm(newlyn, b = b, bias_adjust = "N")

test_that(paste("BB1 vs N, b is OK"), {
  testthat::expect_equal(resBB1$theta_dj["BB2018"],
                             resN$theta_dj["BB2018"], tolerance = my_tol)
})

resBB1 <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "BB1"))
resN <- suppressWarnings(spm(newlyn, b = b_low, bias_adjust = "N"))

test_that(paste("BB1 vs N, b is too low"), {
  testthat::expect_equal(resBB1$theta_dj["BB2018"],
                             resN$theta_dj["BB2018"], tolerance = my_tol)
})
