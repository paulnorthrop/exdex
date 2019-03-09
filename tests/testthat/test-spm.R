context("spm vs spm_slow")

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
                             summary(cbind(res_sl$N2015_data,
                                           res_sl$BB2018_data)),
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_dj"), {
      testthat::expect_equal(summary(res$data_dj),
                             summary(cbind(res_dj$N2015_data,
                                           res_dj$BB2018_data)),
                             tolerance = my_tol)
    })
  }
}

