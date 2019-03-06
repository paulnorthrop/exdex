# Check that the functions disjoint_maxima(), all_disjoint_maxima(),
# sliding_maxima() and all_maxima() agree and give the correct results
# on some simple examples

my_tol <- 1e-10

# ---------------------------- all_disjoint_maxima() ------------------------ #

context("disjoint_maxima")

x <- 1:9
temp <- all_disjoint_maxima(x, b = 3)
index <- c(3, 6, 9)
y_mat <- matrix(x[index], 3, 1)
x_mat <- cbind(1:9)
test_that("x = 1:9, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, y_mat, tol = my_tol)
})
test_that("x = 1:9, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, x_mat, tol = my_tol)
})
temp2 <- disjoint_maxima(x, b = 3)
test_that("x = 1:9, b = 3, disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, as.matrix(temp2$y), tol = my_tol)
})
test_that("x = 1:9, b = 3, disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, as.matrix(temp2$x), tol = my_tol)
})

x <- 1:11
temp <- all_disjoint_maxima(x, b = 3)
index <- c(3, 6, 9, 4, 7, 10, 5, 8, 11)
y_mat <- matrix(x[index], 3, 3)
x_mat <- cbind(1:9, 2:10, 3:11)
test_that("x = 1:11, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, y_mat, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, x_mat, tol = my_tol)
})
temp2 <- disjoint_maxima(x, b = 3)
test_that("x = 1:11, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat[, 1], temp2$y, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat[, 1], temp2$x, tol = my_tol)
})

context("all_maxima()")

# Check that the function all_maxima() gives the same results as the
# functions all_disjoint_maxima(), disjoint_maxima() and sliding_maxima(),
# using the newlyn data

# Sliding maxima and disjoint maxima
a_res <- all_maxima(newlyn, 100)
# All disjoint maxima
d_res <- all_disjoint_maxima(newlyn, 100)
# Sliding maxima
s_res <- sliding_maxima(newlyn, 100)
# Disjoint maxima, starting only from the first observation
d1_res <- disjoint_maxima(newlyn, 100)

test_that("newlyn: disjoint contributing values", {
  testthat::expect_identical(a_res$xd, d_res$x_mat)
})
test_that("newlyn: disjoint maxima", {
  testthat::expect_identical(a_res$yd, d_res$y_mat)
})
test_that("newlyn: sliding contributing values", {
  testthat::expect_identical(a_res$xs, s_res$x)
})
test_that("newlyn: sliding maxima", {
  testthat::expect_identical(a_res$ys, s_res$y)
})
test_that("newlyn: disjoint contributing values vs disjoint_maxima()", {
  testthat::expect_identical(a_res$xd[, 1], d1_res$x)
})
test_that("newlyn: disjoint maxima vs disjoint_maxima()", {
  testthat::expect_identical(a_res$yd[, 1], d1_res$y)
})
