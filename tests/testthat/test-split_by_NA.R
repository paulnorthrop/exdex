#context("Splitting by NAs")

# Check that if there are no missing values then the matrix is not modified

# Create a simple matrix
x <- matrix(1:50, 10, 5)
res <- split_by_NAs(x)

test_that("A matrix with no missings is not modified", {
  testthat::expect_equal(x, res, ignore_attr = TRUE)
})

# Check that splitting a matrix into non missing sequences by column gives
# the correct results

### Example 1

# Create a simple matrix
x <- matrix(1:50, 10, 5)
x[3, 1] <- NA
x[8, 1] <- NA
x[1:2, 3] <- NA
x[5, 3] <- NA
x[10, 3] <- NA
x[1:3, 4] <- NA
x[7:10, 5] <- NA

# Find newx by hand
newx <- matrix(NA, nrow = 10, ncol = 8)
newx[1:2, 1] <- 1:2
newx[1:4, 2] <- 4:7
newx[1:2, 3] <- 9:10
newx[, 4] <- 11:20
newx[1:2, 5] <- 23:24
newx[1:4, 6] <- 26:29
newx[1:7, 7] <- 34:40
newx[1:6, 8] <- 41:46

res <- split_by_NAs(x)

test_that("split_by_NAs is correct, example 1", {
  testthat::expect_equal(newx, res, ignore_attr = TRUE)
})

### Example 2

# Create a simple matrix
x2 <- matrix(1:27, 9, 3)
x2[1:2, 1] <- NA
x2[9, 1] <- NA
x2[1, 2] <- NA
x2[8:9, 2] <- NA
x2[1, 3] <- NA
x2[4, 3] <- NA

# Find newx by hand
newx2 <- matrix(NA, nrow = 6, ncol = 4)
newx2[1:6, 1] <- 3:8
newx2[1:6, 2] <- 11:16
newx2[1:2, 3] <- 20:21
newx2[1:5, 4] <- 23:27

res2 <- split_by_NAs(x2)

test_that("split_by_NAs is correct, example 2", {
  testthat::expect_equal(newx2, res2, ignore_attr = TRUE)
})

### Basic checks for the cheeseboro data

# By hand: the number of sequences of NAs for each year
# 2000: 5
# 2001: 4
# 2002: 10
# 2003: 5
# 2004: 1
# 2005: 2
# 2006: 2
# 2007: 1
# 2008: 1
# 2009: 1
# total: 32

res1 <- split_by_NAs(cheeseboro[, 1])
res2 <- split_by_NAs(cheeseboro[, 2])
res3 <- split_by_NAs(cheeseboro[, 3])
res4 <- split_by_NAs(cheeseboro[, 4])
res5 <- split_by_NAs(cheeseboro[, 5])
res6 <- split_by_NAs(cheeseboro[, 6])
res7 <- split_by_NAs(cheeseboro[, 7])
res8 <- split_by_NAs(cheeseboro[, 8])
res9 <- split_by_NAs(cheeseboro[, 9])
res10 <- split_by_NAs(cheeseboro[, 10])
restotal <- split_by_NAs(cheeseboro)

byhand <- c(5, 4, 10, 5, 1, 2, 2, 1, 1, 1, 32)
res <- c(ncol(res1), ncol(res2), ncol(res3), ncol(res4), ncol(res5),
         ncol(res6), ncol(res7), ncol(res8), ncol(res9), ncol(res10),
         ncol(restotal))

test_that("split_by_NAs is correct, cheeseboro", {
  testthat::expect_equal(byhand, res)
})

# Check that the number of non-missing values has been preserved

no_NA1 <- sum(!is.na(cheeseboro))
no_NA2 <- sum(!is.na(restotal))

test_that("split_by_NAs is correct, cheeseboro non-missings", {
  testthat::expect_equal(no_NA1, no_NA2)
})
