test_that("datatoMatrix creates correct matrix without hyperstate", {
  dat <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    chrom = c(5, 7, 10),
    prob = c(1, 1, 1)
  )
  mat <- datatoMatrix(x = dat, hyper = FALSE)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 3)
  expect_equal(ncol(mat), length(5:10))
  expect_equal(rownames(mat), c("sp1", "sp2", "sp3"))
})

test_that("datatoMatrix creates correct matrix with hyperstate", {
  dat <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    chrom = c(5, 7, 10),
    prob = c(0.8, 0.5, 1.0)
  )
  mat <- datatoMatrix(x = dat, hyper = TRUE)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 3)
  # With hyperstate, columns are doubled

  expect_equal(ncol(mat), 2 * length(5:10))
})

test_that("datatoMatrix respects range argument", {
  dat <- data.frame(
    species = c("sp1", "sp2"),
    chrom = c(5, 10),
    prob = c(1, 1)
  )
  mat <- datatoMatrix(x = dat, range = c(3, 12), hyper = FALSE)
  expect_equal(ncol(mat), length(3:12))
})

test_that("datatoMatrix respects buffer argument", {
  dat <- data.frame(
    species = c("sp1", "sp2"),
    chrom = c(5, 10),
    prob = c(1, 1)
  )
  mat <- datatoMatrix(x = dat, buffer = 2, hyper = FALSE)
  expect_equal(ncol(mat), length(3:12))
})
