test_that("startVals returns correct number of values", {
  vals <- startVals(n = 10, min = 0, max = 1, dist = "unif")
  expect_length(vals, 10)
})

test_that("startVals uniform values are within range", {
  vals <- startVals(n = 100, min = 0.1, max = 0.5, dist = "unif")
  expect_true(all(vals >= 0.1))
  expect_true(all(vals <= 0.5))
})

test_that("startVals normal distribution works", {
  vals <- startVals(n = 100, mean = 0.5, sd = 0.1, dist = "norm")
  expect_length(vals, 100)
  expect_true(is.numeric(vals))
})

test_that("startVals rejects invalid distribution", {
  expect_error(startVals(n = 5, min = 0, max = 1, dist = "gamma"))
})
