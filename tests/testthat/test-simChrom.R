test_that("simChrom works with 2010 model", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  pars <- c(0.1, 0.1, 0.0, 0.01, 5)
  result <- simChrom(tree = tree, pars = pars,
                     limits = c(3, 12), model = "2010")
  expect_true(is.numeric(result))
  expect_equal(length(result), 20)
})

test_that("simChrom works with ChromPlus model", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  pars <- c(0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.01, 0.01, 0.05, 0.05, 5, 0)
  result <- simChrom(tree = tree, pars = pars,
                     limits = c(3, 12), model = "ChromPlus")
  expect_type(result, "list")
  expect_named(result, c("binary.state", "chrom.num"))
})

test_that("simChrom rejects incorrect pars length", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  expect_error(simChrom(tree = tree, pars = c(0.1, 0.1),
                        limits = c(3, 12), model = "2010"))
})
