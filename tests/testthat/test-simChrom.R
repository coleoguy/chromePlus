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

test_that("simChrom state.names labels binary states in ChromPlus", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  pars <- c(0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.01, 0.01, 0.05, 0.05, 5, 0)
  result <- simChrom(tree = tree, pars = pars,
                     limits = c(3, 12), model = "ChromPlus",
                     state.names = c("matched", "mismatched"))
  expect_type(result, "list")
  expect_named(result, c("binary.state", "chrom.num"))
  expect_true(all(result$binary.state %in% c("matched", "mismatched")))
})

test_that("simChrom accepts named pars for ChromPlus", {
  skip_if_not_installed("diversitree")
  set.seed(123)
  tree <- ape::rcoal(20)
  pars_named <- c(asc1 = 0.1, asc2 = 0.1, desc1 = 0.1, desc2 = 0.1,
                   dem1 = 0.0, dem2 = 0.0, pol1 = 0.01, pol2 = 0.01,
                   tran12 = 0.05, tran21 = 0.05,
                   root.chrom = 5, root.state = 0)
  result <- simChrom(tree = tree, pars = pars_named,
                     limits = c(3, 12), model = "ChromPlus")
  expect_type(result, "list")
  expect_named(result, c("binary.state", "chrom.num"))
})

test_that("simChrom accepts labeled named pars with state.names", {
  skip_if_not_installed("diversitree")
  set.seed(456)
  tree <- ape::rcoal(20)
  pars_labeled <- c(asc.matched = 0.1, asc.mismatched = 0.1,
                     desc.matched = 0.1, desc.mismatched = 0.1,
                     dem.matched = 0.0, dem.mismatched = 0.0,
                     pol.matched = 0.01, pol.mismatched = 0.01,
                     tran.matched.to.mismatched = 0.05,
                     tran.mismatched.to.matched = 0.05,
                     root.chrom = 5, root.state = 0)
  result <- simChrom(tree = tree, pars = pars_labeled,
                     limits = c(3, 12), model = "ChromPlus",
                     state.names = c("matched", "mismatched"))
  expect_type(result, "list")
  expect_true(all(result$binary.state %in% c("matched", "mismatched")))
})

test_that("simChrom state.names validates input", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  pars <- c(0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.01, 0.01, 0.05, 0.05, 5, 0)
  expect_error(
    simChrom(tree = tree, pars = pars, limits = c(3, 12),
             model = "ChromPlus", state.names = "only_one"),
    "state.names must be a character vector of length 2"
  )
})

test_that("simChrom named pars for 2010 model", {
  skip_if_not_installed("diversitree")
  tree <- ape::rcoal(20)
  pars_named <- c(asc = 0.1, desc = 0.1, dem = 0.0, pol = 0.01, root.chrom = 5)
  result <- simChrom(tree = tree, pars = pars_named,
                     limits = c(3, 12), model = "2010")
  expect_true(is.numeric(result))
  expect_equal(length(result), 20)
})
