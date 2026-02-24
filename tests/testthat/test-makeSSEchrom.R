test_that("makeSSEchrom returns expected structure", {
  skip_if_not_installed("diversitree")
  set.seed(42)
  result <- makeSSEchrom(
    h = 10, lambda1 = 0.3, lambda2 = 0.3,
    mu1 = 0.1, mu2 = 0.1,
    asc1 = 0.1, asc2 = 0.1,
    desc1 = 0.1, desc2 = 0.1,
    trans1 = 0.05, trans2 = 0.05,
    max.taxa = 20, x0 = 5
  )
  expect_type(result, "list")
  expect_named(result, c("tree", "binary.char", "haploid.num"))
  expect_s3_class(result$tree, "phylo")
})
