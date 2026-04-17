## Sanity checks on the bundled beetleChrom dataset and a mini end-to-end
## workflow so the worked-example vignette can never silently break.

skip_if_not_installed("diversitree")
skip_if_not_installed("ape")

test_that("beetleChrom has expected structure and sizes", {
  data(beetleChrom, package = "chromePlus")
  expect_type(beetleChrom, "list")
  expect_named(beetleChrom, c("tree", "data", "metadata"))
  expect_s3_class(beetleChrom$tree, "phylo")
  expect_equal(ape::Ntip(beetleChrom$tree), 100)
  expect_setequal(beetleChrom$tree$tip.label, beetleChrom$data$species)
  expect_setequal(unique(beetleChrom$data$state), c("standard", "fused"))
  expect_true(all(beetleChrom$data$chrom >= 6))
  expect_true(all(beetleChrom$data$chrom <= 15))
})

test_that("beetleChrom runs the full datatoMatrix -> constrainMkn pipeline", {
  data(beetleChrom, package = "chromePlus")
  dat.mat <- datatoMatrix(
    x = beetleChrom$data[, c("species", "chrom", "prob")],
    hyper = TRUE,
    state.names = c("standard", "fused")
  )
  expect_equal(nrow(dat.mat), 100)
  expect_true(any(grepl("\\.fused$", colnames(dat.mat))))

  suppressMessages({
    lik <- diversitree::make.mkn(beetleChrom$tree,
                                 states = dat.mat,
                                 k = ncol(dat.mat),
                                 strict = FALSE,
                                 control = list(method = "ode"))
    con <- constrainMkn(data = dat.mat, lik = lik, hyper = TRUE,
                        state.names = c("standard", "fused"),
                        constrain = list(drop.demi = TRUE,
                                         drop.poly = TRUE))
  })
  args <- diversitree::argnames(con)
  expect_true(all(c("asc.standard", "asc.fused",
                    "tran.standard.to.fused",
                    "tran.fused.to.standard") %in% args))
  # likelihood evaluable at a reasonable starting point
  sv <- setNames(rep(0.05, length(args)), args)
  ll <- con(sv)
  expect_true(is.finite(ll))
})
