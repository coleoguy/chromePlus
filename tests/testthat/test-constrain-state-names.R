## Regression tests for the state.names -> asc.<name> naming contract.
## These lock in the convention that state.names[1] maps to the first half
## of the data matrix / Q-matrix (the half that gets the probability from
## column 3 of datatoMatrix), and state.names[2] maps to the second half.

skip_if_not_installed("diversitree")
skip_if_not_installed("ape")

make_fit_inputs <- function(state.names = NULL) {
  dat <- data.frame(
    species = paste0("sp", 1:5),
    chrom = c(5, 6, 7, 8, 10),
    prob = c(0.8, 0.6, 0.9, 0.5, 1.0)
  )
  dat.mat <- datatoMatrix(x = dat, hyper = TRUE, state.names = state.names)
  tree <- ape::rcoal(5, tip.label = dat$species)
  list(dat = dat, dat.mat = dat.mat, tree = tree)
}

test_that("datatoMatrix state.names labels the second-half columns", {
  inp <- make_fit_inputs(state.names = c("matched", "mismatched"))
  cn <- colnames(inp$dat.mat)
  split <- ncol(inp$dat.mat) / 2
  # first half is unlabeled chromosome numbers, second half carries the suffix
  expect_equal(cn[1:split], as.character(5:10))
  expect_true(all(grepl("\\.mismatched$", cn[(split + 1):length(cn)])))
  # prob from column 3 lands in the first half (state.names[1] = matched)
  expect_equal(inp$dat.mat[1, "5"], 0.8)
  expect_equal(inp$dat.mat[1, "5.mismatched"], 0.2)
})

test_that("datatoMatrix rejects bad state.names", {
  dat <- data.frame(species = "sp1", chrom = 5, prob = 1)
  expect_error(
    datatoMatrix(x = dat, hyper = TRUE, state.names = "just_one"),
    "state.names must be a character vector of length 2"
  )
})

## mkn needs deterministic tip states; build a matrix with certain tips.
make_mkn_inputs <- function(state.names = NULL) {
  dat <- data.frame(
    species = paste0("sp", 1:5),
    chrom = c(5, 6, 7, 8, 10),
    prob = c(1, 0, 1, 0, 1)  # certain: sp1/3/5 in state1, sp2/4 in state2
  )
  dat.mat <- datatoMatrix(x = dat, hyper = TRUE, state.names = state.names)
  tree <- ape::rcoal(5, tip.label = dat$species)
  list(dat = dat, dat.mat = dat.mat, tree = tree)
}

test_that("constrainMkn emits asc.<name> parameters when state.names given", {
  inp <- make_mkn_inputs(state.names = c("matched", "mismatched"))
  suppressMessages({
    lik <- diversitree::make.mkn(inp$tree, states = inp$dat.mat,
                                 k = ncol(inp$dat.mat),
                                 strict = FALSE,
                                 control = list(method = "ode"))
    con <- constrainMkn(data = inp$dat.mat, lik = lik, hyper = TRUE,
                        state.names = c("matched", "mismatched"),
                        constrain = list(drop.demi = TRUE, drop.poly = TRUE))
  })
  args <- diversitree::argnames(con)
  expect_true("asc.matched"   %in% args)
  expect_true("asc.mismatched" %in% args)
  expect_false("asc1" %in% args)
  expect_false("asc2" %in% args)
  expect_true("tran.matched.to.mismatched" %in% args)
})

test_that("constrainMkn defaults to asc1/asc2 and warns when state.names is NULL", {
  inp <- make_mkn_inputs(state.names = NULL)
  suppressMessages({
    lik <- diversitree::make.mkn(inp$tree, states = inp$dat.mat,
                                 k = ncol(inp$dat.mat),
                                 strict = FALSE,
                                 control = list(method = "ode"))
  })
  expect_message(
    con <- constrainMkn(data = inp$dat.mat, lik = lik, hyper = TRUE,
                        constrain = list(drop.demi = TRUE, drop.poly = TRUE)),
    "state.names not provided"
  )
  args <- diversitree::argnames(con)
  expect_true("asc1" %in% args)
  expect_true("asc2" %in% args)
})

test_that("constrainMuSSE emits asc.<name> parameters when state.names given", {
  inp <- make_fit_inputs(state.names = c("matched", "mismatched"))
  suppressMessages({
    lik <- diversitree::make.musse(inp$tree, states = inp$dat.mat,
                                   k = ncol(inp$dat.mat),
                                   strict = FALSE,
                                   control = list(method = "ode"))
    con <- constrainMuSSE(data = inp$dat.mat, lik = lik,
                          s.lambda = FALSE, s.mu = FALSE,
                          state.names = c("matched", "mismatched"),
                          constrain = list(drop.demi = TRUE, drop.poly = TRUE))
  })
  args <- diversitree::argnames(con)
  expect_true("asc.matched"    %in% args)
  expect_true("asc.mismatched" %in% args)
  expect_true("lambda.matched" %in% args)
  expect_true("mu.mismatched"  %in% args)
  expect_false("asc1" %in% args)
  expect_false("lambda1" %in% args)
})
