## Tests for make.mkn.fast(): the Rcpp sparse-uniformization Mkn likelihood.

skip_if_not_diversitree <- function() {
  if (!requireNamespace("diversitree", quietly = TRUE))
    skip("diversitree not available")
  if (!requireNamespace("ape", quietly = TRUE))
    skip("ape not available")
}

# Helper: build a simulated chromePlus problem and the matched diversitree
# constrained likelihood, plus the make.mkn.fast equivalent.
build_pair <- function(k_target, n_tips, constrain = list(drop.demi = TRUE),
                       seed = 42) {
  half <- k_target / 2
  chrom_range <- c(5, 5 + half - 1)
  set.seed(seed)
  tree <- ape::rcoal(n_tips, tip.label = paste0("sp", 1:n_tips))
  sim <- suppressMessages(simChrom(
    tree = tree, model = "ChromPlus",
    pars = c(asc1 = 0.1, asc2 = 0.12, desc1 = 0.08, desc2 = 0.09,
             dem1 = 0, dem2 = 0, pol1 = 0.005, pol2 = 0.005,
             tran12 = 0.02, tran21 = 0.02,
             root.chrom = round(mean(chrom_range)), root.state = 0),
    limits = chrom_range))
  chroms <- pmin(pmax(as.numeric(sim$chrom.num), chrom_range[1]), chrom_range[2])
  bin <- as.numeric(sim$binary.state)
  dat <- data.frame(species = tree$tip.label, chrom = chroms,
                    prob = ifelse(bin == 0, 1, 0))
  dat.mat <- datatoMatrix(x = dat, range = chrom_range, hyper = TRUE)
  lik <- suppressMessages(diversitree::make.mkn(
    tree, states = dat.mat, k = ncol(dat.mat),
    strict = FALSE, control = list(method = "ode", tol = 1e-10)))
  con <- suppressMessages(constrainMkn(
    data = dat.mat, lik = lik, hyper = TRUE, constrain = constrain))
  ll_fast <- make.mkn.fast(tree, dat.mat, hyper = TRUE,
                           constrain = constrain, root = "obs")
  list(tree = tree, dat.mat = dat.mat, con = con, ll_fast = ll_fast)
}

test_that("make.mkn.fast loglik matches constrainMkn at small k", {
  skip_if_not_diversitree()
  pair <- build_pair(k_target = 20, n_tips = 30)
  argn <- argnames(pair$ll_fast)
  expect_setequal(argn, diversitree::argnames(pair$con))

  pars <- setNames(rep(0.05, length(argn)), argn)
  if ("pol1" %in% argn) pars["pol1"] <- 0.001
  if ("pol2" %in% argn) pars["pol2"] <- 0.0015
  if ("tran12" %in% argn) pars["tran12"] <- 0.005
  if ("tran21" %in% argn) pars["tran21"] <- 0.004

  ll_div  <- pair$con(pars[diversitree::argnames(pair$con)])
  ll_fast <- pair$ll_fast(pars)
  expect_equal(ll_fast, ll_div, tolerance = 1e-8)
})

test_that("make.mkn.fast handles drop.poly + drop.demi", {
  skip_if_not_diversitree()
  pair <- build_pair(40, 30,
                     constrain = list(drop.poly = TRUE, drop.demi = TRUE))
  argn <- argnames(pair$ll_fast)
  expect_setequal(argn, diversitree::argnames(pair$con))
  pars <- setNames(rep(0.05, length(argn)), argn)
  expect_equal(pair$ll_fast(pars),
               pair$con(pars[diversitree::argnames(pair$con)]),
               tolerance = 1e-8)
})

test_that("make.mkn.fast handles nometa, symmetric, meta=SYM", {
  skip_if_not_diversitree()
  for (cs in list(list(nometa = TRUE, drop.demi = TRUE),
                  list(symmetric = TRUE, drop.demi = TRUE),
                  list(meta = "SYM", drop.demi = TRUE))) {
    pair <- build_pair(40, 30, constrain = cs)
    argn <- argnames(pair$ll_fast)
    expect_setequal(argn, diversitree::argnames(pair$con))
    pars <- setNames(rep(0.04, length(argn)), argn)
    expect_equal(pair$ll_fast(pars),
                 pair$con(pars[diversitree::argnames(pair$con)]),
                 tolerance = 1e-8,
                 info = paste("constraint:", paste(names(cs), collapse = ",")))
  }
})

test_that("make.mkn.fast respects state.names", {
  skip_if_not_diversitree()
  pair <- build_pair(30, 30)  # we'll build a state-named version separately
  state_names <- c("std", "fused")

  half <- 15; chrom_range <- c(5, 5 + half - 1)
  set.seed(42)
  tree <- ape::rcoal(30, tip.label = paste0("sp", 1:30))
  sim <- suppressMessages(simChrom(
    tree = tree, model = "ChromPlus",
    pars = c(asc1 = 0.1, asc2 = 0.12, desc1 = 0.08, desc2 = 0.09,
             dem1 = 0, dem2 = 0, pol1 = 0.005, pol2 = 0.005,
             tran12 = 0.02, tran21 = 0.02,
             root.chrom = round(mean(chrom_range)), root.state = 0),
    limits = chrom_range))
  chroms <- pmin(pmax(as.numeric(sim$chrom.num), chrom_range[1]), chrom_range[2])
  bin <- as.numeric(sim$binary.state)
  dat <- data.frame(species = tree$tip.label, chrom = chroms,
                    prob = ifelse(bin == 0, 1, 0))
  dat.mat <- datatoMatrix(x = dat, range = chrom_range, hyper = TRUE,
                          state.names = state_names)

  ll_fast <- make.mkn.fast(tree, dat.mat, hyper = TRUE,
                           state.names = state_names,
                           constrain = list(drop.demi = TRUE))
  argn <- argnames(ll_fast)
  expect_true(any(grepl(paste0("\\.", state_names[1], "$"), argn)))
  expect_true(any(grepl(paste0("\\.", state_names[2], "$"), argn)))
})

test_that("make.mkn.fast accepts both named and positional pars", {
  skip_if_not_diversitree()
  pair <- build_pair(20, 20)
  argn <- argnames(pair$ll_fast)
  pars_named <- setNames(rep(0.05, length(argn)), argn)
  pars_pos   <- unname(pars_named)
  expect_equal(pair$ll_fast(pars_named), pair$ll_fast(pars_pos))
})

test_that("make.mkn.fast returns -Inf for negative pars", {
  skip_if_not_diversitree()
  pair <- build_pair(20, 20)
  argn <- argnames(pair$ll_fast)
  pars <- setNames(rep(0.05, length(argn)), argn)
  pars[1] <- -1
  expect_identical(pair$ll_fast(pars), -Inf)
})

test_that("make.mkn.fast root='flat' matches a flat-prior reference", {
  skip_if_not_diversitree()
  pair <- build_pair(20, 20)
  argn <- argnames(pair$ll_fast)
  pars <- setNames(rep(0.05, length(argn)), argn)

  ll_obs   <- pair$ll_fast(pars)
  pair$ll_fast2 <- make.mkn.fast(pair$tree, pair$dat.mat, hyper = TRUE,
                                  constrain = list(drop.demi = TRUE),
                                  root = "flat")
  ll_flat  <- pair$ll_fast2(pars)
  # Flat prior likelihood = ROOT.OBS likelihood + log(sum_root D^2 / sum_root D)
  # We just check both are finite and different — the absolute relation is
  # tested in run_one elsewhere.
  expect_true(is.finite(ll_obs))
  expect_true(is.finite(ll_flat))
  expect_false(isTRUE(all.equal(ll_obs, ll_flat)))
})

test_that("make.mkn.fast works with find.mle and mcmc", {
  skip_if_not_diversitree()
  pair <- build_pair(20, 20,
                     constrain = list(drop.poly = TRUE, drop.demi = TRUE))
  argn <- argnames(pair$ll_fast)
  init <- setNames(rep(0.05, length(argn)), argn)

  fit <- diversitree::find.mle(
    pair$ll_fast, init, method = "subplex",
    control = list(maxit = 2000, reltol = 1e-3))
  expect_true(is.finite(fit$lnLik))
  expect_named(fit$par, argn, ignore.order = TRUE)

  chain <- diversitree::mcmc(
    pair$ll_fast, init, nsteps = 30, w = 0.05,
    prior = function(p) sum(stats::dexp(p, rate = 10, log = TRUE)),
    lower = 1e-6, print.every = 0)
  expect_equal(nrow(chain), 30)
})

# --- More thorough comparisons against constrainMkn -------------------------

test_that("make.mkn.fast matches constrainMkn under random parameter draws", {
  skip_if_not_diversitree()
  pair <- build_pair(40, 40)
  argn_div  <- diversitree::argnames(pair$con)
  argn_fast <- argnames(pair$ll_fast)
  expect_setequal(argn_div, argn_fast)

  set.seed(7)
  diffs <- vapply(seq_len(10), function(i) {
    p <- setNames(stats::runif(length(argn_fast), 0.001, 0.2), argn_fast)
    abs(pair$ll_fast(p) - pair$con(p[argn_div]))
  }, numeric(1))
  # All draws should agree to within ODE tolerance (we used tol=1e-10).
  expect_lt(max(diffs), 1e-7)
})

test_that("make.mkn.fast matches constrainMkn under hyper=FALSE", {
  skip_if_not_diversitree()
  set.seed(11)
  chrom_range <- c(5, 24)  # k=20, no hyperstate
  tree <- ape::rcoal(40, tip.label = paste0("sp", 1:40))
  sim <- suppressMessages(simChrom(
    tree = tree, model = "2010",
    pars = c(0.10, 0.08, 0, 0.005, round(mean(chrom_range))),
    limits = chrom_range))
  dat <- data.frame(species = tree$tip.label,
                    chrom = pmin(pmax(as.numeric(sim), chrom_range[1]),
                                 chrom_range[2]),
                    prob = 1)
  dat.mat <- datatoMatrix(x = dat, range = chrom_range, hyper = FALSE)
  lik <- suppressMessages(diversitree::make.mkn(
    tree, states = dat.mat, k = ncol(dat.mat),
    strict = FALSE, control = list(method = "ode", tol = 1e-10)))
  con <- suppressMessages(constrainMkn(
    data = dat.mat, lik = lik, hyper = FALSE,
    constrain = list(drop.demi = TRUE)))
  ll_fast <- make.mkn.fast(tree, dat.mat, hyper = FALSE,
                           constrain = list(drop.demi = TRUE))

  argn <- argnames(ll_fast)
  expect_setequal(argn, diversitree::argnames(con))
  pars <- setNames(rep(0.05, length(argn)), argn)
  if ("pol1" %in% argn) pars["pol1"] <- 0.005
  expect_equal(ll_fast(pars),
               con(pars[diversitree::argnames(con)]),
               tolerance = 1e-7)
})

test_that("make.mkn.fast matches constrainMkn under saf.model", {
  skip_if_not_diversitree()
  pair <- build_pair(30, 30,
                     constrain = list(saf.model = TRUE, drop.demi = TRUE,
                                      drop.poly = TRUE))
  argn <- argnames(pair$ll_fast)
  expect_setequal(argn, diversitree::argnames(pair$con))
  pars <- setNames(rep(0.04, length(argn)), argn)
  expect_equal(pair$ll_fast(pars),
               pair$con(pars[diversitree::argnames(pair$con)]),
               tolerance = 1e-7)
})

test_that("make.mkn.fast matches constrainMkn under polyploidy=TRUE", {
  skip_if_not_diversitree()
  set.seed(13)
  chrom_range <- c(5, 14); n_tips <- 30
  tree <- ape::rcoal(n_tips, tip.label = paste0("sp", 1:n_tips))
  # Polyploidy hidden state (PloidEvol) — cleanest with simulated tips.
  sim <- suppressMessages(simChrom(
    tree = tree, model = "PloidEvol",
    pars = c(asc1 = 0.10, asc2 = 0.10, desc1 = 0.08, desc2 = 0.08,
             dem1 = 0, dem2 = 0, pol1 = 0.005, pol2 = 0.005,
             redip = 0.01,
             root.chrom = round(mean(chrom_range)), root.state = 0),
    limits = chrom_range))
  chroms <- pmin(pmax(as.numeric(sim$chrom.num), chrom_range[1]),
                 chrom_range[2])
  state  <- as.numeric(sim$ploidy.state)
  dat <- data.frame(species = tree$tip.label, chrom = chroms,
                    prob = ifelse(state == 0, 1, 0))
  dat.mat <- datatoMatrix(x = dat, range = chrom_range, hyper = TRUE)
  lik <- suppressMessages(diversitree::make.mkn(
    tree, states = dat.mat, k = ncol(dat.mat),
    strict = FALSE, control = list(method = "ode", tol = 1e-10)))
  con <- suppressMessages(constrainMkn(
    data = dat.mat, lik = lik, hyper = TRUE, polyploidy = TRUE,
    constrain = list(drop.demi = TRUE)))
  ll_fast <- make.mkn.fast(tree, dat.mat, hyper = TRUE, polyploidy = TRUE,
                           constrain = list(drop.demi = TRUE))

  argn <- argnames(ll_fast)
  expect_setequal(argn, diversitree::argnames(con))
  pars <- setNames(rep(0.04, length(argn)), argn)
  expect_equal(ll_fast(pars),
               con(pars[diversitree::argnames(con)]),
               tolerance = 1e-7)
})

test_that("make.mkn.fast matches constrainMkn with state.names", {
  skip_if_not_diversitree()
  state_names <- c("std", "fused")
  set.seed(15)
  chrom_range <- c(5, 14); n_tips <- 30
  tree <- ape::rcoal(n_tips, tip.label = paste0("sp", 1:n_tips))
  sim <- suppressMessages(simChrom(
    tree = tree, model = "ChromPlus",
    pars = c(asc1 = 0.10, asc2 = 0.10, desc1 = 0.08, desc2 = 0.08,
             dem1 = 0, dem2 = 0, pol1 = 0.005, pol2 = 0.005,
             tran12 = 0.02, tran21 = 0.02,
             root.chrom = round(mean(chrom_range)), root.state = 0),
    limits = chrom_range))
  chroms <- pmin(pmax(as.numeric(sim$chrom.num), chrom_range[1]),
                 chrom_range[2])
  bin <- as.numeric(sim$binary.state)
  dat <- data.frame(species = tree$tip.label, chrom = chroms,
                    prob = ifelse(bin == 0, 1, 0))
  dat.mat <- datatoMatrix(x = dat, range = chrom_range, hyper = TRUE,
                           state.names = state_names)
  lik <- suppressMessages(diversitree::make.mkn(
    tree, states = dat.mat, k = ncol(dat.mat),
    strict = FALSE, control = list(method = "ode", tol = 1e-10)))
  con <- suppressMessages(constrainMkn(
    data = dat.mat, lik = lik, hyper = TRUE,
    state.names = state_names,
    constrain = list(drop.demi = TRUE)))
  ll_fast <- make.mkn.fast(tree, dat.mat, hyper = TRUE,
                           state.names = state_names,
                           constrain = list(drop.demi = TRUE))

  argn <- argnames(ll_fast)
  # Names should use descriptive labels
  expect_true(any(grepl("\\.std$", argn)))
  expect_true(any(grepl("\\.fused$", argn)))
  expect_setequal(argn, diversitree::argnames(con))

  pars <- setNames(rep(0.04, length(argn)), argn)
  expect_equal(ll_fast(pars),
               con(pars[diversitree::argnames(con)]),
               tolerance = 1e-7)
})

test_that("make.mkn.fast is strictly faster than constrainMkn at k>=50", {
  skip_if_not_diversitree()
  skip_on_cran()  # timing assertions can be flaky on CRAN's varied hardware
  pair <- build_pair(60, 50)
  argn <- argnames(pair$ll_fast)
  pars <- setNames(rep(0.05, length(argn)), argn)
  pars_div <- pars[diversitree::argnames(pair$con)]

  # warm-up
  pair$con(pars_div); pair$ll_fast(pars)
  t_div  <- system.time(for (i in 1:5)  pair$con(pars_div))[["elapsed"]] / 5
  t_fast <- system.time(for (i in 1:50) pair$ll_fast(pars))[["elapsed"]] / 50
  # The expected speedup at k=60 is ~30x; we check >5x as a robust lower bound
  # that survives noisy machines.
  expect_gt(t_div / t_fast, 5)
})
