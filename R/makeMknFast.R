#' Fast Mkn Likelihood for Large Chromosome State Spaces
#'
#' Builds a likelihood function equivalent to `make.mkn()` followed by
#' [constrainMkn()], but evaluates the per-call likelihood with a sparse
#' uniformization-based Felsenstein pruning routine implemented in C++. At
#' large `k` (the number of states), this is typically 30-100x faster than
#' diversitree's ODE pruning, while producing log-likelihoods that agree to
#' near machine precision.
#'
#' The returned function has the same calling convention as a constrained
#' diversitree likelihood: a numeric `pars` vector whose names match
#' `argnames(ll)`. It works directly with `diversitree::find.mle()` and
#' `diversitree::mcmc()`.
#'
#' This function is most useful when the chromosome state space is large
#' (`k >= 100`) and the per-evaluation cost in diversitree is dominated by
#' branch-wise transition-probability computation. For small `k`, the
#' diversitree path is already fast and the speedup margin is modest.
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param data A probability matrix as produced by [datatoMatrix()]. Rows are
#'   species (with `rownames` matching `tree$tip.label`), columns are
#'   chromosome states (and optionally hyperstate columns).
#' @param hyper Logical. If `TRUE` (default), includes a binary hyperstate
#'   allowing different rates of chromosome evolution in each state.
#' @param polyploidy Logical. If `TRUE`, the hyperstate represents ploidy
#'   level. Defaults to `FALSE`.
#' @param oneway Logical. If `TRUE`, the transition rate from state 2 back to
#'   state 1 is set to zero. Defaults to `FALSE`.
#' @param state.names Optional character vector of length 2 giving descriptive
#'   names for the two binary states. When provided, all parameter names use
#'   these labels (e.g. `asc.matched` instead of `asc1`).
#' @param constrain A list of additional model constraints. Same options as
#'   [constrainMkn()]: `drop.poly`, `drop.demi`, `symmetric`, `nometa`,
#'   `meta`, `saf.model`, `sym.hyperstates`.
#' @param root Character. Root prior to use when computing the likelihood.
#'   One of `"obs"` (default; use `D / sum(D)` at root, matching diversitree's
#'   `ROOT.OBS`), `"flat"` (uniform), or `"given"` (use the `root.p` argument
#'   passed to the returned function).
#'
#' @return A function `ll(pars, root.p = NULL)` that computes the
#'   log-likelihood. The returned function has class `"mkn.fast"`,
#'   `"function"` and an `"argnames"` attribute giving the free-parameter
#'   names. Use `argnames(ll)` to retrieve the parameter order.
#'
#' @seealso [constrainMkn()] for the diversitree-backed equivalent.
#'
#' @examples
#' \donttest{
#' library(diversitree); library(ape)
#' set.seed(1)
#' tree <- rcoal(50, tip.label = paste0("sp", 1:50))
#' sim <- simChrom(tree = tree, model = "ChromPlus",
#'   pars = c(asc1=0.1, asc2=0.1, desc1=0.08, desc2=0.08,
#'            dem1=0, dem2=0, pol1=0.005, pol2=0.005,
#'            tran12=0.02, tran21=0.02,
#'            root.chrom = 12, root.state = 0),
#'   limits = c(5, 20))
#' dat <- data.frame(species = tree$tip.label,
#'                   chrom = as.numeric(sim$chrom.num),
#'                   prob = as.numeric(sim$binary.state == 0))
#' dat.mat <- datatoMatrix(x = dat, range = c(5, 20), hyper = TRUE)
#' ll <- make.mkn.fast(tree, dat.mat, hyper = TRUE,
#'                    constrain = list(drop.demi = TRUE))
#' argnames(ll)
#' pars <- setNames(rep(0.05, length(argnames(ll))), argnames(ll))
#' ll(pars)
#' fit <- find.mle(ll, pars, method = "subplex")
#' }
#'
#' @export
make.mkn.fast <- function(tree, data,
                          hyper = TRUE,
                          polyploidy = FALSE,
                          oneway = FALSE,
                          state.names = NULL,
                          constrain = list(),
                          root = c("obs", "flat", "given")) {

  root <- match.arg(root)
  root_mode <- switch(root, flat = 0L, obs = 1L, given = 2L)

  if (!inherits(tree, "phylo"))
    stop("'tree' must be of class \"phylo\"")
  if (!is.matrix(data) || is.null(colnames(data)) || is.null(rownames(data)))
    stop("'data' must be a matrix with row- and colnames (use datatoMatrix())")
  if (!setequal(rownames(data), tree$tip.label))
    stop("rownames(data) must match tree$tip.label")

  # Build the rate-name parMat using the same logic as constrainMkn().
  built <- .chromeplus_build_parMat(
    data = data, hyper = hyper, polyploidy = polyploidy,
    oneway = oneway, state.names = state.names, constrain = constrain,
    emit_messages = FALSE
  )
  parMat <- built$parMat

  spec <- .chromeplus_parse_parMat(parMat)
  k <- spec$k

  # Prepare tree structure (postorder, 0-indexed).
  ord <- stats::reorder(tree, "postorder")
  edge_parent <- as.integer(ord$edge[, 1] - 1L)
  edge_child  <- as.integer(ord$edge[, 2] - 1L)
  edge_length <- as.numeric(ord$edge.length)
  n_tips      <- length(ord$tip.label)
  n_internal  <- as.integer(ord$Nnode)

  # Tip-state matrix in tree's tip order.
  tip_states <- data[ord$tip.label, , drop = FALSE]
  storage.mode(tip_states) <- "double"

  free_pars <- spec$free_pars

  ll <- function(pars, root.p = NULL) {
    if (!is.null(names(pars))) {
      missing <- setdiff(free_pars, names(pars))
      if (length(missing))
        stop("missing pars: ", paste(missing, collapse = ", "))
      pars <- pars[free_pars]
    } else {
      if (length(pars) != length(free_pars))
        stop("pars length must equal length(argnames(ll)) = ",
             length(free_pars))
    }
    if (any(pars < 0)) return(-Inf)

    qx <- as.numeric(pars[spec$par_id]) * spec$coef
    rp <- if (is.null(root.p)) numeric(0) else as.numeric(root.p)
    if (root_mode == 2L && length(rp) != k)
      stop("root.p must have length k = ", k)

    mkn_loglik_sparse_cpp(
      k = k, n_tips = n_tips, n_internal = n_internal,
      edge_parent = edge_parent, edge_child = edge_child,
      edge_length = edge_length,
      qi = spec$qi, qj = spec$qj, qx = qx,
      tip_states = tip_states,
      root_mode = root_mode, root_prior = rp
    )
  }
  attr(ll, "argnames") <- free_pars
  class(ll) <- c("mkn.fast", "function")
  ll
}

#' @export
#' @method argnames mkn.fast
#' @importFrom diversitree argnames
argnames.mkn.fast <- function(x, ...) attr(x, "argnames")

# ---- Internal: parse a parMat with rate-name entries into sparse triples ---

.chromeplus_parse_parMat <- function(parMat) {
  k <- nrow(parMat)
  rmat <- row(parMat); cmat <- col(parMat)
  is_off <- rmat != cmat
  vals   <- as.character(parMat)
  keep   <- is_off & !is.na(vals) & vals != "0" & vals != ""
  qi <- as.integer(rmat[keep] - 1L)
  qj <- as.integer(cmat[keep] - 1L)
  v  <- vals[keep]
  is_half <- startsWith(v, ".5*")
  base    <- ifelse(is_half, substring(v, 4), v)
  coef    <- ifelse(is_half, 0.5, 1.0)
  free_pars <- sort(unique(base))
  par_id <- match(base, free_pars)
  list(qi = qi, qj = qj, par_id = par_id, coef = coef,
       free_pars = free_pars, k = as.integer(k))
}

# ---- Internal: build the parMat + rate.table for chromePlus mkn models -----
#
# Used by both constrainMkn() and make.mkn.fast(). Returns the rate-name
# parMat (k x k character matrix), the rate.table (data frame of i, j, rate
# name), and emits the same informational messages as constrainMkn() about
# which model is being constructed.
#
# `emit_messages = FALSE` suppresses the "Constraining model..." messages
# (used by make.mkn.fast() to keep its build silent).

.chromeplus_build_parMat <- function(data, hyper, polyploidy, oneway,
                                     state.names, constrain,
                                     emit_messages = TRUE) {

  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2)
      stop("state.names must be a character vector of length 2")
  }

  if (is.null(constrain$drop.poly))      constrain$drop.poly      <- FALSE
  if (is.null(constrain$drop.demi))      constrain$drop.demi      <- FALSE
  if (is.null(constrain$symmetric))      constrain$symmetric      <- FALSE
  if (is.null(constrain$nometa))         constrain$nometa         <- FALSE
  if (is.null(constrain$saf.model))      constrain$saf.model      <- FALSE
  if (is.null(constrain$sym.hyperstates))constrain$sym.hyperstates<- FALSE
  if (is.null(constrain$meta))           constrain$meta           <- "ARD"

  pad <- nchar(as.character(ncol(data)))
  parMat <- matrix(0, ncol(data), ncol(data))
  colnames(parMat) <- sprintf(paste0("%0", pad, "d"), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  split <- ncol(parMat) / 2

  if (hyper) chroms <- as.numeric(colnames(data)[1:split])
  else       chroms <- as.numeric(colnames(data))

  if (!hyper) {
    if (emit_messages)
      message("Constraining model to simple chromevol version")
    for (i in 1:(nrow(parMat) - 1)) {
      if ((chroms[i] * 2) <= max(chroms))
        parMat[i, which(chroms == (chroms[i] * 2))] <- 5
      if (!constrain$drop.demi) {
        if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
          x <- chroms[i] * 1.5
          if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10
          if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
        }
      }
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
  }

  if (hyper && polyploidy) {
    if (emit_messages)
      message("Creating rate matrix for chosen chromosome model")
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms))
        parMat[i, (which(chroms[i] * 2 == chroms) + split)] <- 5
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 10
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 11
      }
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms))
        parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13
      }
      parMat[i, (i - split)] <- 7
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 7
      parMat[i, (i + 1)] <- 3
      parMat[(i + 1), i] <- 4
    }
  }

  if (hyper && !polyploidy && !constrain$saf.model) {
    if (emit_messages)
      message("Creating rate matrix for chosen chromosome model")
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms))
        parMat[i, which(chroms == (chroms[i] * 2))] <- 5
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
      }
      parMat[i, (i + split)] <- 8
      if (i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- 8
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms))
        parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13
      }
      parMat[i, (i - split)] <- 9
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 9
      parMat[i, (i + 1)] <- 3
      parMat[(i + 1), i] <- 4
    }
  }

  if (hyper && !polyploidy && constrain$saf.model) {
    if (emit_messages)
      message("Creating rate matrix for chosen chromosome model")
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms))
        parMat[i, which(chroms == (chroms[i] * 2))] <- 5
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
      }
      if (i > 1) parMat[i, (i + split - 1)] <- 14
      if (i == (split - 1)) parMat[(i + 1), (i + split)] <- 14
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms))
        parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13
      }
      parMat[i, (i - split)] <- 15
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 15
      parMat[i, (i + 1)] <- 3
      parMat[(i + 1), i] <- 4
    }
  }

  rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat), 3))
  rate.table[, 1] <- rep(as.character(rownames(parMat)), each = ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[, 2], ]

  rate.table[rate.table[, 3] == 1, 3]  <- "asc1"
  rate.table[rate.table[, 3] == 2, 3]  <- "desc1"
  rate.table[rate.table[, 3] == 3, 3]  <- "asc2"
  rate.table[rate.table[, 3] == 4, 3]  <- "desc2"
  rate.table[rate.table[, 3] == 5, 3]  <- "pol1"
  rate.table[rate.table[, 3] == 6, 3]  <- "pol2"
  rate.table[rate.table[, 3] == 7, 3]  <- "redip"
  rate.table[rate.table[, 3] == 8, 3]  <- "tran12"
  rate.table[rate.table[, 3] == 9, 3]  <- "tran21"
  rate.table[rate.table[, 3] == 10, 3] <- "dem1"
  rate.table[rate.table[, 3] == 11, 3] <- ".5*dem1"
  rate.table[rate.table[, 3] == 12, 3] <- "dem2"
  rate.table[rate.table[, 3] == 13, 3] <- ".5*dem2"
  rate.table[rate.table[, 3] == 14, 3] <- "tranSAF"
  rate.table[rate.table[, 3] == 15, 3] <- "tranRo"

  if (constrain$nometa) {
    rate.table[rate.table[, 3] == "asc2", 3]    <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3]   <- "desc1"
    rate.table[rate.table[, 3] == "pol2", 3]    <- "pol1"
    rate.table[rate.table[, 3] == "dem2", 3]    <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }

  if (constrain$drop.poly) {
    rate.table[rate.table[, 3] == "pol1", 3] <- "0"
    rate.table[rate.table[, 3] == "pol2", 3] <- "0"
  }

  if (constrain$drop.demi) {
    rate.table[rate.table[, 3] == "dem1", 3]    <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3]    <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }

  if (constrain$symmetric) {
    rate.table[rate.table[, 3] == "desc1", 3]   <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3]   <- "asc2"
    rate.table[rate.table[, 3] == "pol2", 3]    <- "pol1"
    rate.table[rate.table[, 3] == "dem2", 3]    <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }

  if (constrain$meta == "SYM") {
    rate.table[rate.table[, 3] == "tran21", 3] <- "tran12"
  }

  if (constrain$saf.model) {
    rate.table[rate.table[, 3] == "pol1", 3]    <- "0"
    rate.table[rate.table[, 3] == "pol2", 3]    <- "0"
    rate.table[rate.table[, 3] == "dem1", 3]    <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3]    <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }

  if (constrain$sym.hyperstates) {
    rate.table[rate.table[, 3] == "asc2", 3]  <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc1"
  }

  if (oneway) {
    rate.table[rate.table[, 3] == "tran21", 3] <- "0"
  }

  if (!is.null(state.names)) {
    s1 <- state.names[1]; s2 <- state.names[2]
    name_map <- c(
      "asc1"     = paste0("asc.", s1),
      "asc2"     = paste0("asc.", s2),
      "desc1"    = paste0("desc.", s1),
      "desc2"    = paste0("desc.", s2),
      "pol1"     = paste0("pol.", s1),
      "pol2"     = paste0("pol.", s2),
      "dem1"     = paste0("dem.", s1),
      "dem2"     = paste0("dem.", s2),
      ".5*dem1"  = paste0(".5*dem.", s1),
      ".5*dem2"  = paste0(".5*dem.", s2),
      "tran12"   = paste0("tran.", s1, ".to.", s2),
      "tran21"   = paste0("tran.", s2, ".to.", s1),
      "tranSAF"  = "tran.SAF",
      "tranRo"   = "tran.Ro"
    )
    for (old_name in names(name_map)) {
      rate.table[rate.table[, 3] == old_name, 3] <- name_map[old_name]
    }
  }

  parMat[as.matrix(rate.table[, 1:2])] <- rate.table[, 3]
  colnames(parMat) <- rownames(parMat) <- colnames(data)
  list(parMat = parMat, rate.table = rate.table)
}
