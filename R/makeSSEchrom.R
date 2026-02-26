#' Simulate Chromosome Data Under a Binary SSE Model
#'
#' Simulates chromosome number evolution alongside a binary trait under a
#' State-dependent Speciation and Extinction (SSE) model using diversitree's
#' `tree.musse()`. This generates a phylogenetic tree together with tip data
#' for both chromosome number and binary character state.
#'
#' The binary state at the root is assigned randomly. Make sure `h` is set
#' high enough that it does not act as a hard barrier on chromosome number
#' evolution.
#'
#' @param h Integer. The maximum haploid chromosome number for the simulation.
#' @param lambda1 Numeric. Speciation rate in binary state 1.
#' @param lambda2 Numeric. Speciation rate in binary state 2.
#' @param mu1 Numeric. Extinction rate in binary state 1.
#' @param mu2 Numeric. Extinction rate in binary state 2.
#' @param asc1 Numeric. Chromosome gain rate in binary state 1.
#' @param asc2 Numeric. Chromosome gain rate in binary state 2.
#' @param desc1 Numeric. Chromosome loss rate in binary state 1.
#' @param desc2 Numeric. Chromosome loss rate in binary state 2.
#' @param trans1 Numeric. Transition rate from binary state 1 to state 2.
#' @param trans2 Numeric. Transition rate from binary state 2 to state 1.
#' @param max.taxa Integer. The simulation runs until this many extant taxa
#'   are present.
#' @param x0 Integer. Starting haploid chromosome number at the root.
#' @param state.names An optional character vector of length 2 giving
#'   descriptive names for the two binary states, e.g. `c("matched",
#'   "mismatched")`. When provided, the `binary.char` output uses these
#'   labels instead of 0/1. Defaults to `NULL`.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{tree}{A phylogenetic tree of class `"phylo"`.}
#'     \item{binary.char}{A named vector of binary character states for each
#'       tip. If `state.names` is provided, values are the descriptive labels;
#'       otherwise integer 0 or 1.}
#'     \item{haploid.num}{A named integer vector of haploid chromosome numbers
#'       for each tip.}
#'   }
#'
#' @seealso [simChrom()] for simulating chromosome evolution on a fixed tree.
#'
#' @references
#' Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
#' drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
#' 511--523.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' result <- makeSSEchrom(
#'   h = 10, lambda1 = 0.3, lambda2 = 0.3,
#'   mu1 = 0.1, mu2 = 0.1,
#'   asc1 = 0.1, asc2 = 0.1,
#'   desc1 = 0.1, desc2 = 0.1,
#'   trans1 = 0.05, trans2 = 0.05,
#'   max.taxa = 30, x0 = 5,
#'   state.names = c("matched", "mismatched")
#' )
#' plot(result$tree)
#' table(result$binary.char)
#' }
#'
#' @export
makeSSEchrom <- function(h,                 # max haploid number
                         lambda1, lambda2,  # speciation in state 1, 2
                         mu1, mu2,          # extinction in state 1, 2
                         asc1, asc2,        # chromosome gain in state 1, 2
                         desc1, desc2,      # chromosome loss in state 1, 2
                         trans1, trans2,     # transitions in binary char
                         max.taxa,          # max taxa
                         x0,                # starting haploid number randomly
                         state.names = NULL  # descriptive state labels
) {
  # Validate state.names
  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2) {
      stop("state.names must be a character vector of length 2")
    }
  }

  # internal function to make parameter string
  make.pars <- function(k, lambda1, lambda2, mu1, mu2,
                        asc1, asc2, desc1, desc2, trans1, trans2) {
    m <- matrix(rep(0, k^2), k, k)
    for (r in 1:k) {
      for (c in 1:k) {
        if (r <= (k / 2)) {        # upper half
          if (c <= (k / 2)) {      # left half
            if (r == (c - 1)) m[r, c] <- "asc1"
            if (r == (c + 1)) m[r, c] <-  "dsc1"
          }
          if (c > (k / 2)) {         # right half
            if ((r + (k / 2)) == c) m[r, c] <- "trans1"
          }
        }
        if (r > (k / 2)) {         # lower half
          if (c <= (k / 2)) {      # left half
            if (r == (c + (k / 2)))  m[r, c] <- "trans2"
          }
          if (c > (k / 2)) {       # right half
            if (r == (c - 1)) m[r, c] <- "asc2"
            if (r == (c + 1)) m[r, c] <-  "dsc2"
          }
        }
      }
    }
    m[m == "asc1"] <- asc1
    m[m == "asc2"] <- asc2
    m[m == "dsc1"] <- desc1
    m[m == "dsc2"] <- desc2
    m[m == "trans1"] <- trans1
    m[m == "trans2"] <- trans2

    diag(m) <- NA
    z <- as.vector(t(m))
    z <- as.numeric(z[!is.na(z)])
    par.vals <- c(rep(lambda1, times = (k / 2)), rep(lambda2, times = (k / 2)),
                  rep(mu1, times = (k / 2)), rep(mu2, times = (k / 2)), z)
    return(as.numeric(par.vals))
  }

  # internal function to convert from musse to chromosomes
  convert.musse <- function(musse, k, state.names) {
    tree <- musse[c(1, 2, 3, 7)]
    class(tree) <- "phylo"
    z <- musse$tip.state
    bin.char <- chrom.char <- musse$tip.state
    bin.char[musse$tip.state <= (k / 2)] <- 0
    bin.char[musse$tip.state > (k / 2)] <- 1
    for (i in 1:length(musse$tip.state)) {
      if (musse$tip.state[i] > (k / 2)) chrom.char[i] <-
          (chrom.char[i] - (k / 2))
    }

    # Apply state names if provided
    if (!is.null(state.names)) {
      bin.char.labeled <- ifelse(bin.char == 0, state.names[1], state.names[2])
      names(bin.char.labeled) <- names(bin.char)
      results <- list(tree, bin.char.labeled, chrom.char)
    } else {
      results <- list(tree, bin.char, chrom.char)
    }
    names(results) <- c("tree", "binary.char", "haploid.num")
    return(results)
  }

  k <- h * 2                    # setup real dimensions of matrix
  x0 <- x0 + sample(c(0, h), 1)  # randomly assign binary state
  t.pars <- make.pars(k = k, lambda1 = lambda1, lambda2 = lambda2,
                      mu1 = mu1, mu2 = mu2,
                      asc1 = asc1, asc2 = asc2,
                      desc1 = desc1, desc2 = desc2,
                      trans1 = trans1, trans2 = trans2)
  results.musse <- tree.musse(pars = t.pars, max.taxa = max.taxa, x0 = x0)
  results.chrom <- convert.musse(results.musse, k = k, state.names = state.names)
  return(results.chrom)
}
