#' Constrain an mkn Likelihood Function for Chromosome Evolution
#'
#' Constrains a diversitree `mkn` likelihood function to match a biologically
#' realistic model of chromosome number evolution. Supports three major model
#' types: a simple chromEvol model (no binary trait), a model where ploidy is
#' the hidden state, and a model where a binary trait affects chromosome
#' evolution rates. Additionally supports a sex chromosome-autosome fusion
#' (SAF) model.
#'
#' @param data A probability matrix as produced by [datatoMatrix()]. Rows are
#'   species, columns are chromosome states (and optionally hyperstate columns).
#' @param lik A likelihood function created by `diversitree::make.mkn()`.
#' @param hyper Logical. If `TRUE` (default), includes a binary hyperstate
#'   allowing different rates of chromosome evolution in each state.
#' @param polyploidy Logical. If `TRUE`, the hyperstate represents ploidy
#'   level (diploid vs. polyploid), and transitions between states also change
#'   chromosome number. Defaults to `FALSE`.
#' @param verbose Logical. If `TRUE`, returns a list containing the
#'   constrained likelihood function, the parameter identity matrix, and the
#'   rate table. Defaults to `FALSE`.
#' @param oneway Logical. If `TRUE`, the transition rate from state 2 back to
#'   state 1 is set to zero. Defaults to `FALSE`.
#' @param state.names An optional character vector of length 2 giving
#'   descriptive names for the two binary states, e.g. `c("matched",
#'   "mismatched")` or `c("diploid", "polyploid")`. When provided, all
#'   parameter names use these labels (e.g. `asc.matched` instead of `asc1`,
#'   `tran.matched.to.mismatched` instead of `tran12`). Defaults to `NULL`
#'   (numeric suffixes).
#' @param constrain A list of additional model constraints. Can include:
#'   \describe{
#'     \item{drop.poly}{Logical. If `TRUE`, polyploidy rate is set to zero.}
#'     \item{drop.demi}{Logical. If `TRUE`, demiploidy rate is set to zero.}
#'     \item{symmetric}{Logical. If `TRUE`, chromosome change rates are equal
#'       across binary states.}
#'     \item{nometa}{Logical. If `TRUE`, chromosome rates are constrained to
#'       be equal across hyperstates.}
#'     \item{meta}{Character. Either `"ARD"` (all rates different, default) or
#'       `"SYM"` (symmetric transitions between hyperstates).}
#'     \item{saf.model}{Logical. If `TRUE`, transitions between hyperstates
#'       follow the sex chromosome-autosome fusion model.}
#'     \item{sym.hyperstates}{Logical. If `TRUE`, ascending and descending
#'       rates are equal across hyperstates.}
#'   }
#'
#' @return If `verbose = FALSE` (default), returns a constrained likelihood
#'   function compatible with `diversitree::find.mle()` and
#'   `diversitree::mcmc()`. If `verbose = TRUE`, returns a list with elements:
#'   \describe{
#'     \item{`likelihood function`}{The constrained likelihood function.}
#'     \item{`parameter matrix`}{A matrix showing which rate category each
#'       transition belongs to.}
#'     \item{`ratetable`}{A data frame mapping transitions to rate names.}
#'   }
#'
#' @details
#' The rates in the model are (shown with default numeric suffixes; when
#' `state.names` is provided, the suffixes are replaced with descriptive
#' labels):
#' \describe{
#'   \item{asc1/asc2}{Ascending aneuploidy in state 1/2}
#'   \item{desc1/desc2}{Descending aneuploidy in state 1/2}
#'   \item{pol1/pol2}{Polyploidization in state 1/2}
#'   \item{dem1/dem2}{Demipolyploidy in state 1/2}
#'   \item{redip}{Rediploidization}
#'   \item{tran12/tran21}{Transitions between hyperstates}
#'   \item{tranSAF/tranRo}{SAF transitions (sex chromosome-autosome fusions
#'     and Robertsonian transitions)}
#' }
#'
#' @seealso [constrainMuSSE()] for the MuSSE (state-dependent diversification)
#'   version, [datatoMatrix()] for preparing input data.
#'
#' @references
#' Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
#' drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
#' 511--523.
#'
#' @examples
#' \donttest{
#' library(diversitree)
#' # Create example data
#' dat <- data.frame(
#'   species = paste0("sp", 1:5),
#'   chrom = c(5, 6, 7, 8, 10),
#'   prob = c(1, 1, 1, 1, 1)
#' )
#' dat.mat <- datatoMatrix(x = dat, hyper = FALSE)
#' # Create a random tree
#' tree <- ape::rcoal(5, tip.label = dat$species)
#' lik <- make.mkn(tree, states = dat.mat, k = ncol(dat.mat),
#'                 strict = FALSE, control = list(method = "ode"))
#' con.lik <- constrainMkn(data = dat.mat, lik = lik, hyper = FALSE,
#'                          constrain = list(drop.demi = TRUE))
#' argnames(con.lik)
#'
#' # With descriptive state names
#' dat2 <- data.frame(
#'   species = paste0("sp", 1:5),
#'   chrom = c(5, 6, 7, 8, 10),
#'   prob = c(0.8, 0.6, 0.9, 0.5, 1.0)
#' )
#' dat.mat2 <- datatoMatrix(x = dat2, hyper = TRUE)
#' tree2 <- ape::rcoal(5, tip.label = dat2$species)
#' lik2 <- make.mkn(tree2, states = dat.mat2, k = ncol(dat.mat2),
#'                  strict = FALSE, control = list(method = "ode"))
#' con.lik2 <- constrainMkn(data = dat.mat2, lik = lik2, hyper = TRUE,
#'                           state.names = c("matched", "mismatched"),
#'                           constrain = list(drop.demi = TRUE,
#'                                            drop.poly = TRUE))
#' argnames(con.lik2)
#' # Returns: "asc.matched" "desc.matched" "asc.mismatched" ...
#' }
#'
#' @export
constrainMkn <- function(data,
                         lik,
                         hyper = TRUE,
                         polyploidy = FALSE,
                         verbose = FALSE,
                         oneway = FALSE,
                         state.names = NULL,
                         constrain = list(drop.poly = FALSE,
                                        drop.demi = FALSE,
                                        symmetric = FALSE,
                                        nometa = FALSE,
                                        saf.model = FALSE,
                                        sym.hyperstates = FALSE,
                                        meta = "ARD")){

  # Validate state.names
  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2) {
      stop("state.names must be a character vector of length 2")
    }
  }

  # This fills out the list of constraints the default are no constraints
  if (is.null(constrain$drop.poly)) constrain$drop.poly <- FALSE
  if (is.null(constrain$drop.demi)) constrain$drop.demi <- FALSE
  if (is.null(constrain$symmetric)) constrain$symmetric <- FALSE
  if (is.null(constrain$nometa)) constrain$nometa <- FALSE
  if (is.null(constrain$saf.model)) constrain$saf.model <- FALSE
  if (is.null(constrain$sym.hyperstates)) constrain$sym.hyperstates <- FALSE
  if (is.null(constrain$meta)) constrain$meta <- "ARD"

  ## BUILD AN EMPTY MATRIX MATCHING OUR MODEL
  # create and store variable for padding rate names
  if (ncol(data) < 100) pad <- 2
  if (ncol(data) >= 100) pad <- 3
  if (ncol(data) < 10) pad <- 1

  # make the matrix of rates
  parMat <- matrix(0, ncol(data), ncol(data))
  # make the components of the rate names the column and row
  # names this will allow for easy creation of constraints later
  colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  # now we have a matrix with all zeros but the right state names
  # in the column and row names

  # we need to know where our duplication
  # of the matrix begins so here is that
  split <- ncol(parMat) / 2

  # we also need the actual chromosome numbers
  if (hyper == TRUE) chroms <- as.numeric(colnames(data)[1:split])
  if (hyper == FALSE) chroms <- as.numeric(colnames(data))

  ## NOW WE HAVE A SERIES OF LOOPS THAT FILL IN OUR parMAT
  ## MATRIX WITH NUMBERS 1:9 INDICATIVE OF THE DIFFERENT POSSIBLE
  ## RATES WE WISH TO INCLUDE IN OUR MODEL.  EACH OF THESE LOOPS
  ## REPRESENT A DIFFERENT MODEL OF CHROMOSOME EVOLUTION

  ## OLD CHROMEVOL MODEL
  if (hyper == FALSE) {
    message("Constraining model to simple chromevol version")
    for (i in 1:(nrow(parMat) - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- 5 #polyploidy
      if (constrain$drop.demi == FALSE) {
        if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
          x <- chroms[i] * 1.5
          if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10 #demiploidy state1 even
          if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
        }
      }
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
    }
  }

  # MODEL 1 PLOIDY IS HIDDEN STATE
  if (hyper == TRUE & polyploidy == TRUE) {
    message("Creating rate matrix for chosen chromosome model")
    # diploid rates
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, (which(chroms[i] * 2 == chroms) + split)] <- 5 #polyploidy-1
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 10 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 11 #demiploidy state 1 odd
      }
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - diploids
      parMat[(i + 1), i] <- 2 #descending aneuploidy - diploids
    }
    # polyploid rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms)) parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      # demiploidy
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13 #demiploidy state 1 odd
      }
      parMat[i, (i - split)] <- 7 #rediploidization
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 7 #rediploidization
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - polyploids
      parMat[(i + 1), i] <- 4 #descending aneuploidy - polyploids
    }
  }

  # MODEL 2 PLOIDY IS NOT THE HYPER STATE
  if (hyper == TRUE & polyploidy == FALSE & constrain$saf.model == FALSE) {
    message("Creating rate matrix for chosen chromosome model")
    # state 1 rates
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- 5 #polyploidy - 1
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
      }

      parMat[i, (i + split)] <- 8 # transitions state 1->2
      # special case for last row
      if (i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- 8 # transitions state 1->2
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - 1
      parMat[(i + 1), i] <- 2 #descending aneuploidy - 1
    }
    # state 2 rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms)) parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      # demiploidy
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13 #demiploidy state 2 odd
      }
      parMat[i, (i - split)] <- 9 #transition state 2->1
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 9 # transitions state 2->1
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - 2
      parMat[(i + 1), i] <- 4 #descending aneuploidy - 2
    }
  }

  # MODEL 3 SAF IS HYPER STATE
  if (hyper == TRUE & polyploidy == FALSE & constrain$saf.model == TRUE) {
    message("Creating rate matrix for chosen chromosome model")
    # state 1 rates
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- 5 #polyploidy - 1
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
      }

      if (i > 1) {
        parMat[i, (i + split - 1)] <- 14 # transitions state 1->2 (SAF)
      }
      # special case for last row
      if (i == (split - 1)) parMat[(i + 1), (i + split)] <- 14 # transitions state 1->2 (SAF)
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - 1
      parMat[(i + 1), i] <- 2 #descending aneuploidy - 1
    }
    # state 2 rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms)) parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      # demiploidy
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- 12 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13 #demiploidy state 2 odd
      }
      parMat[i, (i - split)] <- 15 #transition state 2->1
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 15 # transitions state 2->1
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - 2
      parMat[(i + 1), i] <- 4 #descending aneuploidy - 2
    }
  }

  rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat), 3))
  rate.table[, 1] <- rep(as.character(row.names(parMat)), each = ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[, 2], ]

  rate.table[rate.table[, 3] == 1, 3] <- "asc1"
  rate.table[rate.table[, 3] == 2, 3] <- "desc1"
  rate.table[rate.table[, 3] == 3, 3] <- "asc2"
  rate.table[rate.table[, 3] == 4, 3] <- "desc2"
  rate.table[rate.table[, 3] == 5, 3] <- "pol1"
  rate.table[rate.table[, 3] == 6, 3] <- "pol2"
  rate.table[rate.table[, 3] == 7, 3] <- "redip"
  rate.table[rate.table[, 3] == 8, 3] <- "tran12"
  rate.table[rate.table[, 3] == 9, 3] <- "tran21"
  rate.table[rate.table[, 3] == 10, 3] <- "dem1"
  rate.table[rate.table[, 3] == 11, 3] <- ".5*dem1"
  rate.table[rate.table[, 3] == 12, 3] <- "dem2"
  rate.table[rate.table[, 3] == 13, 3] <- ".5*dem2"
  rate.table[rate.table[, 3] == 14, 3] <- "tranSAF"
  rate.table[rate.table[, 3] == 15, 3] <- "tranRo"


  if (constrain$nometa == TRUE) {
    rate.table[rate.table[, 3] == "asc2", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc1"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }

  if (constrain$drop.poly == TRUE) {
    rate.table[rate.table[, 3] == "pol1", 3] <- "0"
    rate.table[rate.table[, 3] == "pol2", 3] <- "0"
  }

  if (constrain$drop.demi == TRUE) {
    rate.table[rate.table[, 3] == "dem1", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }

  if (constrain$symmetric == TRUE) {
    rate.table[rate.table[, 3] == "desc1", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "asc2"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == 7, 3] <- "redip"
    rate.table[rate.table[, 3] == 8, 3] <- "tran12"
    rate.table[rate.table[, 3] == 9, 3] <- "tran21"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }

  if (constrain$meta == "SYM") {
    rate.table[rate.table[, 3] == "tran21", 3] <- "tran12"
  }

  if (constrain$saf.model == TRUE) {
    rate.table[rate.table[, 3] == "asc1", 3] <- "asc1"
    rate.table[rate.table[, 3] == "asc2", 3] <- "asc2"
    rate.table[rate.table[, 3] == "desc1", 3] <- "desc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc2"
    rate.table[rate.table[, 3] == "pol1", 3] <- "0"
    rate.table[rate.table[, 3] == "pol2", 3] <- "0"
    rate.table[rate.table[, 3] == "dem1", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }

  if (constrain$sym.hyperstates == TRUE) {
    rate.table[rate.table[, 3] == "asc2", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc1"
  }

  if (oneway == TRUE) {
    rate.table[rate.table[, 3] == "tran21", 3] <- "0"
  }

  # Apply descriptive state names to rate parameters
  if (!is.null(state.names)) {
    s1 <- state.names[1]
    s2 <- state.names[2]
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

  formulae <- vector(mode = "character", length = nrow(rate.table))
  for (i in 1:nrow(rate.table)) {
    formulae[i] <- paste("q",
                         rate.table[i, 1],
                         rate.table[i, 2],
                         " ~ ",
                         rate.table[i, 3],
                         collapse = "", sep = "")
  }

  # Build extras vector with appropriate names
  if (!is.null(state.names)) {
    s1 <- state.names[1]
    s2 <- state.names[2]
    extras <- c(paste0("asc.", s1), paste0("desc.", s1),
                paste0("asc.", s2), paste0("desc.", s2),
                paste0("pol.", s1), paste0("pol.", s2),
                "redip",
                paste0("tran.", s1, ".to.", s2),
                paste0("tran.", s2, ".to.", s1),
                paste0("dem.", s1), paste0("dem.", s2),
                "tran.SAF", "tran.Ro")
  } else {
    extras <- c("asc1", "desc1",
                "asc2", "desc2",
                "pol1", "pol2",
                "redip", "tran12", "tran21",
                "dem1", "dem2",
                "tranSAF", "tranRo")
  }

  # Update parameter matrix to match rate table
  for (i in 1:nrow(rate.table)) {
    parMat[paste0("", rate.table[i, 1], ""), paste0("", rate.table[i, 2], "")] <-
      rate.table[i, 3]
  }

  lik.con <- constrain(lik, formulae = formulae, extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(data)

  if (verbose == TRUE) {
    result.list <- list(lik.con, parMat, rate.table)
    names(result.list) <- c("likelihood function", "parameter matrix", "ratetable")
    return(result.list)
  }
  if (verbose == FALSE) return(lik.con)
}
