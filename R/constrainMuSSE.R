#' Constrain a MuSSE Likelihood Function for Chromosome Evolution
#'
#' Constrains a diversitree `MuSSE` (Multi-State Speciation and Extinction)
#' likelihood function to model chromosome number evolution with
#' state-dependent diversification. This extends [constrainMkn()] by allowing
#' speciation and extinction rates to differ between binary states.
#'
#' @param data A probability matrix as produced by [datatoMatrix()]. Rows are
#'   species, columns are chromosome states (and optionally hyperstate columns).
#' @param lik A likelihood function created by `diversitree::make.musse()`.
#' @param hyper Logical. If `TRUE` (default), includes a binary hyperstate
#'   allowing different rates of chromosome evolution in each state.
#' @param polyploidy Logical. If `TRUE`, the hyperstate represents ploidy
#'   level (diploid vs. polyploid). Defaults to `FALSE`.
#' @param s.lambda Logical. If `TRUE` (default), a single speciation rate is
#'   estimated across all states. If `FALSE`, separate speciation rates are
#'   estimated for each hyperstate.
#' @param s.mu Logical. If `TRUE` (default), a single extinction rate is
#'   estimated across all states. If `FALSE`, separate extinction rates are
#'   estimated for each hyperstate.
#' @param verbose Logical. If `TRUE`, returns a list containing the
#'   constrained likelihood function and the parameter identity matrix.
#'   Defaults to `FALSE`.
#' @param state.names An optional character vector of length 2 giving
#'   descriptive names for the two binary states, e.g. `c("matched",
#'   "mismatched")`. When provided, all parameter names use these labels
#'   (e.g. `asc.matched` instead of `asc1`, `lambda.matched` instead of
#'   `lambda1`). Defaults to `NULL` (numeric suffixes).
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
#'   }
#'
#' @return If `verbose = FALSE` (default), returns a constrained likelihood
#'   function compatible with `diversitree::find.mle()` and
#'   `diversitree::mcmc()`. If `verbose = TRUE`, returns a list with the
#'   constrained likelihood function and parameter identity matrix.
#'
#' @details
#' The chromosome evolution rates are the same as in [constrainMkn()]:
#' ascending/descending aneuploidy, polyploidy, demiploidy, and transitions
#' between hyperstates. In addition, this function constrains speciation
#' (`lambda`) and extinction (`mu`) rates, which can either be shared across
#' states or estimated separately depending on `s.lambda` and `s.mu`.
#'
#' When `state.names` is provided, speciation and extinction parameters
#' are also labeled descriptively (e.g. `lambda.matched`, `mu.mismatched`).
#'
#' @seealso [constrainMkn()] for the mkn (no diversification) version,
#'   [datatoMatrix()] for preparing input data.
#'
#' @references
#' Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
#' drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
#' 511--523.
#'
#' @examples
#' \donttest{
#' library(diversitree)
#' # Create example data with uncertainty in binary trait
#' dat <- data.frame(
#'   species = paste0("sp", 1:5),
#'   chrom = c(5, 6, 7, 8, 10),
#'   prob = c(0.8, 0.6, 0.9, 0.5, 1.0)
#' )
#' dat.mat <- datatoMatrix(x = dat, hyper = TRUE)
#' tree <- ape::rcoal(5, tip.label = dat$species)
#' lik <- make.musse(tree, states = dat.mat, k = ncol(dat.mat),
#'                   strict = FALSE, control = list(method = "ode"))
#' con.lik <- constrainMuSSE(data = dat.mat, lik = lik,
#'                            s.lambda = FALSE, s.mu = FALSE,
#'                            state.names = c("matched", "mismatched"),
#'                            constrain = list(drop.demi = TRUE))
#' argnames(con.lik)
#' }
#'
#' @export
constrainMuSSE <- function(data,
                           lik,
                           hyper = TRUE,
                           polyploidy = FALSE,
                           s.lambda = TRUE,
                           s.mu = TRUE,
                           verbose = FALSE,
                           state.names = NULL,
                           constrain = list(drop.poly = FALSE,
                                          drop.demi = FALSE,
                                          symmetric = FALSE,
                                          nometa = FALSE,
                                          meta = "ARD")){

  # Validate state.names
  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2) {
      stop("state.names must be a character vector of length 2")
    }
  } else if (isTRUE(hyper)) {
    message("constrainMuSSE: state.names not provided. Output rates will use ",
            "'1'/'2' suffixes, where suffix '1' refers to the FIRST half of ",
            "the columns of `data` (the half without the 'h' suffix from ",
            "datatoMatrix) and '2' refers to the second half. Pass ",
            "state.names = c(\"name1\", \"name2\") to get unambiguous rate ",
            "names like asc.name1 / asc.name2.")
  }

  # This fills out the list of constraints the default are no constraints
  if (is.null(constrain$drop.poly)) constrain$drop.poly <- FALSE
  if (is.null(constrain$drop.demi)) constrain$drop.demi <- FALSE
  if (is.null(constrain$symmetric)) constrain$symmetric <- FALSE
  if (is.null(constrain$nometa)) constrain$nometa <- FALSE
  if (is.null(constrain$meta)) constrain$meta <- "ARD"

  ## BUILD AN EMPTY MATRIX MATCHING OUR MODEL
  # padding rate names so column labels sort as strings
  pad <- nchar(as.character(ncol(data)))
  # make the matrix of rates
  parMat <- matrix(0, ncol(data), ncol(data))
  # make the components of the rate names the column and row
  # names this will allow for easy creation of constraints later
  colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)

  # we need to know where our duplication of the matrix begins so here is that
  split <- ncol(parMat) / 2

  # we also need the actual chromosome numbers
  if (hyper == TRUE) chroms <- as.numeric(colnames(data)[1:split])
  if (hyper == FALSE) chroms <- as.numeric(colnames(data))

  ## NOW WE HAVE A SERIES OF LOOPS THAT FILL IN OUR parMat
  ## MATRIX WITH NUMBERS 1:13 INDICATIVE OF THE DIFFERENT POSSIBLE
  ## RATES WE WISH TO INCLUDE IN OUR MODEL.  EACH OF THESE LOOPS
  ## REPRESENT A DIFFERENT MODEL OF CHROMOSOME EVOLUTION

  ## OLD CHROMEVOL MODEL
  if (hyper == FALSE) {
    message("Constraining model to simple chromevol version")
    for (i in 1:(nrow(parMat) - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- 5 #polyploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10 #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
      }
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
    }
  }

  # BiSCE MODEL 1 PLOIDY IS hyper STATE
  if (hyper == TRUE & polyploidy == TRUE) {
    message("Constraining model where ploidy is a meta state and different rates of chromosome evolution are possible based on being polyploid or diploid")
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

  # BiSPCE 2 PLOIDY IS NOT THE HYPER STATE
  if (hyper == TRUE & polyploidy == FALSE) {
    message("Constraining model with a hyper state that may have different rates of chromosome number evolution")
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
  # we now have a matrix with a number 1-13 that matches the rates present
  # under one of our models we will use this to build our
  # arguments for the standard diversitree constrain function
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
      "tran21"   = paste0("tran.", s2, ".to.", s1)
    )
    for (old_name in names(name_map)) {
      rate.table[rate.table[, 3] == old_name, 3] <- name_map[old_name]
    }
  }

  formulae <- paste0("q", rate.table[, 1], rate.table[, 2], " ~ ",
                     rate.table[, 3])

  ## Lambda and Mu
  lambda <- mu <- vector()
  for (i in 1:nrow(parMat)) {
    # Determine lambda/mu names based on state.names
    lam1_name <- if (!is.null(state.names)) paste0("lambda.", state.names[1]) else "lambda1"
    lam2_name <- if (!is.null(state.names)) paste0("lambda.", state.names[2]) else "lambda2"
    mu1_name <- if (!is.null(state.names)) paste0("mu.", state.names[1]) else "mu1"
    mu2_name <- if (!is.null(state.names)) paste0("mu.", state.names[2]) else "mu2"

    # no hyper state
    if (hyper == FALSE) {
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ ", lam1_name, sep = ""))
      mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ ", mu1_name, sep = ""))
    }
    # hyper model with single lambda
    if (hyper == TRUE & s.lambda == TRUE) {
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ ", lam1_name, sep = ""))
    }
    # hyper model with single mu
    if (hyper == TRUE & s.mu == TRUE) {
      mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ ", mu1_name, sep = ""))
    }
    # hyper model with two lambdas
    if (hyper == TRUE & s.lambda == FALSE) {
      if (i <= split) {
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ ", lam1_name, sep = ""))
      }
      if (i > split) {
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ ", lam2_name, sep = ""))
      }
    }
    # hyper model with two mus
    if (hyper == TRUE & s.mu == FALSE) {
      if (i <= split) {
        mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ ", mu1_name, sep = ""))
      }
      if (i > split) {
        mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ ", mu2_name, sep = ""))
      }
    }
  }

  # Build extras vector with appropriate names
  if (!is.null(state.names)) {
    s1 <- state.names[1]
    s2 <- state.names[2]
    extras <- c(paste0("asc.", s1), paste0("desc.", s1),
                paste0("asc.", s2), paste0("desc.", s2),
                paste0("pol.", s1), paste0("pol.", s2),
                paste0("dem.", s1), paste0("dem.", s2),
                "redip",
                paste0("tran.", s1, ".to.", s2),
                paste0("tran.", s2, ".to.", s1),
                paste0("lambda.", s1), paste0("mu.", s1),
                paste0("lambda.", s2), paste0("mu.", s2))
  } else {
    extras <- c("asc1", "desc1", "asc2", "desc2",
                "pol1", "pol2", "dem1", "dem2", "redip", "tran12", "tran21",
                "lambda1", "mu1", "lambda2", "mu2")
  }

  lik.con <- constrain(lik, formulae = c(formulae, lambda, mu), extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(data)
  if (verbose == TRUE) return(list(lik.con, parMat))
  if (verbose == FALSE) return(lik.con)
}
