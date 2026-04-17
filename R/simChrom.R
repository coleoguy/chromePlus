#' Simulate Chromosome Number Evolution on a Phylogeny
#'
#' Simulates chromosome number evolution along a fixed phylogenetic tree using
#' one of four built-in models or a user-provided Q-matrix. Internally uses
#' `diversitree::sim.character()` with an `mkn` model.
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param pars Parameter values for the model. Can be either:
#'   \describe{
#'     \item{Unnamed numeric vector}{Positional parameters (original behavior).
#'       Length depends on the model:
#'       \describe{
#'         \item{`"2010"` (length 5)}{gain, loss, demiploidy, polyploidy, root
#'           chromosome number}
#'         \item{`"ChromPlus"` or `"SAF"` (length 12)}{gain1, gain2, loss1,
#'           loss2, demiploidy1, demiploidy2, polyploidy1, polyploidy2,
#'           transition 1->2, transition 2->1, root chromosome number, root
#'           binary state (0 or 1)}
#'         \item{`"PloidEvol"` (length 11)}{gain-diploid, gain-polyploid,
#'           loss-diploid, loss-polyploid, demiploidy-diploid,
#'           demiploidy-polyploid, polyploidy-diploid, polyploidy-polyploid,
#'           rediploidization, root chromosome number, root ploidy state
#'           (0=diploid, 1=polyploid)}
#'         \item{`NULL` with Qmat (length 1)}{root chromosome number}
#'       }
#'     }
#'     \item{Named numeric vector}{Parameters identified by name rather than
#'       position. Names can use numeric suffixes (e.g. `asc1`, `desc2`) or
#'       descriptive labels matching `state.names` (e.g. `asc.matched`,
#'       `desc.mismatched`). For `"2010"` model: `asc`, `desc`, `dem`, `pol`,
#'       `root.chrom`. For `"ChromPlus"`/`"SAF"`: `asc1`, `asc2`, `desc1`,
#'       `desc2`, `dem1`, `dem2`, `pol1`, `pol2`, `tran12`, `tran21`,
#'       `root.chrom`, `root.state`. For `"PloidEvol"`: `asc1`, `asc2`,
#'       `desc1`, `desc2`, `dem1`, `dem2`, `pol1`, `pol2`, `redip`,
#'       `root.chrom`, `root.state`.}
#'   }
#' @param limits Numeric vector of length 2 with the minimum and maximum
#'   chromosome numbers, e.g. `c(3, 20)`. Required for built-in models.
#'   Ignored when using a custom Q-matrix.
#' @param model Character string specifying the model: `"2010"`, `"ChromPlus"`,
#'   `"PloidEvol"`, or `"SAF"`. Defaults to `NULL` (user-provided Q-matrix).
#' @param Qmat A user-provided square Q-matrix describing chromosome evolution
#'   rates. Row and column names should be chromosome state labels. Diagonals
#'   are automatically set. Only used when `model = NULL`.
#' @param verbose Logical. If `TRUE`, the parameter matrix used to build the
#'   model is also returned. Defaults to `FALSE`.
#' @param state.names An optional character vector of length 2 giving
#'   descriptive names for the two binary states, e.g. `c("matched",
#'   "mismatched")`. When provided, output uses these labels instead of
#'   0/1 for binary states. Also allows named `pars` to use these labels.
#'   Defaults to `NULL`.
#'
#' @return Depends on the model:
#'   \describe{
#'     \item{`"2010"` or `NULL`}{A named numeric vector of chromosome numbers
#'       at tree tips. If `verbose = TRUE`, a list with `chrom.num` and
#'       `parameter.matrix`.}
#'     \item{`"ChromPlus"`}{A list with `binary.state` (labeled with
#'       `state.names` if provided, otherwise 0/1) and `chrom.num`. If
#'       `verbose = TRUE`, also includes `parameter.matrix`.}
#'     \item{`"PloidEvol"`}{A list with `ploidy.state` and `chrom.num`.}
#'     \item{`"SAF"`}{A list with `fusion.state` and `chrom.num`.}
#'   }
#'
#' @seealso [makeSSEchrom()] for simulating both tree and chromosome data
#'   under an SSE model, [constrainMkn()] and [constrainMuSSE()] for fitting
#'   models.
#'
#' @references
#' Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
#' drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
#' 511--523.
#'
#' @examples
#' \donttest{
#' tree <- ape::rcoal(20)
#' # Simple 2010 model
#' result <- simChrom(tree = tree,
#'                    pars = c(0.1, 0.1, 0.0, 0.01, 5),
#'                    limits = c(3, 12), model = "2010")
#' head(result)
#'
#' # ChromPlus model with named parameters and state labels
#' result2 <- simChrom(tree = tree,
#'   pars = c(asc.matched = 0.1, asc.mismatched = 0.1,
#'            desc.matched = 0.1, desc.mismatched = 0.1,
#'            dem.matched = 0, dem.mismatched = 0,
#'            pol.matched = 0.01, pol.mismatched = 0.01,
#'            tran.matched.to.mismatched = 0.05,
#'            tran.mismatched.to.matched = 0.05,
#'            root.chrom = 5, root.state = 0),
#'   limits = c(3, 12), model = "ChromPlus",
#'   state.names = c("matched", "mismatched"))
#' table(result2$binary.state)
#' }
#'
#' @export
simChrom <- function(tree, pars, limits = NULL, model = NULL, Qmat = NULL,
                     verbose = FALSE, state.names = NULL) {

  # Validate state.names
  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2) {
      stop("state.names must be a character vector of length 2")
    }
  } else if (!is.null(model) && model %in% c("ChromPlus", "SAF", "PloidEvol")) {
    message("simChrom: state.names not provided. In the positional pars ",
            "vector, the '1' suffix refers to the FIRST half of the ",
            "Q-matrix (state 0 in the returned binary.state / ploidy.state / ",
            "fusion.state vector) and '2' to the second half. Pass ",
            "state.names = c(\"name1\", \"name2\") for unambiguous output ",
            "labels.")
  }

  # Resolve named pars to positional vector
  if (!is.null(names(pars))) {
    pars <- .resolve_named_pars(pars, model, state.names)
  }

  if (is.null(model) == FALSE && model == "2010") {
    message("building q-matrix")
    root <- pars[5] - limits[1] + 1
    if (length(pars) != 5) stop("pars should have length of 5")
    # set up an empty matrix
    q <- matrix(0, length(limits[1]:limits[2]), length(limits[1]:limits[2]))
    rownames(q) <- colnames(q) <- chroms <- limits[1]:limits[2]
    # fill in the matrix
    for (i in 1:(nrow(q) - 1)) {
      q[i, (i + 1)] <- pars[1] # gain
      q[(i + 1), i] <- pars[2] # loss
      if ((chroms[i] * 2) <= max(chroms)) q[i, which(chroms == (chroms[i] * 2))] <- pars[4] #poly
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  q[i, which(chroms == x)] <- q[i, which(chroms == x)] + pars[3] #demi even
        if (x %% 1 != 0)  q[i, which(chroms %in% c(floor(x), ceiling(x)))] <- q[i, which(chroms %in% c(floor(x), ceiling(x)))] + pars[3] / 2 #demi odd
      }
      # special fix for chromosome num = 1
      diag(q) <- 0
    }
  }

  ## ChromPlus MODEL Q MATRIX
  if (is.null(model) == FALSE && model == "ChromPlus") {
    message("building q-matrix")
    if (length(pars) != 12) stop("pars should have length of 12")
    # set up an empty matrix
    parMat <- matrix(0, 2 * length(limits[1]:limits[2]), 2 * length(limits[1]:limits[2]))
    pad <- nchar(as.character(ncol(parMat)))
    colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep = ""), 1:ncol(parMat))
    rownames(parMat) <- colnames(parMat)
    split <- ncol(parMat) / 2
    chroms <- limits[1]:limits[2]
    # state 1 rates
    for (i in 1:(split - 1)) {
      parMat[i, (i + 1)] <- pars[1] #ascending aneuploidy - 1
      parMat[(i + 1), i] <- pars[3] #descending aneuploidy - 1
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- parMat[i, which(chroms == (chroms[i] * 2))] + pars[7] #polyploidy - 1
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- parMat[i, which(chroms == x)] + pars[5] #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] + (pars[5] / 2) #demiploidy state 1 odd
      }
      parMat[i, (i + split)] <- pars[9] # transitions state 1->2
      # special case for last row
      if (i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- pars[9] # transitions state 1->2

    }
    # state 2 rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      parMat[i, (i + 1)] <- pars[2] #ascending aneuploidy - 2
      parMat[(i + 1), i] <- pars[4] #descending aneuploidy - 2
      if ((chroms[i - split] * 2) <= max(chroms)) parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] + pars[8] #polyploidy-2
      # demiploidy
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- parMat[i, (which(chroms == x) + split)] + pars[6] #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] + (pars[6] / 2) #demiploidy state 2 odd
      }
      parMat[i, (i - split)] <- pars[10] #transition state 2->1
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- pars[10] # transitions state 2->1
    }
    q <- parMat
    if (pars[12] == 0) {
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[11])])
    } else if (pars[12] == 1) {
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral binary state not recognized")
    }
  }

  ## PloidEvol MODEL Q MATRIX
  if (is.null(model) == FALSE && model == "PloidEvol") {
    message("building q-matrix")
    if (length(pars) != 11) stop("pars should have length of 11")
    # set up an empty matrix
    parMat <- matrix(0, 2 * length(limits[1]:limits[2]), 2 * length(limits[1]:limits[2]))
    pad <- nchar(as.character(ncol(parMat)))
    colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep = ""), 1:ncol(parMat))
    rownames(parMat) <- colnames(parMat)
    split <- ncol(parMat) / 2
    chroms <- limits[1]:limits[2]
    # diploidy rates
    for (i in 1:(split - 1)) {
      parMat[i, (i + 1)] <- pars[1] #ascending aneuploidy - 1
      parMat[(i + 1), i] <- pars[3] #descending aneuploidy - 1
      #polyploidy - 1
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, (split + which(chroms == (chroms[i] * 2)))] <-
        pars[7] + parMat[i, which(chroms == (chroms[i] * 2))]
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        #demiploidy state1 even
        if (x %% 1 == 0)  parMat[i, (split + which(chroms == x))] <-
            pars[5] + parMat[i, (split + which(chroms == x))]
        #demiploidy state 1 odd
        if (x %% 1 != 0)  parMat[i, (split + which(chroms %in% c(floor(x), ceiling(x))))] <-
            (pars[5] / 2) + parMat[i, (split + which(chroms %in% c(floor(x), ceiling(x))))]
      }
    }
    # polyploidy rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      parMat[i, (i + 1)] <- pars[2] + parMat[i, (i + 1)] #ascending aneuploidy - 2
      parMat[(i + 1), i] <- pars[4] + parMat[(i + 1), i] #descending aneuploidy - 2
      #polyploidy-2
      if ((chroms[i - split] * 2) <= max(chroms)) {
        parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <-
          pars[8] + parMat[i, (which(chroms[i - split] * 2 == chroms) + split)]
      }
      # demiploidy-2
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        #demiploidy-2 even
        if (x %% 1 == 0) {
          parMat[i, (which(chroms == x) + split)] <-
            pars[6] + parMat[i, (which(chroms == x) + split)]
        }
        #demiploidy-2 odd
        if (x %% 1 != 0) {
          parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <-
            (pars[6] / 2) + parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)]
        }
      }
      # rediploidization
      parMat[i, (i - split)] <- pars[9]
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- pars[9]
    }
    q <- parMat
    if (pars[11] == 0) {
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[10])])
    } else if (pars[11] == 1) {
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[10])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral ploidy state not recognized")
    }
  }

  ## SAF MODEL Q MATRIX
  if (is.null(model) == FALSE && model == "SAF") {
    message("building q-matrix")
    if (length(pars) != 12) stop("pars should have length of 12")
    # set up an empty matrix
    parMat <- matrix(0, 2 * length(limits[1]:limits[2]), 2 * length(limits[1]:limits[2]))
    pad <- nchar(as.character(ncol(parMat)))
    colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep = ""), 1:ncol(parMat))
    rownames(parMat) <- colnames(parMat)
    split <- ncol(parMat) / 2
    chroms <- limits[1]:limits[2]
    # state 1 rates
    for (i in 1:(split - 1)) {
      parMat[i, (i + 1)] <- pars[1] #ascending aneuploidy - 1
      parMat[(i + 1), i] <- pars[3] #descending aneuploidy - 1
      if ((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- parMat[i, which(chroms == (chroms[i] * 2))] + pars[7] #polyploidy - 1
      # demiploidy
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x %% 1 == 0)  parMat[i, which(chroms == x)] <- parMat[i, which(chroms == x)] + pars[5] #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] + (pars[5] / 2) #demiploidy state 1 odd
      }
      # SAF transitions reduce chromosome count by 1 on fusion, so the
      # off-diagonal target is (i + split - 1). Column (i + split - 1) only
      # makes sense for i > 1 (otherwise it points back into the state-1
      # half). This matches constrainMkn's SAF model.
      if (i > 1) parMat[i, (i + split - 1)] <- pars[9] # transitions state 1->2
      # special case for last row
      if (i == (split - 1)) parMat[(i + 1), (i + split)] <- pars[9] # transitions state 1->2

    }
    # state 2 rates
    for (i in (split + 1):(nrow(parMat) - 1)) {
      parMat[i, (i + 1)] <- pars[2] #ascending aneuploidy - 2
      parMat[(i + 1), i] <- pars[4] #descending aneuploidy - 2
      if ((chroms[i - split] * 2) <= max(chroms)) parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] <- parMat[i, (which(chroms[i - split] * 2 == chroms) + split)] + pars[8] #polyploidy-2
      # demiploidy
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x %% 1 == 0)  parMat[i, (which(chroms == x) + split)] <- parMat[i, (which(chroms == x) + split)] + pars[6] #demiploidy state1 even
        if (x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] + (pars[6] / 2) #demiploidy state 2 odd
      }
      parMat[i, (i - split)] <- pars[10] #transition state 2->1
      # special case for last row
      if (i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- pars[10] # transitions state 2->1
    }
    q <- parMat
    if (pars[12] == 0) {
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[11])])
    } else if (pars[12] == 1) {
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral fusion state not recognized")
    }
  }

  ## USER PROVIDED Q MATRIX
  if (is.null(model) == TRUE) {
    if (is.null(Qmat) == TRUE) stop("no model or q-matrix provided")

    if (is.null(limits) == FALSE) {
      message("ignoring provided limits")
    }

    if (is.matrix(Qmat) == FALSE) stop("q-matrix provided is not matrix")

    if (nrow(Qmat) != ncol(Qmat)) stop("q-matrix provided has unequal number of rows and columns")

    if (setequal(rownames(Qmat), colnames(Qmat)) == FALSE) stop("q-matrix provided has nonmatching column and row names")

    if (length(pars) != 1) stop("pars should have length of 1")

    #Assign user supplied Q-matrix to q
    message("using user provided q-matrix")
    q <- Qmat

    if (pars[1] %in% colnames(q) == FALSE) stop("root provided does not fall within range of q-matrix")

    #Get root state
    root <- pars[1] - as.numeric(colnames(q)[1]) + 1

  }

  diag(q) <- 0
  diag(q) <- -rowSums(q)

  # simulate the chromosome numbers
  message("performing simulation")
  dsims <- sim.character(tree, pars = q, x0 = root, model = "mkn")
  attr(dsims, "node.state") <- NULL
  # save the names for various uses below
  tips <- names(dsims)
  # in the case of the 2010 model we have column names = to the
  # chromosome numbers so we can just use them
  if (model == "2010" || is.null(model) == TRUE) {
    dsims[] <- as.numeric(colnames(q)[dsims])
    #Check if parameter matrix is returned
    if (verbose == TRUE) {
      result.list <- list(dsims, parMat)
      names(result.list) <- c("chrom.num", "parameter.matrix")
      return(result.list)
    } else {
      return(dsims)
    }
  }

  # under the ChromPlus model things are bit more complex and have
  # to be converted back to chromosome number and binary state
  if (model == "ChromPlus") {
    # for chromRate we need to return two vectors
    # 1) binary state
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1

    # Apply state names if provided
    if (!is.null(state.names)) {
      b.state.labeled <- ifelse(b.state == 0, state.names[1], state.names[2])
      names(b.state.labeled) <- tips
    }

    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]

    #Check if parameter matrix is returned
    if (verbose == TRUE) {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims, parMat)
      names(result) <- c("binary.state", "chrom.num", "parameter.matrix")
      return(result)
    } else {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims)
      names(result) <- c("binary.state", "chrom.num")
      return(result)
    }
  }

  # under the PloidEvol model things are bit more complex and have
  # to be converted back to chromosome number and ploidy state
  if (model == "PloidEvol") {
    # for chromRate we need to return two vectors
    # 1) binary state
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1

    # Apply state names if provided
    if (!is.null(state.names)) {
      b.state.labeled <- ifelse(b.state == 0, state.names[1], state.names[2])
      names(b.state.labeled) <- tips
    }

    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]

    #Check if parameter matrix is returned
    if (verbose == TRUE) {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims, parMat)
      names(result) <- c("ploidy.state", "chrom.num", "parameter.matrix")
      return(result)
    } else {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims)
      names(result) <- c("ploidy.state", "chrom.num")
      return(result)
    }
  }

  # under the SAF model things are bit more complex and have
  # to be converted back to chromosome number and fusion state
  if (model == "SAF") {
    # for SAF we need to return two vectors
    # 1) binary state
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1

    # Apply state names if provided
    if (!is.null(state.names)) {
      b.state.labeled <- ifelse(b.state == 0, state.names[1], state.names[2])
      names(b.state.labeled) <- tips
    }

    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]

    if (verbose == TRUE) {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims, parMat)
      names(result) <- c("fusion.state", "chrom.num", "parameter.matrix")
      return(result)
    } else {
      result <- list(if (!is.null(state.names)) b.state.labeled else b.state,
                     dsims)
      names(result) <- c("fusion.state", "chrom.num")
      return(result)
    }
  }
}


#' Resolve Named Parameters to Positional Vector
#'
#' Internal helper function that maps named parameter vectors to the
#' positional order expected by each model in [simChrom()].
#'
#' @param pars Named numeric vector of parameters.
#' @param model Character string specifying the model.
#' @param state.names Optional character vector of length 2 for state labels.
#'
#' @return Unnamed numeric vector in the correct positional order.
#'
#' @keywords internal
.resolve_named_pars <- function(pars, model, state.names) {
  if (is.null(model)) {
    # User-provided Q-matrix: just need root.chrom
    if ("root.chrom" %in% names(pars)) return(unname(pars["root.chrom"]))
    return(unname(pars))
  }

  if (model == "2010") {
    canonical <- c("asc", "desc", "dem", "pol", "root.chrom")
    if (all(canonical %in% names(pars))) {
      return(unname(pars[canonical]))
    }
    return(unname(pars))
  }

  if (model %in% c("ChromPlus", "SAF")) {
    # Canonical names (numeric suffix)
    canonical <- c("asc1", "asc2", "desc1", "desc2",
                    "dem1", "dem2", "pol1", "pol2",
                    "tran12", "tran21", "root.chrom", "root.state")

    # Try canonical names first
    if (all(canonical %in% names(pars))) {
      return(unname(pars[canonical]))
    }

    # Try state.names-based labels
    if (!is.null(state.names)) {
      s1 <- state.names[1]
      s2 <- state.names[2]
      labeled <- c(paste0("asc.", s1), paste0("asc.", s2),
                    paste0("desc.", s1), paste0("desc.", s2),
                    paste0("dem.", s1), paste0("dem.", s2),
                    paste0("pol.", s1), paste0("pol.", s2),
                    paste0("tran.", s1, ".to.", s2),
                    paste0("tran.", s2, ".to.", s1),
                    "root.chrom", "root.state")
      if (all(labeled %in% names(pars))) {
        return(unname(pars[labeled]))
      }
    }

    # Fall through to positional
    return(unname(pars))
  }

  if (model == "PloidEvol") {
    canonical <- c("asc1", "asc2", "desc1", "desc2",
                    "dem1", "dem2", "pol1", "pol2",
                    "redip", "root.chrom", "root.state")

    if (all(canonical %in% names(pars))) {
      return(unname(pars[canonical]))
    }

    if (!is.null(state.names)) {
      s1 <- state.names[1]
      s2 <- state.names[2]
      labeled <- c(paste0("asc.", s1), paste0("asc.", s2),
                    paste0("desc.", s1), paste0("desc.", s2),
                    paste0("dem.", s1), paste0("dem.", s2),
                    paste0("pol.", s1), paste0("pol.", s2),
                    "redip", "root.chrom", "root.state")
      if (all(labeled %in% names(pars))) {
        return(unname(pars[labeled]))
      }
    }

    return(unname(pars))
  }

  # Unknown model, return as-is
  return(unname(pars))
}
