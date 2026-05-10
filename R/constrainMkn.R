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
#' # With descriptive state names. Note: make.mkn requires deterministic
#' # tip-state assignments (prob must be 0 or 1). For probabilistic tip
#' # states use make.musse / constrainMuSSE instead.
#' dat2 <- data.frame(
#'   species = paste0("sp", 1:5),
#'   chrom = c(5, 6, 7, 8, 10),
#'   prob = c(1, 0, 1, 0, 1)
#' )
#' dat.mat2 <- datatoMatrix(x = dat2, hyper = TRUE,
#'                          state.names = c("matched", "mismatched"))
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
  } else if (isTRUE(hyper)) {
    message("constrainMkn: state.names not provided. Output rates will use ",
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
  if (is.null(constrain$saf.model)) constrain$saf.model <- FALSE
  if (is.null(constrain$sym.hyperstates)) constrain$sym.hyperstates <- FALSE
  if (is.null(constrain$meta)) constrain$meta <- "ARD"

  # Build the rate-name parMat and rate.table via the shared helper. The
  # same helper is used by make.mkn.fast(), so the model-construction logic
  # lives in one place (R/makeMknFast.R, .chromeplus_build_parMat).
  built <- .chromeplus_build_parMat(
    data = data, hyper = hyper, polyploidy = polyploidy,
    oneway = oneway, state.names = state.names, constrain = constrain,
    emit_messages = TRUE
  )
  parMat     <- built$parMat
  rate.table <- built$rate.table

  formulae <- paste0("q", rate.table[, 1], rate.table[, 2], " ~ ",
                     rate.table[, 3])

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

  lik.con <- constrain(lik, formulae = formulae, extra = extras)

  if (verbose == TRUE) {
    result.list <- list(lik.con, parMat, rate.table)
    names(result.list) <- c("likelihood function", "parameter matrix", "ratetable")
    return(result.list)
  }
  if (verbose == FALSE) return(lik.con)
}
