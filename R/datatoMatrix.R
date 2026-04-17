#' Convert Tip Data to a Probability Matrix
#'
#' Converts a data frame of species names, chromosome counts, and (optionally)
#' binary trait probabilities into a probability matrix suitable for use with
#' diversitree's `make.mkn()` or `make.musse()` functions.
#'
#' When `hyper = TRUE`, the matrix is doubled to accommodate two hidden states
#' (hyperstates). Column 3 of the input data is used to assign the probability
#' of being in hyperstate 1 vs. hyperstate 2 for each species. When
#' `hyper = FALSE`, each species is assigned probability 1 at its observed
#' chromosome number.
#'
#' @param x A data frame with at least 2 columns:
#'   \describe{
#'     \item{Column 1}{Character. Species names.}
#'     \item{Column 2}{Numeric. Chromosome counts (will be rounded to whole
#'       numbers if non-integer).}
#'     \item{Column 3}{Numeric. Probability of being in the FIRST binary
#'       state (i.e. the state whose columns appear in the left half of the
#'       returned matrix, and which downstream chromePlus functions refer to
#'       with a `1` suffix or with `state.names[1]`). Required when
#'       `hyper = TRUE`. Values should be between 0 and 1; the second-state
#'       probability is computed as `1 - prob`.}
#'   }
#' @param range Numeric vector of length 2 giving the minimum and maximum
#'   chromosome numbers to include in the model, e.g. `c(5, 20)`. Defaults to
#'   `NULL`, in which case the range is determined from the observed data plus
#'   any `buffer`.
#' @param hyper Logical. If `TRUE` (default), a binary hyperstate is included
#'   and the matrix is doubled in width.
#' @param buffer Integer. Number of chromosome states to add above and below the
#'   observed range. Only used when `range = NULL`. Defaults to `0`.
#' @param state.names Optional character vector of length 2 giving descriptive
#'   names for the two binary states. When provided, columns for the second
#'   hyperstate are suffixed with `state.names[2]` (joined by `.`) instead of
#'   `"h"`, making the output self-documenting. `state.names[1]` corresponds
#'   to the left half / `prob` column from `x`. Defaults to `NULL`.
#'
#' @return A numeric matrix with species as rows and chromosome states as
#'   columns. Column names are the chromosome numbers (and, when
#'   `hyper = TRUE`, chromosome numbers suffixed with `"h"` for the second
#'   hyperstate, or with `state.names[2]` if `state.names` is provided). Row
#'   names are the species names from `x`.
#'
#' @seealso [constrainMkn()], [constrainMuSSE()] for constraining likelihood
#'   functions built from this matrix.
#'
#' @examples
#' # Simple example without hyperstate
#' dat <- data.frame(
#'   species = c("sp1", "sp2", "sp3"),
#'   chrom = c(5, 7, 10),
#'   prob = c(1, 1, 1)
#' )
#' datatoMatrix(x = dat, hyper = FALSE)
#'
#' # With hyperstate and uncertainty in binary trait
#' dat2 <- data.frame(
#'   species = c("sp1", "sp2", "sp3"),
#'   chrom = c(5, 7, 10),
#'   prob = c(0.8, 0.5, 1.0)
#' )
#' datatoMatrix(x = dat2, hyper = TRUE)
#'
#' @export
datatoMatrix <- function(x, range = NULL, hyper = TRUE, buffer = 0,
                         state.names = NULL){
  if (!is.null(state.names)) {
    if (!is.character(state.names) || length(state.names) != 2) {
      stop("state.names must be a character vector of length 2")
    }
  }

  if (sum(x[, 2] != round(x[, 2])) > 0) {
    message("Chromosome numbers are being rounded to nearest whole number")
    x[, 2] <- round(x[, 2])
  }

  # automate the range argument based on input data
  if (is.null(range)) {
    low <- min(x[, 2]) - buffer
    high <- max(x[, 2]) + buffer
    if (low < 1) low <- 1
    range <- c(low, high)
  }

  matsize <- range[2] - range[1] + 1
  # Single match() per species, then vectorised cbind-index assignment
  col_idx <- match(as.character(x[, 2]), as.character(range[1]:range[2]))
  if (anyNA(col_idx)) {
    bad <- x[is.na(col_idx), 2]
    stop("chromosome number(s) ", paste(unique(bad), collapse = ", "),
         " fall outside range [", range[1], ", ", range[2], "]")
  }
  row_idx <- seq_len(nrow(x))

  if (hyper == TRUE) {
    dmat <- matrix(0, nrow(x), matsize * 2)
    hyper.suffix <- if (is.null(state.names)) "h" else paste0(".", state.names[2])
    states <- c(as.character(range[1]:range[2]),
                paste(as.character(range[1]:range[2]), hyper.suffix, sep = ""))
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    dmat[cbind(row_idx, col_idx)]           <- x[, 3]
    dmat[cbind(row_idx, col_idx + matsize)] <- 1 - x[, 3]
  }
  if (hyper == FALSE) {
    dmat <- matrix(0, nrow(x), matsize)
    states <- as.character(range[1]:range[2])
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    dmat[cbind(row_idx, col_idx)] <- 1
  }
  return(dmat)
}
