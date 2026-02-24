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
#'     \item{Column 3}{Numeric. Probability of being in binary state 1
#'       (required when `hyper = TRUE`). Values should be between 0 and 1.}
#'   }
#' @param range Numeric vector of length 2 giving the minimum and maximum
#'   chromosome numbers to include in the model, e.g. `c(5, 20)`. Defaults to
#'   `NULL`, in which case the range is determined from the observed data plus
#'   any `buffer`.
#' @param hyper Logical. If `TRUE` (default), a binary hyperstate is included
#'   and the matrix is doubled in width.
#' @param buffer Integer. Number of chromosome states to add above and below the
#'   observed range. Only used when `range = NULL`. Defaults to `0`.
#'
#' @return A numeric matrix with species as rows and chromosome states as
#'   columns. Column names are the chromosome numbers (and, when
#'   `hyper = TRUE`, chromosome numbers suffixed with `"h"` for the second
#'   hyperstate). Row names are the species names from `x`.
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
datatoMatrix <- function(x, range=NULL, hyper = T, buffer = 0){
  if(sum(x[,2] != round(x[,2])) > 0){
    print("Chromosome numbers are being rounded to nearest whole number")
    x[,2] <- round(x[,2])
  }

  # automate the range argument based on input data
  if(is.null(range)){
    low <- min(x[,2]) - buffer
    high <- max(x[,2]) + buffer
    if(low < 1) low <- 1
    range <- c(low, high)
  }


  matsize <- range[2]-range[1]+1
  if(hyper == T){
    dmat <- matrix(0, nrow(x), matsize * 2)
    states <- c(as.character(range[1]:range[2]),
                paste(as.character(range[1]:range[2]), "h", sep = ""))
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    for(i in 1:nrow(x)){
      dmat[i, which(colnames(dmat) == x[i, 2])] <- x[i,3]
      dmat[i, which(colnames(dmat) == x[i, 2]) + matsize] <- 1-x[i,3]
    }
  }
  if(hyper == F){
    dmat <- matrix(0, nrow(x), matsize)
    states <- as.character(range[1]:range[2])
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    for(i in 1:nrow(x)){
      dmat[i, which(colnames(dmat) == x[i, 2])] <- 1
    }
  }
  return(dmat)
}
