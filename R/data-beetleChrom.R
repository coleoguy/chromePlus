#' Simulated Beetle Karyotype Dataset
#'
#' A non-trivial example dataset for illustrating end-to-end chromePlus
#' workflows. It represents a simulated 100-species beetle clade with joint
#' evolution of haploid chromosome number and a binary sex-chromosome-system
#' trait ("standard" vs "fused"). Parameter values were chosen to produce a
#' realistic beetle-like karyotype distribution (haploid numbers 7-13, modal
#' around 9-10) while keeping model fits tractable on a laptop.
#'
#' The data were generated with [makeSSEchrom()]; see
#' `beetleChrom$metadata$makeSSEchrom.call` for the exact call. Tip labels
#' combine real beetle genus names with random species epithets and should
#' not be read as real taxa.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{tree}{A phylogenetic tree of class `"phylo"` with 100 tips.}
#'   \item{data}{A data frame with 100 rows and four columns:
#'     \describe{
#'       \item{species}{Character. Tip label matching `tree$tip.label`.}
#'       \item{chrom}{Integer. Haploid chromosome number (7-13).}
#'       \item{state}{Character. Sex-chromosome-system state: `"standard"`
#'         or `"fused"`.}
#'       \item{prob}{Numeric. Probability of being in the `"standard"` state
#'         (1 or 0 here because the simulation gives deterministic tip
#'         states; real datasets would use values in `[0, 1]`).}
#'     }}
#'   \item{metadata}{A list documenting provenance: description, source,
#'     random seed, and the full `makeSSEchrom()` call used to generate the
#'     data.}
#' }
#'
#' @source Simulated with [makeSSEchrom()]. Not real biological data.
#'
#' @examples
#' data(beetleChrom)
#' str(beetleChrom, max.level = 2)
#' table(beetleChrom$data$state)
#' table(beetleChrom$data$chrom)
#' \donttest{
#' # Minimal end-to-end fit on your laptop (~a few seconds)
#' library(diversitree)
#' dat.mat <- datatoMatrix(beetleChrom$data[, c("species", "chrom", "prob")],
#'                         hyper = TRUE,
#'                         state.names = c("standard", "fused"))
#' lik <- make.mkn(beetleChrom$tree, states = dat.mat, k = ncol(dat.mat),
#'                 strict = FALSE, control = list(method = "ode"))
#' con <- constrainMkn(data = dat.mat, lik = lik, hyper = TRUE,
#'                     state.names = c("standard", "fused"),
#'                     constrain = list(drop.demi = TRUE, drop.poly = TRUE))
#' diversitree::argnames(con)
#' }
"beetleChrom"
