#' Generate Starting Values from a Distribution
#'
#' Samples starting values from a uniform or normal distribution for use as
#' initial parameter values in likelihood optimization or MCMC analysis.
#'
#' @param n Integer. The number of starting values to generate.
#' @param min Numeric. The minimum value for a uniform distribution. Required
#'   when `dist = "unif"`.
#' @param max Numeric. The maximum value for a uniform distribution. Required
#'   when `dist = "unif"`.
#' @param mean Numeric. The mean for a normal distribution. Required when
#'   `dist = "norm"`.
#' @param sd Numeric. The standard deviation for a normal distribution. Required
#'   when `dist = "norm"`.
#' @param dist Character string specifying the distribution. Either `"unif"`
#'   (default) for uniform or `"norm"` for normal.
#'
#' @return A numeric vector of length `n` containing randomly sampled values.
#'
#' @seealso [constrainMkn()], [constrainMuSSE()] for functions that use
#'   starting values in model fitting.
#'
#' @examples
#' # Generate 5 uniform starting values between 0 and 1
#' startVals(n = 5, min = 0, max = 1)
#'
#' # Generate 5 normal starting values
#' startVals(n = 5, mean = 0.5, sd = 0.1, dist = "norm")
#'
#' @export
startVals <- function(n, min, max, mean, sd, dist="unif"){
  if(!dist %in% c("unif", "norm")) stop("Specified distribution not implemented")
  if(dist=="unif") start.vals <- runif(n, min=min, max=max)
  if(dist=="norm") start.vals <- rnorm(n, mean=mean, sd=sd)
  return(start.vals)
}
