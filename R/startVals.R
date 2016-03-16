# Heath Blackmon
# coleoguy@gmail.com
# this contains helper functions for using 
# diversitree to implement chromevol models

# takes n for the number of values needed and 
# samples from a uniform distribution between l and h
startVals <- function(n, min, max, mean, sd, dist="unif"){
  if(!dist %in% c("unif", "norm")) stop("Specified distribution not implemented")
  if(dist=="unif") start.vals <- runif(n, min=min, max=max)
  if(dist=="norm") start.vals <- rnorm(n, mean=mean, sd=sd)
  return(start.vals)
}

