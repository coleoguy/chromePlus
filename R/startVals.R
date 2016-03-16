# Heath Blackmon
# coleoguy@gmail.com
# this contains helper functions for using 
# diversitree to implement chromevol models

# takes n for the number of values needed and 
# samples from a uniform distribution between l and h
startVals <- function(n, l, h){
  start.vals <- runif(n, min=l, max=h)
  return(start.vals)
}

