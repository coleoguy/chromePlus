simChrom <- function(tree, pars, limits){
  # args: tree, pars, limits
  # tree: phylo object
  # pars: c(ascending, descending, polyploidy, root)
  # limits: c(low, high)
  
  
  # set up an empty matrix
  q <- matrix(0, length(limits), length(limits))
  rownames(q) <- colnames(q) <- limits
  
  # fill in the matrix
  # the rows
  for(i in 1:nrow(q)){
    # the cols
    for(j in 1:ncol(q)){
      if(i - j == -1) q[i,j] <- pars[1] #increase
      if(i - j == 1) q[i,j] <- pars[2] #decrease
      if(as.numeric(colnames(q)[j]) / 
         as.numeric(rownames(q)[i]) == 2) q[i,j] <- pars[3] #polyploid
    }
  }
  diag(q) <- -rowSums(q)
  
  # simulate the chromosome numbers
  # return polyploidy state
  dsims <- geiger::sim.char(tree, q, model="discrete", root=pars[4])[,,1]
  tips <- names(dsims)
  dsims <- as.numeric(colnames(q)[dsims])
  names(dsims) <- tips
  return(dsims)
}