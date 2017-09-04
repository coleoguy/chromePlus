## TODO FIX third column checks so ignored when 
## hyper = F

# this makes the data matrix for chromevol analysis
# x: is a dataframe col 1 is names col 2 is chromosome
# counts column 3 is probability in metastate 1.
# excess: this how far over highest count we should
# include in the model for instance we allow that
# the highes number on phylogeny may be some point
# in the past rather than an extant lineage
# polyploids: if T then an unobserved state for polypoidy
# is created if F then basically you have the older style 
# chromevol.
datatoMatrix <- function(x, range=NULL, hyper = T, buffer = 4){
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