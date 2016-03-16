# this makes the data matrix for chromevol analysis
# x: is a dataframe col 1 is names col 2 is chromosome
# counts.
# excess: this how far over highest count we should
# include in the model for instance we allow that
# the highes number on phylogeny may be some point
# in the past rather than an extant lineage
# polyploids: if T then an unobserved state for polypoidy
# is created if F then basically you have the older style 
# chromevol.
datatoMatrix <- function(x, excess = .1, hidden = T){
  matsize <- round(max(as.numeric(x[, 2])) * (1 + excess))
  if(hidden == T){
    dmat <- matrix(0, nrow(x), matsize * 2)
    states <- c(as.character(1:matsize),
                paste(as.character(1:matsize), "p", sep = ""))
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    for(i in 1:nrow(x)){
      dmat[i, which(colnames(dmat) == x[i, 2])] <- .5
      dmat[i, which(colnames(dmat) == x[i, 2]) + matsize] <- .5
    }
  }
  if(hidden == F){
    dmat <- matrix(0, nrow(x), matsize)
    states <- as.character(1:matsize)
    colnames(dmat) <- states
    row.names(dmat) <- x[, 1]
    for(i in 1:nrow(x)){
      dmat[i, which(colnames(dmat) == x[i, 2])] <- 1
    }
  }
  return(dmat)
}