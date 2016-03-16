
# returns a constrained diversitree likelihood 
# equation
#
# arguments 
# data: this is the data matrix in the format for
# for make.mkn in diversitree
# polyploidy is present
#
# rates that are implemented include
# rate1 ascending aneuploidy - diploid      ascdip
# rate2 descending aneuploidy - diploid     descdip
# rate5 polyploidization                    polypl
#
# rates 3,4,6 are only present when hidden = T

# can additional constraints can be added after this
# function by using the normal constrain approach
constrainMkn <- function(data, lik){
  # create and store variable for padding rate names
  if(ncol(data) < 100) pad <- 2
  if(ncol(data) >= 100) pad <- 3
  if(ncol(data) < 10) pad <- 1
  # make the matrix of rates
  parMat <- matrix(0,ncol(data),ncol(data))
  # make the components of the rate names the column and row
  # names this will allow for easy creation of constraints later
  colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep=""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  # now we have a matrix with all zeros but the right state names
  # in the column and row names
  # this will allow us to fit the old fashioned chromevol 
  # model just as easily with very little overhead
    for(i in 1:(nrow(parMat) - 1)){
      if((i * 2) <= ncol(parMat)) parMat[i, i*2] <- 5 #polyploidy
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
    }
    split <- ncol(parMat)/2
    for(i in 1:(split - 1)){
      if((i * 2) <= split) parMat[i, (i * 2 + split)] <- 5 #polyploidy
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - diploids
      parMat[(i + 1), i] <- 2 #descending aneuploidy - diploids
    }
  # each of these vectors will hold the formulae for that class of
  # parameters (described up at the top)
  restricted <- ascdip <- descdip <- ascpol <- descpol <- polypl <- redip <- vector()
  for(i in 1:nrow(parMat)){ # by rows then
    for(j in 1:ncol(parMat)){ # by cols
      if(parMat[i, j] == 0 & i != j){
        restricted <- c(restricted, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 0", sep="" ))
      }
      if(parMat[i, j] == 1){
        ascdip <- c(ascdip, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ ascdip", sep="" ))
      }
      if(parMat[i, j] == 2){
        descdip <- c(descdip, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ descdip", sep="" ))
      }
      if(parMat[i, j] == 5){
        polypl <- c(polypl, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ polypl", sep="" ))
      }
      if(parMat[i, j] == 9){
        spec <- c(spec, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 'polypl'+'ascdip'", sep="" ))
      }
    }
  }
  # lets store these in realy obvious names
  formulae <- c(restricted, ascdip, descdip, polypl)
  extras <- c("restricted", "ascdip", "descdip", "polypl")
  lik.con <- constrain(lik, formulae=formulae, extra=extras)
  return(lik.con)
}