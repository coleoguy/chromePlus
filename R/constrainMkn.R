
# returns a constrained diversitree likelihood 
# equation
#
# arguments 
# data: this is the data matrix in the format for
# for make.mkn in diversitree
#
# rates that are implemented include
# rate1 ascending aneuploidy - diploid      ascdip
# rate2 descending aneuploidy - diploid     descdip
# rate5 polyploidization                    polypl
#
constrainMkn <- function(data, lik, model="single"){
  
  ### setup
  
  #diversitree needs states 1->n but we need to keep track
  #of actual numbers so that is the first step
  if(model=="single"){
    chrom.numbs <- colnames(data)
  }
  if(model=="hyper"){
    split <- ncol(data)/2
    chrom.numbs <- c(colnames(data)[1:split], colnames(data)[1:split])
  }
  # create and store variable for padding rate names
  
  if(ncol(data) < 100) pad <- 2
  if(ncol(data) >= 100) pad <- 3
  if(ncol(data) < 10) pad <- 1
  
  # make a matrix of rates
  parMat <- matrix(0,ncol(data),ncol(data))
  # make the components of the rate names the column and row
  # names this will allow for easy creation of constraints later
  colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep=""), 1:ncol(data))
  rownames(parMat) <- colnames(parMat)
  # now we have a matrix with all zeros but the right state names
  # in the column and row names
  ### Now we use our chrom.numbs to lookup values and fill in the matrix
  if(model=="single"){
    for(i in 1:(nrow(parMat) - 1)){
      if((as.numeric(chrom.numbs)[i] * 2) <= max(as.numeric(chrom.numbs))){
        parMat[i, i+as.numeric(chrom.numbs[i])] <- 5 #polyploidy
      } 
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
    }
    # special case for 1->2
    if(as.numeric(chrom.numbs[1]) == 1) parMat[1,2] <- 9
  }
  
  if(model=="hyper"){
    for(i in 1:(nrow(parMat) - 1)){
      if(i < split){
        if((as.numeric(chrom.numbs)[i] * 2) <= max(as.numeric(chrom.numbs))){
          parMat[i, i+as.numeric(chrom.numbs[i])] <- 5 #polyploidy
        } 
        parMat[i, (i + 1)] <- 1 #ascending aneuploidy
        parMat[(i + 1), i] <- 2 #descending aneuploidy
      }
      if(i >= split){
        if((as.numeric(chrom.numbs)[i] * 2) <= max(as.numeric(chrom.numbs))){
          parMat[i, i+as.numeric(chrom.numbs[i])] <- 6 #polyploidy
        } 
        parMat[i, (i + 1)] <- 3 #ascending aneuploidy
        parMat[(i + 1), i] <- 4 #descending aneuploidy
      }
      # special case for 1->2
      if(as.numeric(chrom.numbs[1]) == 1){
        parMat[1,2] <- 9
        parmat[(1+split),(2+split)] <- 10
      }
    }
  }
  if(model=="ploidy"){
    stop("Ploidy model not finished for mkn function")
  }
  
  
  # each of these vectors will hold the formulae for that class of
  # parameters (described up at the top)
  
  # create arguments for constrain
  restricted <- 
    asc   <- desc   <- polypl   <- spec <- 
    asc.h <- desc.h <- polypl.h <- spec.h <- vector()
  
  for(i in 1:nrow(parMat)){ # by rows then
    for(j in 1:ncol(parMat)){ # by cols
      if(parMat[i, j] == 0 & i != j){
        restricted <- c(restricted, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 0", sep="" ))
      }
      if(parMat[i, j] == 1){
        asc <- c(asc, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ asc", sep="" ))
      }
      if(parMat[i, j] == 2){
        desc <- c(desc, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ desc", sep="" ))
      }
      if(parMat[i, j] == 5){
        polypl <- c(polypl, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ polypl", sep="" ))
      }
      if(parMat[i, j] == 9){
        spec <- c(spec, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 'polypl'+'ascdip'", sep="" ))
      }
      
      if(parMat[i, j] == 1){
        asc.h <- c(asc.h, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ asc.h", sep="" ))
      }
      if(parMat[i, j] == 2){
        desc.h <- c(desc.h, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ desc.h", sep="" ))
      }
      if(parMat[i, j] == 5){
        polypl.h <- c(polypl.h, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ polypl.h", sep="" ))
      }
      if(parMat[i, j] == 9){
        spec.h <- c(spec.h, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 'polypl.h'+'asc.h'", sep="" ))
      }
    }
  }
  # lets store these in realy obvious names
  formulae <- c(restricted, asc,   desc,   polypl,   spec, 
                            asc.h, desc.h, polypl.h, spec.h)
  extras <- c("restricted", "asc",   "desc",   "polypl",   "spec", 
                            "asc.h", "desc.h", "polypl.h", "spec.h")
  lik.con <- constrain(lik, formulae=formulae, extra=extras)
  return(lik.con)
}