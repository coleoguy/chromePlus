
# rates that are implemented include
# rate1 ascending aneuploidy - diploid      asc1
# rate2 descending aneuploidy - diploid     desc1
# rate3 ascending aneuploidy - polyploid    asc2
# rate4 descending aneuploidy - polyploid   desc2
# rate5 polyploidization of diploid          poly1
# rate6 polploidization of polyploid         poly2
# rate7 rediploidization                    redip
# rate8 rediploidization                    tran12
# rate9 rediploidization                    tran21
#

# can additional constraints can be added after this
# function by using the normal constrain approach

constrainMuSSE <- function(data, lik, hidden = T, 
                           s.lambda = T, s.mu = T, polyploidy = T){
  
  ## BUILD AN EMPTY MATRIX MATCHING OUR MODEL
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
  
  # we need to know where our duplication 
  # of the matrix begins so here is that
  split <- ncol(parMat)/2
  
  # we also need the actual chromosome numbers
  if(hidden==T) chroms <- as.numeric(colnames(chrom.matrix)[1:split])
  if(hidden==F) chroms <- as.numeric(colnames(chrom.matrix))
  
  ## NOW WE HAVE A SERIES OF LOOPS THAT FILL IN OUR parMAT 
  ## MATRIX WITH NUMBERS 1:9 INDICATIVE OF THE DIFFERENT POSSIBLE
  ## RATES WE WISH TO INCLUDE IN OUR MODEL.  EACH OF THESE LOOPS
  ## REPRESENT A DIFFERENT MODEL OF CHROMOSOME EVOLUTION
  
  ## OLD CRHOMEVOL MODEL
  if(hidden==F){
    for(i in 1:(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
      if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms==(chroms[i]*2))] <- 5 #polyploidy
    }
    # currently this has the issue of missing polyploidy for q12
    # this transition should be = ascending + polyploidy this should
  }
  
  # MODEL 1 PLOIDY IS HIDDEN STATE
  if(hidden==T & polyploidy == T){
    # diploid rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - diploids
      parMat[(i + 1), i] <- 2 #descending aneuploidy - diploids
      ###THIS IS WHERE WE NEED TO START FIXING THINGS USING THIS NEW CHROMS VECTOR
      if((i * 2) <= split) parMat[i, (i * 2 + split)] <- 5 #polyploidy
    }
    # polyploid rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - polyploids
      parMat[(i + 1), i] <- 4 #descending aneuploidy - polyploids
      
      if(  (split + ((i-split)*2)) <= ncol(parMat)     ) parMat[i, (i * 2 + split)] <- 5 #polyploidy
      
      
      parMat[i, (i - split)] <- 7 #rediploidization
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 7 #rediploidization
    }
  }
  
  # this is the new model with hidden state not being ploidy
  if(hidden==T & polyploidy == F){
    split <- ncol(parMat)/2
    for(i in 1:(split - 1)){
      if((i * 2) <= split) parMat[i, (i * 2)] <- 5 #polyploidy state 1
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - state 1
      parMat[(i + 1), i] <- 2 #descending aneuploidy - state 1
    }
    for(i in (split + 1):(nrow(parMat) - 1)){
      if(split + (i-split)*2 <= ncol(parMat)) parMat[i, (split+2*(i-split))] <- 5 #polyploidy state 2
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - state 2
      parMat[(i + 1), i] <- 4 #descending aneuploidy - state 2
      parMat[i, (i - split)] <- 6 # state 2 to state 1
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 6 #rediploidization
    }
  }
  # we now have a matrix with a number 1-6 that matches the rates present
  # under either of our models hidden or not we will use this to build our 
  # arguments for the standard diversitree constrain function
  #
  # each of these vectors will hold the formulae for that class of
  # parameters (described up at the top)
  restricted <- ascdip <- descdip <- ascpol <- descpol <- polypl <- redip <- lambda <- mu <- vector()
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
      if(parMat[i, j] == 3){
        ascpol <- c(ascpol, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ ascpol", sep="" ))
      }
      if(parMat[i, j] == 4){
        descpol <- c(descpol, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ descpol", sep="" ))
      }
      if(parMat[i, j] == 5){
        polypl <- c(polypl, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ polypl", sep="" ))
      }
      if(parMat[i, j] == 6){
        redip <- c(redip, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ redip", sep="" ))
      }
      if(parMat[i, j] == 9){
        spec <- c(spec, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 'polypl'+'ascdip'", sep="" ))
      }
    }
  }
  
  # since we are using musse now we have to build in the 
  # speciation and extinction parameters.  It seems that 
  # the only sensible way to deal with these is to either
  # have a single rate for everyone or two have seperate 
  # rates for diploids and polyploids.
  for(i in 1:nrow(parMat)){
    # no hidden state
    if(hidden==F){
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ lambda1", sep = ""))
      mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ mu1", sep = ""))
    }
    # hidden model with single lambda
    if(hidden == T & s.lambda == T){
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ lambda1", sep = ""))
    }
    # hidden model with single mu
    if(hidden == T & s.mu == T){
      mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ mu1", sep = ""))
    }
    # hidden model with two lambdas
    if(hidden == T & s.lambda == F){
      if(i <= split){
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ lambda1", sep = ""))
      }
      if(i > split){
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], " ~ lambda2", sep = ""))
      }
    }
    # hidden model with two mus
    if(hidden == T & s.mu == F){
      if(i <= split){
        mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ mu1", sep = ""))
      }
      if(i > split){
        mu <- c(mu, paste("mu", colnames(parMat)[i], " ~ mu2", sep = ""))
      }
    }
  }
  
  # lets store these in realy obvious names
  formulae <- c(restricted, ascdip, descdip, ascpol, descpol, polypl, redip, lambda, mu)
  extras <- c("restricted", "ascdip", "descdip", "ascpol", "descpol", "polypl", "redip", "lambda1", "mu1", "lambda2", "mu2")
  lik.con <- constrain(lik, formulae=formulae, extra=extras)
  return(lik.con)
}


