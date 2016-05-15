
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

constrainMkn <- function(data, lik, hidden = T, 
                           s.lambda = T, s.mu = T, polyploidy = T, verbose=F){
  
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
    print("Constraining model to simple chromevol version")
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
    print("Constraining model where ploidy is a meta state and different rates of chromosome evolution are possible based on being polyploid or diploid")
    # diploid rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - diploids
      parMat[(i + 1), i] <- 2 #descending aneuploidy - diploids
      if((chroms[i] * 2) <= max(chroms)) parMat[i, (which(chroms[i] * 2 == chroms) + split)] <- 5 #polyploidy-1
    }
    # polyploid rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - polyploids
      parMat[(i + 1), i] <- 4 #descending aneuploidy - polyploids
      if((chroms[i-split] * 2) <= max(chroms)) parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      parMat[i, (i - split)] <- 7 #rediploidization
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 7 #rediploidization
    }
  }
  
  # MODEL 2 PLOIDY IS NOT THE HYPER STATE
  if(hidden==T & polyploidy == F){
    print("Constraining model with a hyper state that may have different rates of chromsome number evolution")
    # state 1 rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - 1
      parMat[(i + 1), i] <- 2 #descending aneuploidy - 1
      if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms==(chroms[i]*2))] <- 5 #polyploidy - 1
      parMat[i, (i+split)] <- 8 # transitions state 1->2
      # special case for last row
      if(i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- 8 # transitions state 2->1
      
    }
    # state 2 rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - 2
      parMat[(i + 1), i] <- 4 #descending aneuploidy - 2
      if((chroms[i-split] * 2) <= max(chroms)) parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      parMat[i, (i - split)] <- 9 #transition state 2->1
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 9 # transitions state 2->1
    }
  }
  
  
  
  
  
  
  # we now have a matrix with a number 1-9 that matches the rates present
  # under one of our models we will use this to build our 
  # arguments for the standard diversitree constrain function
  #
  # each of these vectors will hold the formulae for that class of
  # parameters (described up at the top)
  restricted <- asc1 <- desc1 <- asc2 <- desc2 <- pol1 <- pol2 <- redip <- tran12 <- tran21 <- vector()
  for(i in 1:nrow(parMat)){ # by rows then
    for(j in 1:ncol(parMat)){ # by cols
      if(parMat[i, j] == 0 & i != j){
        restricted <- c(restricted, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ 0", sep="" ))
      }
      if(parMat[i, j] == 1){
        asc1 <- c(asc1, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ asc1", sep="" ))
      }
      if(parMat[i, j] == 2){
        desc1 <- c(desc1, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ desc1", sep="" ))
      }
      if(parMat[i, j] == 3){
        asc2 <- c(asc2, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ asc2", sep="" ))
      }
      if(parMat[i, j] == 4){
        desc2 <- c(desc2, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ desc2", sep="" ))
      }
      if(parMat[i, j] == 5){
        pol1 <- c(pol1, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ pol1", sep="" ))
      }
      if(parMat[i, j] == 6){
        pol2 <- c(pol2, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ pol2", sep="" ))
      }
      if(parMat[i, j] == 7){
        redip <- c(redip, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ redip", sep="" ))
      }
      if(parMat[i, j] == 8){
        tran12 <- c(tran12, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ tran12", sep=""))
      }
      if(parMat[i, j] == 9){
        tran21 <- c(tran21, paste("q", row.names(parMat)[i], colnames(parMat)[j], " ~ tran21", sep="" ))
      }
    }
  }
  
  
  
  # lets store these in realy obvious names
  formulae <- c(restricted, asc1, desc1, asc2, desc2, pol1, pol2, redip, tran12, tran21)
  extras <- c("restricted", "asc1", "desc1", "asc2", "desc2", 
              "pol1", "pol2", "redip", "tran12", "tran21")
  lik.con <- constrain(lik, formulae=formulae, extra=extras)
  if(verbose==T) return(list(lik.con, parMat))
  if(verbose==F) return(lik.con)
}


