## TODO finish filling out the constrain argument
## I want it to be a list that could have the following elements

## drop.poly=T: polyoploidy is dropped from model
## drop.demi=T: demiploidy is dropped from model
## singlerate=T: gain and loss rate are equal for each meta state
## nometa=T: gain1 and gain2 are equal and loss1 and loss2 are equal
## meta="ARD": assumes q01 and q10 are different if "SYM" then single rate

# rates that are implemented include
# rate1 ascending aneuploidy - diploid      asc1
# rate2 descending aneuploidy - diploid     desc1
# rate3 ascending aneuploidy - polyploid    asc2
# rate4 descending aneuploidy - polyploid   desc2
# rate5 polyploidization of diploid          poly1
# rate6 polploidization of polyploid         poly2
# rate7 rediploidization of a polyploid
# rate8 transitions from 1 to 2 for hyperstate
# rate9 transitions from 2 to 1 for hyperstate
# rate10 demipolyploidy for state1            dem1 even
# rate11 demipolyploidy for state1            dem1 odd
# rate12 demipolyploidy for state2            dem2 even
# rate13 demipolyploidy for state2            dem2 odd

constrainMuSSE <- function(data, lik, hidden = T, s.lambda = T, s.mu = T, 
                           polyploidy = T, verbose=F, 
                           constrain=list(drop.poly=F, 
                                          drop.demi=F, 
                                          symmetric=F, 
                                          nometa=F, 
                                          meta="ARD")){
  # This fills out the list of constraints the default are no constraints
  if(length(constrain) < 5){
    if(is.null(constrain$drop.pol)) constrain$drop.poly=F
    if(is.null(constrain$drop.demi)) constrain$drop.demi=F
    if(is.null(constrain$symmetric)) constrain$symmetric=F
    if(is.null(constrain$nometa)) constrain$nometa=F
    if(is.null(constrain$meta)) constrain$meta="ARD"
  }

  ## BUILD AN EMPTY MATRIX MATCHING OUR MODEL
  # padding rate names
  if(ncol(data) < 100) pad <- 2
  if(ncol(data) >= 100) pad <- 3
  if(ncol(data) < 10) pad <- 1
  # make the matrix of rates
  parMat <- matrix(0,ncol(data),ncol(data))
  # make the components of the rate names the column and row
  # names this will allow for easy creation of constraints later
  colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep=""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)

  # we need to know where our duplication of the matrix begins so here is that
  split <- ncol(parMat)/2
  
  ## TODO CHANGE HIDDEN TO HYPER
  # we also need the actual chromosome numbers
  if(hidden==T) chroms <- as.numeric(colnames(data)[1:split])
  if(hidden==F) chroms <- as.numeric(colnames(data))
  
  ## NOW WE HAVE A SERIES OF LOOPS THAT FILL IN OUR parMat
  ## MATRIX WITH NUMBERS 1:13 INDICATIVE OF THE DIFFERENT POSSIBLE
  ## RATES WE WISH TO INCLUDE IN OUR MODEL.  EACH OF THESE LOOPS
  ## REPRESENT A DIFFERENT MODEL OF CHROMOSOME EVOLUTION
  
  ## OLD CRHOMEVOL MODEL
  if(hidden==F){
    print("Constraining model to simple chromevol version")
    for(i in 1:(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
      if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms==(chroms[i]*2))] <- 5 #polyploidy
      if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
        x <- chroms[i] * 1.5
        if(x %% 1 == 0)  parMat[i, which(chroms==x)] <- 10 #demiploidy state1 even
        if(x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
      }
    }
    # currently this has the issue of missing polyploidy for q12 should only be an issue when low chrom number is 1
    # this transition should be = ascending + polyploidy this should
  }
  
  # BiSCE MODEL 1 PLOIDY IS HIDDEN STATE
  if(hidden==T & polyploidy == T){
    print("Constraining model where ploidy is a meta state and different rates of chromosome evolution are possible based on being polyploid or diploid")
    # diploid rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - diploids
      parMat[(i + 1), i] <- 2 #descending aneuploidy - diploids
      if((chroms[i] * 2) <= max(chroms)) parMat[i, (which(chroms[i] * 2 == chroms) + split)] <- 5 #polyploidy-1
      # demiploidy
      if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
        x <- chroms[i] * 1.5
        if(x %% 1 == 0)  parMat[i, (which(chroms==x) + split)] <- 10 #demiploidy state1 even
        if(x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 11 #demiploidy state 1 odd
      }
    }
    # polyploid rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - polyploids
      parMat[(i + 1), i] <- 4 #descending aneuploidy - polyploids
      if((chroms[i-split] * 2) <= max(chroms)) parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      # demiploidy
      if((ceiling(chroms[i-split] * 1.5)) <= max(chroms)){
        x <- chroms[i-split] * 1.5
        if(x %% 1 == 0)  parMat[i, (which(chroms==x) + split)] <- 12 #demiploidy state1 even
        if(x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13 #demiploidy state 1 odd
      }
      parMat[i, (i - split)] <- 7 #rediploidization
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 7 #rediploidization
    }
  }
  
  # BiSPCE 2 PLOIDY IS NOT THE HYPER STATE
  if(hidden==T & polyploidy == F){
    print("Constraining model with a hyper state that may have different rates of chromsome number evolution")
    # state 1 rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy - 1
      parMat[(i + 1), i] <- 2 #descending aneuploidy - 1
      if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms==(chroms[i]*2))] <- 5 #polyploidy - 1
      # demiploidy
      if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
        x <- chroms[i] * 1.5
        if(x %% 1 == 0)  parMat[i, which(chroms==x)] <- 10 #demiploidy state1 even
        if(x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
      }
      parMat[i, (i+split)] <- 8 # transitions state 1->2
      # special case for last row
      if(i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- 8 # transitions state 2->1
      
    }
    # state 2 rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- 3 #ascending aneuploidy - 2
      parMat[(i + 1), i] <- 4 #descending aneuploidy - 2
      if((chroms[i-split] * 2) <= max(chroms)) parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- 6 #polyploidy-2
      # demiploidy
      if((ceiling(chroms[i-split] * 1.5)) <= max(chroms)){
        x <- chroms[i-split] * 1.5
        if(x %% 1 == 0)  parMat[i, (which(chroms==x) + split)] <- 12 #demiploidy state1 even
        if(x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 13 #demiploidy state 2 odd
      }
      parMat[i, (i - split)] <- 9 #transition state 2->1
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- 9 # transitions state 2->1
    }
  }
  # we now have a matrix with a number 1-13 that matches the rates present
  # under one of our models we will use this to build our 
  # arguments for the standard diversitree constrain function
  rate.table <- as.data.frame(matrix(,nrow(parMat) * ncol(parMat), 3))
  rate.table[, 1] <- rep(as.character(row.names(parMat)), each=ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[2], ]
  rate.table[rate.table[, 3] == 1, 3] <- "asc1"
  rate.table[rate.table[, 3] == 2, 3] <- "desc1"
  rate.table[rate.table[, 3] == 3, 3] <- "asc2"
  rate.table[rate.table[, 3] == 4, 3] <- "desc2"
  rate.table[rate.table[, 3] == 5, 3] <- "pol1"
  rate.table[rate.table[, 3] == 6, 3] <- "pol2"
  rate.table[rate.table[, 3] == 7, 3] <- "redip"
  rate.table[rate.table[, 3] == 8, 3] <- "tran12"
  rate.table[rate.table[, 3] == 9, 3] <- "tran21"
  rate.table[rate.table[, 3] == 10, 3] <- "dem1"
  rate.table[rate.table[, 3] == 11, 3] <- ".5*dem1"
  rate.table[rate.table[, 3] == 12, 3] <- "dem2"
  rate.table[rate.table[, 3] == 13, 3] <- ".5*dem2"
  
  
  if(constrain$nometa == T){
    rate.table[rate.table[, 3] == "asc2", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc1"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }
  
  if(constrain$drop.poly == T){
    rate.table[rate.table[, 3] == "pol1", 3] <- "0"
    rate.table[rate.table[, 3] == "pol2", 3] <- "0"
  }
  
  if(constrain$drop.demi == T){
    rate.table[rate.table[, 3] == "dem1", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }
  
  if(constrain$symmetric == T){
    rate.table[rate.table[, 3] == "desc1", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "asc2"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == 7, 3] <- "redip"
    rate.table[rate.table[, 3] == 8, 3] <- "tran12"
    rate.table[rate.table[, 3] == 9, 3] <- "tran21"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }

  if(constrain$meta == "SYM"){
    rate.table[rate.table[, 3] == "tran21", 3] <- "tran12"
  }
  
  

  formulae <- vector(mode="character", length=nrow(rate.table))
  for(i in 1:nrow(rate.table)){
    formulae[i] <- paste("q",
                         rate.table[i, 1], 
                         rate.table[i, 2],
                         " ~ ",
                         rate.table[i, 3],
                         collapse="", sep="")
  }
  
  ## This section could be sped up by getting rid of loop
  
  lambda <- mu <- vector()
  # Lambda and Mu
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
    extras <- c("asc1", "desc1", "asc2", "desc2", 
              "pol1", "pol2", "dem1", "dem2", "redip", "tran12", "tran21", 
              "lambda1", "mu1", "lambda2", "mu2")
  lik.con <- constrain(lik, formulae=c(formulae, lambda, mu), extra=extras)
  colnames(parMat) <- rownames(parMat) <- colnames(data)
  if(verbose==T) return(list(lik.con, parMat))
  if(verbose==F) return(lik.con)
}


