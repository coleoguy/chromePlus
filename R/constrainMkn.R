
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

# can additional constraints can be added after this
# function by using the normal constrain approach

constrainMkn <- function(data, 
                         lik, 
                         hyper = T, 
                         polyploidy = F, 
                         verbose=F, 
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
  if(hyper == T) chroms <- as.numeric(colnames(data)[1:split])
  if(hyper == F) chroms <- as.numeric(colnames(data))
  
  ## NOW WE HAVE A SERIES OF LOOPS THAT FILL IN OUR parMAT 
  ## MATRIX WITH NUMBERS 1:9 INDICATIVE OF THE DIFFERENT POSSIBLE
  ## RATES WE WISH TO INCLUDE IN OUR MODEL.  EACH OF THESE LOOPS
  ## REPRESENT A DIFFERENT MODEL OF CHROMOSOME EVOLUTION
  
  ## OLD CRHOMEVOL MODEL
  if(hyper == F){
    print("Constraining model to simple chromevol version")
    for(i in 1:(nrow(parMat) - 1)){
      if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms == (chroms[i] * 2))] <- 5 #polyploidy
      if(constrain$drop.demi == F){
        if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
          x <- chroms[i] * 1.5
          if(x %% 1 == 0)  parMat[i, which(chroms == x)] <- 10 #demiploidy state1 even
          if(x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11 #demiploidy state 1 odd
        }
      }
      parMat[i, (i + 1)] <- 1 #ascending aneuploidy
      parMat[(i + 1), i] <- 2 #descending aneuploidy
    }
  }
  
  # MODEL 1 PLOIDY IS HIDDEN STATE
  if(hyper==T & polyploidy == T){
    print("Creating rate matrix for chosen chromosome model")
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
  
  # MODEL 2 PLOIDY IS NOT THE HYPER STATE
  if(hyper==T & polyploidy == F){
    print("Creating rate matrix for chosen chromosome model")
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

  # lets store these in realy obvious names
  extras <- c("asc1", "desc1", 
              "asc2", "desc2", 
              "pol1", "pol2", 
              "redip", "tran12", "tran21", 
              "dem1", "dem2")
  
  lik.con <- constrain(lik, formulae=formulae, extra=extras)
  colnames(parMat) <- rownames(parMat) <- colnames(data)
  if(verbose==T){
    result.list <- list(lik.con, parMat)
    names(result.list) <- c("likelihood function", "parameter matrix")
    return(result.list)
  } 
  if(verbose==F) return(lik.con)
}


