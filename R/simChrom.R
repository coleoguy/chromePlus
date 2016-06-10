simChrom <- function(tree, pars, limits, model){
  # args: tree, pars, limits
  # tree: phylo object
  # pars: c(gain[1:2], loss[1:2], demi[1:2], poly[1:2], root)
  # limits: c(low, high)
  # model: "2010", "chromRate", "ploidEvol"
  
  
  # having a demi rate of 0 creates error
  #pars[pars == 0] <- .000000000000001
  # the simulation process is the same for all models so we seperate out
  # the process of generating the qmat but that it
  if(model == "2010"){
    root <- pars[5]
    if(length(pars) != 5) stop("pars should have length of 5")
    # set up an empty matrix
    q <- matrix(0, length(limits[1]:limits[2]), length(limits[1]:limits[2]))
    rownames(q) <- colnames(q) <- chroms <- limits[1]:limits[2] 
    # fill in the matrix
    for(i in 1:(nrow(q) - 1)){
      q[i, (i + 1)] <- pars[1] # gain
      q[(i + 1), i] <- pars[2] # loss
      if((chroms[i] * 2) <= max(chroms)) q[i, which(chroms==(chroms[i]*2))] <- pars[4] #poly
      if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
        x <- chroms[i] * 1.5
        if(x %% 1 == 0)  q[i, which(chroms==x)] <- q[i, which(chroms==x)] + pars[3] #demi even
        if(x %% 1 != 0)  q[i, which(chroms %in% c(floor(x), ceiling(x)))] <- q[i, which(chroms %in% c(floor(x), ceiling(x)))] + pars[3]/2 #demi odd
      }
    }
  }

  if(model == "chromRate"){
    root <- pars[11]
    if(length(pars) != 11) stop("pars should have length of 11")
    # set up an empty matrix
    if(limits[2] < 100) pad <- 2
    if(limits[2] >= 100) pad <- 3
    if(limits[2] < 10) pad <- 1
    
    parMat <- matrix(0, 2*length(limits[1]:limits[2]), 2*length(limits[1]:limits[2]))
    colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep=""), 1:ncol(parMat))
    rownames(parMat) <- colnames(parMat)
    split <- ncol(parMat)/2
    chroms <- limits[1]:limits[2]
      # state 1 rates
      for(i in 1:(split - 1)){
        parMat[i, (i + 1)] <- pars[1] #ascending aneuploidy - 1
        parMat[(i + 1), i] <- pars[3] #descending aneuploidy - 1
        if((chroms[i] * 2) <= max(chroms)) parMat[i, which(chroms==(chroms[i]*2))] <- parMat[i, which(chroms==(chroms[i]*2))] + pars[7] #polyploidy - 1
        # demiploidy
        if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
          x <- chroms[i] * 1.5
          if(x %% 1 == 0)  parMat[i, which(chroms==x)] <- parMat[i, which(chroms==x)] + pars[5] #demiploidy state1 even
          if(x %% 1 != 0)  parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] + (pars[5]/2) #demiploidy state 1 odd
        }
        parMat[i, (i+split)] <- pars[9] # transitions state 1->2
        # special case for last row
        if(i == (split - 1)) parMat[(i + 1), (i + 1 + split)] <- pars[9] # transitions state 1->2
        
      }
      # state 2 rates
      for(i in (split + 1):(nrow(parMat) - 1)){
        parMat[i, (i + 1)] <- pars[2] #ascending aneuploidy - 2
        parMat[(i + 1), i] <- pars[4] #descending aneuploidy - 2
        if((chroms[i-split] * 2) <= max(chroms)) parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] + pars[8] #polyploidy-2
        # demiploidy
        if((ceiling(chroms[i-split] * 1.5)) <= max(chroms)){
          x <- chroms[i-split] * 1.5
          if(x %% 1 == 0)  parMat[i, (which(chroms==x) + split)] <- parMat[i, (which(chroms==x) + split)] + pars[6] #demiploidy state1 even
          if(x %% 1 != 0)  parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] + (pars[6]/2) #demiploidy state 2 odd
        }
        parMat[i, (i - split)] <- pars[10] #transition state 2->1
        # special case for last row
        if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- pars[10] # transitions state 2->1
      }
    q <- parMat
    }

  
  
    diag(q) <- -rowSums(q)
    # simulate the chromosome numbers
    dsims <- geiger::sim.char(tree, q, model="discrete", root=root)[,,1]
    # save the names for various uses below
    tips <- names(dsims)
    # in the case of the 2010 model we have column names = to the
    # chromosome numbers so we can just use them
    if(model == "2010"){
      dsims[] <- as.numeric(colnames(q)[dsims])
      return(dsims)
    } 

    # under the chromRate model things are bit more complex and have
    # to be converted back to chromosome number and binary state
    if(model == "chromRate"){
      # for chromRate we need to return two vectors
      # 1) binary state 
      b.state <- rep(0, length(dsims))
      names(b.state) <- tips
      # if we are in binary state 1 add that in
      b.state[dsims > ncol(q) / 2] <- 1
      
      # 2) chromosome number
      dsims[] <- c(chroms, chroms)[dsims]
      result <- list(b.state, dsims)
      return(result)
    } 
}