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
    parMat <- matrix(0, 2*length(limits[1]:limits[2]), 2*length(limits[1]:limits[2]))
    colnames(parMat) <- sprintf(paste('%0', pad, 'd', sep=""), 1:ncol(parMat))
    rownames(parMat) <- colnames(parMat)
    # now we have a matrix with all zeros but the right state names
    # in the column and row names
    
    # we need to know where our duplication 
    # of the matrix begins so here is that
    split <- ncol(parMat)/2
    
    
    
    
    q <- matrix(0, length(limits[1]:limits[2]), length(limits[1]:limits[2]))
    rownames(q) <- colnames(q) <- chroms <- limits[1]:limits[2] 
    
    
    
    
    # MODEL 2 PLOIDY IS NOT THE HYPER STATE
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
  }

  
  
    diag(q) <- -rowSums(q)
    # simulate the chromosome numbers
    # return polyploidy state
    dsims <- geiger::sim.char(tree, q, model="discrete", root=root)[,,1]
    tips <- names(dsims)
    dsims <- as.numeric(colnames(q)[dsims])
    names(dsims) <- tips
    return(dsims)
}