simChrom <- function(tree, pars, limits = NULL, model = NULL, Qmat = NULL, verbose = F){
  # args: tree, pars, limits
  # tree: phylo object
  # pars: c(gain[1:2], loss[1:2], demi[1:2], poly[1:2], root)
  # limits: c(low, high), NULL
  # model: "2010", "ChromTrait", "PloidEvol", "SAF", NULL
  # Qmat: NULL (default) or user supplied Q-matrix
  # verbose: option to return parameter matrix

  if(is.null(model)==F && model == "2010"){
    print("building q-matrix")
    root <- pars[5] - limits[1] +1
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
      # special fix for chromosome num = 1
      diag(q) <- 0
    }
  }
  
  ## ChromPlus MODEL Q MATRIX
  if(is.null(model)==F && model == "ChromPlus"){
    print("building q-matrix")
    if(length(pars) != 12) stop("pars should have length of 12")
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
    if(pars[12] == 0){
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[11])])
    }
    if(pars[12] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    }
  }

  ## PloidEvol MODEL Q MATRIX
  if(is.null(model)==F && model == "PloidEvol"){
    print("building q-matrix")
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
    # diploidy rates
    for(i in 1:(split - 1)){
      parMat[i, (i + 1)] <- pars[1] #ascending aneuploidy - 1
      parMat[(i + 1), i] <- pars[3] #descending aneuploidy - 1
      #polyploidy - 1
      if((chroms[i] * 2) <= max(chroms)) parMat[i, (split + which(chroms==(chroms[i]*2)))] <- 
        pars[7] + parMat[i, which(chroms==(chroms[i]*2))]
      # demiploidy
      if((ceiling(chroms[i] * 1.5)) <= max(chroms)){
        x <- chroms[i] * 1.5
        #demiploidy state1 even
        if(x %% 1 == 0)  parMat[i, (split + which(chroms==x))] <- 
            pars[5] + parMat[i, (split + which(chroms==x))] 
        #demiploidy state 1 odd
        if(x %% 1 != 0)  parMat[i, (split + which(chroms %in% c(floor(x), ceiling(x))))] <- 
            (pars[5]/2) + parMat[i, (split + which(chroms %in% c(floor(x), ceiling(x))))]
      }
    }
    # polyploidy rates
    for(i in (split + 1):(nrow(parMat) - 1)){
      parMat[i, (i + 1)] <- pars[2] + parMat[i, (i + 1)] #ascending aneuploidy - 2
      parMat[(i + 1), i] <- pars[4] + parMat[(i + 1), i] #descending aneuploidy - 2
    }
    #polyploidy-2
    if((chroms[i-split] * 2) <= max(chroms)){
      parMat[i, (which(chroms[i-split] * 2 == chroms) + split)] <- 
        pars[8] + parMat[i, (which(chroms[i-split] * 2 == chroms) + split)]
    }
    # demiploidy-2
    if((ceiling(chroms[i-split] * 1.5)) <= max(chroms)){
      x <- chroms[i-split] * 1.5
      #demiploidy-2 even
      if(x %% 1 == 0){
        parMat[i, (which(chroms==x) + split)] <- 
          pars[6] + parMat[i, (which(chroms==x) + split)]
      }
      #demiploidy-2 odd
      if(x %% 1 != 0){
        parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)] <- 
          (pars[6]/2) + parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + split)]
      }
      # rediploidization
      parMat[i, (i - split)] <- pars[9] 
      # special case for last row
      if(i == (nrow(parMat) - 1)) parMat[(i + 1), (i + 1 - split)] <- pars[9]
    }
    q <- parMat 
    if(pars[11] == 0){
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[10])])
    }
    if(pars[11] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[10])
      root <- as.numeric(colnames(q)[x])
    }
  }
  
  ## SAF MODEL Q MATRIX
  if(is.null(model)==F && model == "SAF"){
    print("building q-matrix")
    if(length(pars) != 12) stop("pars should have length of 12")
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
      parMat[i, (i+split-1)] <- pars[9] # transitions state 1->2
      # special case for last row
      if(i == (split - 1)) parMat[(i + 1), (i + split)] <- pars[9] # transitions state 1->2
      
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
    if(pars[12] == 0){
      root <- as.numeric(colnames(q)[which(limits[1]:limits[2] == pars[11])])
    }
    if(pars[12] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    }
  }
  
  ## USER PROVIDED Q MATRIX
  if(is.null(model) == T){
    if(is.null(Qmat) == T) stop ("no model or q-matrix provided")
    
    if(is.null(limits) == F){
      print("ignoring provided limits")
    }
    
    if(is.matrix(Qmat) == F) stop ("q-matrix provided is not matrix")
    
    if(nrow(Qmat) != ncol(Qmat)) stop ("q-matrix provided has unequal number of rows and columns")
    
    if(setequal(rownames(q),colnames(q)) == F) stop ("q-matrix provided has nonmatching column and row names")
    
    if(length(pars) != 1) stop("pars should have length of 1")
    
    #Assign user supplied Q-matrix to q
    print("using user provided q-matrix")
    q <- Qmat
    
    if(pars[1] %in% colnames(q) == F) stop ("root provided does not fall within range of q-matrix")
    
    #Get root state
    root <- pars[1] - as.numeric(colnames(q)[1]) + 1
    
  }

  diag(q) <- 0
  diag(q) <- -rowSums(q)
  
  # simulate the chromosome numbers
  print("performing simulation")
  dsims <- sim.character(tree.test, pars=q, x0=root, model="mkn")
  attr(dsims, "node.state") <- NULL
  # save the names for various uses below
  tips <- names(dsims)
  # in the case of the 2010 model we have column names = to the
  # chromosome numbers so we can just use them
  if(model == "2010" || is.null(model)==T){
    dsims[] <- as.numeric(colnames(q)[dsims])
    #Check if parameter matrix is returned
    if(verbose == T){
      result.list <- list(dsims, parMat)
      names(result.list) <- c("chrom.num", "parameter.matrix")
      return(result.list)
    } else {
      return(dsims)
    }
  } 
  
  # under the ChromPlus model things are bit more complex and have
  # to be converted back to chromosome number and binary state
  if(model == "ChromPlus"){
    # for chromRate we need to return two vectors
    # 1) binary state 
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1
    
    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]
    
    #Check if parameter matrix is returned
    if(verbose == T){
      result <- list(b.state,dsims,parMat)
      names(result) <- c("binary.state", "chrom.num","parameter.matrix")
      return(result)
    } else {
      result <- list(b.state,dsims)
      names(result)<- c("binary.state","chrom.num")
      return(result)
    }
  } 
  
  # under the PloidEvol model things are bit more complex and have
  # to be converted back to chromosome number and ploidy state
  if(model == "PloidEvol"){
    # for chromRate we need to return two vectors
    # 1) binary state 
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1
    
    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]
    
    #Check if parameter matrix is returned
    if(verbose == T){
      result <- list(b.state,dsims,parMat)
      names(result) <- c("ploidy.state", "chrom.num","parameter.matrix")
      return(result)
    } else {
      result <- list(b.state,dsims)
      names(result)<- c("ploidy.state","chrom.num")
      return(result)
    }
  } 
  
  # under the SAF model things are bit more complex and have
  # to be converted back to chromosome number and fusion state
  if(model == "SAF"){
    # for SAF we need to return two vectors
    # 1) binary state 
    b.state <- rep(0, length(dsims))
    names(b.state) <- tips
    # if we are in binary state 1 add that in
    b.state[dsims > ncol(q) / 2] <- 1
    
    # 2) chromosome number
    dsims[] <- c(chroms, chroms)[dsims]
    
    
    if(verbose == T){
      result <- list(b.state,dsims,parMat)
      names(result) <- c("fusion.state", "chrom.num","parameter.matrix")
      return(result)
    } else {
      result <- list(b.state,dsims)
      names(result)<- c("fusion.state","chrom.num")
      return(result)
    }
  } 
}