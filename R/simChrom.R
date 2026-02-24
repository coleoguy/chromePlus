#' Simulate Chromosome Number Evolution on a Phylogeny
#'
#' Simulates chromosome number evolution along a fixed phylogenetic tree using
#' one of four built-in models or a user-provided Q-matrix. Internally uses
#' `diversitree::sim.character()` with an `mkn` model.
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param pars Numeric vector of model parameters. Length depends on the model:
#'   \describe{
#'     \item{`"2010"` (length 5)}{gain, loss, demiploidy, polyploidy, root
#'       chromosome number}
#'     \item{`"ChromPlus"` or `"SAF"` (length 12)}{gain1, gain2, loss1, loss2,
#'       demiploidy1, demiploidy2, polyploidy1, polyploidy2, transition 1->2,
#'       transition 2->1, root chromosome number, root binary state (0 or 1)}
#'     \item{`"PloidEvol"` (length 11)}{gain-diploid, gain-polyploid,
#'       loss-diploid, loss-polyploid, demiploidy-diploid, demiploidy-polyploid,
#'       polyploidy-diploid, polyploidy-polyploid, rediploidization, root
#'       chromosome number, root ploidy state (0=diploid, 1=polyploid)}
#'     \item{`NULL` with Qmat (length 1)}{root chromosome number}
#'   }
#' @param limits Numeric vector of length 2 with the minimum and maximum
#'   chromosome numbers, e.g. `c(3, 20)`. Required for built-in models.
#'   Ignored when using a custom Q-matrix.
#' @param model Character string specifying the model: `"2010"`, `"ChromPlus"`,
#'   `"PloidEvol"`, or `"SAF"`. Defaults to `NULL` (user-provided Q-matrix).
#' @param Qmat A user-provided square Q-matrix describing chromosome evolution
#'   rates. Row and column names should be chromosome state labels. Diagonals
#'   are automatically set. Only used when `model = NULL`.
#' @param verbose Logical. If `TRUE`, the parameter matrix used to build the
#'   model is also returned. Defaults to `FALSE`.
#'
#' @return Depends on the model:
#'   \describe{
#'     \item{`"2010"` or `NULL`}{A named numeric vector of chromosome numbers
#'       at tree tips. If `verbose = TRUE`, a list with `chrom.num` and
#'       `parameter.matrix`.}
#'     \item{`"ChromPlus"`}{A list with `binary.state` (0/1 vector) and
#'       `chrom.num`. If `verbose = TRUE`, also includes `parameter.matrix`.}
#'     \item{`"PloidEvol"`}{A list with `ploidy.state` (0/1 vector) and
#'       `chrom.num`. If `verbose = TRUE`, also includes `parameter.matrix`.}
#'     \item{`"SAF"`}{A list with `fusion.state` (0/1 vector) and
#'       `chrom.num`. If `verbose = TRUE`, also includes `parameter.matrix`.}
#'   }
#'
#' @seealso [makeSSEchrom()] for simulating both tree and chromosome data
#'   under an SSE model, [constrainMkn()] and [constrainMuSSE()] for fitting
#'   models.
#'
#' @references
#' Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
#' drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
#' 511--523.
#'
#' @examples
#' \donttest{
#' tree <- ape::rcoal(20)
#' # Simple 2010 model
#' result <- simChrom(tree = tree,
#'                    pars = c(0.1, 0.1, 0.0, 0.01, 5),
#'                    limits = c(3, 12), model = "2010")
#' head(result)
#'
#' # ChromPlus model with binary trait
#' result2 <- simChrom(tree = tree,
#'   pars = c(0.1, 0.1, 0.1, 0.1, 0, 0, 0.01, 0.01, 0.05, 0.05, 5, 0),
#'   limits = c(3, 12), model = "ChromPlus")
#' table(result2$binary.state)
#' }
#'
#' @export
simChrom <- function(tree, pars, limits = NULL, model = NULL, Qmat = NULL, verbose = F){

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
    } else if(pars[12] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral binary state not recognied")
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
    } else if(pars[11] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[10])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral ploidy state not recognied")
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
    } else if(pars[12] == 1){
      x <- length(limits[1]:limits[2]) + which(limits[1]:limits[2] == pars[11])
      root <- as.numeric(colnames(q)[x])
    } else {
      stop("ancestral fusion state not recognied")
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
  dsims <- sim.character(tree, pars=q, x0=root, model="mkn")
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