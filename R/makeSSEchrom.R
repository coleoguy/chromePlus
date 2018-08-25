
makeSSEchrom <- function(h,                 # max haploid number
                         lambda1, lambda2,  # speciation in state 1, 2
                         mu1, mu2,          # extinction in state 1, 2
                         asc1, asc2,        # chromosome gain in state 1, 2
                         desc1, desc2,      # chromosome loss in state 1, 2
                         trans1, trans2,    # transitions in binary char
                         max.taxa,          # max taxa
                         x0                 # starting haploid number randomly
){                                          # picks the binary state to be 
                                            # assigned to
  # internal function to make parameter string
  make.pars <- function(k, lambda1, lambda2, mu1, mu2, 
                        asc1, asc2, desc1, desc2, trans1, trans2){
    m <- matrix(rep(0, k^2), k, k)
    for(r in 1:k){  
      for(c in 1:k){      
        if(r <= (k / 2)){        # upper half
          if(c <= (k / 2)){      # left half
            if(r == (c - 1)) m[r, c] <- "asc1"
            if(r == (c + 1)) m[r, c] <-  "dsc1"
          }
          if(c > (k / 2)){         # right half
            if((r + (k / 2)) == c) m[r, c] <- "trans1"
          }
        }
        if(r > (k / 2)){         # lower half
          if(c <= (k / 2)){      # left half
            if(r == (c + (k / 2)))  m[r, c] <- "trans2"
          }
          if(c > (k / 2)){       # right half
            if(r == (c - 1)) m[r, c] <- "asc2"
            if(r == (c + 1)) m[r, c] <-  "dsc2"
          }
        }
      }
    }
    m[m == "asc1"] <- asc1
    m[m == "asc2"] <- asc2
    m[m == "dsc1"] <- desc1
    m[m == "dsc2"] <- desc2
    m[m == "trans1"] <- trans1
    m[m == "trans2"] <- trans2
    
    diag(m) <- NA
    z <- as.vector(t(m))
    z <- as.numeric(z[!is.na(z)])
    par.vals <- c(rep(lambda1, times = (k / 2)), rep(lambda2, times = (k / 2)),
                  rep(mu1, times = (k / 2)), rep(mu2, times = (k / 2)), z)
    return(as.numeric(par.vals))
  }
  
  # internal function to convert from musse to chromosomes
  convert.musse <- function(musse, k){
    tree <- musse[c(1, 2, 3, 7)]
    class(tree) <- "phylo"
    z <- musse$tip.state
    bin.char <- chrom.char <- musse$tip.state
    bin.char[musse$tip.state <= (k / 2)] <- 0
    bin.char[musse$tip.state > (k / 2)] <- 1
    for(i in 1:length(musse$tip.state)){
      if(musse$tip.state[i] > (k / 2)) chrom.char[i] <- 
          (chrom.char[i] - (k / 2))
    }
    results <- list(tree, bin.char, chrom.char)
    names(results) <- c("tree", "binary.char", "haploid.num")
    return(results)
  }
  k <- h * 2                    # setup real dimensions of matrix
  x0 <- x0 + sample(c(0, h), 1)  # randomly assign binary state
  t.pars <- make.pars(k = k, lambda1 = lambda1, lambda2 = lambda2, 
                      mu1 = mu1, mu2 = mu2, 
                      asc1 = asc1, asc2 = asc2, 
                      desc1 = desc1, desc2 = desc2, 
                      trans1 = trans1, trans2 = trans2)
  results.musse <- tree.musse(pars = t.pars, max.taxa = max.taxa, x0 = x0)
  results.chrom <- convert.musse(results.musse, k = k)
  return(results.chrom)
}
