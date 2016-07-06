library(chromevolR)
library(diversitree)
### Loading data
data("tree")
data("chrom")
# convert chromosome number to format for diversitree with hyperstate
d.data <- cbind(names(chrom), chrom)
hp.mat <- datatoMatrix(x=d.data, range=c(8,17), hyper=T)
# Now we make the full mkn likelihood function (w hyper state)
h.lik <- make.mkn(tree, states=hp.mat, k=20, strict=F)
####### This illustrates the outstanding problem of 
####### how we deal with all unknowns in the mkn framework
# For our purposes we can manually fix one
hp.mat[1, 6]<-0
h.lik <- make.mkn(tree, states=hp.mat, k=20, strict=F)
# Constrain to chromevol (w hyperstate)
h.lik.con <- constrainMkn(hp.mat, h.lik, model="hyper")
# figure this out later
foo <- find.mle(h.lik.con, x.init = startVals(6, 0, 1))

  