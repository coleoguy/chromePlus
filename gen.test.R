library(devtools)
install_github('coleoguy/chromevolR')
library(chromevolR)
library(geiger)
library(diversitree)


####
#### Data simulation section
####
# Make a birth-death tree
set.seed(2)
tree <- trees(pars = c(.5, .1), 
              type = "bd", n = 1, 
              max.taxa = 100, 
              include.extinct = FALSE)[[1]]
# Simplify tree & rescale this tree to unit length
tree <- tree[c(1, 3, 5, 2)]
class(tree) <- "phylo"
tree <- geiger::rescale(tree, model = "depth", 1)
plot(tree)
# Evolve chromosome dataset using our simChrom function
set.seed(2)
data <- simChrom(tree, pars=c(2, 1, 0, 10), limits=2:100)
hist(data, breaks=range(data)[2]-range(data)[1]) 
###
###
###

###
### Loading data
###
data("tree")
data("chrom")
###
###
###

# convert chromosome number to format for diversitree
d.data <- cbind(names(chrom), chrom)
p.mat <- datatoMatrix(x=d.data, range=c(8,17), hyper=F)

# convert chromosome number to format for diversitree with hyperstate
d.data <- cbind(names(chrom), chrom)
hp.mat <- datatoMatrix(x=d.data, range=c(8,17), hyper=T)

# Now we make the full mkn likelihood function (w/o hyper state)
lik <- make.mkn(tree, states=p.mat, k=10, strict=F)

# Now we make the full mkn likelihood function (w hyper state)
h.lik <- make.mkn(tree, states=hp.mat, k=20, strict=F)
#######
#######
####### This illustrates the outstanding problem of 
####### how we deal with all unknowns in the mkn framework
#######
#######

# For our purposes we can manually fix one
hp.mat[1, 6]<-0
h.lik <- make.mkn(tree, states=hp.mat, k=20, strict=F)

# Constrain to chromevol (w/o hyperstate)
lik.con <- constrainMkn(p.mat, lik, model="single")

# Constrain to chromevol (w hyperstate)
h.lik.con <- constrainMkn(hp.mat, h.lik, model="hyper")

# find MLE
foo <- find.mle(lik.con, x.init = startVals(3, 0, 1))

# figure this out later
foo <- find.mle(h.lik.con, x.init = startVals(6, 0, 1))

  