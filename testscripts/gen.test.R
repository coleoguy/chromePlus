library(devtools)
install_github('coleoguy/chromevolR')
library(chromevolR)
library(geiger)
library(diversitree)


########################
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
########################

###
### Loading data
###
data("tree")
data("chrom")
###
###
###

range <- c(min(chrom)-2,max(chrom)+2)

# convert chromosome number to format for diversitree
d.data <- cbind(names(chrom), chrom)
p.mat <- datatoMatrix(x=d.data, range=range, hyper=F)

# convert chromosome number to format for diversitree with hyperstate
d.data <- cbind(names(chrom), chrom)
hp.mat <- datatoMatrix(x=d.data, range=range, hyper=T)

# Now we make the full mkn likelihood function (w/o hyper state)
lik <- make.mkn(tree, states=p.mat, k=ncol(p.mat), strict=F)
# Constrain to chromevol (w/o hyperstate)
lik.con <- constrainMkn(p.mat, lik, model="single")
# find MLE
foo <- find.mle(lik.con, x.init = startVals(length(argnames(lik.con)), 0, 1))

# Now we make the full mkn likelihood function (w hyper state)
h.lik <- make.mkn(tree, states=hp.mat, k=ncol(hp.mat), strict=F)

####### This illustrates the outstanding problem of 
####### how we deal with all unknowns in the diversitree framework

# For our purposes we can manually fix one
hp.mat[1, 6] <- 0 

# this will allow us to run make.mkn without error
h.lik <- make.mkn(tree, states=hp.mat, k=ncol(hp.mat), strict=F)

# but we can't get a likelihood
h.lik(startVals(380,0,1))

# we could check that nothing odd is going on by getting rid of all uncertainty
hp.mat[,] <- 0
for(i in 1:nrow(hp.mat)){
  hp.mat[i, sample(1:20, 1)] <- 1
}

# now it works
h.lik <- make.mkn(tree, states=hp.mat, k=ncol(hp.mat), strict=F)
h.lik(startVals(380,0,1))

# add one taxa uncertain
hp.mat[1,] <- 0
hp.mat[1, c(1,11)]<-.5

# mkn now fails
h.lik <- make.mkn(tree, states=hp.mat, k=ncol(hp.mat), strict=F)
h.lik(startVals(380,0,1))
