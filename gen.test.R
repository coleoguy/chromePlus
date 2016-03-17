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

###
###
###
d.data <- cbind(names(data), data)
d.data <- datatoMatrix(x=d.data, excess=.25, hidden=F)
d.data <- d.data[, 5:18]

# Now we make the full likelihood functions
lik <- make.mkn(tree, states=d.data,k=14,strict=F)

# Constrain to chromevol
lik.con <- constrainMkn(d.data, lik)

# find MLE
foo <- find.mle(lik.con,x.init=runif(min = 0, max = 3, 3))

  