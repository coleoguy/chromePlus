library(devtools)
install_github('coleoguy/chromevolR')
library(chromevolR)
library(geiger)
library(diversitree)
# Make a birth-death tree
set.seed(1)
tree <- trees(pars = c(.5, .1), 
              type = "bd", n = 1, 
              max.taxa = 200, 
              include.extinct = FALSE)[[1]]
# Simplify tree
tree <- tree[c(1, 3, 5, 2)]
class(tree) <- "phylo"
# Lets rescale this tree to unit length
tree <- geiger::rescale(tree, model = "depth", 1)
plot(tree)

# Evolve chromosome dataset using our simChrom function
set.seed(1)
data <- simChrom(tree, pars=c(2, 1, .05, 40), limits=2:100)
hist(data, breaks=max(data)) 

# now lets analyze this with old chromevol
setwd("~/Desktop/Dropbox/gitrepos/development/chromevol/PloSSE/testing")
# save tree
write.tree(tree, file = "tree.new")
# save simulated data
chrom.data <- vector()
for(i in 1:length(data)){
  chrom.data <- paste(chrom.data, ">", names(data)[i], "\n", data[i], "\n", sep = "")
}
write(chrom.data, file="chrom.txt")

system(command = "./chromEvol params.txt")

# For optimization the tree branches were multiplied by 0.219758
# Final Model parameters
# LOSS_CONST	2.3666 = 0.5201
# GAIN_CONST	1.98355 = 0.4359
# DUPL	1.67884e-10 = 0ish

## THIS RESULT FROM CHROMEVOL IS NOT 
## STOCHASTIC SAME RESULT ON MULTIPLE RUNS (AT LEAST TO REPORTED ACCURACY)

# True paramters 
# LOSS_CONST	.5
# GAIN_CONST	.5
# DUPL	0

# Now lets repeat this with the diversitree version
# First we get our chromosome numbers into a probability matrix
d.data <- cbind(names(data), data)
d.data <- datatoMatrix(x=d.data, excess=.25, hidden=F)[,2:20]

# Now we make the full likelihood function
lik <- make.musse(tree,states=d.data,k=19,strict=F)
# Make sure it is a working function
lik(rep(.5, 380))
# Constrain to chromevol
lik.con <- constrainChrom(d.data, lik, hidden=F, s.lambda=T, s.mu=T, polyploidy=T)
# still works
lik.con(rep(.5,5))
# find MLE
result <- find.mle(lik.con,x.init=runif(min = 0, max = 1, 5))
# LOSS_CONST	0.5393
# GAIN_CONST	0.4577
# DUPL	1.4791e-06 = 0ish

# so these results are within 1% of true values



# I built a new version that should work with MK lets test it:
# Now lets repeat this with the diversitree version
# First we get our chromosome numbers into a probability matrix
d.data <- cbind(names(data), data)
d.data <- datatoMatrix(x=d.data, excess=.25, hidden=F)[,2:20]

# Now we make the full likelihood function
lik <- make.mkn(tree,states=d.data,k=19,strict=F)
# Make sure it is a working function
lik(rep(.5, 342))
# Constrain to chromevol
lik.con <- makeMKChrom(d.data, lik)
# still works
lik.con(rep(.5,3))
# find MLE
result <- find.mle(lik.con,x.init=runif(min = 0, max = 1, 3))

# ascdip  0.4576690021
# descdip 0.5392637317
# polypl  0.0000995091