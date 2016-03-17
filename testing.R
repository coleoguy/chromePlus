library(devtools)
install_github('coleoguy/chromevolR')
library(chromevolR)
library(geiger)
library(diversitree)
# Make a birth-death tree
set.seed(2)
tree <- trees(pars = c(.5, .1), 
              type = "bd", n = 1, 
              max.taxa = 100, 
              include.extinct = FALSE)[[1]]
# Simplify tree
tree <- tree[c(1, 3, 5, 2)]
class(tree) <- "phylo"
# Lets rescale this tree to unit length
tree <- geiger::rescale(tree, model = "depth", 1)
plot(tree)



# Evolve chromosome dataset using our simChrom function
set.seed(2)
data <- simChrom(tree, pars=c(2, 1, 0, 10), limits=2:100)
hist(data, breaks=range(data)[2]-range(data)[1]) 



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

#For optimization issues the tree branches were multiplied by 0.229043
#To preserve the original time unit the model parameters should be multiplied by the same factor !!!
#Final Model parameters
# LOSS_CONST	4.57825
# GAIN_CONST	7.70226
# DUPL	1.67884e-10
# LogLikelihood = -120.253
# Likelihood = 5.95139e-53
# AIC (Akaike information criterion) = 246.507





# I built a new version that should work with MK lets test it:
# Now lets repeat this with the diversitree version
# First we get our chromosome numbers into a probability matrix
d.data <- cbind(names(data), data)
d.data <- datatoMatrix(x=d.data, excess=.25, hidden=F)
d.data <- d.data[, 5:18]
# Now we make the full likelihood function
lik <- make.mkn(tree,states=d.data,k=14,strict=F)
# Make sure it is a working function
lik(rep(.5, 182))
# Constrain to chromevol
lik.con <- constrainMkn(d.data, lik)
# still works
lik.con(rep(.5,3))
# find MLE
result <- find.mle(lik.con,x.init=runif(min = 0, max = 1, 3))

# LOSS_CONST	1.072065
# GAIN_CONST	1.702783
# DUPL	5.792263e-08 = 0ish




# Now lets repeat this with the diversitree version
# First we get our chromosome numbers into a probability matrix
d.data <- cbind(names(data), data)
d.data <- datatoMatrix(x=d.data, excess=.25, hidden=F)
d.data <- d.data[, 5:18]
# Now we make the full likelihood function
lik <- make.musse(tree,states=d.data,k=16,strict=F)
# Make sure it is a working function
lik(rep(.5, 272))
# Constrain to chromevol
lik.con <- constrainMuSSE(d.data, lik, hidden=F, s.lambda=T, s.mu=T, polyploidy=T)
# still works
lik.con(rep(.5,67))
# find MLE
result <- find.mle(lik.con,x.init=runif(min = 0, max = 1, 5))
# LOSS_CONST	1.072065
# GAIN_CONST	1.702783
# DUPL	5.792263e-08 = 0ish

# so these results are within 1% of true values
