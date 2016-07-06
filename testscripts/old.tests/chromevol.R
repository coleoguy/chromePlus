library(devtools)
install_github('coleoguy/chromevolR')
library(chromevolR)
library(diversitree)
### Loading data
data("tree")
data("chrom")
# convert chromosome number to format for diversitree
d.data <- cbind(names(chrom), chrom)
p.mat <- datatoMatrix(x=d.data, range=c(8,17), hyper=F)
# Now we make the full mkn likelihood function (w/o hyper state)
lik <- make.mkn(tree, states=p.mat, k=10, strict=F)
# Constrain to chromevol (w/o hyperstate)
lik.con <- constrainMkn(p.mat, lik, model="single")
# find MLE
foo <- find.mle(lik.con, x.init = startVals(3, 0, 1))
