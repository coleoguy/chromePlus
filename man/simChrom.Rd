\name{simChrom}
\alias{simChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulates chromosome number data
}
\description{
Provides a simulation of chromosome number evolution based on tree and one of 4 models or a user provided q-matrix.  This function uses diversitrees sim.character to perform simulation.
}
\usage{
simChrom(tree, pars, limits, model, Qmat, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
a phylogenetic tree of the class phylo
}
  \item{pars}{
numeric vector of paramters the length of this vector depends on the model specified see details below
}
  \item{limits}{
numeric vector of length 2 with upper and lower limit for chromosome number. Defaults to NULL
}
  \item{model}{
A text string indicating the model for the simulation.  Options are "2010", "ChromPlus", "PloidEvol", or "SAF". Defaults to NULL 
}
  \item{Qmat}{
A user provided q-matrix which describes the desired model of chromosome evolution. Transitions should be filled in with desired rates. Diagonals and rowsums will automatically be set to equal 0. Function currently doesn't support user provided q-matrices with hyperstates. Defaults to NULL
}
  \item{verbose}{
Binary argument giving whether or not to return parameter matrix used to build model. Defaults to FALSE  
}
}
\details{
The model chosen will determine the correct form of the numeric vector supplied to the pars argument.  
\cr \cr 
If model="2010" the pars vector should be length 5: 
\cr chromosome gain 
\cr chromosome loss 
\cr demipolyploidy 
\cr polyploidy 
\cr chromosome number at the root of the tree.  
\cr \cr \cr 
If model="ChromPlus" or "SAF" then the numeric vector supplied to pars should be length 12: 
\cr chromosome gain in hyperstate/fusion state 1  
\cr chromosome gain in hyperstate/fusion state 2 
\cr chromosome loss in hyperstate/fusion state 1 
\cr chromosome loss in hyperstate/fusion state 2 
\cr demiploidy in hyperstate/fusion state 1 
\cr demiploidy in hyperstate/fusion state 2 
\cr polyploidy in hyperstate/fusion state 1 
\cr polyploidy in hyperstate/fusion state 2 
\cr transition from hyperstate/fusion state 1 to 2 
\cr transition from hyperstate/fusion state 2 to 1 
\cr chromosome number at the root of the tree
\cr hyperstate/fusion state at the root of the tree (0=hyperstate 1/unfused, 1 =hyperstate 2/fused)
\cr \cr \cr 
If model="PloidEvol" then the numeric vector supplied to pars should be length 11: 
\cr chromosome gain in diploids 
\cr chromosome gain in polyploids
\cr chromosome loss in diploids
\cr chromosome loss in polyploids
\cr demiploidy in diploids
\cr demiploidy in polyploids
\cr polyploidy in diploids
\cr polyploidy in polyploids
\cr rediploidization
\cr chromosome number at the root of the tree. 
\cr ploidy level at the root of the tree (0=diploid, 1=polyploid)
\cr \cr \cr
If model=NULL then the numeric vector supplied to pars should be length 1:
\cr chromosome number at the root of the tree.
}
\value{
a named numeric vector representing chromosome numbers at the tips of the tree is returned for model="2010" and model=NULL. For model="ChromPlus", model="PloidEvol", and model="SAF" a list of length 2 is returned containing a vector of chromosome numbers for tips of the tree as well as a vector containing the hyperstate or ploidy for each tip. If verbose=T, the parameter matrix used to simulate the data is also returned.
}

\author{
Heath Blackmon}

