\name{simChrom}
\alias{simChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulates chromosome number data
}
\description{
Provides a simulation of chromosome number evolution based on tree and one of 3 models.
}
\usage{
simChrom(tree, pars, limits, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
a phylogenetic tree of the class phylo
}
  \item{pars}{
numeric vector paramters the length of this vector depends on the model specified see details below
}
  \item{limits}{
numeric vector of length 2 with upper and lower limit for chromosome number
}
  \item{model}{
A text string indicating the model for the simulation.  Options are "2010", "chromRate", or "ploidEvol"
}
}
\details{
The model chosen will determine the correct form of the numeric vector supplied to the pars argument.  

If model="2010" the pars vector should be length 5: rates of chromosome gain, chromosome loss, demipolyploidy, polyploidy, and the final value is the chromosome number at the root of the tree.  

If model="chromRate" then the numeric vector supplied to pars should be length 11: chromosome gain in hyperstate 1 and hyperstate 2, chromosome loss in hyperstate 1 and hyperstate 2, demiploidy in hyperstate 1 and hyperstate 2, polyploidy in hyperstate 1 and hyperstate 2, transition from hyperstate 1 to 2, transition from hyperstate 2 to 1, and the chromosome number at the root of the tree.  

If model="ploidEvol" then the numeric vector supplied to pars should be length 10: chromosome gain in diploids and polyploids, chromosome loss in diploids and polyploids, demiploidy in diploids and polyploids, polyploidy in diploids and polyploids, rediploidization, and the chromosome number at the root of the tree. 
}
\value{
a named numeric vector representing chromosome numbers at the tips of the tree is returned
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Heath Blackmon}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}