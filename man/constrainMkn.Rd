\name{constrainMkn}
\alias{constrainMkn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constrain a diversitree mkn model
}
\description{
Returns a likelihood function that has been constrained to match an mkn model of chromosome evolution yielding essentially the original chromevol model.
}
\usage{
constrainMkn(data, lik, hyper, polyploidy, verbose, oneway, constrain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
  \item{lik}{
This is a likelihood function created with make.mkn
}
  \item{hyper}{
\code{TRUE} indicates that a hyperstate col 3 of data is presesnt 
}
  \item{polyploidy}{
\code{TRUE} indicates that hyperstate is ploidy 
}
  \item{verbose}{
If \code{TRUE} then a list is returned with the first element containing the likelihood function and the second containing the rate identity matrix.
}
  \item{oneway}{
Default is \code{FALSE} if \code{TRUE} then transition rate from state 2 back to 1 is zero.
}
  \item{constrain}{
Takes a list that can include up to 6 items to create different types of constrained models.  Options are explained below in the details section.
}

}
\details{
If the rate identity matrix is returned rates are indicated by the numbers:
\cr rate1 ascending aneuploidy - state1
\cr rate2 descending aneuploidy - state1
\cr rate3 ascending aneuploidy - state2
\cr rate4 descending aneuploidy - state2
\cr rate5 polyploidization of state1
\cr rate6 polploidization of state2
\cr rate7 rediploidization of a state2
\cr rate8 transitions from 1 to 2 for hyperstate
\cr rate9 transitions from 2 to 1 for hyperstate
\cr rate10 demipolyploidy for state1 - even
\cr rate11 demipolyploidy for state1 - odd
\cr rate12 demipolyploidy for state2 - even
\cr rate13 demipolyploidy for state2 - odd
\cr rate14 sex chromosome - autosome fusions (transitions from hyperstate 1 to 2)
\cr rate15 Ro transitions (transitions from hyperstate 2 to 1)
\cr 
\cr 
\cr 
The argument \code{constrain} can contain any of the items listed below.  The default state is shown.
\cr 
\cr \code{drop.poly=F}
\cr Sets polyploidy rate to zero
\cr 
\cr \code{drop.demi=F}
\cr Sets demiploidy rate to zero
\cr 
\cr \code{symmetric=F}
\cr Sets all chromosome change rates to be symmetric
\cr 
\cr \code{nometa=F}
\cr Drops metastate from analysis
\cr 
\cr \code{meta="ARD"}
\cr Can be set as either \code{"ARD"} or \code{"SYM"} these have the same meanings as \code{ace} in the APE package.
\cr
\cr \code{saf.model=F}
\cr Models transitions between hyperstates as sex chromosome-autosome fusions (1 -> 2) and Ro transitions (2 -> 1). Follows model described in Anderson et al. 2020.
\cr
\cr \code{sym.hyperstates=F}
\cr If \code{TRUE}, sets one ascendeing and descending transition rate for hyperstates one and two.
\cr

}
\value{
constrained likelihood function is returned
}
\author{
Heath Blackmon}
