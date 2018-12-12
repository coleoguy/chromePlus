\name{constrainMuSSE}
\alias{constrainMuSSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constrain a diversitree MuSSE model
}
\description{
Returns a likelihood function that has been constrained to match an MuSSE model of chromosome evolution yielding essentially the original chromevol model.
}
\usage{
constrainMuSSE <- function(data, lik, hyper, s.lambda, s.mu, polyploidy, verbose, constrain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
  \item{lik}{
This is a likelihood function created with make.MuSSE
}
  \item{hyper}{
If \code{TRUE} then a binary metastate is assumed
}
  \item{s.lambda}{
If \code{TRUE} a single value of lambda is assumed across metastates
}
  \item{s.mu}{
If \code{TRUE} a single value of mu is assumed across metastates
}
  \item{polyploidy}{
If \code{TRUE} then polyploidy is assumed to be the metastate
}
  \item{verbose}{
If \code{TRUE} then a list is returned with the first element containing the likelihood function and the second containing the rate identity matrix.
}
  \item{constrain}{
a list capable of carrying additional model constraints these include see details for specifications.

}
}
\details{
If the rate identity matrix is returned rates are indicated by the numbers:

rate1 ascending aneuploidy - diploid or state 1 of hypertrait

rate2 descending aneuploidy - diploid or state 1 of hypertrait

rate3 ascending aneuploidy - polyploid or state 2 of hypertrait

rate4 descending aneuploidy - polyploid or state 2 of hypertrait

rate5 polyploidization of a diploid or state 1 of hypertrait

rate6 polploidization of a polyploid or state 2 of hypertrait

rate7 rediploidization of a polyploid

rate8 transitions from state 1 to 2 of hypertrait

rate9 transitions from state 2 to 1 of hypertrait

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



}
\value{
constrained likelihood function is returned
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