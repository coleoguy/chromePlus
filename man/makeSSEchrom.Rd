\name{makeSSEchrom}
\alias{makeSSEchrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulates chromosome number data under a binary SSE model
}
\description{
Provides a simulation of chromosome number evolution based on a 
model where an evolving binary trait effects can effect the rate 
of speciation, extinction, chromosome gain, and chromosome loss.
}
\usage{
makeSSEchrom(h, lambda1, lambda2, mu1, mu2, 
             asc1, asc2, desc1, desc2, trans1, trans2, 
             max.taxa, max taxa, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
the maximum haploid number for the simulation
}
  \item{lambda1}{
speciation in binary state 1
}
  \item{lambda2}{
speciation in binary state 2
}
  \item{mu1}{
extinction in state 1
}
  \item{mu2}{
extinction in state 2
}
  \item{asc1}{
chromosome gain in state 1
}
  \item{asc2}{
chromosome gain in state 2
}
  \item{desc1}{
chromosome loss in state 1
}
  \item{desc2}{
chromosome loss in state 2
}
  \item{trans1}{
transition from state 1 to state 2
}
  \item{trans2}{
transition from state 2 to state 1
}
  \item{max.taxa}{
simulation will run till this many taxa are present
}
  \item{x0}{
starting haploid number; the binary state is assigned randomly
}
}
\details{
When using this function you should insure that h is set high enough that it does not act as a hard barrier on chromosme number.
}
\value{
This function returns a list with three elements.  The first is the resulting phylogeny of class phylo. The second is a named vector of the binary states of the tips and the third is a named vector of haploid chromosome numbers.  

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