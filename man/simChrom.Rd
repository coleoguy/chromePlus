\name{simChrom}
\alias{simChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulates chromosome number data
}
\description{
Provides a simulation of chromosome number evolution based on tree and a 5 paramter model.
}
\usage{
simChrom(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
a phylogenetic tree of the class phylo
}
  \item{pars}{
numeric vector of length 4 in the order ascending rate, descending rate, polyploidy, root value
}
  \item{limits}{
numeric vector of length 2 with upper and lower limit for chromosome number
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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