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
constrainMuSSE <- function(data, lik, hidden, s.lambda, s.mu, polyploidy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
  \item{lik}{
This is a likelihood function created with make.MuSSE
}
  \item{hidden}{
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