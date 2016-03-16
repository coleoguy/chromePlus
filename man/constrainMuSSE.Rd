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
constrainMuSSE(x)
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
If true then a binary metastate is assumed
}
  \item{polyploidy}{
If true then polyploidy is assumed to be the metastate
}
  \item{s.lambda}{
If true a single value of lambda is assumed
}
  \item{s.mu}{
If true a single value of mu is assumed
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