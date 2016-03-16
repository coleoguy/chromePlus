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
constrainMkn(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
  \item{lik}{
This is a likelihood function created with make.mkn
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