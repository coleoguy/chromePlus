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
constrainMkn(data, lik, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is the data matrix created by running \code{dataMatrx}
}
  \item{lik}{
This is a likelihood function created with make.mkn
}
  \item{model}{
This is text string specifying the model being used. Options are: "single" for a single chromosome number transition matrix. "hyper" for a model where the state of a binary trait determines which transition matrix is used.  This allows for each binary state to have idiosyncratic rates of chromosome number evolution. "ploidy" in this model the binary trait is ploidy such that diploids and polyploids can have their own rate of evolution.
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