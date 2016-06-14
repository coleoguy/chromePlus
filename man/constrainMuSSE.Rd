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
constrainMuSSE <- function(data, lik, hidden, s.lambda, s.mu, polyploidy, verbose, constrain)
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
  \item{verbose}{
If \code{TRUE} then a list is returned with the first element containing the likelihood function and the second containing the rate identity matrix.
}
  \item{constrain}{
a list capable of carrying additional model constraints these include: 

drop.poly if T polyoploidy is dropped from model

drop.demi if T demiploidy is dropped from model

singlerate if T model is constrained so gain1=loss1 and gain2=loss2

nometa if T model is constrained so gain1=gain2 and loss1=loss2

meta model specificatin for binary character options are ARD and SYM the same as in APE

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