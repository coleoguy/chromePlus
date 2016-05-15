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
constrainMkn(data, lik, hyper, polyploidy, verbose)
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