\name{datatoMatrix}
\alias{datatoMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a typical data table to probability matrix 
}
\description{
Returns a probability matrix where each column represents a discrete state and each row is a species.  Values in the matrix describe the probability of being in any of these states.
}
\usage{
datatoMatrix(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a dataframe where the first column is the species name and the second column is the state observed in each tip
}
  \item{excess}{
how far beyond the observed states do you want to allow in your analysis
}
  \item{hidden}{
if true model will include a binary meta state that will can have distinct rates of chromosome number evolution
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