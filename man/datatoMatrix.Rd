\name{datatoMatrix}
\alias{datatoMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a typical data table to probability matrix 
}
\description{
Returns a probability matrix where each column represents a discrete state and each row is a species.  Values in the matrix describe the probability of being in any of these states. If a hyper state is specified each species is assumed to be present with equal probability in either hyper state.
}
\usage{
datatoMatrix(x, range, hyper = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a dataframe where the first column is the species name, the second column is the number of chromosomes, and the third is the probability in being in state 1 of trait 2
}
  \item{range}{
a vector of length two with the first element equal to the lowest chromosome number to be included in the model and the second element equal to the highest chromosome number to be included in the model.
}
  \item{hyper}{
If `TRUE` model will include a binary hyper state (col 3 of input data) that can have distinct rates of chromosome number evolution
}
}
\value{
a probability matrix is returned
}

\author{
Heath Blackmon}
