\name{startVals}
\alias{startVals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provides starting values from various distributions
}
\description{
Returns a vector of starting values for use in functions that require starting values.
}
\usage{
startVals(n, min, max, mean, sd, dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of values to return
}
  \item{min}{
the minimum value
}
  \item{max}{
the maximum value
}
  \item{mean}{
the mean value
}
  \item{sd}{
the standard deviation of the distribution
}
  \item{dist}{
the distribution to be used currently just "norm" or "unif"
}

  \value{
constrained likelihood function is returned
}

  \author{
Heath Blackmon
}

