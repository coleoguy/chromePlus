\name{plot.chromeplus}
\alias{plot.chromeplus}
\title{
Create a density plot from Chromeplus results
}
\description{
Function makes a density plot of Chroeplus results.  
}
\usage{
plot.chromeplus(data, poly, colors, x_title, y_title, main_title, 
legend_title, alpha_geom, alpha_line)}
\arguments{
  \item{data}{A data frame containing results returned from ChromePlus. The first column should be the MCMC iteration number, the last column is the probability of each iteration and the }
  \item{poly}{A logical vector of length one. TRUE indicates the model includes 
  polyploidy. Default is FALSE and indicates the model does not include polyploidy}
  \item{colors}{A character vector containing hex codes for plotting. The length 
  of the vector is dependent on the number of parameters included in the model.}
  \item{x_title}{The text for the x-axis title.}
  \item{y_title}{The text for the y-axis title. Default is "density" as this is 
  making a density plot.}
  \item{main_title}{The text for the main title for the plot.}
  \item{legend_title}{he text for the legend title for the plot. Default is 
  "parameters".}
  \item{alpha_geom}{A numeric value indicating the alpha parameter for the 
  geometric shapes within the plot.}
  \item{alpha_line}{A numeric value indicating the alpha parameter for the 
  geometric lines for the HPD intervals within the plot.}
}
\value{Plots a density map. A vector with highest posterior density values for each parameter is returned invisibly.}
\author{Michelle M. Jonika and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, plotting}
