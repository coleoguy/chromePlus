\name{plotChromeplus}
\alias{plotChromeplus}
\title{
Create a density plot from Chromeplus results
}
\description{
Function makes a density plot of Chromeplus results.  
}
\usage{
plotChromeplus(data, colors, x_title, y_title = "density", main_title, 
legend_title = "parameters", alpha_geom = 0.75, alpha_line = 0.75)}
\arguments{
  \item{data}{A data frame containing results returned from ChromePlus. The first column should be the MCMC iteration number, the last column is the probability of each iteration and the intermediate columns are the parameters included in the model.}
  \item{colors}{A character vector containing hex codes for plotting. The length 
  of the vector should be the same as the number of parameters included in the model.}
  \item{x_title}{The text for the x-axis title.}
  \item{y_title}{The text for the y-axis title. Default is "density".}
  \item{main_title}{The text for the main title for the plot.}
  \item{legend_title}{The text for the legend title. Default is "parameters".}
  \item{alpha_geom}{A numeric value indicating the alpha parameter for the 
  geometric shapes within the plot. Defaults to 0.75}
  \item{alpha_line}{A numeric value indicating the alpha parameter for the 
  geometric lines for the HPD intervals within the plot. Defaults to 0.75}
}
\value{Plots a density map. A vector with highest posterior density values for each parameter is returned invisibly.}
\author{Michelle M. Jonika and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, plotting}
