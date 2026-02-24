#' Density Plot of chromePlus MCMC Results
#'
#' Creates a density plot of parameter estimates from a chromePlus MCMC
#' analysis, with highest posterior density (HPD) intervals displayed as
#' horizontal lines below the density curves.
#'
#' @param data A data frame where each column represents a parameter from an
#'   MCMC run. Column names are used as parameter labels in the legend.
#' @param colors A character vector of hex color codes, one per parameter
#'   (column) in `data`. If `NULL`, colors are generated automatically using
#'   [viridis::viridis()].
#' @param x_title Character string for the x-axis label.
#' @param y_title Character string for the y-axis label. Defaults to
#'   `"density"`.
#' @param main_title Character string for the plot title.
#' @param legend_title Character string for the legend title. Defaults to
#'   `"parameters"`.
#' @param alpha_geom Numeric value between 0 and 1 controlling the
#'   transparency of the density fill. Defaults to `0.75`.
#' @param alpha_line Numeric value between 0 and 1 controlling the
#'   transparency of the HPD interval lines. Defaults to `0.75`.
#'
#' @return A `ggplot` object displaying the density plot.
#'
#' @seealso [constrainMkn()], [constrainMuSSE()] for generating the MCMC
#'   results that can be plotted with this function.
#'
#' @examples
#' # Create example MCMC-like data
#' set.seed(42)
#' mcmc_data <- data.frame(
#'   asc1 = rnorm(500, 0.1, 0.02),
#'   desc1 = rnorm(500, 0.15, 0.03)
#' )
#' plotChromeplus(
#'   data = mcmc_data,
#'   colors = c("#440154", "#21918c"),
#'   x_title = "rate",
#'   main_title = "Parameter estimates"
#' )
#'
#' @export
plotChromeplus <- function(data,
                           colors,
                           x_title,
                           y_title = "density",
                           main_title,
                           legend_title = "parameters",
                           alpha_geom = 0.75,
                           alpha_line = 0.75){

  ### --- perform checks --- ###
  if(is.null(colors)){
    colnum <- ncol(data)
    colors <- viridis(colnum)
  }

  ### --- set up plotting theme --- ###
  ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border=element_blank(),
                                axis.line = element_line(colour="grey30"),
                                axis.title = element_text(colour="grey20"),
                                axis.text = element_text(colour="grey30"),
                                legend.title = element_text(colour="grey20"),
                                legend.text = element_text(colour="grey30"))
  # get the range of values that we will be plotting on our Y axis
  getMax <- function(x){
    max(density(x)$y)
  }
  maxy <- max(apply(data, MARGIN=2,FUN=getMax))


    ### --- create long data for ggplot --- ###
    long_data <- data.frame(c(as.numeric(unlist(data))),
                            rep(colnames(data), each = nrow(data)))
    colnames(long_data) <- c("rate", "type")

    #check that the approporiate number of colors is supplied
    if(length(colors) != ncol(data)){
      stop("\n Please check the number of supplied colors for the given number
           of parameters.")
    }

    ### --- calculate HPD intervals --- ###
    myHPD <- function(x){
      HPDinterval(as.mcmc(x))
    }
    HPDs <- apply(data, MARGIN = 2, FUN=myHPD)

    HPDs <- data.frame(as.numeric(unlist(as.data.frame(HPDs))),
                       rep(-1:(-ncol(HPDs)), each=2),
                       rep(colnames(HPDs), each=2))
    colnames(HPDs) <-c("X","Y","rate")
    HPDs$Y <- HPDs$Y * .03*maxy
    ### --- plot --- ###
    pplot <- ggplot(long_data, aes(x=long_data$rate)) +
    geom_density(aes(fill=as.factor(long_data$type),
                     y=after_stat(density)),
                 stat="density",
                 position="identity",
                 alpha = c(alpha_geom)) +
    geom_line(data=HPDs,
              aes(x=HPDs$X, y=HPDs$Y, color=as.factor(HPDs$rate)),
              alpha = c(alpha_line),
              linewidth=1.6,
              lineend="round") +
    geom_vline(xintercept=0,
               linetype="dashed",
               color="grey40") +
    scale_color_manual(values = c(colors))+
    scale_fill_manual(values = c(colors)) +
    guides(fill=guide_legend(title="parameter"),
           color="none") +
    labs(title = c(main_title),
         x = c(x_title),
         y = c(y_title)) +
    ggtheme
    return(pplot)
  }
