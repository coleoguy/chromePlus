# Michelle Jonika
# 19 December
# This code is a function that creates a density plot of ChromPlus data 


plotChromeplus <- function(data,
                           colors,
                           x_title,
                           y_title = "density",
                           main_title,
                           legend_title = "parameters",
                           alpha_geom = 0.75,
                           alpha_line = 0.75){
  ### --- define inputs --- ###
  # data: A data frame containing results returned from ChromePlus.
  # colors: A character vector containing hex codes for plotting. The length of
  #the vector is dependent on the numbers of parameters included in the model
  # x_title: The text for the x-axis title.
  #y_title: The text for the y-axis title. Default is "density" as this is 
  #making a density plot.
  # main_title: The text for the main title for the plot.
  # legend_title: The text for the legend title for the plot. Default is 
  #"parameters".
  # alpha_geom: A numeric value indicating the alpha parameter for the geometric
  # shapes within the plot.
  #alpha_line: A numeric value indicating the alpha parameter for the geometric
  #lines for the HPD intervals within the plot.

  
  ### --- perform checks --- ###
  #check the number of colors supplied by the user dependent on the number of
  #parameters
  if(is.null(colors)){
    colnum <- ncol(data)
    colors <- viridis(colnum)
  }
  
  ### --- set up plotting theme --- ###
 #set up ggplot theme for the plot
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
    pplot <- ggplot(long_data, aes(x=rate)) +
    geom_density(aes(fill=as.factor(type),
                     y=after_stat(density)),
                 stat="density", 
                 position="identity", 
                 alpha = c(alpha_geom)) +
    geom_line(data=HPDs, 
              aes(x=X, y=Y, color=as.factor(rate)),
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
  
