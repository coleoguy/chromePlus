# Michelle Jonika
# 19 December
# This code is a function that creates a density plot of ChromPlus data 

library(ggplot2)
library(coda)

plot.chromeplus <- function(data = NULL,
                            poly = FALSE,
                            colors = NULL,
                            x_title = NULL,
                            y_title = "density",
                            main_title = NULL,
                            legend_title = "parameters",
                            alpha_geom = 0.75,
                            alpha_line = 0.75){
  ### --- define inputs --- ###
  # data: A data frame containing results returned from ChromePlus.
  # poly: A logical vector of length one. TRUE indicates the model includes 
  #polyploidy. Default is FALSE and indicates the model does not include polyploidy
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
  
  #store the number of colors supplied 
  check <- length(colors)
  
  #loop to check that if there is not polyploidy, there should be two colors
  if(poly == FALSE){
    if(check == 2){
      print("Correct number of supplied colors for the given number of parameters.")
    } else{
      stop("\n Please check the number of supplied colors for the given number 
           of parameters.")
    } 
  }
  
  #loop to check that if there is polyploidy, there should be three colors
  if(poly == TRUE){
    if(check == 3){
      print("Correct number of supplied colors for the given number of parameters.")
    } else{
      stop("\n Please check the number of supplied colors for the given number 
           of parameters.")
    } 
  }
  
  
  
  
  #set up ggplot theme for the plot
  ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border=element_blank(),
                                axis.line = element_line(colour="grey30"),
                                axis.title = element_text(colour="grey20"),
                                axis.text = (element_text(colour="grey30")),
                                legend.title = element_text(colour="grey20"),
                                legend.text = element_text(colour="grey30"))
  
  if(poly == FALSE){
    ### --- create long data for ggplot --- ###
    long_data <- data.frame(c(data$asc2 - data$asc1,
                              data$desc2 - data$desc1),
                            rep(c("fission", "fusion"), 
                                each = nrow(data)))
    colnames(long_data) <- c("rate", "type")
    
    ### --- calculate HPD intervals --- ###
    hpd_wopoly <- data.frame(X = c(HPDinterval(as.mcmc(long_data$rate[long_data$type=="fission"]))[1,],
                                   HPDinterval(as.mcmc(long_data$rate[long_data$type=="fusion"]))[1,]),
                             Y = c(-1, -1, -2, -2),
                             types = rep(c("fission", "fusion"), each = 2))
    
    ### --- plot --- ###
    ggplot(long_data, aes(x=rate)) +
    geom_density(aes(fill=as.factor(type),
                     y=..density..),
                 stat="density", 
                 position="identity", 
                 alpha = c(alpha_geom)) +
    geom_line(data=hpd_wopoly, 
              aes(x=X, y=Y, color=as.factor(types)),
              alpha = c(alpha_line), 
              size=1.4, 
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
  }
  
  if(poly == TRUE){
    ### --- create long data for ggplot --- ###
    long_data <- data.frame(c(data$asc2 - data$asc1,
                              data$desc2 - data$desc1,
                              data$pol2 - data$pol1),
                            rep(c("fission", "fusion", "wgd"), 
                                each = nrow(data)))
    colnames(long_data) <- c("rate", "type")
    
    ### --- calculate HPD intervals --- ###
    hpd_poly <- data.frame(X = c(HPDinterval(as.mcmc(long_data$rate[long_data$type=="fission"]))[1,],
                                 HPDinterval(as.mcmc(long_data$rate[long_data$type=="fusion"]))[1,],
                                 HPDinterval(as.mcmc(long_data$rate[long_data$type=="wgd"]))[1,]),
                           Y = c(-2, -2, -4, -4, -6, -6),
                           types = rep(c("fission", "fusion", "wgd"), each = 2))
    
    ### --- plot --- ###
    ggplot(long_data, aes(x=rate)) +
      geom_density(aes(fill=as.factor(type),
                       y=..density..),
                   stat="density", 
                   position="identity", 
                   alpha= c(alpha_geom)) +
      geom_line(data=hpd_poly, 
                aes(x=X, y=Y, color=as.factor(types)),
                alpha= c(alpha_line), 
                size=1.4, 
                lineend="round") +
      geom_vline(xintercept=0, 
                 linetype="dashed", 
                 color="grey40", 
                 size=.5) +
      scale_color_manual(values = c(colors)) +
      scale_fill_manual(values = c(colors)) +
      guides(fill=guide_legend(title="parameter"),
             color="none") +
      labs(title = c(main_title),
           x = c(x_title),
           y = c(y_title)) +
      ggtheme
  }
}