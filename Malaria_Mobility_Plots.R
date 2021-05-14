## Mobility and Malaria transmission model plotting functions
## Created for the Animal and Parasitism Book Chapter by Hannah Meredith and Amy Wesolowski
## Last updated May 13, 2021

##### Plotting time course for cases and movement for different scenarios, compared to baseline (no movement)

plot.baseline.comparison <- function(model1, model2, model1_name, model2_name){
  
  model.compare <- rbind.data.frame(model1, model2)
  models <- unique(model.compare$model)
  
  infections <- ggplot(subset(model.compare, SEIR %in% c("all.infected_h")), aes(t.months, prop.pop, linetype = as.factor(subpop), color = model))+
    geom_line(size = 1.5)+
    scale_color_manual(name = "Model", 
                       limits = models,
                       labels = c(model1_name, model2_name),
                       values=c("black", "blue"))+
    scale_linetype(name = "Sub-population")+
    # facet_wrap(~model)
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'right') +
    labs(x = "Time (months)", y = "Infections (proportion)")+
    lims(y = c(0,0.5))
  
  movement <- ggplot(subset(model.compare, hum.or.moz == "h" ))+
    geom_line(aes(x = t.months, y = prop.leaving.tot, linetype = as.factor(subpop), color = model), size = 1.5)+
    scale_color_manual(name = "Model", 
                       limits = models,
                       labels = c(model1_name, model2_name),
                       values=c("black", "blue"))+
    scale_linetype(name = "Sub-population")+
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'right') +
    labs(x = "Time (months)", y = "Travelers (proportion)")
  # lims(y = c(0,0.15))
  
  time.courses <- ggarrange(infections, movement, ncol = 1, nrow = 2,
                            common.legend = T,
                            legend = "right")
  
  return(time.courses)
}

plot.baseline.comparison.diff.I <- function(model1, model2, model1_name, model2_name){
  
  model.compare <- rbind.data.frame(model1, model2)
  models <- unique(model.compare$model)
  
  infections <- ggplot(subset(model.compare, SEIR %in% c("all.infected_h")), aes(t.months, prop.pop, linetype = as.factor(subpop), color = model))+
    geom_line(size = 1.5)+
    scale_color_manual(name = "Model", 
                       limits = models,
                       labels = c(model1_name, model2_name),
                       values=c("black", "blue"))+
    scale_linetype(name = "Sub-population")+
    # facet_wrap(~model)
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'right') +
    labs(x = "Time (months)", y = "Infections (proportion)")+
    lims(y = c(0,0.5))
  
  movement <- ggplot(subset(model.compare, hum.or.moz == "h" ))+
    geom_line(aes(x = t.months, y = prop.leaving.tot, linetype = as.factor(subpop), color = model), size = 1.5)+
    scale_color_manual(name = "Model", 
                       limits = models,
                       labels = c(model1_name, model2_name),
                       values=c("black", "blue"))+
    scale_linetype(name = "Sub-population")+
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'right') +
    labs(x = "Time (months)", y = "Travelers (proportion)")
  # lims(y = c(0,0.15))
  
  time.courses <- ggarrange(infections, movement, ncol = 1, nrow = 2,
                            common.legend = T,
                            legend = "right")
  
  return(time.courses)
}
