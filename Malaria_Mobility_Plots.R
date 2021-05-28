## Mobility and Malaria transmission model plotting functions
## Created for the Animal and Parasitism Book Chapter by Hannah Meredith and Amy Wesolowski
## Last updated May 25, 2021

##### Plotting time course for cases and movement for different scenarios, compared to baseline (no movement)

plot.baseline.comparison <- function(model1, model2, model1_name, model2_name, title){
  
  model.compare <- rbind.data.frame(model1, model2)
  models <- unique(model.compare$model)
  alpha = ifelse(model.compare$model == model1_name, 0.5, 1)
  
  infections <- ggplot(subset(model.compare, SEIR %in% c("all.infected_h")), aes(t.months, prop.pop, linetype = model, color = as.factor(subpop), alpha = model))+
    geom_line(size = 1.5)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_linetype(name = "Movement",
                   labels = c(model1_name, model2_name))+
    scale_alpha_discrete(range = c(0.5, 1))+
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'none') +
    labs(x = "Time (months)", y = "Infections (proportion)")+
    lims(y = c(0,0.5))+
    ggtitle(title)
  
  movement <- ggplot(subset(model.compare, SEIR %in% c("all.infected_h")))+
    geom_line(aes(x = t.months, y = leaving.prop.total.i, linetype = model, color = as.factor(subpop), alpha = model), size = 2)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_linetype(name = "Movement",
                   labels = c(model1_name, model2_name))+
    scale_x_date(labels = date_format("%b"))+
    scale_alpha_discrete(range = c(0.5, 1))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'none') +
    labs(x = "Time (months)", y = "Travelers (proportion)")+
    lims(y = c(0,0.06))

  time.courses <- ggarrange(infections, movement, ncol = 1, nrow = 2,
                            common.legend = T,
                            legend = "none")
  
  return(time.courses)
}

plot.baseline.comparison.diff.I <- function(model1, model2, model1_name, model2_name, title){
  
  model.compare <- rbind.data.frame(model1, model2)
  models <- unique(model.compare$model)
  
  infections <- ggplot(subset(model.compare, SEIR %in% c("all.infected_h")), aes(t.months, prop.pop, linetype = model, color = as.factor(subpop)))+
    geom_line(size = 1.5)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_linetype(name = "Movement",
                   labels = c(model1_name, model2_name))+
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'none') +
    labs(x = "Time (months)", y = "Infections (proportion)")+
    lims(y = c(0,0.5))+
    ggtitle(title)
  
  movement <- ggplot(subset(model.compare, hum.or.moz == "h" ))+
    geom_line(aes(x = t.months, y = leaving.prop.tot.i, linetype = as.factor(subpop), color = model), size = 1.5)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_linetype(name = "Movement",
                   labels = c(model1_name, model2_name))+
    scale_x_date(labels = date_format("%b"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 14),
          axis.title = element_text(size = 16, hjust = 0.5, vjust = 0.5),
          legend.position = 'none') +
    labs(x = "Time (months)", y = "Travelers (proportion)")+
  lims(y = c(0,0.06))
  
  time.courses <- ggarrange(infections, movement, ncol = 1, nrow = 2,
                            common.legend = T,
                            legend = "none")
  
  return(time.courses)
}

## Plotting time courses of cases and travelers as well as % change in annual cases
plot.cases.movement.impact <- function(movement.scenarios, base.movement){

  movement.totals <- subset(movement.scenarios, SEIR == "all.infected_h") %>%
    group_by(model.name, subpop) %>%
    summarise(total.infections = sum(count))
  baseline.totals <- subset(base.movement, SEIR == "all.infected_h") %>%
    group_by(subpop) %>%
    summarise(total.infections = sum(count))
  movement.totals <- left_join(movement.totals, baseline.totals, by = "subpop")
  movement.totals <- movement.totals %>%
    mutate(change.in.infection = round((total.infections.x - total.infections.y)/total.infections.y*100,2))
  
  fig.infections <- ggplot(subset(movement.scenarios, SEIR %in% c("all.infected_h")), aes(t.months, prop.pop.x, color = as.factor(subpop)))+
    geom_line(size = 1, linetype = 1)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_x_date(breaks = seq(as.Date("2020-01-01"), 
                              as.Date("2020-12-01"), by = "3 months"), date_labels ="%b")+
    scale_y_continuous(breaks = c(0, 0.25, 0.5), limits = c(0, 0.6))+
    facet_wrap(~model.name, nrow = 1)+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_blank(),
          legend.position = 'none',
          strip.text = element_blank(),
          strip.background = element_blank()) +
    labs(x = "", y = "Infections (proportion)")+ 
    geom_line(aes(t.months, prop.pop.y, color = as.factor(subpop)),size = 1.5, linetype = 1, alpha = 0.3)
  
  fig.movement <- ggplot(subset(movement.scenarios, SEIR %in% c("all.infected_h")), aes(t.months, leaving.prop.total.i.x, color = as.factor(subpop)))+
    geom_line(size = 1, linetype = 1)+
    scale_color_manual(name = "Population", 
                       values=c("black", "blue"))+
    scale_x_date(breaks = seq(as.Date("2020-01-01"), 
                              as.Date("2020-12-01"), by = "3 months"), date_labels ="%b")+
    scale_y_continuous(breaks = c(0, 0.025, 0.05), limits = c(0,0.06))+ 
    facet_wrap(~model.name, nrow = 1)+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          legend.position = 'none',
          strip.text = element_blank(),
          strip.background = element_blank()) +
    labs(x = "Time (months)", y = "Travelers (proportion)")+
    geom_line(aes(t.months, leaving.prop.total.i.y, color = as.factor(subpop)),size = 1.5, linetype = 1, alpha = 0.3)
  
  
  fig.case.impact <- ggplot(movement.totals, aes(change.in.infection, model.name, fill = as.factor(subpop)))+
    geom_bar(stat = "identity")+
    labs(y = "Movement scenario", x = "Annual % change in cases,\n relative to baseline movement", fill = "Population")+
    scale_fill_manual(values = c("black", "blue"))+
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 12),
          axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          legend.position = 'none',
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  plots <- list(fig.infections, fig.movement, fig.case.impact)
  return(plots)

}

