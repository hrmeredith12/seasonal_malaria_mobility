## Mobility and Malaria transmission model functions
## Created for the Animal and Parasitism Book Chapter by Hannah Meredith and Amy Wesolowski
## Last updated May 26, 2021


make.constant.travel.rate.matrix <- function(travel.dates, prob.travel.base, prob.travel.holiday){
  
  travel.days <- lubridate::yday(travel.days) + P[["t_ss"]]
  days <- seq(0, P[["t_ss"]] + P[["treatment_pd"]], 1)
  
  origins = destinations = seq(1,subpop, 1)
  OD <- expand.grid(origins, destinations)
  colnames(OD) <- c("origin", "destination")
  
  sim.days <- as.data.frame(matrix(prob.travel.base, nrow = nrow(OD), ncol = length(days)))
  colnames(sim.days) <- days
  sim.days[colnames(sim.days) %in% travel.days] <- prob.travel.holiday
  
  travel.rate.matrix <- cbind(OD, sim.days)
  
  return(travel.rate.matrix)
}
make.specific.date.travel.rate.matrix <- function(travel.out.days, travel.back.days, prob.travel.base, prob.travel.holiday){
  
  travel.out.days <- lubridate::yday(travel.out.days) + P[["t_ss"]]
  travel.back.days <- lubridate::yday(travel.back.days) + P[["t_ss"]]
  days <- seq(0, P[["t_ss"]] + P[["treatment_pd"]], 1)
  
  origins = destinations = seq(1,subpop, 1)
  OD <- expand.grid(origins, destinations)
  colnames(OD) <- c("origin", "destination")
  
  sim.days <- as.data.frame(matrix(prob.travel.base, nrow = nrow(OD), ncol = length(days)))
  colnames(sim.days) <- days  
  
  sim.days[colnames(sim.days) %in% travel.out.days] <- prob.travel.holiday[,1]
  sim.days[colnames(sim.days) %in% travel.back.days] <- prob.travel.holiday[,2]
  
  travel.rate.matrix <- cbind(OD, sim.days)
  
  return(travel.rate.matrix)
  
}

## mal.***.mob.models run the simulation "biting". 
## Inputs are vector of initial variables y0 and vector of parameters P
## Outputs are dataframe of human and mosquitoes in each subpopulation at different malaria stages (SEIRA/SEI) over time and a dataframe of the total numbers for each subpopulation

mal.constant.mob.model <- function(y0, P, model) {
  
  step = 1                                             # the step size of the ODE solver (set to 1 day)
  t = seq(0, P[["t_ss"]] + P[["treatment_pd"]], by = step) #duration of simulation = t_ss (time to reach steady state) + treatment period
  
  model.output <-
    rk(
      y = y0,
      times = t,
      func = biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  model.output.cleaned <- clean.model.output(model.output, model, P)
  return(model.output.cleaned)
}
mal.dynamic.mob.model <- function(y0, P, model){
  step = 1                                             # the step size of the ODE solver (set to 1 day)
  t = seq(0, P[["t_ss"]] + P[["treatment_pd"]], by = step) #duration of simulation = t_ss (time to reach steady state) + treatment period
  
  model.output <-
    rk(
      y = y0,
      times = t,
      func = biting.travel.matrix,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  model.output.cleaned <- clean.model.output.OD.matrix(model.output, model, P)
  return(model.output.cleaned)
}
mal.constant.SEIRA.mob.model <- function(y0, P, model){
  step = 1                                             # the step size of the ODE solver (set to 1 day)
  t = seq(0, P[["t_ss"]] + P[["treatment_pd"]], by = step) #duration of simulation = t_ss (time to reach steady state) + treatment period
  
  model.output <-
    rk(
      y = y0,
      times = t,
      func = biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  model.output.cleaned <- clean.model.output.diff.I(model.output, model, P)
  return(model.output.cleaned)
}
mal.threshold.mob.model <- function(y0, P, model){
  step = 1                                             # the step size of the ODE solver (set to 1 day)
  t = seq(0, P[["t_ss"]] + P[["treatment_pd"]], by = step) #duration of simulation = t_ss (time to reach steady state) + treatment period
  
  model.output <-
    rk(
      y = y0,
      times = t,
      func = biting.threshold.mob,
      parms = P,
      method = "ode45",
      atol = 1e-7,#1e-10,
      rtol =  1e-7#1e-10
    )
  
  model.output.cleaned <- clean.model.output.threshold(model.output, model, P)
  return(model.output.cleaned)
}



clean.model.output <- function(model.output, model, P){
  
  model.output <- as.data.frame(model.output)
  
  # Remove period when system is reaching steady state and wrangle dates 
  model.output$time.idx <- (model.output$time - P[["t_ss"]]) 
  model.output <- model.output[model.output$time.idx > 0, ][-1]
  
  # Determine the proportion of each subpopulation at a given time and also form a "Total infected" category
  model.output_long <- model.output %>%
    pivot_longer(cols = c(1:ncol(model.output)-1),
                 names_to = "SEIR",
                 values_to = "count")%>%
    separate(subpop, into = c('SEIR', 'subpop'), sep = -1, convert = TRUE) %>%
    mutate(hum.or.moz = sub(".*_", "", SEIR))
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz)%>%
    mutate(total.pop = sum(count),
           prop.pop = count/total.pop)
  
  model.output_long$t.months <- as.Date(as.numeric(model.output_long$time.idx), origin ="2020-01-01")
  
  # determine count and proportion of population that traveled at the time point
  model.output_long$leaving.SEIRA.i <- ifelse(model.output_long$hum.or.moz == "h", rowSums(P[["travelers"]])/model.output_long$total.pop * model.output_long$count, 0)
  model.output_long$leaving.prop.total.i <- ifelse(model.output_long$hum.or.moz == "h", rowSums(P[["travelers"]])/model.output_long$total.pop, 0)
  model.output_long$arriving.SEIRA.i <- ifelse(model.output_long$hum.or.moz == "h", colSums(P[["travelers"]])/model.output_long$total.pop * model.output_long$count, 0)
  model.output_long$arriving.prop.total.i <- ifelse(model.output_long$hum.or.moz == "h", colSums(P[["travelers"]])/model.output_long$total.pop, 0)
  
  # create category for all infectious people (I and A)
  infected <- subset(model.output_long, SEIR %in% c("I_h", "A_h"))
  infected.1 <- infected %>%
    group_by(time.idx, subpop)%>%
    summarise(SEIR = "all.infected_h",
              count = sum(count),
              hum.or.moz = "h", 
              total.pop = total.pop,
              prop.pop = count/total.pop,
              t.months = t.months,
              leaving.SEIRA.i = sum(leaving.SEIRA.i),
              leaving.prop.total.i = leaving.prop.total.i, 
              arriving.SEIRA.i = sum(arriving.SEIRA.i),
              arriving.prop.total.i = arriving.prop.total.i
              )%>%
    distinct(time.idx, subpop, .keep_all = T)
  
  model.output_long <- rbind.data.frame(model.output_long, infected.1)
  model.output_long$model <- rep(model, nrow(model.output_long))
  
  return(model.output_long)  
}
clean.model.output.OD.matrix <- function(model.output, model, P){
  
  model.output <- as.data.frame(model.output)
  
  # Remove period when system is reaching steady state and wrangle dates 
  model.output$time.idx <- (model.output$time - P[["t_ss"]]) 
  model.output <- model.output[model.output$time.idx > 0, ][-1]
  
  travelers <-P[["travelers"]]
  travelers <- travelers[ ,!colnames(travelers) %in% seq(0, P[["t_ss"]],1)] 
  colnames(travelers) <- c("origin", "destination", seq(as.numeric(colnames(travelers[3]))-P[["t_ss"]],
                                                        as.numeric(colnames(travelers[ncol(travelers)]))-P[["t_ss"]],
                                                        1)) 
  # Calculate total inbound and outbound travelers for each origin
  travelers.long <- travelers %>%
    pivot_longer((cols = c(3:ncol(travelers))),
                 names_to = "time.idx",
                 values_to = "travelers")
  travelers.long$time.idx <- as.numeric(travelers.long$time.idx)
  
  travelers.out <- travelers.long %>%
    group_by(time.idx, origin)%>%
    summarise(outbound = sum(travelers))
  
  travelers.in <- travelers.long %>%
    group_by(time.idx, destination)%>%
    summarise(inbound = sum(travelers))
  
  travelers.in.out <- left_join(travelers.in, travelers.out , by = c("time.idx", "destination" = "origin"))
  
  # Determine the proportion of each subpopulation at a given time and also form a "Total infected" category
  model.output_long <- model.output %>%
    pivot_longer(cols = c(1:ncol(model.output)-1),
                 names_to = "SEIR",
                 values_to = "count")%>%
    separate(subpop, into = c('SEIR', 'subpop'), sep = -1, convert = TRUE) %>%
    mutate(hum.or.moz = sub(".*_", "", SEIR))
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz)%>%
    mutate(total.pop = sum(count),
           prop.pop = count/total.pop)
  
  model.output_long <- left_join(model.output_long, travelers.in.out, by = c("time.idx", "subpop" = "destination"))
  
  model.output_long$t.months <- as.Date(as.numeric(model.output_long$time.idx), origin ="2020-01-01")
  
  # determine count and proportion of population that traveled at the time point
  model.output_long$leaving.SEIRA.i <- ifelse(model.output_long$hum.or.moz == "h", model.output_long$outbound/model.output_long$total.pop * model.output_long$count, 0) # proportion of specific category that was outbound
  model.output_long$leaving.prop.total.i <- ifelse(model.output_long$hum.or.moz == "h", model.output_long$outbound/model.output_long$total.pop, 0) #proportion of population that was outbound
  model.output_long$arriving.SEIRA.i <- ifelse(model.output_long$hum.or.moz == "h", model.output_long$inbound/model.output_long$total.pop * model.output_long$count, 0) # proportion of specific category that was inbound
  model.output_long$arriving.prop.total.i <- ifelse(model.output_long$hum.or.moz == "h", model.output_long$inbound/model.output_long$total.pop, 0) # proportion of population that was inbound
  
  # create category for all infectious people (I and A)
  infected <- subset(model.output_long, SEIR %in% c("I_h", "A_h"))
  infected.1 <- infected %>%
    group_by(time.idx, subpop)%>%
    summarise(SEIR = "all.infected_h",
              count = sum(count),
              hum.or.moz = "h", 
              total.pop = total.pop,
              prop.pop = count/total.pop,
              inbound = NA,
              outbound = NA,
              t.months = t.months,
              leaving.SEIRA.i = sum(leaving.SEIRA.i),
              leaving.prop.total.i = leaving.prop.total.i, 
              arriving.SEIRA.i = sum(arriving.SEIRA.i),
              arriving.prop.total.i = arriving.prop.total.i
              )%>%
    distinct(time.idx, subpop, .keep_all = T)
  
  model.output_long <- rbind.data.frame(model.output_long, infected.1)
  model.output_long$model <- rep(model, nrow(model.output_long))
  
  return(model.output_long)  
}
clean.model.output.diff.I <- function(model.output, model, P){
  
  model.output <- as.data.frame(model.output)
  
  # Remove period when system is reaching steady state and wrangle dates 
  model.output$time.idx <- (model.output$time - P[["t_ss"]]) 
  model.output <- model.output[model.output$time.idx > 0, ][-1]
  
  # Determine the proportion of each subpopulation at a given time and also form a "Total infected" category
  model.output_long <- model.output %>%
    pivot_longer(cols = c(1:ncol(model.output)-1),
                 names_to = "SEIR",
                 values_to = "count")%>%
    separate(subpop, into = c('SEIR', 'subpop'), sep = -1, convert = TRUE) %>%
    mutate(hum.or.moz = sub(".*_", "", SEIR))
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz)%>%
    mutate(total.pop = sum(count),
           prop.pop = count/total.pop)
  
  model.output_long$t.months <- as.Date(as.numeric(model.output_long$time.idx), origin ="2020-01-01")
  
  # determine count and proportion of population that traveled at the time point
  
  model.output_long$leaving.SEIRA.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h"), rowSums(P[["travelers"]])/model.output_long$total.pop * model.output_long$count,
                                              ifelse(model.output_long$SEIR %in% c("I_h"), rowSums(P[['adjust.I']] * P[["travelers"]])/model.output_long$total.pop * model.output_long$count, 0))
  
  model.output_long$leaving.prop.total.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h"), rowSums(P[["travelers"]])/model.output_long$total.pop,
                                                   ifelse(model.output_long$SEIR %in% c("I_h"), rowSums(P[['adjust.I']] * P[["travelers"]])/model.output_long$total.pop, 0))
  
  model.output_long$arriving.SEIRA.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h"), colSums(P[["travelers"]])/model.output_long$total.pop * model.output_long$count,
                                               ifelse(model.output_long$SEIR %in% c("I_h"), colSums(P[['adjust.I']] * P[["travelers"]])/model.output_long$total.pop * model.output_long$count, 0))
  
  model.output_long$arriving.prop.total.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h"), colSums(P[["travelers"]])/model.output_long$total.pop,
                                                    ifelse(model.output_long$SEIR %in% c("I_h"), colSums(P[['adjust.I']] * P[["travelers"]])/model.output_long$total.pop, 0))
  
  
  # create category for all infectious people (I and A)
  infected <- subset(model.output_long, SEIR %in% c("I_h", "A_h"))
  infected.1 <- infected %>%
    group_by(time.idx, subpop)%>%
    summarise(SEIR = "all.infected_h",
              count = sum(count),
              hum.or.moz = "h", 
              total.pop = total.pop,
              prop.pop = count/total.pop,
              t.months = t.months,
              leaving.SEIRA.i = sum(leaving.SEIRA.i),
              leaving.prop.total.i = leaving.prop.total.i, 
              arriving.SEIRA.i = sum(arriving.SEIRA.i),
              arriving.prop.total.i = arriving.prop.total.i
    )%>%
    distinct(time.idx, subpop, .keep_all = T)
  
  model.output_long <- rbind.data.frame(model.output_long, infected.1)
  model.output_long$model <- rep(model, nrow(model.output_long))
  
  return(model.output_long)  
}
clean.model.output.threshold <- function(model.output, model, P){
  
  model.output.orig <-model.output
  model.output <- as.data.frame(model.output)
  
  # Remove period when system is reaching steady state and wrangle dates
  model.output$time.idx <- (model.output$time - P[["t_ss"]])
  model.output <- model.output[model.output$time.idx > 0, ][-1]
  
  # Determine the proportion of each subpopulation at a given time and also form a "Total infected" category
  model.output_long <- model.output %>%
    pivot_longer(cols = c(1:ncol(model.output)-1),
                 names_to = "SEIR",
                 values_to = "count")%>%
    separate(subpop, into = c('SEIR', 'subpop'), sep = -1, convert = TRUE) %>%
    mutate(hum.or.moz = sub(".*_", "", SEIR))
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz)%>%
    mutate(total.pop = sum(count),
           prop.pop = count/total.pop)
  
  model.output_long$t.months <- as.Date(as.numeric(model.output_long$time.idx), origin ="2020-01-01")
  
  # Create scaling matrix
  
  
  
  # determine count and proportion of population that traveled at the time point
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz) %>%
    mutate(all.infected = sum(count[SEIR %in% c("I_h", "A_h")]),
           all.infected.prop = all.infected/total.pop)
  
  # ggplot(subset(model.output_long, SEIR %in% c("I_h")), aes(time.idx, prop.pop, color = as.factor(subpop)))+
  #   geom_line()+
  #   geom_line(aes(time.idx, all.infected.prop))
  
  
  model.output_long$leaving.i <- model.output_long$arriving.i <-rep(0, nrow = nrow(model.output_long))
  
  for(i in 1:length(unique(model.output_long$time.idx))){
    infected.df <- subset(model.output_long, time.idx == i & hum.or.moz == "h") %>%
      distinct(time.idx, subpop, all.infected.prop)
    
    ii <- infected.df$subpop[which(infected.df$all.infected.prop > P[["threshold"]])]
    # scale.factor <- matrix( 1, nrow = subpop, ncol = subpop)
    # scale.factor[ , ii] <- P[["outbound.scaler"]]  # scale outbound travel from places with prevalence > threshold
    # scale.factor[ii , ] <- P[["inbound.scaler"]]
    #
    scale.factor.SERA <- scale.factor.I <- matrix( 1, nrow = subpop, ncol = subpop)
    scale.factor.SERA[ , ii] <- P[["outbound.scaler"]]  # scale outbound travel from places with prevalence > threshold
    scale.factor.SERA[ii, ] <- P[["inbound.scaler"]]  # scale inbound travel to places with prevalence > threshold
    scale.factor.I[ , ii] <- P[["outbound.scaler.I"]]
    scale.factor.I[ii, ] <- P[["inbound.scaler.I"]]
    #
    # print(infected.df$all.infected.prop)
    # print(scale.factor)
    
    model.output_long$leaving.i[model.output_long$time.idx == i] <- ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("S_h"), model.output_long$count[model.output_long$time.idx == i] * rowSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                           ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("E_h"), model.output_long$count[model.output_long$time.idx == i] * rowSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                  ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("R_h"), model.output_long$count[model.output_long$time.idx == i] * rowSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                         ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("A_h"), model.output_long$count[model.output_long$time.idx == i] * rowSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                                ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("I_h"), model.output_long$count[model.output_long$time.idx == i] * rowSums(P[["prob.travel.I"]] * scale.factor.I),0)))))
    
    model.output_long$arriving.i[model.output_long$time.idx == i] <- ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("S_h"), model.output_long$count[model.output_long$time.idx == i] * colSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                            ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("E_h"), model.output_long$count[model.output_long$time.idx == i] * colSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                   ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("R_h"), model.output_long$count[model.output_long$time.idx == i] * colSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                          ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("A_h"), model.output_long$count[model.output_long$time.idx == i] * colSums(P[["prob.travel.SERA"]] * scale.factor.SERA),
                                                                                                 ifelse(model.output_long$SEIR[model.output_long$time.idx == i] %in% c("I_h"), model.output_long$count[model.output_long$time.idx == i] * colSums(P[["prob.travel.I"]] * scale.factor.I), 0)))))
    
  }
  
  # model.output_long$leaving.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h") & model.output_long$all.infected.prop > P[["threshold"]], model.output_long$count * rowSums(P[["prob.travel.SERA"]]) * P[["threshold.scaler"]],
  #                                       ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h") & model.output_long$all.infected.prop > P[["threshold"]], model.output_long$count * rowSums(P[["prob.travel.SERA"]]) ,
  #                                       ifelse(model.output_long$SEIR %in% c("I_h") & model.output_long$all.infected.prop > P[["threshold"]], model.output_long$count * rowSums(P[["prob.travel.I"]]) * P[["threshold.scaler"]],
  #                                       ifelse(model.output_long$SEIR %in% c("I_h") & model.output_long$all.infected.prop > P[["threshold"]], model.output_long$count * rowSums(P[["prob.travel.I"]]),0))))
  #
  #
  # model.output_long$arriving.i <- ifelse(model.output_long$SEIR %in% c("S_h", "E_h", "R_h", "A_h"), model.output_long$count * colSums(P[["prob.travel.SERA"]]),
  #                                        ifelse(model.output_long$SEIR %in% c("I_h"), model.output_long$count * colSums(P[["prob.travel.I"]]), 0))
  #
  
  # ggplot(subset(model.output_long, SEIR %in% c("I_h") & subpop == "1"), aes(time.idx, leaving.i))+
  #   geom_line()
  
  
  model.output_long <- model.output_long %>%
    group_by(time.idx, subpop, hum.or.moz)%>%
    mutate(leaving.total = sum(leaving.i),
           arriving.total = sum(arriving.i),
           prop.leaving.i = leaving.i/total.pop,
           prop.arriving.i = arriving.i/total.pop,
           prop.leaving.tot = leaving.total/total.pop,
           prop.arriving.tot = arriving.total/total.pop)
  
  # ggplot(subset(model.output_long, SEIR %in% c("I_h") & subpop == "1"), aes(time.idx, prop.arriving.tot))+
  #   geom_line()
  
  # model.output_long <- model.output_long[ , c("time.idx", "subpop", "SEIR", "count", "hum.or.moz", "total.pop",
  #                                             "prop.pop", "t.months", "leaving.i", "arriving.i", "leaving.total",
  #                                             "arriving.total", "prop.leaving.i", "prop.arriving.i", "prop.leaving.tot",
  #                                             "prop.arriving.tot")]
  
  # create category for all infectious people (I and A)
  infected <- subset(model.output_long, SEIR %in% c("I_h", "A_h"))
  infected.1 <- infected %>%
    group_by(time.idx, subpop)%>%
    summarise(SEIR = "all.infected_h",
              count = sum(count),
              hum.or.moz = "h",
              total.pop = total.pop,
              prop.pop = count/total.pop,
              t.months = t.months,
              leaving.i = sum(leaving.i),
              arriving.i = sum(arriving.i),
              leaving.total = leaving.total,
              arriving.total = arriving.total,
              prop.leaving.i = leaving.i/total.pop,
              prop.arriving.i = arriving.i/total.pop,
              prop.leaving.tot = leaving.total/total.pop,
              prop.arriving.tot = arriving.total/total.pop
    )%>%
    distinct(time.idx, subpop, .keep_all = T)
  
  # ggplot(subset(infected.1, subpop == "1"), aes(time.idx, prop.arriving.tot))+
  #   geom_line()
  
  model.output_long <- model.output_long[ , colnames(model.output_long) %in% colnames(infected.1)]
  
  model.output_long <- rbind.data.frame(model.output_long, infected.1)
  model.output_long$model <- rep(model, nrow(model.output_long))
  
  return(model.output_long)
}



## System of equations for malaria tranmission with mobility between multiple subpopulations

biting <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    
    O <- unname(y[names(y) %like% 'O'])
    L <- unname(y[names(y) %like% 'L'])
    P <- unname(y[names(y) %like% 'P'])
    S_m <- unname(y[names(y) %like% 'S_m'])
    E_m <- unname(y[names(y) %like% 'E_m'])
    I_m <- unname(y[names(y) %like% 'I_m'])
    S_h <- unname(y[names(y) %like% 'S_h'])
    E_h <- unname(y[names(y) %like% 'E_h'])
    I_h <- unname(y[names(y) %like% 'I_h'])
    R_h <- unname(y[names(y) %like% 'R_h'])
    A_h <- unname(y[names(y) %like% 'A_h'])
    dO <- rep(0,subpop)
    dL <- rep(0,subpop)
    dP <- rep(0,subpop)
    dS_m <- rep(0,subpop)
    dE_m <- rep(0,subpop)
    dI_m <- rep(0,subpop)
    dS_h <- rep(0,subpop)
    dE_h <- rep(0,subpop)
    dI_h <- rep(0,subpop)
    dR_h <- rep(0,subpop)
    dA_h <- rep(0,subpop)
   
    # EQs and ODEs
    V <- S_m + E_m + I_m                # total mosquito population
    H <- S_h + E_h + I_h + R_h + A_h    # total human population
    
    b_h <-  a * b    #  proportion of bites resulting in human infections 
    b_m <-  a * c     #  proportion of bites resulting in mosquito infections    
    
    K = K_0 * (1 + delta * cos(2 * pi * (t/365 + omega)))              # seasonal carrying capacity of environment (sinusoidal seasonal pattern, need to convert t is in years)
    xi_m = xi_m_0 / ( 1 - delta * cos(2 * pi * (t/365 + omega)))       # seasonal maturation rate of P. falciparum maturation
    mu_m = mu_m_0 / ( 1 + delta * cos(2 * pi * (t/365 + omega)))       # seasonal death rate of mosquitoes
    # seasonal.travel = (1 + delta_t * cos(2 * pi * (t/365 + omega_t)))  # seasonal fluctuation in travel
    travel.rate <- travelers/H
    
    # humans
    dS_h <- r * (E_h + I_h + A_h) + q2 * R_h - b_h * m * S_h * I_m  + colSums(travel.rate * S_h) - rowSums(travel.rate * S_h)      #+ colSums(prob.travel * seasonal.travel * S_h) - rowSums(prob.travel * seasonal.travel * S_h)      # susceptible humans
    dE_h <- b_h * m * S_h * I_m  - (xi_h + r) * E_h                 + colSums(travel.rate * E_h) - rowSums(travel.rate * E_h)      #+ colSums(prob.travel * seasonal.travel * E_h) - rowSums(prob.travel * seasonal.travel * E_h)      # exposed humans
    dI_h <- xi_h * E_h - (q1 + r) * I_h                             + adjust.I * colSums(travel.rate * I_h) - adjust.I * rowSums(travel.rate * I_h)      #+ colSums(prob.travel * seasonal.travel * I_h) - rowSums(prob.travel * seasonal.travel * I_h)      # infected humans
    dR_h <- q1 * (I_h + A_h) - (theta * b_h * m * I_m + q2) * R_h   + colSums(travel.rate * R_h) - rowSums(travel.rate * R_h)      #+ colSums(prob.travel * seasonal.travel * R_h) - rowSums(prob.travel * seasonal.travel * R_h)      # recovered/immune humans
    dA_h <- theta * b_h * m * R_h * I_m - (q1 + r) * A_h            + colSums(travel.rate * A_h) - rowSums(travel.rate * A_h)      # + colSums(prob.travel * seasonal.travel * A_h) - rowSums(prob.travel * seasonal.travel * A_h)      # asymptomatically infected humans
    
    # mosquitoes
    dO <- beta_m * V - d_o * O - mu_o * (1 + (O + L)/ K) * O                 # eggs and early larval instars
    dL <- d_o * O - d_l * L - mu_l * (1 + gam * (O + L)/ K) * L              # late larval instars
    dP <- d_l * L - d_p * P - mu_p * P                                       # pupae
    dS_m <- 1/2 * d_p * P - b_m * (I_h + sigma * A_h) * S_m - mu_m * S_m     # susceptible mosquitoes
    dE_m <- b_m * (I_h + sigma * A_h) * S_m - (xi_m + mu_m) * E_m            # exposed mosquitoes  
    dI_m <- xi_m * E_m - mu_m * I_m                                          # infected mosquitoes
    
   # output
    res <- c(dS_h, dE_h, dI_h, dR_h, dA_h, dO, dL, dP, dS_m, dE_m, dI_m)
    list(res)
  })
}
biting.travel.matrix <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    
    O <- unname(y[names(y) %like% 'O'])
    L <- unname(y[names(y) %like% 'L'])
    P <- unname(y[names(y) %like% 'P'])
    S_m <- unname(y[names(y) %like% 'S_m'])
    E_m <- unname(y[names(y) %like% 'E_m'])
    I_m <- unname(y[names(y) %like% 'I_m'])
    S_h <- unname(y[names(y) %like% 'S_h'])
    E_h <- unname(y[names(y) %like% 'E_h'])
    I_h <- unname(y[names(y) %like% 'I_h'])
    R_h <- unname(y[names(y) %like% 'R_h'])
    A_h <- unname(y[names(y) %like% 'A_h'])
    dO <- rep(0,subpop)
    dL <- rep(0,subpop)
    dP <- rep(0,subpop)
    dS_m <- rep(0,subpop)
    dE_m <- rep(0,subpop)
    dI_m <- rep(0,subpop)
    dS_h <- rep(0,subpop)
    dE_h <- rep(0,subpop)
    dI_h <- rep(0,subpop)
    dR_h <- rep(0,subpop)
    dA_h <- rep(0,subpop)
    
    # EQs and ODEs
    V <- S_m + E_m + I_m                # total mosquito population
    H <- S_h + E_h + I_h + R_h + A_h    # total human population
    
    b_h <-  a * b    #  proportion of bites resulting in human infections 
    b_m <-  a * c     #  proportion of bites resulting in mosquito infections    
    
    K = K_0 * (1 + delta * cos(2 * pi * (t/365 + omega)))              # seasonal carrying capacity of environment (sinusoidal seasonal pattern, need to convert t is in years)
    xi_m = xi_m_0 / ( 1 - delta * cos(2 * pi * (t/365 + omega)))       # seasonal maturation rate of P. falciparum maturation
    mu_m = mu_m_0 / ( 1 + delta * cos(2 * pi * (t/365 + omega)))       # seasonal death rate of mosquitoes
    
    # prob.travel matrix
    t.idx <- floor(t)  # index time such that probability of travel for incremental time steps is associated with the correct day
    travel.rate <- (matrix(travelers[ ,colnames(travelers) == t.idx], ncol = subpop, nrow = subpop))/H
    
    # humans
    dS_h <- r * (E_h + I_h + A_h) + q2 * R_h - b_h * m * S_h * I_m  + colSums(travel.rate * S_h) - rowSums(travel.rate * S_h)      # susceptible humans
    dE_h <- b_h * m * S_h * I_m  - (xi_h + r) * E_h                 + colSums(travel.rate * E_h) - rowSums(travel.rate * E_h)      # exposed humans
    dI_h <- xi_h * E_h - (q1 + r) * I_h                             + colSums(travel.rate * I_h) - rowSums(travel.rate * I_h)      # infected humans
    dR_h <- q1 * (I_h + A_h) - (theta * b_h * m * I_m + q2) * R_h   + colSums(travel.rate * R_h) - rowSums(travel.rate * R_h)      # recovered/immune humans
    dA_h <- theta * b_h * m * R_h * I_m - (q1 + r) * A_h            + colSums(travel.rate * A_h) - rowSums(travel.rate * A_h)      # asymptomatically infected humans
    
    # mosquitoes
    dO <- beta_m * V - d_o * O - mu_o * (1 + (O + L)/ K) * O                 # eggs and early larval instars
    dL <- d_o * O - d_l * L - mu_l * (1 + gam * (O + L)/ K) * L              # late larval instars
    dP <- d_l * L - d_p * P - mu_p * P                                       # pupae
    dS_m <- 1/2 * d_p * P - b_m * (I_h + sigma * A_h) * S_m - mu_m * S_m     # susceptible mosquitoes
    dE_m <- b_m * (I_h + sigma * A_h) * S_m - (xi_m + mu_m) * E_m            # exposed mosquitoes  
    dI_m <- xi_m * E_m - mu_m * I_m                                          # infected mosquitoes
    
    # output
    res <- c(dS_h, dE_h, dI_h, dR_h, dA_h, dO, dL, dP, dS_m, dE_m, dI_m)
    list(res)
  })
}
biting.threshold.mob <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    
    O <- unname(y[names(y) %like% 'O'])
    L <- unname(y[names(y) %like% 'L'])
    P <- unname(y[names(y) %like% 'P'])
    S_m <- unname(y[names(y) %like% 'S_m'])
    E_m <- unname(y[names(y) %like% 'E_m'])
    I_m <- unname(y[names(y) %like% 'I_m'])
    S_h <- unname(y[names(y) %like% 'S_h'])
    E_h <- unname(y[names(y) %like% 'E_h'])
    I_h <- unname(y[names(y) %like% 'I_h'])
    R_h <- unname(y[names(y) %like% 'R_h'])
    A_h <- unname(y[names(y) %like% 'A_h'])
    dO <- rep(0,subpop)
    dL <- rep(0,subpop)
    dP <- rep(0,subpop)
    dS_m <- rep(0,subpop)
    dE_m <- rep(0,subpop)
    dI_m <- rep(0,subpop)
    dS_h <- rep(0,subpop)
    dE_h <- rep(0,subpop)
    dI_h <- rep(0,subpop)
    dR_h <- rep(0,subpop)
    dA_h <- rep(0,subpop)
    
    # EQs and ODEs
    V <- S_m + E_m + I_m                # total mosquito population
    H <- S_h + E_h + I_h + R_h + A_h    # total human population
    A.I <- (I_h + A_h)/H                    # total infections
    
    # determine if threshold has been surpassed & modify scaler as necessary
    ii <- which(A.I > threshold)
    scale.factor.SERA <- scale.factor.I <- matrix( 1, nrow = subpop, ncol = subpop)
    scale.factor.SERA[ , ii] <- outbound.scaler  # scale outbound travel from places with prevalence > threshold
    scale.factor.SERA[ii, ] <- inbound.scaler  # scale inbound travel to places with prevalence > threshold
    scale.factor.I[ , ii] <- outbound.scaler.I
    scale.factor.I[ii, ] <- inbound.scaler.I
    
    b_h <-  a * b    #  proportion of bites resulting in human infections
    b_m <-  a * c     #  proportion of bites resulting in mosquito infections
    
    K = K_0 * (1 + delta * cos(2 * pi * (t/365 + omega)))              # seasonal carrying capacity of environment (sinusoidal seasonal pattern, need to convert t is in years)
    xi_m = xi_m_0 / ( 1 - delta * cos(2 * pi * (t/365 + omega)))       # seasonal maturation rate of P. falciparum maturation
    mu_m = mu_m_0 / ( 1 + delta * cos(2 * pi * (t/365 + omega)))       # seasonal death rate of mosquitoes
    
    # humans
    dS_h <- r * (E_h + I_h + A_h) + q2 * R_h - b_h * m * S_h * I_m  + colSums(prob.travel.SERA * scale.factor.SERA * S_h) - rowSums(prob.travel.SERA * scale.factor.SERA * S_h)      # susceptible humans
    dE_h <- b_h * m * S_h * I_m  - (xi_h + r) * E_h                 + colSums(prob.travel.SERA * scale.factor.SERA * E_h) - rowSums(prob.travel.SERA * scale.factor.SERA * E_h)      # exposed humans
    dI_h <- xi_h * E_h - (q1 + r) * I_h                             + colSums(prob.travel.I * scale.factor.I * I_h) - rowSums(prob.travel.I * scale.factor.I * I_h)      # infected humans
    dR_h <- q1 * (I_h + A_h) - (theta * b_h * m * I_m + q2) * R_h   + colSums(prob.travel.SERA * scale.factor.SERA * R_h) - rowSums(prob.travel.SERA * scale.factor.SERA * R_h)      # recovered/immune humans
    dA_h <- theta * b_h * m * R_h * I_m - (q1 + r) * A_h            + colSums(prob.travel.SERA * scale.factor.SERA * A_h) - rowSums(prob.travel.SERA * scale.factor.SERA * A_h)      # asymptomatically infected humans
    
    # mosquitoes
    dO <- beta_m * V - d_o * O - mu_o * (1 + (O + L)/ K) * O                 # eggs and early larval instars
    dL <- d_o * O - d_l * L - mu_l * (1 + gam * (O + L)/ K) * L              # late larval instars
    dP <- d_l * L - d_p * P - mu_p * P                                       # pupae
    dS_m <- 1/2 * d_p * P - b_m * (I_h + sigma * A_h) * S_m - mu_m * S_m     # susceptible mosquitoes
    dE_m <- b_m * (I_h + sigma * A_h) * S_m - (xi_m + mu_m) * E_m            # exposed mosquitoes
    dI_m <- xi_m * E_m - mu_m * I_m                                          # infected mosquitoes
    
    ##
    # print(t)
    # print(A.I)
    # print(prob.travel.SERA*scale.factor)
    
    # output
    res <- c(dS_h, dE_h, dI_h, dR_h, dA_h, dO, dL, dP, dS_m, dE_m, dI_m)
    list(res)
  })
}



