# This ODE model simulates the transmission of malaria between mosquito and human populations. 
# No interventions are present. 
# Seasonality impacts 1) carrying capacity for mosquito population (K)
#                     2) incubation period of malaria (xi_m)
#                     3) death rate of mosquito (mu_m)
# Outputs include plots of seasonal dynamics of malaria cases and mosquito population

# by Hannah Meredith
# last updated May 25, 2021


# libraries

library("pracma")
library("deSolve")
library("ggplot2")
library("readxl")
library("grid")
library("gridExtra")
library("reshape2")
library("ggpubr")
library("data.table")
library("tidyr")
library("dplyr")
library("scales")

source("Malaria_Mobility_Functions.R")
source("Malaria_Mobility_Plots.R")

# define parameters and  ranges-----

subpop = 2

# Parameters
P <- list(
  beta_m = 21.19, # number of eggs a mosquito lays per day [1/day]
  d_o = 0.15, # rate at which early larval instars mature into late larval instars [1/day]
  d_l = 0.27, # rate at which late larval instars mature into pupae [1/day]
  d_p = 1.56, # rate at which pupae matures into mosquitoes [1/day]
  a = 0.01, # 0.2,    # biting frequency [1/day]
  b = 0.2, #0.5,    # proportion of bites that produce infection in humans [unitless]
  c = 0.5,    # proportion of bites that produce infection in mosquitoes [unitless]
  sigma = 0.25, # adjustment factor for asymptomatic infection transmissibility to vector [unitless]
  xi_m_0 = 1/10, # average rate of P. falciparum maturation in mosquito [1/day]
  xi_h = 1/21, # rate of P. falciparum maturation in human [1/day]
  q1 = 1/200, # rate of immunity acquisition [1/day]
  q2 = 1/1000, # rate of immunity loss [1/day]
  theta = 0.5, # level of reduced susceptibility to secondary infection [unitless]
  r = 0.01, # rate of recovery [1/day]
  mu_o = 0.034, # death rate of early larval instars [1/day]
  mu_l = 0.035, # death rate of late larval instars [1/day]
  K_0 = c(2.5,5),#100, # average carrying capacity of environment assuming 50 mm water over the past 4.5 days and 150 mosquitoes captured
  gam = 13.25, # adjustment factor to correct for different density dependence of late vs early instars
  mu_p = 0.25, # death rate of pupae [1/day]
  mu_m_0 = 0.12, # average death rate of mosquitoes [1/day]
  m = c(5,10), # ratio of mosquitoes to humans
  delta = 0.15, # amplitude of seasonal forcing
  omega = 0, # phase of seasonal forcing
  t_ss = 200,   # time allowed to reach steady state [days] - this is done to let the model stabilize before you introduce interventions
  treatment_pd = 1 * 365,   # treatment period [days]
  subpop = subpop, # number of subpopulations to model
  delta_t = matrix(c(0, 0, 0 ,0), nrow = 2, ncol = 2),   # amplitude of seasonal mobility : 0.5
  omega_t = matrix(c(0 ,0 ,0, 0), nrow = 2, ncol = 2),  # phase of seasonal mobility
  adjust.I = 1) # scaling factor to adjust number of Infected travelers, relative to other groups

# define initial conditions
y0 <- c(
  S_h = c(1000, 10000), #rep(100000, subpop),    # susceptible humans
  E_h = rep(0, subpop),   # exposed humans
  I_h = rep(0, subpop),   # infected humans
  R_h = rep(0, subpop),   # recovered humans
  A_h = rep(0, subpop),   # asymptomatic humans
  O = rep(100,  subpop),    # eggs
  L = rep(0, subpop),     # larvae
  P = rep(0,  subpop),    # pupae
  S_m = rep(100, subpop),   # susceptible mosquitoes
  E_m = rep(0, subpop),   # exposed mosquitoes
  I_m = rep(50,  subpop))  # infected mosquitoes


##1. Simulate transmission without any movement (baseline) ----
travelers = matrix(c(0, 0, 0, 0), nrow = subpop, ncol = subpop) # Assuming symmetrical mobility patterns
P[['travelers']] <- travelers
no.movement <- as.data.frame(mal.constant.mob.model(y0, P, "no.movement"))  

## 2. Simulate transmission with constant movement (constant.movement) ----
travelers = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop) # Assuming symmetrical mobility patterns
P[['travelers']] <- travelers
base.movement <- as.data.frame(mal.constant.mob.model(y0, P, "constant.movement"))

# ## 3. Simulate transmission with bi-directional movement peaking in sync with malaria (peak.movement) ----
# travelers.base = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop)
# travelers.holiday = matrix(c(0, 50, 50, 0), nrow = subpop, ncol = subpop)
# travel.days = seq(as.Date("2020-06-01"), as.Date("2020-09-01"),1)  # enter start and end date
# travel.days = c(travel.days)
# 
# travelers <- make.constant.travel.rate.matrix(travel.dates, travelers.base, travelers.holiday)
# P[['travelers']] <- travelers
# peak.movement <- as.data.frame(mal.dynamic.mob.model(y0, P, "peak.movement"))  
# peak.movement.1 <- peak.movement[ ,!colnames(peak.movement) %in% c("inbound", "outbound")]
# 
# ## 4. Simulate transmission with set days of movement peaking in sync with malaria (peak.movement.stay) ----
# travelers.base = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop)
# travelers.holiday <- cbind.data.frame(c(0, 50, 5, 0),
#                                         c(0, 5, 50, 0))
# travel.out.days <- seq(as.Date("2020-06-01"), as.Date("2020-06-10"), 1) # enter start and end date for outbound trips
# travel.back.days <- seq(as.Date("2020-08-26"), as.Date("2020-09-04"), 1)  # enter start and end date for return trips
# travelers <- make.specific.date.travel.rate.matrix(travel.out.days, travel.back.days, travelers.base, travelers.holiday)
# P[['travelers']] <- travelers
# peak.movement.2 <- as.data.frame(mal.dynamic.mob.model(y0, P, "peak.and.stay.movement"))  
# peak.movement.2 <- peak.movement.2[ ,!colnames(peak.movement.2) %in% c("inbound", "outbound")]
# 
# ## 5. Simulate transmission with bi-directional movement peaking out of sync with malaria (off.peak.movement) ----
# travelers.base = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop)
# travelers.holiday = matrix(c(0, 50, 50, 0), nrow = subpop, ncol = subpop)
# travel.days = seq(as.Date("2020-01-01"), as.Date("2020-04-01"),1)  # enter start and end date
# travel.days = c(travel.days)
# 
# travelers <- make.constant.travel.rate.matrix(travel.dates, travelers.base, travelers.holiday)
# P[['travelers']] <- travelers
# off.peak.movement <- as.data.frame(mal.dynamic.mob.model(y0, P, "off.peak.movement"))  
# off.peak.movement <- off.peak.movement[ ,!colnames(off.peak.movement) %in% c("inbound", "outbound")]
# 
# ## 6. Simulate transmission with set days of movement peaking in sync with malaria (off.peak.movement.stay) ----
# travelers.base = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop)
# travelers.holiday <- cbind.data.frame(c(0, 50, 5, 0),
#                                       c(0, 5, 50, 0))
# travel.out.days <- seq(as.Date("2020-01-01"), as.Date("2020-01-10"), 1) # enter start and end date for outbound trips
# travel.back.days <- seq(as.Date("2020-04-26"), as.Date("2020-05-05"), 1)  # enter start and end date for return trips
# 
# # travel.days <- c(as.Date("2020-01-01"), as.Date("2020-04-01"))
# travelers <- make.specific.date.travel.rate.matrix(travel.out.days, travel.back.days, travelers.base, travelers.holiday)
# P[['travelers']] <- travelers
# off.peak.movement.2 <- as.data.frame(mal.dynamic.mob.model(y0, P, "off.peak.and.stay.movement"))  
# off.peak.movement.2 <- off.peak.movement.2[ ,!colnames(off.peak.movement.2) %in% c("inbound", "outbound")]
# 
# 
# ## FIGURE 1: Different forms of seasonal mobility and its impact on transmission
# 
# # Comobine movement scenarios for plotting
# movement.scenarios <- rbind(peak.movement.1, peak.movement.2, off.peak.movement,off.peak.movement.2)
# movement.scenarios <- left_join(movement.scenarios, base.movement[, c("time.idx", "SEIR", "subpop", "prop.pop", "leaving.prop.total.i")], by = c("time.idx", "SEIR", "subpop"))
# movement.scenarios$model.name <- factor(movement.scenarios$model, 
#                                         levels = c("off.peak.and.stay.movement", "peak.and.stay.movement", "off.peak.movement", "peak.movement"),
#                                         labels = c("Relocate for Off-peak", "Relocate for Peak", "Increase for Off-peak", "Increase for Peak"))
# # Generate plots
# figure1 <- plot.cases.movement.impact(movement.scenarios, base.movement)
# 
# # # print figures/.save as png # uncommment to save figures
# # png("Fig1_cases.png", width = 6, height = 2.25, units = 'in', res = 600)
# # figure1[[1]]
# # dev.off()
# # png("Fig1_travelers.png", width = 6, height = 2.25, units = 'in', res = 600)
# # figure1[[2]]
# # dev.off()
# # png("Fig1_case_impact.png", width = 4, height = 4, units = 'in', res = 600)
# # figure1[[3]]
# # dev.off()


## 7. Simulate transmission with reduced mobility rates for infected people ------
travelers = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop) 
P[['travelers']] <- travelers
P[['adjust.I']] <- 0.1 # infected, symptomatic people move at 1/10th the rate
reduced.I.constant.movement <- as.data.frame(mal.constant.SEIRA.mob.model(y0, P, "reduced.I.movement"))  
base.reduced.I.constant.movement.plot <- plot.baseline.comparison(base.movement, reduced.I.constant.movement, "Constant movement", "Constant movement \n Reduced I rate", "X")

## 8. Simulate transmission with reduced mobility into/out of location when cases are above threshold------
travelers = matrix(c(0, 5, 5, 0), nrow = subpop, ncol = subpop) 
P[['travelers']] <- travelers
P[['adjust.I']] <- 0.1 # infected, symptomatic people move at 1/10th the rate
P[["inbound.scaler"]] = 0.1  # once threshold is exceeded, travel rates from all groups into/out of location are reduced to 10% of normal
P[["outbound.scaler"]] = 0.1  
P[["inbound.scaler.I"]] = 0.1  
P[["outbound.scaler.I"]] = 0.1  
P[["threshold"]] = 0.35        

threshold.reduced.all.movement <- as.data.frame(mal.threshold.mob.model(y0, P, "threshold.movement"))  # this is where the first function, baseline, is called. Baseline then calls the ODE solver function.
threshold.reduced.all.movement <- threshold.reduced.all.movement[ ,!colnames(threshold.reduced.all.movement) %in% c("inbound", "outbound")]
threshold.reduced.all.movement.plot <- plot.baseline.comparison(baseline, threshold.reduced.all.movement, "No movement", "Constant movement \n Threshold Reduce In/Out", "test")

## 9. Simulate transmission with reduced mobility into and increased mobility out of location when cases are above threshold------
prob.travel.SERA = matrix(c(0, 0.01, 0.01, 0), nrow = subpop, ncol = subpop) # Assuming 1 % of the population from SERA categories moves per day
prob.travel.I = matrix(c(0, 0.001, 0.001, 0), nrow = subpop, ncol = subpop) # Assuming 0.1 % of the population from I category moves per day
P[['prob.travel.SERA']] <- prob.travel.SERA
P[['prob.travel.I']] <- prob.travel.I
P[["inbound.scaler"]] = 0.1  # once threshold is exceeded, travel rates into location are reduced to 10% of normal
P[["outbound.scaler"]] = 1.1  # once threshold is exceeded, travel rates out of location increased by 50% normal
P[["inbound.scaler.I"]] = 0.1  # once threshold is exceeded, travel rates into location are reduced to 10% of normal
P[["outbound.scaler.I"]] = 0.1  # once threshold is exceeded, travel rates out of location are reduced to 10% of normal
P[["threshold"]] = 0.35        # once prevalence exceeds 35%, scale back travel to/from location with outbreak

threshold.reduced.in.increased.out.movement <- as.data.frame(mal.threshold.mob.model(y0, P, "threshold.movement2"))  # this is where the first function, baseline, is called. Baseline then calls the ODE solver function.
threshold.reduced.in.increased.out.movement <- threshold.reduced.in.increased.out.movement[ ,!colnames(threshold.reduced.in.increased.out.movement) %in% c("inbound", "outbound")]
threshold.reduced.in.increased.out.movement.plot <- plot.baseline.comparison(baseline, threshold.reduced.in.increased.out.movement, "No movement", "Constant movement \n Threshold Reduce In/Increase Out", "test")

test.no.mob <- subset(baseline, SEIR == "all.infected_h")%>%
  group_by(subpop, model)%>%
  summarise(yearly.infected = sum(count))

test.mob <-subset(threshold.reduced.in.increased.out.movement, SEIR == "all.infected_h")%>%
  group_by(subpop, model)%>%
  summarise(yearly.infected = sum(count))



