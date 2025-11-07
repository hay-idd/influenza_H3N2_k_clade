setwd("~/Documents/GitHub/influenza_H3N2_k_clade")
## Functions we need
source("auxiliary_funcs.R")
source("sir_functions.R")

## Libraries we need
library(socialmixr)
library(pracma)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

## Demography data
## Population of England
N_tot <- 60000000

## POLYMOD contact matrix
data("polymod")
polymod_c <- contact_matrix(polymod,
                            countries="United Kingdom",
                            age.limits = c(0,5,15,25,35,45,55,65,75,85),
                            symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)
image(C)

## Solve base model for comparison
N_props <- polymod_c$participants$proportion
N_age_classes <- length(N_props)
N_immunity_classes <- 2

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)

## Relative susceptibility of each age group
## Vary R0, prop immune and relative susceptibility to give a reasonable epidemic size and duration
## Fix prop susc in 0-5 at 0
## Vary 5-15 
## Set all older ages to be the same
## 5 parameters: R0, prop immune in 5-15, 15-25, prop immune in older ages, relative susceptibility in immune group
R0 <- 1.2
gamma <- 3
prop_immune_younger <- 0.75
prop_immune_older <- 0.75
prop_immune <- c(0,prop_immune_younger,prop_immune_younger, rep(prop_immune_older,7))
beta_scales <- rep(1, N_age_classes)
alphas <- c(1.2,1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

beta_scales <- rep(1,N_age_classes)

C_use <- setup_C(C, N, beta_scales)
beta_par <- get_beta(C,polymod_c$participants$proportion,gamma, R0)
y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                       alphas=alphas, age_seed=4,immunity_seed=1,return_compartments=TRUE)
y_base$t <- 1:nrow(y_base)
y_base <- y_base %>% pivot_longer(-t)
y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
                            age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
                            immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
total_inf <- y_base %>% filter(t == max(t),compartment == "inc")  %>% pull(value) %>% sum()

inc <- y_base %>% filter(compartment == "inc")
inc <- inc %>% group_by(age, immunity) %>% mutate(value = value - lag(value, 1))



inc <- inc %>% left_join(polymod_c$participants %>% mutate(age = 1:nrow(polymod_c$participants)))
inc$age.group <- factor(inc$age.group,levels=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                               "[55,65)", "[65,75)", "75+"))

inc <- inc %>% group_by(age.group,t) %>% summarize(value = sum(value)) %>% group_by(age.group) %>% mutate(gr = log(value/lag(value,1)))

ggplot(inc%>% filter(t < 75) %>% ungroup() %>% mutate(value = value/max(value,na.rm=TRUE))) + geom_line(aes(x=t,y=value,col=age.group)) + scale_color_viridis_d() 
#ggplot(inc %>% filter(t < 200)) + geom_line(aes(x=t,y=gr,col=age.group)) + scale_color_viridis_d() 

print(total_inf/sum(N))
