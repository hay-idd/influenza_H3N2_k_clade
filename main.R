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
## Distribution of cases
china_cases <- read_csv("data/case_data_2020.02.11.csv")

## Age distribution of China
age_distribution <- read_csv("data/population_demography_10year.csv")
age_distribution <- age_distribution %>% mutate(propn = number/sum(number))
age_distribution$age_group <- factor(age_distribution$age_group, levels=unique(age_distribution$age_group))
age_distribution <- age_distribution %>% mutate(wuhan=N_tot*propn)

## Population of Wuhan
N_tot <- 11080000

## POLYMOD contact matrix
data("polymod")
N_props1 <- c(0.118574496482027, 0.115752336341461, 0.128634781768929, 0.15898446700289, 
              0.150148368364849, 0.154367713154031, 0.105371761585134, 0.0496726586040798, 
              0.0184934166965997)
age_dat <- data.frame(lower.age.limit=seq(0,80,by=10),population=N_props1)
polymod_c <- contact_matrix(polymod,survey.pop=age_dat,age.limits = seq(0,80,by=10),symmetric=TRUE,
                            missing.contact.age = "sample",
                            missing.participant.age = "remove")
C <- polymod_c$matrix
row.names(C) <- colnames(C)


## Solve base model for comparison
N_props <- age_distribution %>% pull(propn)
N_age_classes <- length(N_props1)
N_immunity_classes <- 1

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)

## Relative susceptibility of each age group
alphas <- c(1)

R0 <- 1.4
gamma <- 5

prop_immune <- rep(0, N_age_classes)
beta_scales <- rep(1, N_age_classes)
alphas <- c(1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

beta_scales <- rep(1,9)

C_use <- setup_C(C, N, beta_scales)
beta_par <- get_beta(C,age_dat$population,gamma, R0)
y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
                       alphas=alphas, age_seed=4,immunity_seed=1,return_compartments=TRUE)
head(y_base)
