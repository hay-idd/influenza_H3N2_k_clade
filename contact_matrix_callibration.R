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
library(data.table)
library(lazymcmc)


prop_work_contacts_in_hols <- 0.4
prop_home_contacts_in_hols <- 2

## Demography data
## Population of England
N_tot <- 60000000

########################################################
## CREATE CONTACT MATRICES
########################################################

## POLYMOD contact matrix
data("polymod")

contacts <- polymod$contacts
## Create new contact survey with no school contacts, decreased work contacts, and increased home contacts
contacts_no_schools <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home==0)# & cnt_leisure == 0)
## Resample to give frac of work contacts
contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(prop_work_contacts_in_hols)
## Resample with replacement to give more home contacts
contacts_home <- contacts %>% filter(cnt_home == 1)
contacts_home <- contacts_home[sample(1:nrow(contacts_home), size = nrow(contacts_home)*prop_home_contacts_in_hols, replace=TRUE), ]
contacts_holidays <- bind_rows(contacts_no_schools,contacts_work,contacts_home)

## Create holiday contact matrix
polymod1 <- polymod
polymod1$contacts <- contacts_holidays
polymod_c_holidays <- contact_matrix(polymod1,
                                     countries="United Kingdom",
                                     age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                     symmetric=TRUE,
                                     missing.contact.age = "sample",
                                     missing.participant.age = "remove")
C_holidays <- polymod_c_holidays$matrix
row.names(C_holidays) <- colnames(C_holidays)

## Create term time contact matrix
polymod_c_term <- contact_matrix(polymod,
                                 countries="United Kingdom",
                                 age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                 symmetric=TRUE,
                                 missing.contact.age = "sample",
                                 missing.participant.age = "remove")
C_term <- polymod_c_term$matrix
row.names(C_term) <- colnames(C_term)

########################################################
## READ IN INCIDENCE DATA FOR FITTING
########################################################
start_date <- as.Date("2023-09-01")
end_date <- as.Date("2024-05-01")
start_day <- as.numeric(start_date)
end_day <- as.numeric(end_date)
date_seq <- seq(start_date,end_date,by="1 day")
ts <- seq(start_day,end_day,by=1) - start_day + 1

## From https://www.kennedycollegeoxford.org.uk/term-dates
## Holiday dates

school_days_weighted <- setup_holiday_tibble(start_date,end_date, half_term_start,half_term_end,
                                             winter_holiday_start,winter_holiday_end, smooth_time = 7)
prop_immune <- c(prop_immune_younger,prop_immune_younger2, rep(prop_immune_older,5), prop_immune_oldest,prop_immune_oldest)

## Setup population proportions
N_props <- polymod_c_term$participants$proportion
N_age_classes <- length(N_props)
N_immunity_classes <- 2

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)
beta_scales <- rep(1, N_age_classes)
N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
beta_par <- get_beta(C_term,polymod_c_term$participants$proportion,gamma, R0)

beta_scales <- rep(1,N_age_classes)
C_use_holiday <- setup_C(C_holidays, N, beta_scales)
C_use_term <- setup_C(C_term, N, beta_scales)

C_list <- NULL
## Create list of interpolated contact matrices
for(i in 1:nrow(school_days_weighted)){
  C_list[[i]] <- C_use_term*school_days_weighted$weight[i] + C_use_holiday*(1 - school_days_weighted$weight[i])
}