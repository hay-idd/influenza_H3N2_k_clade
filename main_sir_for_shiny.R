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

start_date <- as.Date("2023-09-01")
end_date <- as.Date("2024-05-01")

half_term_start <- as.Date("2023-10-21")
half_term_end <- as.Date("2023-10-30")
winter_holiday_start <- as.Date("2023-12-15")
winter_holiday_end <- as.Date("2024-01-04")

## Proportion of work contacts kept during the school holidays
prop_work_contacts_in_hols <- 0.4
## Increase in at-home contacts during the holidays
prop_home_contacts_in_hols <- 2

## Proportion of 0-4 year olds in the immune category
prop_immune_younger <- 0.5
## Proportion of 5-14 year olds in the immune category
prop_immune_younger2 <- 0.5
## Proportion of 15-64 year olds in the immune category
prop_immune_older <- 0.7
## Proportion of 65+ year olds in the immune category
prop_immune_oldest <- 0.7

## Relative susceptability of the two immune classes
alphas <- c(1.2,1)

## Basic reproductive number
R0 <- 1.2
## Infectious period mean
gamma <- 3
## Number of seed viruses as of 2023-09-01
seed_size <- 200

## Proportion of infections leading to symptoms by age group
## https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(21)00141-8/fulltext
symp_frac <- c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.75,0.75)

## Proportion of symptomatic cases reported
reporting_rate <- 0.01


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

start_day <- as.numeric(start_date)
end_day <- as.numeric(end_date)
date_seq <- seq(start_date,end_date,by="1 day")
ts <- seq(start_day,end_day,by=1) - start_day + 1

########################################################
## SETUP MODEL FOR 2023/24 SEASON
########################################################

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

## Trial run
y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                       alphas=alphas, age_seed=4,immunity_seed=1,seed_size=seed_size,return_compartments=TRUE)
use_cols <- which(colnames(y_base) %like% "inc")

#rt <- y_base %>% select(time,Rt)
#ggplot(rt) + geom_line(aes(x=time,y=Rt))

#y_base$t <- 1:nrow(y_base)
y_base <- y_base %>% pivot_longer(-time)
y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
                            age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
                            immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
total_inf <- y_base %>% filter(time == max(time),compartment == "inc")  %>% pull(value) %>% sum()

print(total_inf/sum(N))

inc <- y_base %>% filter(compartment == "inc")
inc <- inc %>% group_by(age, immunity) %>% mutate(value = value - lag(value, 1))



inc <- inc %>% left_join(polymod_c_term$participants %>% mutate(age = 1:nrow(polymod_c_term$participants)))
inc$age.group <- factor(inc$age.group,levels=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                               "[55,65)", "[65,75)", "75+"))

inc <- inc %>% group_by(age.group,time) %>% summarize(value = sum(value)) %>% group_by(age.group) %>% mutate(gr = log(value/lag(value,1)))

symp_frac_rep <- rep(symp_frac,each=2)
symp_frac_dat <- data.frame(age.group=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                        "[55,65)", "[65,75)", "75+"),
                            symp_frac=symp_frac)


inc <- inc %>% left_join(symp_frac_dat)

date_key <- data.frame(time=ts,date=as.Date("2023-09-01") + ts-1)


ggplot(inc %>% ungroup() %>% 
         #mutate(value = value/max(value,na.rm=TRUE)) %>%
         mutate(value = value*symp_frac*reporting_rate) %>%
         
         left_join(date_key)) + geom_line(aes(x=date,y=value,col=age.group)) + scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month") +
  theme_bw() +
  xlab("Date") +
  ylab("Reported symptomatic cases (daily)")

print(total_inf/sum(N))
