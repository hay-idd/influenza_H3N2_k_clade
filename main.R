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

contacts <- polymod$contacts

## Create new contact survey with no school contacts, decreased work contacts, and increased home contacts
contacts_no_schools <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home==0)# & cnt_leisure == 0)

## Resample to give frac of work contacts
contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(0.4)

## Resample with replacement to give more home contacts
contacts_home <- contacts %>% filter(cnt_home == 1)
contacts_home <- contacts_home[sample(1:nrow(contacts_home), size = nrow(contacts_home)*2, replace=TRUE), ]

contacts_holidays <- bind_rows(contacts_no_schools,contacts_work,contacts_home)

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
#image(C_holidays)


polymod_c_term <- contact_matrix(polymod,
                                     countries="United Kingdom",
                                     age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                     symmetric=TRUE,
                                     missing.contact.age = "sample",
                                     missing.participant.age = "remove")
C_term <- polymod_c_term$matrix
row.names(C_term) <- colnames(C_term)
#image(C_term)

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
prop_immune_younger <- 0.5
prop_immune_older <- 0.7
prop_immune <- c(0,prop_immune_younger,prop_immune_younger, rep(prop_immune_older,6))
beta_scales <- rep(1, N_age_classes)
alphas <- c(1.2,1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

beta_scales <- rep(1,N_age_classes)

C_use_holiday <- setup_C(C_holidays, N, beta_scales)
C_use_term <- setup_C(C_term, N, beta_scales)

beta_par <- get_beta(C_term,polymod_c$participants$proportion,gamma, R0)
C_list <- NULL
ts=seq(1,240,by=1)
holiday_start <- 110
smooth_time <- 7
add_holiday <- TRUE
for(i in seq_along(ts)){
  C_list[[i]] <- C_use_term
  if(add_holiday){
  if(ts[i] > (holiday_start - smooth_time) & ts[i] < (holiday_start + smooth_time)){
    C_list[[i]] <- C_use_holiday*(1/(2*smooth_time))*(ts[i] - (holiday_start - smooth_time)) + 
      C_use_term*(1 - (1/(2*smooth_time))*(ts[i] - (holiday_start - smooth_time)))
  } else if(ts[i] >= holiday_start){
    C_list[[i]] <- C_use_holiday
  }
  }
}

#y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
 #                      alphas=alphas, age_seed=4,immunity_seed=1,return_compartments=TRUE)
library(data.table)
use_cols <- which(colnames(y_base) %like% "inc")

## https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(21)00141-8/fulltext
symp_frac <- rep(c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.55,0.55),each=2)
symp_frac_dat <- data.frame(age.group=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                                                "[55,65)", "[65,75)", "75+"),
                            symp_frac=c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.55,0.55))
cost_function <- function(pars){
  prop_immune_younger <- pars[1]
  prop_immune_older <- pars[2]
  prop_immune <- c(0,prop_immune_younger,prop_immune_younger, rep(prop_immune_older,6))
  alphas <- c(pars[3],1)
  reporting_rate <- pars[4]
  
  y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                         alphas=alphas, age_seed=4,immunity_seed=1,seed_size=200,return_compartments=TRUE)
  return(y_base[,use_cols] * symp_frac * reporting_rate)
}
cost_function(c(0.5,0.7,1.2,0.2))


y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                       alphas=alphas, age_seed=4,immunity_seed=1,seed_size=200,return_compartments=TRUE)

#rt <- y_base %>% select(time,Rt)
#ggplot(rt) + geom_line(aes(x=time,y=Rt))

#y_base$t <- 1:nrow(y_base)
y_base <- y_base %>% pivot_longer(-time)
y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
                            age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
                            immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
total_inf <- y_base %>% filter(time == max(time),compartment == "inc")  %>% pull(value) %>% sum()

inc <- y_base %>% filter(compartment == "inc")
inc <- inc %>% group_by(age, immunity) %>% mutate(value = value - lag(value, 1))



inc <- inc %>% left_join(polymod_c$participants %>% mutate(age = 1:nrow(polymod_c$participants)))
inc$age.group <- factor(inc$age.group,levels=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                               "[55,65)", "[65,75)", "75+"))

inc <- inc %>% group_by(age.group,time) %>% summarize(value = sum(value)) %>% group_by(age.group) %>% mutate(gr = log(value/lag(value,1)))

inc <- inc %>% left_join(symp_frac_dat)

date_key <- data.frame(time=ts,date=as.Date("2025-09-01") + ts-1)
reporting_rate <- 0.01
ggplot(inc %>% ungroup() %>% 
         #mutate(value = value/max(value,na.rm=TRUE)) %>%
         mutate(value = value*symp_frac*reporting_rate) %>%
         left_join(date_key)) + geom_line(aes(x=date,y=value,col=age.group)) + scale_color_viridis_d() +
  geom_vline(xintercept= as.Date("2025-09-01") + holiday_start, linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month")
#ggplot(inc %>% filter(t < 200)) + geom_line(aes(x=t,y=gr,col=age.group)) + scale_color_viridis_d() 

inc %>% ungroup() %>% 
  #mutate(value = value/max(value,na.rm=TRUE)) %>%
  mutate(value = value*symp_frac*reporting_rate) %>% group_by(age.group) %>% summarize(total_inf=sum(value,na.rm=TRUE))

print(total_inf/sum(N))
