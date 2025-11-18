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

inc_dat <- read_csv("data/final/ili_plus_datasets_by_age.csv")
inc_dat <- read_csv("data/final/flunet_h3_cases_historic.csv")
inc_dat_subset <- inc_dat %>% filter(group != "All",Season == "2023 to 2024") %>% select(date,group, ILIplus) %>% filter(date >= start_date, date <=end_date)
ggplot(inc_dat_subset) + geom_line(aes(x=date,y=ILIplus,col=group))
inc_dat_wide <- inc_dat_subset %>% pivot_wider(names_from=group, values_from=ILIplus) %>% select(-date)
########################################################
## SETUP MODEL FOR 2023/24 SEASON
########################################################
R0 <- 1.8
gamma <- 3

## https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(21)00141-8/fulltext
symp_frac <- rep(c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.75,0.75),each=2)
symp_frac_dat <- data.frame(age.group=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
                                        "[55,65)", "[65,75)", "75+"),
                            symp_frac=c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.75,0.75))

## From https://www.kennedycollegeoxford.org.uk/term-dates
## Holiday dates
half_term_start <- as.Date("2023-10-21")
half_term_end <- as.Date("2023-10-30")
winter_holiday_start <- as.Date("2023-12-15")
winter_holiday_end <- as.Date("2024-01-04")
school_days_weighted <- setup_holiday_tibble(start_date,end_date, half_term_start,half_term_end,
                                             winter_holiday_start,winter_holiday_end, smooth_time = 7)

## Setup population proportions
N_props <- polymod_c_term$participants$proportion
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
prop_immune_younger <- 0.5
prop_immune_younger2 <- 0.5
prop_immune_older <- 0.7
prop_immune <- c(prop_immune_younger,prop_immune_younger2, rep(prop_immune_older,7))
beta_scales <- rep(1, N_age_classes)

## First index is susceptible, second index is immune
alphas <- c(1.2,1)

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
beta_scales <- rep(1,N_age_classes)
beta_par <- get_beta(C_term,polymod_c_term$participants$proportion,gamma, R0)

C_use_holiday <- setup_C(C_holidays, N, beta_scales)
C_use_term <- setup_C(C_term, N, beta_scales)


C_list <- NULL
#for(i in seq_along(ts)){
#  C_list[[i]] <- C_use_term
#  if(add_holiday){
#  if(ts[i] > (holiday_start - smooth_time) & ts[i] < (holiday_start + smooth_time)){
#    C_list[[i]] <- C_use_holiday*(1/(2*smooth_time))*(ts[i] - (holiday_start - smooth_time)) + 
#      C_use_term*(1 - (1/(2*smooth_time))*(ts[i] - (holiday_start - smooth_time)))
#  } else if(ts[i] >= holiday_start){
#    C_list[[i]] <- C_use_holiday
#  }
#  }
#}
## Create list of interpolated contact matrices
for(i in 1:nrow(school_days_weighted)){
  C_list[[i]] <- C_use_term*school_days_weighted$weight[i] + C_use_holiday*(1 - school_days_weighted$weight[i])
}

## Trial run
y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                       alphas=alphas, age_seed=4,immunity_seed=1,seed_size=200,return_compartments=TRUE)
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

inc <- inc %>% left_join(symp_frac_dat)

date_key <- data.frame(time=ts,date=as.Date("2023-09-01") + ts-1)
reporting_rate <- 0.01
ggplot(inc %>% ungroup() %>% 
         #mutate(value = value/max(value,na.rm=TRUE)) %>%
         mutate(value = value*symp_frac*reporting_rate) %>%
         left_join(date_key)) + geom_line(aes(x=date,y=value,col=age.group)) + scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month")
#ggplot(inc %>% filter(t < 200)) + geom_line(aes(x=t,y=gr,col=age.group)) + scale_color_viridis_d() 

inc %>% ungroup() %>% 
  #mutate(value = value/max(value,na.rm=TRUE)) %>%
  mutate(value = value*symp_frac*reporting_rate) %>% group_by(age.group) %>% summarize(total_inf=sum(value,na.rm=TRUE))

print(total_inf/sum(N))
break

#y_base <- epi_ode_size(C_use, beta_par, gamma, N, ts=seq(1,365,by=1),
 #                      alphas=alphas, age_seed=4,immunity_seed=1,return_compartments=TRUE)

parTab <- data.frame(values=c(0.7,0.7,0.7,0.7,0.5,0.0005,log(20),3,2),
                     names=c("prop_immune_younger","prop_immune_younger2","prop_immune_older","prop_immune_oldest","alpha1","reporting_rate","seed_size","obs_sd","R0"),
                     fixed=c(0,0,0,0,0,0,0,1,1),
                     lower_bound=c(0,0,0,0,0,0,0,0,1),
                     upper_bound=c(1,1,1,1,1,0.01,log(10000),10,10),
                     steps=c(0.1,0.1,0.1,0.1,0.1,0.0001,0.1,0.1,0.1),
                     stringsAsFactors=FALSE)

my_prior <- function(pars){
  a <- dbeta(pars["prop_immune_younger"],1,5,1)
  b <- dbeta(pars["prop_immune_younger2"],2,1.5,1)
  c <- dbeta(pars["prop_immune_older"],10,6,1)
  d <- dbeta(pars["prop_immune_oldest"],10,6,1)
  e <- dbeta(pars["alpha1"],1,1,1)
  f <- dunif(pars["reporting_rate"], 0, 0.01) #dbeta(pars["reporting_rate"],2,1000,1)
  g <- dnorm(pars["seed_size"],log(50),2,1)
  tot_prior <- a + b + c + d + e + f + g
  return(tot_prior)
}


posterior_func <- function(parTab, data, PRIOR_FUNC, return_dat=FALSE,aggregate=FALSE,...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################
  
  ## Somewhat redundant example
  parameter_names <- parTab$names
  
  cost_function <- function(pars, return_dat=FALSE,aggregate=FALSE){
    #lik <- -10000
    names(pars) <- parameter_names
    #if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    #return(lik)
    
    prop_immune_younger <- pars[1]
    prop_immune_younger2 <- pars[2]
    prop_immune_older <- pars[3]
    prop_immune_oldest <- pars[4]
    prop_immune <- c(prop_immune_younger,prop_immune_younger2, rep(prop_immune_older,5),prop_immune_oldest,prop_immune_oldest)
    alphas <- c(1,pars[5])
    reporting_rate <- pars[6]
    seed_size <- exp(pars[7])
    
    N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
    #N_props_long <- N_props
    N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
    
    y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                           alphas=alphas, age_seed=4,immunity_seed=1,seed_size=seed_size,return_compartments=TRUE)
    ret <- as.matrix(y_base[,use_cols])
    ret <- apply(ret, 2, function(x) c(0, diff(x)))
    
    ret <- t((reporting_rate*symp_frac) * t(ret))
    # sum every 2 columns
    ret2 <- ret[, seq(1, ncol(ret), 2)] + ret[, seq(2, ncol(ret), 2)]
    
    # then sum every 7 rows
    row_groups <- gl(nrow(ret2) %/% 7 + (nrow(ret2) %% 7 > 0), 7, nrow(ret2))
    ret_sum <- aggregate(ret2, by = list(row_groups), FUN = sum)[, -1]
    
    ## Final size 
    y_base <- y_base %>% pivot_longer(-time)
    y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
                                age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
                                immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
    total_inf <- y_base %>% filter(time == max(time),compartment == "inc")  %>% pull(value) %>% sum()
    
    final_size <- total_inf/sum(N)
    
    if(return_dat){
      if(aggregate){
        ret_sum <- aggregate_to_groups(ret_sum)
      }
      print(final_size)
      return(ret_sum)
    } else {
      ret_test <- aggregate_to_groups(ret_sum)
      lik <- sum(dnorm(c(as.matrix(inc_dat_wide)),c(as.matrix(ret_test)),sd=pars[8],log=TRUE)) + dnorm(final_size,0.15,0.01,log=TRUE)
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      
      return(lik)
    }
  }
  return(cost_function)
}
model_func <- posterior_func(parTab,NULL,NULL,TRUE,TRUE)
model_func(parTab$values)
mcmcPars <- c("iterations"=500,"popt"=0.44,"opt_freq"=100,
              "thin"=1,"adaptive_period"=500,"save_block"=100)


output <- run_MCMC(parTab=parTab, data=NULL, mcmcPars=mcmcPars, filename="test", 
                   CREATE_POSTERIOR_FUNC=posterior_func, mvrPars=NULL, 
                   PRIOR_FUNC = my_prior, OPT_TUNING=0.2)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

bestPars <- get_best_pars(chain)

## Calculate the covariance matrix of the model parameters from the 
## previous MCMC chain. Need to exclude the first and last column
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)

## List of multivariate proposal parameters.
## Specifying 3 things; 1) the nxn covariance matrix of all model parameters (including fixed ones);
## 2) the starting step scaler for the proposals. The equation below should work well as an
## initial guess, but it will get updated as the sampler progresses; 
## 3) the sampler updates the covariance matrix (by opt_freq) but gives some weight to
## the old covariance matrix to preserver markov properties. A weighting of 
## 0.8 means that the new covariance matrix after an adaptive step is 
## `covMat <- (1-w)*oldCovMat + w*newCovMat`
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

startTab <- parTab
startTab$values <- bestPars
output2 <- run_MCMC(startTab, NULL, mcmcPars, filename="test", posterior_func, mvrPars, PRIOR_FUNC = my_prior  ,0.2)
chain <- read.csv(output2$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
bestPars <- get_best_pars(chain)

pred <- as.data.frame(model_func(bestPars,TRUE,TRUE))
pred$date <- date_seq[seq(1,length(date_seq),by=7)]
plot_dat <- inc_dat_wide
plot_dat$date <- date_seq[seq(1,length(date_seq),by=7)]
plot_dat_long <- plot_dat %>% pivot_longer(-date)
plot_dat_long$name <- factor(plot_dat_long$name, levels=c("1-4","5-14", "15-44","45-64","65+"))
pred_long <- pred %>% pivot_longer(-date)
pred_long$name <- factor(pred_long$name, levels=c("1-4","5-14", "15-44","45-64","65+"))

ggplot(pred_long) + 
  geom_line(aes(x=date,y=value,col="Pred")) +
  #geom_line(data=plot_dat_long, aes(x=date,y=value,col="Data")) +
  facet_wrap(~name)


create_prediction_intervals <- function(chain,nsamp=100){
  samps <- sample(1:nrow(chain), nsamp)
  tmp_preds <- NULL
  for(i in seq_along(samps)){
    tmp_pars <- get_index_pars(chain, samps[i])
    tmp_pred <- as.data.frame(model_func(tmp_pars,TRUE,TRUE))
    tmp_pred$date <- date_seq[seq(1,length(date_seq),by=7)]
    tmp_pred <- tmp_pred %>% pivot_longer(-date)
    tmp_pred$name <- factor(tmp_pred$name, levels=c("1-4","5-14", "15-44","45-64","65+"))
    tmp_pred$value <- rnorm(nrow(tmp_pred),tmp_pred$value,tmp_pars["obs_sd"])
    tmp_pred$samp_no <- samps[i]
    tmp_preds[[i]] <- tmp_pred
  }
  tmp_preds <- do.call("bind_rows",tmp_preds)
  tmp_preds
}
pred_intervals <- create_prediction_intervals(chain, 100)
pred_intervals <- pred_intervals %>%
  group_by(date,name) %>%
  summarize(lower=quantile(value,0.025),
            upper=quantile(value,0.975),
            median=quantile(value,0.5))


ggplot(pred_intervals) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.5) +
 # geom_line(data=pred_long,aes(x=date,y=value,col="Pred")) +
  #geom_line(data=plot_dat_long, aes(x=date,y=value,col="Data")) +
  facet_wrap(~name)


## Plot scenarios
bestPars <- get_best_pars(chain)
#bestPars[7] <- log(5000)
pred <- as.data.frame(model_func(bestPars,TRUE,TRUE))
pred$date <- date_seq[seq(1,length(date_seq),by=7)]
pred_long <- pred %>% pivot_longer(-date)
pred_long$name <- factor(pred_long$name, levels=c("1-4","5-14", "15-44","45-64","65+"))
pred_long %>% group_by(name) %>% summarise(total_inf = sum(value/bestPars["reporting_rate"]))
p_scenario <- ggplot(pred_long) + geom_line(aes(x=date,y=value/bestPars["reporting_rate"],col=name)) + scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month") +
  theme_bw()

p_base/p_scenario
