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


start_date <- as.Date("2022-09-01")
end_date <- as.Date("2023-05-01")

half_term_start <- as.Date("2022-10-21")
half_term_end <- as.Date("2022-10-30")
winter_holiday_start <- as.Date("2022-12-15")
winter_holiday_end <- as.Date("2023-01-04")


## Proportion of work contacts kept during the school holidays
prop_work_contacts_in_hols <- 0.5
## Increase in at-home contacts during the holidays
prop_home_contacts_in_hols <- 1.5

parTab <- data.frame(values=c(0.7,0.7,0.8,0.7,
                              0.5,0.001,
                              log(20),3,2,
                              0.5,0.4,0.5,0.6),
                     names=c("prop_immune_youngest","prop_immune_younger","prop_immune_older","prop_immune_oldest",
                             "alpha1","reporting_rate",
                             "seed_size","obs_sd","R0",
                             "symp1","symp2","symp3","symp4"),
                     fixed=c(0,0,1,0,
                             0,0,
                             0,0,0, 
                             0,0,1,0),
                     lower_bound=c(0,0,0,0,
                                   0,0,0,
                                   0,1,
                                   0,0,0,0),
                     upper_bound=c(1,1,1,1,
                                   1,0.01,log(10000),
                                   50,10,
                                   1,1,1,1),
                     steps=c(0.1,0.1,0.1,0.1,
                             0.1,0.0001,0.1,
                             0.1,0.1,
                             0.1,0.1,0.1,0.1),
                     stringsAsFactors=FALSE)

symp_frac <- c(0.5,0.4,0.5,0.6)

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
age_breaks <- c(0,5,18,65)
polymod1 <- polymod
polymod1$contacts <- contacts_holidays
polymod_c_holidays <- contact_matrix(polymod1,
                                     countries="United Kingdom",
                                     #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                     age.limits=age_breaks,
                                     symmetric=TRUE,
                                     missing.contact.age = "sample",
                                     missing.participant.age = "remove")
C_holidays <- polymod_c_holidays$matrix
row.names(C_holidays) <- colnames(C_holidays)

## Create term time contact matrix
polymod_c_term <- contact_matrix(polymod,
                                 countries="United Kingdom",
                                 #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                 age.limits=age_breaks,
                                 symmetric=TRUE,
                                 missing.contact.age = "sample",
                                 missing.participant.age = "remove")
C_term <- polymod_c_term$matrix
row.names(C_term) <- colnames(C_term)

########################################################
## READ IN INCIDENCE DATA FOR FITTING
########################################################
flu_dat <- read_csv("data/final/flu_2022_2023.csv")
start_day <- as.numeric(start_date)
end_day <- as.numeric(end_date)
date_seq <- seq(start_date,end_date,by="1 day")
ts <- seq(start_day,end_day,by=1) - start_day + 1
flu_dat <-  flu_dat %>% filter(date >= start_date,date<=end_date)
flu_dat_wide <- flu_dat %>% filter(date >= start_date,date<=end_date) %>%
  select(date, age_group,ILI_flu) %>% pivot_wider(names_from=age_group,values_from=ILI_flu)
  
########################################################
## SETUP MODEL FOR 2023/24 SEASON
########################################################

## From https://www.kennedycollegeoxford.org.uk/term-dates
## Holiday dates
## Setup contact matrix list
school_days_weighted <- setup_holiday_tibble(start_date,end_date, half_term_start,half_term_end,
                                             winter_holiday_start,winter_holiday_end, smooth_time = 7)
prop_immune <- c(prop_immune_younger,prop_immune_younger2,prop_immune_older,prop_immune_oldest)
## Setup population proportions
N_props <- polymod_c_term$demography$proportion
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


f <- posterior_func(parTab, data=flu_dat_wide, PRIOR_FUNC=my_prior,
                    #symp_frac=symp_frac, 
                    N_props=N_props,C_term=C_term, C_list=C_list,return_dat=FALSE)

f(parTab$values)
model_f <- posterior_func(parTab, data=flu_dat_wide, PRIOR_FUNC=NULL,
                          #symp_frac=symp_frac, 
                          N_props=N_props,C_term=C_term, C_list=C_list,return_dat=TRUE)
tmp_dat <- model_f(parTab$values)
tmp_dat$t <- 1:nrow(tmp_dat)
p1 <- tmp_dat %>% pivot_longer(-t) %>% ggplot() + geom_line(aes(x=t,y=value,col=name))
p2 <- flu_dat %>% ggplot() + geom_line(aes(x=floor(as.numeric(date-min(date))/7),y=ILI_flu,col=age_group))
p1/p2

##### RUN MCMC UNIVARIATE CHAIN
mcmcPars <- c("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=50000,"save_block"=1000)


output <- run_MCMC(parTab, data=flu_dat_wide,
                  #symp_frac=symp_frac, 
                   N_props=N_props,C_term=C_term, C_list=C_list,return_dat=FALSE,
                   mcmcPars=mcmcPars, filename="test", 
                   CREATE_POSTERIOR_FUNC=posterior_func, mvrPars=NULL, 
                   PRIOR_FUNC = my_prior, OPT_TUNING=0.2)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))

bestPars <- get_best_pars(chain)
## Calculate the covariance matrix of the model parameters from the 
## previous MCMC chain. Need to exclude the first and last column
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.2, w=0.8)#2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

startTab <- parTab
startTab$values <- bestPars

mcmcPars <- c("iterations"=5000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=5000,"save_block"=1000)
output2 <- run_MCMC(parTab, data=flu_dat_wide,
         #symp_frac=symp_frac, 
         N_props=N_props,C_term=C_term, C_list=C_list,return_dat=FALSE,
         mcmcPars=mcmcPars, filename="test", 
         CREATE_POSTERIOR_FUNC=posterior_func, mvrPars=mvrPars, 
         PRIOR_FUNC = my_prior, OPT_TUNING=0.2)
chain <- read.csv(output2$file)

create_prediction_intervals <- function(chain,model_f,nsamp=100){
  samps <- sample(1:nrow(chain), nsamp)
  tmp_preds <- NULL
  for(i in seq_along(samps)){
    tmp_pars <- get_index_pars(chain, samps[i])
    tmp_pred <- as.data.frame(model_f(tmp_pars))
    tmp_pred$date <- date_seq[seq(1,length(date_seq),by=7)]
    tmp_pred <- tmp_pred %>% pivot_longer(-date)
    #tmp_pred$name <- factor(tmp_pred$name, levels=c("1-4","5-14", "15-44","45-64","65+"))
    tmp_pred$value <- rpois(nrow(tmp_pred),tmp_pred$value)
    tmp_pred$samp_no <- samps[i]
    tmp_preds[[i]] <- tmp_pred
  }
  tmp_preds <- do.call("bind_rows",tmp_preds)
  tmp_preds
}
pred_intervals <- create_prediction_intervals(chain,model_f, 100)
pred_intervals <- pred_intervals %>%
  group_by(date,name) %>%
  summarize(lower=quantile(value,0.025),
            upper=quantile(value,0.975),
            median=quantile(value,0.5))

age_key <- c("inc_1_1"="0-4",
             "inc_2_1"="5-18",
             "inc_3_1"="19-64",
             "inc_4_1"="65+")
pred_intervals$age_group <- age_key[pred_intervals$name]

ggplot(pred_intervals) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.5) +
  geom_line(data=flu_dat,aes(x=date,y=ILI_flu)) +
  # geom_line(data=pred_long,aes(x=date,y=value,col="Pred")) +
  #geom_line(data=plot_dat_long, aes(x=date,y=value,col="Data")) +
  facet_wrap(~age_group)


p <- ggplot(inc %>% ungroup() %>% 
              #mutate(value = value/max(value,na.rm=TRUE)) %>%
              mutate(value = value*symp_frac*reporting_rate) %>%
              
              left_join(date_key)) + geom_line(aes(x=date,y=value,col=age.group)) + scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month") +
  ggtitle(paste0("Final size: ",signif(percent_infected,3)))+
  theme_bw() +
  xlab("Date") +
  ylab("Reported symptomatic cases (daily)")







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



inc <- inc %>% left_join(polymod_c_term$participants %>% mutate(age = 1:nrow(polymod_c_term$participants)))
#inc$age.group <- factor(inc$age.group,levels=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", 
#                                               "[55,65)", "[65,75)", "75+"))
inc$age.group <- factor(inc$age.group,levels=c("[0,5)","[5,18)","[18,65)","65+"))

inc <- inc %>% group_by(age.group,time) %>% summarize(value = sum(value)) %>% group_by(age.group) %>% mutate(gr = log(value/lag(value,1)))

symp_frac_rep <- rep(symp_frac,each=2)
symp_frac_dat <- data.frame(#age.group=c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)","[55,65)", "[65,75)", "75+"),
  age.group = c("[0,5)","[5,18)","[18,65)","65+"),
  symp_frac=symp_frac)


inc <- inc %>% left_join(symp_frac_dat)

date_key <- data.frame(time=ts,date=as.Date("2022-09-01") + ts-1)

percent_infected <- total_inf/sum(N)
inc$age.group <- factor(inc$age.group,levels=c("[0,5)","[5,18)","[18,65)","65+"))


p <- ggplot(inc %>% ungroup() %>% 
              #mutate(value = value/max(value,na.rm=TRUE)) %>%
              mutate(value = value*symp_frac*reporting_rate) %>%
              
              left_join(date_key)) + geom_line(aes(x=date,y=value,col=age.group)) + scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month") +
  ggtitle(paste0("Final size: ",signif(percent_infected,3)))+
  theme_bw() +
  xlab("Date") +
  ylab("Reported symptomatic cases (daily)")

flu_dat <- flu_dat %>% mutate(age_group = case_when(
  age_group == "0-4" ~ "[0,5)",
  age_group == "5-18" ~ "[5,18)",
  age_group == "19-64" ~ "[18,65)",
  age_group == "65+" ~ "65+"
))
flu_dat$age_group <- factor(flu_dat$age_group, levels=c("[0,5)","[5,18)","[18,65)","65+"))
p_dat <-  ggplot(flu_dat %>% filter(date >= start_date,date <= end_date),aes(x=date,y=ILI_flu,col=age_group)) + 
  geom_line() +
  scale_color_viridis_d() +
  geom_vline(xintercept= c(half_term_start,half_term_end,winter_holiday_start,winter_holiday_end), linetype="dashed", color = "red") +
  scale_x_date(breaks="1 month") +
  theme_bw() +
  xlab("Date") +
  ylab("Reported symptomatic cases (daily)")

print(total_inf/sum(N))
print(p/p_dat)
