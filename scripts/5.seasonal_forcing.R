## Test if there is a relationship between growth rate, cumulative cases, school holidays (+/- 1 week) and aboslute humidity
library(tidyverse)
library(terra)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(marginaleffects)

source("R/funcs.R")

#################################################
## Read in climate data
#################################################
# --- open NetCDF files ---
tfile <- "~/Downloads/tas_hadukgrid_uk_country_mon_188401-202412.nc"   # mean temperature (°C)
pvfile <- "~/Downloads/pv_hadukgrid_uk_country_mon_196101-202412.nc"  # vapour pressure (hPa)

tas <- rast(tfile, subds = "tas")   # try variable name if needed
pv  <- rast(pvfile,  subds = "pv")

country_index <- 8
tas_values <-as.numeric(tas[,country_index,][[1]])
pv_values <- as.numeric(pv[,country_index,][[1]])
tas_values <- tas_values[(length(tas_values)-length(pv_values)+1):length(tas_values)]

ah <- 216.7 * (pv_values / (tas_values + 273.15))
dates <- seq(as.Date("1961-01-01"),as.Date("2024-12-01"), by = "1 month")

comb <- data.frame(ah=ah,date=dates)

## Look at estimated absolute humidity over time
ggplot(comb %>% filter(date >= "2012-01-01")) + geom_line(aes(x=date,y=ah))
comb1 <- comb %>% filter(date >= "2012-01-01")

## Enumerate into daily data
comb1 <- comb1 %>%
  mutate(next_date = lead(date, default = max(date) + 30)) %>%
  rowwise() %>%
  mutate(date = list(seq(date, next_date - 1, by = "day"))) %>%
  unnest(cols = c(date)) %>%
  select(-next_date)

## Use the 2024 data for 2025 as a placeholder
tmp_dat_2024 <- comb1 %>% filter(date >= "2024-01-01") %>% mutate(date = date + 365)
comb1 <- bind_rows(comb1, tmp_dat_2024)

## Smooth daily AH using 30-day rolling average
comb1 <- comb1 %>% mutate(ah_smooth = zoo::rollmean(ah, k = 30, fill = NA, align = "center"))

ggplot(comb1 %>% filter(date >= "2012-01-01")) + geom_line(aes(x=date,y=ah_smooth))

## Created absolute humidity dataset

#################################################
## Read in estimate growth rates
#################################################
## Read in previously calculated growth rates and incidence data, first try WHO FluNet
grs_est <- read_csv("results/flunet_all_growth_rates.csv")
inc_dat <- read_csv("results/fit_gams_all.csv") %>% filter(Source=="FluNet")
school_days_full2 <- school_days_full2 %>% mutate(time_since_preweek = ceiling(time_since_preweek/7))
fit_gam_model <- function(climate_dat, gr_dat, inc_dat, school_holidays, formula = y ~ s(qlogis(prop_of_cumu_cases)) + s(ah_smooth_lag) + time_since_preweek + Season,
                          ah_lag=0){
  ## Get cumulative cases in incidence data
  inc_dat <- inc_dat %>% group_by(Season) %>% mutate(cumu_cases = cumsum(cases_total))
  ## Calculate proportion of cumulative cases for the season and put on -Inf to Inf scale
  inc_dat <- inc_dat %>% group_by(Season) %>% mutate(prop_of_cumu_cases =(cumu_cases+1)/(max(cumu_cases)+2)) %>% mutate(qlogis_prop_cases = qlogis(prop_of_cumu_cases))
  
  ## Add in school holidays
  inc_dat <- inc_dat %>% mutate(date = ymd(date)) %>% left_join(school_holidays %>% mutate(date = ymd(date)) %>% select(date,type,time_since_preweek))
  
  ## Combine climate data, growth rate data and incidence data
  comb_for_regression <- left_join(gr_dat,climate_dat) %>% left_join(inc_dat)%>% mutate(date = ymd(date))
  comb_for_regression1 <- comb_for_regression %>% select(date,day_of_year,gr_smooth,y,cumu_cases,ah_smooth,Season,type,time_since_preweek,prop_of_cumu_cases) %>% drop_na()
  
  ## Optionally add lag on absolute humidity
  comb_for_regression1 <- comb_for_regression1 %>% mutate(ah_smooth_lag = lag(ah_smooth,ah_lag)) %>% drop_na()
  
  ## Remove latest season and first partial season
  predict_data <- comb_for_regression1 %>% filter(Season == "2025/26")
  comb_for_regression1 <- comb_for_regression1 %>% filter(Season != "2025/26") %>% filter(Season != "2011/12")
  comb_for_regression1 <- comb_for_regression1 %>% group_by(Season) %>% mutate(time = as.numeric(date-min(date)))
  #fit <- mgcv::gam(y ~ s(log(cumu_cases + 1)) + s(ah_smooth_lag)  + type,data=comb_for_regression1)
  fit <- gam(formula,data=comb_for_regression1,select=TRUE)#,  correlation = nlme::corAR1(form = ~ time | Season))
  dat <- comb_for_regression1
  dat$pred <- predict(fit)
  p_draw <- draw(fit)
  p_appraise <- appraise(fit)
  p_pred_week <- plot_predictions(fit,condition="time_since_preweek")
  p_pred_season <- plot_predictions(fit,condition="Season")
  ## Compare variance explained
  check_var_explained <- function(fit, vars=c("s(qlogis(prop_of_cumu_cases))",
                                              "s(ah_smooth_lag)",
                                              "s(time_since_preweek)",
                                              "Season")){
    full_dev <- summary(fit)$dev.expl
    var_contrib <- sapply(vars, function(v) {
      fml <- update(formula(fit), paste(". ~ . -", v))
      fit_red <- update(fit, formula = fml)
      full_dev - summary(fit_red)$dev.expl
    })
    
    var_contrib / full_dev
  }
  var_explained_by_variable <- check_var_explained(fit)
  
  
  list(model_fit=fit,data=comb_for_regression1,var_explained_by_variable=var_explained_by_variable,predictions=dat,
       p_draw=p_draw,p_appraise=p_appraise,p_pred_week=p_pred_week,p_pred_season=p_pred_season,
       data_2025=predict_data)
}

## Could I use the difference between the predicted growth rate and the true growth rate, and say, if the difference is ascribed to under-predicted cumulative proportion infected, how much lower would the cumu proportion infected need to be to explain the difference?

fit <- fit_gam_model(climate_dat=comb1, gr_dat=grs_est, inc_dat=inc_dat, school_holidays=school_days_full2,
                            formula = y ~ s(qlogis(prop_of_cumu_cases)) + s(ah_smooth_lag) + time_since_preweek + Season)

generate_all_outputs <- function(fit){
  
  ## Get the change in growth rate prediction if AH is changed from Nov 2025 to Dec 2025 and school term set to 0
  ## Predict values for this year, assuming same as 2023/24 for total cases, AH and intercept
  pred_2023_to_2024 <- fit$data %>% filter(Season == "2023/24")#fit_flunet$data_2025 
  ## Get total cases from 2023/24
  total_cases_2023 <- fit$data %>% filter(Season == "2023/24") %>% filter(cumu_cases==max(cumu_cases)) %>% pull(cumu_cases)
  pred_2023_to_2024$Season <- "2023/24"
  add_confint <- function(dat, fit){
    dat$pred <- predict(fit,newdata=dat)
    dat$se <- predict(fit,newdata=dat,se=TRUE)$se.fit
    dat$lower <- dat$pred - 1.96*dat$se
    dat$upper <- dat$pred + 1.96*dat$se
    dat
  }
  pred_2023_to_2024 <- add_confint(pred_2023_to_2024, fit$model_fit)
  ## Prediction vs. real for 2023/24
  p_pred_2023_to_2024 <- ggplot(pred_2023_to_2024) + 
    geom_hline(yintercept=0,linetype="dashed",col="grey") +
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill="Predicted"),alpha=0.2) +
    geom_line(aes(x=date,y=pred,col="Predicted")) +
    geom_line(aes(x=date,y=y,col="Observed",fill="Observed")) +
    scale_fill_brewer("",palette="Set1") +
    scale_color_brewer("",palette="Set1") +
    theme_bw() +
    xlab("Date") +
    ylab("Growth rate (weekly)") +
    coord_cartesian(ylim=c(-0.75,0.75)) +
    ggtitle("Comparison of predicted and observed growth rates in 2023/24")
  
  ## Prediction vs. real for 2025/26
  pred_2025_to_2026 <- fit$data_2025 
  ## Get total cases from 2023/24
  total_cases_2023 <- fit$data %>% filter(Season == "2023/24") %>% filter(cumu_cases==max(cumu_cases)) %>% pull(cumu_cases)
  total_cases_2023 <- total_cases_2023
  
  pred_2025_to_2026$prop_of_cumu_cases <- pred_2025_to_2026$cumu_cases / (total_cases_2023 + 1)
  pred_2025_to_2026$Season <- "2023/24"
  pred_2025_to_2026 <- add_confint(pred_2025_to_2026, fit$model_fit)
  
  p_pred_2025_to_2026 <- ggplot(pred_2025_to_2026) + 
    geom_hline(yintercept=0,linetype="dashed",col="grey") +
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill="Predicted"),alpha=0.2) +
    geom_line(aes(x=date,y=pred,col="Predicted")) +
    geom_line(aes(x=date,y=y,col="Observed",fill="Observed")) +
    scale_fill_brewer("",palette="Set1") +
    scale_color_brewer("",palette="Set1") +
    theme_bw() +
    xlab("Date") +
    ylab("Growth rate (weekly)") +
    coord_cartesian(ylim=c(-0.75,0.75)) +
    ggtitle("Comparison of predicted and observed growth rates in 2025/26")
  
  
  
  ## What should AH be in December?
  ah_smooth_later <- pred_2025_to_2026 %>% mutate(date = date + 30) %>% select(date) %>% left_join(comb1 %>% select(date,ah_smooth) %>% rename(ah_smooth_lag=ah_smooth))
  pred_2025_to_2026_later <- pred_2025_to_2026
  pred_2025_to_2026_later$ah_smooth_lag <- ah_smooth_later$ah_smooth_lag
  pred_2025_to_2026_later <- add_confint(pred_2025_to_2026_later, fit$model_fit)
  pred_2025_to_2026_later$Scenario <- "Outbreak started 1 month later"
  pred_2025_to_2026$Scenario <- "Base prediction"
  
  ## Fewer total cases
  pred_2025_to_2026_fewer <- pred_2025_to_2026
  pred_2025_to_2026_fewer$prop_of_cumu_cases <- pred_2025_to_2026_fewer$cumu_cases / (total_cases_2023*0.5 + 1)
  pred_2025_to_2026_fewer <- add_confint(pred_2025_to_2026_fewer, fit$model_fit)
  pred_2025_to_2026_fewer$Scenario <- "Half total cases of 2023/24"
  
  ## More total cases
  pred_2025_to_2026_more <- pred_2025_to_2026
  pred_2025_to_2026_more$prop_of_cumu_cases <- pred_2025_to_2026_more$cumu_cases / (total_cases_2023*2 + 1)
  pred_2025_to_2026_more <- add_confint(pred_2025_to_2026_more, fit$model_fit)
  pred_2025_to_2026_more$Scenario <- "Double total cases of 2023/24"
  
  pred_2025_to_2026_scenarios <- bind_rows(pred_2025_to_2026, pred_2025_to_2026_later, pred_2025_to_2026_fewer, pred_2025_to_2026_more)
  
  p_pred_2025_to_2026_scenarios <- ggplot(pred_2025_to_2026_scenarios) + 
    geom_hline(yintercept=0,linetype="dashed",col="grey") +
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=Scenario),alpha=0.2) +
    geom_line(aes(x=date,y=pred,col=Scenario)) +
    geom_line(aes(x=date,y=y,col="Observed",fill="Observed")) +
    scale_fill_brewer("",palette="Set1") +
    scale_color_brewer("",palette="Set1") +
    theme_bw() +
    xlab("Date") +
    ylab("Growth rate (weekly)") +
    coord_cartesian(ylim=c(-0.75,0.75)) +
    ggtitle("Comparison of predicted and observed growth rates in 2025/26")
  
  return(list(
  p_pred_2023_to_2024/p_pred_2025_to_2026,
  p_pred_2025_to_2026_scenarios,
  summary(fit$model_fit),
  fit$var_explained_by_variable,
  fit$p_draw,
  fit$p_appraise,
  fit$p_pred_season +  theme_bw() + xlab("Season") + ylab("Intercept (additive on growth rate)"),
  fit$p_pred_week + theme_bw() + coord_cartesian(xlim=c(0,6),ylim=c(-0.2,0.2)) + ylab("Additive effect on growth rate") +
    xlab("Weeks since 1 week prior to start of school holiday \n (i.e., 0 = week before school holiday)")
  ))
}
all_outputs <- generate_all_outputs(fit)


ggsave("figures/climate_models/flunet_predicted_vs_observed_growth_rates_2023_2024.png", all_outputs[[1]], width=8, height=6)
ggsave("figures/climate_models/flunet_predicted_vs_observed_growth_rates_2025_2026.png", all_outputs[[2]], width=8, height=4)
ggsave("figures/climate_models/spline_terms_model_flunet.png", all_outputs[[5]], width=8, height=4)
ggsave("figures/climate_models/diagnostics_model_flunet.png", all_outputs[[6]], width=8, height=7)
ggsave("figures/climate_models/season_terms_model_flunet.png", all_outputs[[7]], width=8, height=3)
ggsave("figures/climate_models/holiday_terms_model_flunet.png", all_outputs[[8]], width=5, height=3)

## Save all_outputs[[3]] to text file
sink(file="figures/climate_models/gam_model_summary_flunet.txt")
print(all_outputs[[3]])
print(all_outputs[[4]])
sink(NULL)

grs_est <- read_csv("results/ukhsa_all_flu_growth_rates.csv")%>% mutate(day_of_year = yday(date))
inc_dat <- read_csv("results/fit_gams_all.csv") %>% filter(Source=="UKHSA") %>% mutate(day_of_year = yday(date))
fit <- fit_gam_model(climate_dat=comb1, gr_dat=grs_est, inc_dat=inc_dat, school_holidays=school_days_full2,                                                                                                 formula = y ~ s(qlogis(prop_of_cumu_cases)) + s(ah_smooth_lag) + time_since_preweek + Season)
all_outputs <- generate_all_outputs(fit)

ggsave("figures/climate_models/ukhsa_predicted_vs_observed_growth_rates_2023_2024.png", all_outputs[[1]], width=8, height=6)
ggsave("figures/climate_models/ukhsa_predicted_vs_observed_growth_rates_2025_2026.png", all_outputs[[2]], width=8, height=4)
ggsave("figures/climate_models/spline_terms_model_ukhsa.png", all_outputs[[5]], width=8, height=4)
ggsave("figures/climate_models/diagnostics_model_ukhsa.png", all_outputs[[6]], width=8, height=7)
ggsave("figures/climate_models/season_terms_model_ukhsa.png", all_outputs[[7]], width=8, height=3)
ggsave("figures/climate_models/holiday_terms_model_ukhsa.png", all_outputs[[8]], width=5, height=3)

## Save all_outputs[[3]] to text file
sink(file="figures/climate_models/gam_model_summary_ukhsa.txt")
print(all_outputs[[3]])
print(all_outputs[[4]])
sink(NULL)

## Save all of these outputs to file
