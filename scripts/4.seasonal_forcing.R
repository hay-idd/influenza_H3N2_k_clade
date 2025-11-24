## Test if there is a relationship between growth rate, cumulative cases, school holidays (+/- 1 week) and aboslute humidity
library(tidyverse)
library(terra)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(marginaleffects)

source("R/funcs.R")
source("R/school_holidays.R")

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

comb1$year <- year(comb1$date)
comb1$day <- yday(comb1$date)

## Created absolute humidity dataset

#################################################
## Read in estimate growth rates
#################################################
## Read in previously calculated growth rates and incidence data, first try WHO FluNet
grs_est <- read_csv("results/flunet_all_growth_rates.csv")
inc_dat <- read_csv("results/fit_gams_all.csv") %>% filter(Source=="FluNet")
school_days_full2 <- school_days_full2 %>% mutate(time_since_preweek = ceiling(time_since_preweek/7))
fit_gam_model <- function(climate_dat, gr_dat, inc_dat, school_holidays, formula = y ~ s(prop_of_cumu_cases) + s(ah_smooth_lag) + time_since_preweek + Season,
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
  comb_for_regression1$prop_of_cumu_cases_mod <- qlogis(comb_for_regression1$prop_of_cumu_cases)
  #fit <- mgcv::gam(y ~ s(log(cumu_cases + 1)) + s(ah_smooth_lag)  + type,data=comb_for_regression1)
  fit <- mgcv::gam(formula,data=comb_for_regression1,select=TRUE)#,  correlation = nlme::corAR1(form = ~ time | Season))
  p_pred_week <- plot_predictions(fit,newdata=comb_for_regression1, condition="time_since_preweek")
  p_pred_season <- plot_predictions(fit,condition="Season")
  
  dat <- comb_for_regression1
  dat$pred <- predict(fit)
  p_draw <- draw(fit)
  p_appraise <- appraise(fit)
  ## Compare variance explained
  check_var_explained <- function(fit, vars=c("s(prop_of_cumu_cases)",
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
  
  ## Get prediction for the 2025/26 season
  ## If cumulative cases is very similar stage to the 2022/23 season
  ## Most recent growth rate -- what % increase over the following month?
  ## Equivalent would be Nov 11th 27
  ## Build new predictor dataset

  
  ## Get data to base our alternative realities on
  new_pred_data <- comb_for_regression1 %>% 
    select(date, day_of_year, ah_smooth_lag, time_since_preweek, Season, prop_of_cumu_cases,cumu_cases)
  
  ## Need to get the same dates onward
  new_pred_data <- new_pred_data %>% filter(day_of_year > yday(as.Date("2025-11-01")))
  
  prop_of_cumu_cases <- new_pred_data %>% filter(Season == "2022/23") %>% mutate(diff_cases = abs(cumu_cases - max(predict_data$cumu_cases))) %>% filter(diff_cases == min(diff_cases)) %>% pull(prop_of_cumu_cases)
  
  
  new_pred_data$prop_of_cumu_cases <-prop_of_cumu_cases
  new_pred_data$pred <- predict(fit,newdata=new_pred_data %>% mutate(Season = "2022/23"))
  
  ## Get difference to growth rate on first day
  new_pred_data_subset <- new_pred_data %>% mutate(days_elapsed = day_of_year - min(day_of_year)) %>% filter(days_elapsed <= 30)
  
  base_pred <- new_pred_data_subset %>% filter(days_elapsed==0) %>% select(Season,pred) %>% rename(base_pred=pred)
  new_pred_data_subset %>% left_join(base_pred) %>% mutate(pred1 = pred - base_pred) %>% ggplot() +
    geom_line(aes(x=days_elapsed,y=pred1,group=Season,col=Season)) +
    geom_line(data = .%>% group_by(days_elapsed) %>% summarize(mean_pred=mean(pred1)), aes(x=days_elapsed,y=mean_pred),col="red",linewidth=1)
  
  list(model_fit=fit,data=comb_for_regression1,var_explained_by_variable=var_explained_by_variable,predictions=dat,
       p_draw=p_draw,p_appraise=p_appraise,p_pred_week=p_pred_week,p_pred_season=p_pred_season,
       data_2025=predict_data)
}

## Could I use the difference between the predicted growth rate and the true growth rate, and say, if the difference is ascribed to under-predicted cumulative proportion infected, how much lower would the cumu proportion infected need to be to explain the difference?

fit <- fit_gam_model(climate_dat=comb1, gr_dat=grs_est, inc_dat=inc_dat, school_holidays=school_days_full2,
                            formula = y ~ s(ah_smooth_lag) + time_since_preweek + Season)

all_outputs <- generate_all_outputs(fit)
p_tmp <- all_outputs[[5]]
fit_complex <- fit_gam_model(climate_dat=comb1, gr_dat=grs_est, inc_dat=inc_dat, school_holidays=school_days_full2,formula = y ~ s(ah_smooth_lag) + s(prop_of_cumu_cases) + time_since_preweek + Season)

all_outputs_complex <- generate_all_outputs(fit_complex)
## Plot contribution of climate to predicted growth rate
preds <- fit$predictions
preds$yday1 <- yday(preds$date)
preds$year <- year(preds$date)

preds_complex <- fit_complex$predictions
preds_complex$yday1 <- yday(preds_complex$date)
preds_complex$year <- year(preds_complex$date)

# Subset the dates you want to label
vlines_df <- comb1 %>%
  filter(date %in% as.Date(c("2024-10-01","2024-11-01","2024-12-01"))) %>%
  mutate(day_of_year = yday(date)) %>%
  select(day_of_year, date, ah_smooth)

p_cor <- ggplot(preds) +
  geom_point(aes(x = ah_smooth_lag, y = y),size=0.75,alpha=0.5) +
  geom_smooth(aes(x = ah_smooth_lag, y = y, fill="LOESS", col="LOESS")) +
  geom_smooth(aes(x = ah_smooth_lag, y = y, fill="LM", col="LM"), method="lm") +
  
  # dashed vertical lines
  geom_vline(data = vlines_df, aes(xintercept = ah_smooth), linetype = "dashed") +
  
  # ---- NEW: date labels next to dashed lines ----
geom_text(
  data = vlines_df,
  aes(x = ah_smooth - 0.2, y = -1, label = format(date, "%Y-%m-%d")),
  col="orange",
  angle = 90,                # rotated text for readability
  hjust = -0.2,              # move slightly to the right of the line
  size = 3
) +
  
  xlab("Absolute humidity") +
  ylab("Observed growth rate") + 
  scale_color_manual("", values=c("LOESS"="purple", "LM"="darkgreen")) +
  scale_fill_manual("", values=c("LOESS"="purple", "LM"="darkgreen")) +
  theme_bw() +
  labs(tag = "B")


p_ah_data <- ggplot(comb1 %>% filter(year < 2025)) + 
  geom_line(aes(x=day,y=ah_smooth,group=year,col="Individual year"),alpha=0.5)+ 
  geom_line(data=.%>% filter(year == 2024), aes(x=day,y=ah_smooth,group=year,col="2024"),alpha=0.5)+ 
  geom_line(data=.%>% group_by(day) %>% summarize(mean_ah=mean(ah_smooth)), aes(x=day,y=mean_ah,col="Mean (2012-2024)"),linewidth=1)+
  #geom_vline(xintercept=yday(as.Date(c("2025-11-01","2025-12-01"))),linetype="dashed",col="black") +
  theme_bw() +
  theme(legend.position=c(0.8,0.8)) +
  scale_color_manual("",values=c("Individual year"="grey","Mean (2012-2024)"="blue","2024"="red")) +
  xlab("Day of year (starting 1st January)") +
  ylab("Absolute humidity") + labs(tag="A") +
  # dashed vertical lines
  geom_vline(data = vlines_df, aes(xintercept = day_of_year + 5), linetype = "dashed") +
  
  # ---- NEW: date labels next to dashed lines ----
geom_text(
  data = vlines_df,
  aes(x = day_of_year, y = 4, label = format(date, "%Y-%m-%d")),
  col="orange",
  angle = 90,                # rotated text for readability
  hjust = -0.2,              # move slightly to the right of the line
  size = 3
)


p_ah_data/p_cor


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
