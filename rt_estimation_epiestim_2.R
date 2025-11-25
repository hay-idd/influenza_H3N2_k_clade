
library(EpiEstim)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(patchwork)

setwd("~/Documents/GitHub/recent_influenza/Untitled")




#read the data:fluenet data 


flu_net_data<-read.csv("data/final/flunet_h3_cases_historic.csv")

#get the total sums:
incidence_weekly<-flu_net_data$flu_pos_sum

flu_net_data$date=as.Date(flu_net_data$date)



#interpolate the values to daily 
interpolate_daily_h3_sum=round(approx(x=flu_net_data$date, y=flu_net_data$H3_sum,xout=all_dates_data$date,method="linear",rule=2)$y)
interpolate_daily_positive_sum=round(approx(x=flu_net_data$date, y=flu_net_data$flu_pos_sum,xout=all_dates_data$date,method="linear",rule=2)$y)
         

#scale them so that daily values are equal to weekly totals 
daily_df <- data.frame(date = all_dates_data, interp_h3_sum = interpolate_daily_h3_sum, 
                       interp_positive_sum=interpolate_daily_positive_sum) %>%
  mutate(
    week_start = flu_net_data$date[findInterval(date, flu_net_data$date)]
  ) %>%
  group_by(week_start) %>%
  mutate(
    sum_interp_h3_sum = sum(interp_h3_sum),
    sum_interp_all_post=sum(interp_positive_sum),
    cases_interp_h3_sum = round(ifelse(sum_interp_h3_sum == 0, 0, (interp_h3_sum * flu_net_data$H3_sum[match(week_start, flu_net_data$date)]) / sum_interp_h3_sum)),
    cases_all_pos=round(ifelse(sum_interp_all_post == 0, 0, (interp_positive_sum * flu_net_data$flu_pos_sum[match(week_start, flu_net_data$date)]) / sum_interp_all_post))
  ) %>%
  ungroup() %>%
  select(date, cases_interp_h3_sum,cases_all_pos)



#enter the dataframe with date and cases as columns:
scale_daily_to_weekly <- function(weekly_data_frame) {
  
  # Create daily sequence
  all_dates <- seq(min(weekly_data_frame$date), max(weekly_data_frame$date) + days(6), by = "day")
  
  # Interpolate daily cases
  interpolate_daily_cases <- round(approx(
    x = weekly_data_frame$date,
    y = weekly_data_frame$cases,
    xout = all_dates,
    method = "linear",
    rule = 2
  )$y)
  
  # Create daily data frame and scale to weekly totals
  daily_df <- data.frame(
    date = all_dates,
    interp_cases = interpolate_daily_cases
  ) %>%
    mutate(
      week_start = weekly_data_frame$date[findInterval(date, weekly_data_frame$date)]
    ) %>%
    group_by(week_start) %>%
    mutate(
      sum_interp_cases = sum(interp_cases),
      cases_interp = round(ifelse(
        sum_interp_cases == 0, 0,
        (interp_cases * weekly_data_frame$cases[match(week_start, weekly_data_frame$date)]) / sum_interp_cases
      ))
    ) %>%
    ungroup() %>%
    select(date, cases_interp)
  
  return(daily_df)
}


df<-data.frame("date"=flu_net_data$date, "cases"=flu_net_data$flu_pos_sum)
xx<-scale_daily_to_weekly(df)

ggplot(data=xx,aes(x=date,y = cases_interp))+
  geom_line()+
  geom_line(data=df,aes(x=date,y=cases), color="red")

#define seriel interval, Mean SI = 3.6 days
#SD SI = 1.6 days
# Cowling et al., Epidemiology, 2009

mean_si <- 3.6 
std_si  <- 1.6 

method <- "parametric_si"
config <- make_config(list(mean_si = mean_si,
                           std_si = std_si))


dat<- data.frame("date"=xx$date,I = xx$cases_interp)

dat <- dat[!is.na(dat$I), ]

res_total <- estimate_R(
  incid = dat,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si
  ))
)



window_size <- 7
t_start <- seq(2, nrow(dat) - window_size)       # start of each window
t_end <- t_start + window_size                   # end of each window

res_total <- estimate_R(
  incid = dat,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si,
    t_start = t_start,
    t_end = t_end
  ))
)

XX=res_total$R
plot(res_total)





# Plot daily incidence
# -----------------------------
p_incidence <- ggplot(daily_df, aes(x = date, y = cases_interp)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "Daily Interpolated Cases", x = "Date", y = "Cases") +
  theme_minimal()
p_incidence



daily_df<-daily_df[4541:5054,]

# Plot Rt estimates

Rt_df <- data.frame(
  date = daily_df$date[res_total$R$t_end],
  Rt_mean = res_total$R$`Mean(R)`,
  Rt_lower = res_total$R$`Quantile.0.025(R)`,
  Rt_upper = res_total$R$`Quantile.0.975(R)`
)

Rt_df<-Rt_df[4541:5054,]

p_Rt <- ggplot(Rt_df, aes(x = date, y = Rt_mean)) +
  geom_line(color = "darkred", size = 1) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
  labs(title = "Estimated Reproduction Number (Rt)", x = "Date", y = "Rt") +
  theme_minimal()
p_Rt

p_incidence/p_Rt









# Create a ggplot
ggplot(res_total$R, aes(x = t_start, y = `Mean(R)`)) +
  geom_line(color = "blue") +
  ylim(c(0,10))+
  geom_ribbon(aes(ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), alpha = 0.2) +
  labs(x = "Time", y = "R", title = "Estimated R") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # <-- horizontal line
  theme_minimal() +
  theme(legend.position = "bottom")



