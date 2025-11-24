library(EpiEstim)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(patchwork)
library(tidyr)
library(ggthemes)

#set the working directory 
setwd("~/Documents/GitHub/recent_influenza/Untitled")


#read the data:fluenet data 
flu_net_data<-read.csv("data/final/flunet_h3_cases_historic.csv")

#covet the dates to a Date class in R 
flu_net_data$date=as.Date(flu_net_data$date)

#rename the columns to make it easy to run the analysis and functions
df<-data.frame("date"=flu_net_data$date, "cases"=flu_net_data$flu_pos_sum)

####function to generate daily incidence from weekly incidence data 
#enter the dataframe with date and cases as columns:
scale_daily_to_weekly <- function(weekly_data_frame,roll_d) {
  #assign the average form the weekly case data to each day, and take the rolling mean 
  
  daily_df <- weekly_data_frame %>%
    rowwise() %>%
    mutate(
      daily_dates = list(seq(date, by = "day", length.out = 7)),
      daily_cases = list(rep(cases/7, 7))
    ) %>%
    unnest(cols = c(daily_dates, daily_cases)) %>%
    select(daily_dates, daily_cases) %>%  # Keep only daily columns
    rename(date = daily_dates, cases = daily_cases)
  
  
  # Compute 7-day rolling mean of cases
  daily_df$rolling_mean <- rollmean(daily_df$cases, k = roll_d, fill = NA, align = "right")
  
  # Compute growth rate: ratio of current rolling mean to previous rolling mean
  daily_df$growth_rate <- daily_df$rolling_mean / lag(daily_df$rolling_mean)
  
  #Remove NA rows 
  daily_df <- na.omit(daily_df)
  
  return(daily_df)
}

#use the function to 
#convert the weekly data to daily data: 
daily_data<-scale_daily_to_weekly(df,14)

#to compare the daily cases vs. weekly cases: 
#ggplot(data=daily_data,aes(x=date,y = rolling_mean))+
# geom_line()+
# geom_line(data=df,aes(x=date,y=cases), color="red")


#####estimate the Rt:
dat<- data.frame("date"=daily_data$date,I = round(daily_data$rolling_mean))
dat <- dat[!is.na(dat$I), ]


#define seriel interval, Mean SI = 3.6 days
#SD SI = 1.6 days
# Cowling et al., Epidemiology, 2009

mean_si <- 3.6 
std_si  <- 1.6 

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



###visulaisation of the data:
#
#daily data label seasons and days 
daily_case_data <- daily_data %>%
  mutate(
    season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
    season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )




daily_case_data <- daily_data %>%
  mutate(
    season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
    season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01"))),
    subtype = case_when(
      season_label == "2020-2021" ~ "NA (COVID-19 pandmeic)",  # COVID year
      season_label %in% c("2015-2016",,"2018-2019", "2023-2024", "2024-2025") ~ "A/H1N1pdm09",
      TRUE ~ "A/H3N2"
    )
  )


p_incidence <- ggplot(daily_case_data, aes(x = day_in_season, y = rolling_mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  facet_wrap(~ season_label, scales = "free") +
  labs(title = "Daily Interpolated Cases", x = "Date", y = "Cases") +
  theme_minimal()
p_incidence



# Plot Rt estimates
Rt_df <- data.frame(
  date = daily_data$date[res_total$R$t_end],
  Rt_mean = res_total$R$`Mean(R)`,
  Rt_lower = res_total$R$`Quantile.0.025(R)`,
  Rt_upper = res_total$R$`Quantile.0.975(R)`
)



#Rt_df <- Rt_df %>%
#  mutate(
  #  season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
  #  season_label = paste0(season_year, "-", season_year + 1),
 #   day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01"))) # day 0 = September 1
 # )







Rt_df <- Rt_df %>%
  mutate(
    season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
    season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01"))), # day 0 = September 1
    subtype = case_when(
      season_label == "2020-2021" ~"NA (COVID-19 pandmeic)",
      season_label %in% c("2015-2016",,"2018-2019", "2023-2024", "2024-2025") ~ "A/H1N1pdm09",
      TRUE ~ "A/H3N2"
    )
  )




# 2. Plot Rt by season using facet_wrap
p_Rt <- ggplot(Rt_df, aes(x = day_in_season, y = Rt_mean)) +
  geom_line(color = "darkred", size = 1) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
  facet_wrap(~ season_label, scales = "free") +
  labs(title = "Estimated reproduction number (Rt) by influenza season",
       x = "Date", y = "Rt") +
  theme_minimal()

p_Rt


xx=unique(daily_case_data$season_label)

for (i in 1:length(xx)) {
  
  p_incidence <- ggplot(subset(daily_case_data,daily_case_data$season_label==xx[i]), aes(x = day_in_season, y = rolling_mean)) +
    geom_line(color = "steelblue", linewidth = 1) +
    facet_wrap(~ season_label, scales = "free") +
    ylab("Cases")+
    xlab("Date")+
    # labs(title = "Daily Interpolated Cases", x = "Date", y = "Cases") +
    theme_minimal()
  
  
  p_Rt <- ggplot(subset(Rt_df,Rt_df$season_label==xx[i]), aes(x = day_in_season, y = Rt_mean)) +
    geom_line(color = "darkred", size = 1) +
    geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
    facet_wrap(~ season_label, scales = "free") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    ylab("Estimated Rt")+
    xlab("Date")+
    theme_minimal()
  
  
  pl<- p_incidence/p_Rt
  name<-paste0("rt_estimates_plots/plot_",xx[i],".png") 
  
  ggsave(name,pl,height =6 ,width = 4)
}





# pandemic_phase: "Pre-pandemic", "Pandemic", "Post-pandemic"
Rt_df <- Rt_df %>%
  mutate(day_of_year = yday(date))  # Align by day of year


Rt_df <- Rt_df %>%
  mutate(pandemic_phase = case_when(
    date < as.Date("2020-03-01") ~ "Pre-pandemic",
    date >= as.Date("2020-03-01") & date <= as.Date("2023-05-31") ~ "Pandemic",
    date > as.Date("2023-05-31") ~ "Post-pandemic"
  ))

#also do the same for the case data: 
daily_case_data <- daily_case_data %>%
  mutate(pandemic_phase = case_when(
    date < as.Date("2020-03-01") ~ "Pre-pandemic",
    date >= as.Date("2020-03-01") & date <= as.Date("2023-05-31") ~ "Pandemic",
    date > as.Date("2023-05-31") ~ "Post-pandemic"
  ))


#force the entire 2022-2023 to be pandemic 
Rt_df$pandemic_phase[Rt_df$season_label=="2022-2023"]="Pandemic"
daily_case_data$pandemic_phase[Rt_df$season_label=="2022-2023"]="Pandemic"


#to plot in the order:
level_f<-c("Pre-pandemic","Pandemic","Post-pandemic")
Rt_df$pandemic_phase=as.factor(Rt_df$pandemic_phase)
Rt_df$pandemic_phase=factor(Rt_df$pandemic_phase,levels = level_f)
daily_case_data$pandemic_phase=as.factor(daily_case_data$pandemic_phase)
daily_case_data$pandemic_phase=factor(daily_case_data$pandemic_phase,levels = level_f)

# Plot Rt estimates aligned by day in season of year

#remove the rt estimates that don't have good estimates for the Rt 
years_to_remove<-c("2011-2012","2012-2013", "2013-2014","2014-2015", "2019-2020","2020-2021",
                   "2021-2022")
Rt_df <- subset(Rt_df, !(season_label %in% years_to_remove))
#do the same for the daily case data: 
daily_case_data<-subset(daily_case_data, !(season_label %in% years_to_remove))


case_pl_1<-ggplot(daily_case_data, aes(x = day_in_season, y = rolling_mean, color = subtype)) +
  geom_line(size = 1) +
  #coord_cartesian(ylim=c(0,1500))+
  theme_bw()+
  facet_wrap(~season_label, ncol = 2,scales="free_y") +  # vertical alignment
  scale_color_brewer(palette = "Dark2")+
  #scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  #scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Daily incidence (interpolated) by influenza season",
       x = "Day in season (season starts September 1)", y = "Daily incidence") +
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")

case_pl_1


p_Rt_1 <- ggplot(Rt_df, aes(x = day_in_season, y = Rt_mean, color = subtype)) +
  geom_line(size = 1) +
  coord_cartesian(ylim = c(0, 2.5)) +
  geom_hline(yintercept = 1.2, linetype = "dotted", color = "darkred", size = .5) +
  geom_hline(yintercept = 1.4, linetype = "dotted", color = "darkred", size = .5) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = subtype), alpha = 0.4, color = NA) +
  facet_wrap(~season_label, ncol = 2) +
  theme_bw()+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
 # scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  #scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Estimated effective reproduction number (Rt) by influenza season",
       x = "Day in season (season starts September 1)", y = expression("Estimated " * R[t])) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

p_Rt_1

pl<-case_pl_1/p_Rt_1
pl

ggsave("rt_estimates_plots/rt_estimates_by_pandemic_phase_and_season_x1.png",pl,height =14 ,width = 10)




# Find peak Rt date for each season (pre-Christmas)
peak_dates <- Rt_df %>%
  group_by(season_label) %>%
  filter(month(date) > 6  & (month(date) < 12 | day(date)<25 )) %>%  
  slice_max(order_by = Rt_mean, n = 1) %>%
  select(season_label, peak_date = date, peak_Rt = Rt_mean, Rt_lower=Rt_lower,   Rt_upper= Rt_upper)


Rt_df_shifted <- Rt_df %>%
  left_join(peak_dates, by = "season_label", suffix = c("", "_peak")) %>%
  mutate(days_since_peak = as.numeric(date - peak_date))

p_Rt_shifted_1 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean, color = subtype)) +
  geom_line(size = 1) +
  coord_cartesian(ylim=c(0,2.5))+
  xlim(c(-50,50))+
  # xlim(c(-100,100))+
  theme_bw()+
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkred") +
  # Add text labels for peak Rt  
  geom_text(data = peak_dates, aes(x = -30, y = 2, 
                                   label = paste0(  "Peak on ", peak_dates$peak_date, " \n with Rt= ",
                                                    round(peak_dates$peak_Rt, 2)," (",round(peak_dates$Rt_lower, 2),",",round(peak_dates$Rt_upper, 2),")", "") ), 
            inherit.aes = FALSE,    color = "darkred",     vjust = -0.5,   size = 3  ) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = subtype), alpha = 0.4, color = NA) +
  facet_wrap(~ season_label, ncol = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  #scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
 # scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Rt estimates aligned to peak pre-Christmas Rt",
       x = "Days since peak Rt", y = expression("Estimated " * R[t])) +
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())

p_Rt_shifted_1



ggsave("rt_estimates_plots/pre-Christmas_rt_estimates_by_season_1_before_chrismas_day_x1.png",p_Rt_shifted_1,height =12 ,width = 10)


# Peak Rt per season 
peak_table <- Rt_df %>%
  group_by(season_label) %>%
  slice_max(order_by = Rt_mean, n = 1) %>%
  ungroup() %>%
  select(season_label, Peak_Rt = Rt_mean, Peak_Date = date)

#  Peak Rt per season ( pre-Christmas)
peak_pre_xmas_table <- Rt_df %>%
  filter(month(date) > 6  & (month(date) < 12 | day(date)<25 )) %>%  
  #filter(month(date) < 12 | (month(date) == 12 & day(date) < 24)) %>%  # pre-Christmas (before Dec 24)
  group_by(season_label) %>%
  slice_max(order_by = Rt_mean, n = 1) %>%
  ungroup() %>%
  select(season_label, Peak_Rt_PreChristmas = Rt_mean, Peak_Date_PreChristmas = date)


# Save as CSV
write.csv(peak_table, "rt_estimates_plots/peak_Rt_per_season.csv", row.names = FALSE)
write.csv(peak_pre_xmas_table, "rt_estimates_plots/peak_Rt_preChristmas.csv", row.names = FALSE)

