library(EpiEstim)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(patchwork)
library(tidyr)
setwd("~/Documents/GitHub/recent_influenza/Untitled")



#read the data:fluenet data 
flu_net_data<-read.csv("data/final/flunet_h3_cases_historic.csv")

flu_net_data$date=as.Date(flu_net_data$date)

#enter the dataframe with date and cases as columns:
scale_daily_to_weekly <- function(weekly_data_frame,roll_d) {
  
  # Create daily sequence
  all_dates <- seq(min(weekly_data_frame$date), max(weekly_data_frame$date) + days(6), by = "day")
ln<-length(all_dates)


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

#: Remove NA rows 
daily_df <- na.omit(daily_df)


  return(daily_df)
}



df<-data.frame("date"=flu_net_data$date, "cases"=flu_net_data$flu_pos_sum)
daily_data<-scale_daily_to_weekly(df,14)

#ggplot(data=daily_data,aes(x=date,y = rolling_mean))+
 # geom_line()+
 # geom_line(data=df,aes(x=date,y=cases), color="red")

#define seriel interval, Mean SI = 3.6 days
#SD SI = 1.6 days
# Cowling et al., Epidemiology, 2009

mean_si <- 3.6 
std_si  <- 1.6 



dat<- data.frame("date"=daily_data$date,I = round(daily_data$rolling_mean))
dat <- dat[!is.na(dat$I), ]




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



#daily data label seasons and days 
daily_data <- daily_data %>%
  mutate(
    season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
    season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )

p_incidence <- ggplot(daily_data, aes(x = date, y = rolling_mean)) +
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



Rt_df <- Rt_df %>%
  mutate(
   season_year = if_else(month(date) >= 9, year(date), year(date) - 1),
   season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01"))) # day 0 = July 1
 )

# 2. Plot Rt by season using facet_wrap
p_Rt <- ggplot(Rt_df, aes(x = date, y = Rt_mean)) +
  geom_line(color = "darkred", size = 1) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
  facet_wrap(~ season_label, scales = "free") +
  labs(title = "Estimated Reproduction Number (Rt) by Influenza Season",
       x = "Date", y = "Rt") +
  theme_minimal()

p_Rt


xx=unique(daily_data$season_label)

for (i in 1:length(xx)) {
  
  p_incidence <- ggplot(subset(daily_data,daily_data$season_label==xx[i]), aes(x = date, y = rolling_mean)) +
    geom_line(color = "steelblue", linewidth = 1) +
    facet_wrap(~ season_label, scales = "free") +
    ylab("Cases")+
    xlab("Date")+
    # labs(title = "Daily Interpolated Cases", x = "Date", y = "Cases") +
    theme_minimal()

  
  p_Rt <- ggplot(subset(Rt_df,Rt_df$season_label==xx[i]), aes(x = date, y = Rt_mean)) +
    geom_line(color = "darkred", size = 1) +
    geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
    facet_wrap(~ season_label, scales = "free") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # <-- horizontal line
    # labs(title = "Estimated Reproduction Number (Rt) by Influenza Season",
    #    x = "Date", y = "Rt") +
    ylab("Estimated Rt")+
    xlab("Date")+
    theme_minimal()
  
  
 pl<- p_incidence/p_Rt
 name<-paste0("rt_estimates_plots/plot_",xx[i],".png") 
 
  ggsave(name,pl,height =6 ,width = 4)
}




# Assume Rt_df has columns: date, Rt_mean, Rt_lower, Rt_upper, season_label, pandemic_phase
# pandemic_phase: "Pre-pandemic", "Pandemic", "Post-pandemic"

Rt_df <- Rt_df %>%
  mutate(day_of_year = yday(date))  # Align by day of year



Rt_df <- Rt_df %>%
  mutate(pandemic_phase = case_when(
    date < as.Date("2020-03-01") ~ "Pre-pandemic",
    date >= as.Date("2020-03-01") & date <= as.Date("2023-05-31") ~ "Pandemic",
    date > as.Date("2023-05-31") ~ "Post-pandemic"
  ))


#force the entire 2022-2023 to be pandemic 
Rt_df$pandemic_phase[Rt_df$season_label=="2022-2023"]="Pandemic"

level_f<-c("Pre-pandemic","Pandemic","Post-pandemic")
Rt_df$pandemic_phase=as.factor(Rt_df$pandemic_phase)
Rt_df$pandemic_phase=factor(Rt_df$pandemic_phase,levels = level_f)


# Plot Rt estimates aligned by day of year

#remove the rt estiates that don't look weird 
Rt_df <- subset(Rt_df, !(season_label %in% c("2011-2012","2012-2013", "2013-2014","2014-2015", "2019-2020","2020-2021")))

p_Rt_1 <- ggplot(Rt_df, aes(x = day_in_season, y = Rt_mean, color = pandemic_phase)) +
  geom_line(size = 1) +
  coord_cartesian(ylim=c(0,2))+
  geom_hline(yintercept = 1.2, linetype = "dashed", color = "darkred",size=.25) + 
  geom_hline(yintercept = 1.4, linetype = "dashed", color = "darkred",size=.25) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = pandemic_phase), alpha = 0.2, color = NA) +
 facet_wrap(~season_label, ncol = 3) +  # vertical alignment
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Estimated reproduction number (Rt) by influenza season",
       x = "Day in season (season starts September 1)", y = "Estimated Rt") +
  theme_minimal()+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")

p_Rt_1
 
ggsave("rt_estimates_plots/rt_estimates_by_pandemic_phase_and_season.png",p_Rt_1,height =8 ,width = 10)

p_Rt_2 <- ggplot(Rt_df, aes(x = day_in_season, y = Rt_mean, color = season_label)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1.2, linetype = "dashed", color = "darkred",size=.25) + 
  geom_hline(yintercept = 1.4, linetype = "dashed", color = "darkred",size=.25) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = season_label), alpha = 0.2, color = NA) +
  facet_wrap(~ season_label, ncol = 2) +  # vertical alignment
  coord_cartesian(ylim=c(0,2))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
 # scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
 # scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Estimated reproduction number (Rt) by influenza season",
       x = "Day in season (season starts September 1)", y = "Estimated Rt") +
  theme_minimal()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_Rt_2
ggsave("rt_estimates_plots/rt_estimates_by_season.png",p_Rt_2,height =8 ,width = 10)




# Find peak Rt date for each season (pre-Christmas)
peak_dates <- Rt_df %>%
  filter(month(date) < 12) %>%  # restrict to pre-Christmas
  group_by(season_label) %>%
  slice_max(order_by = Rt_mean, n = 1) %>%
  select(season_label, peak_date = date)


Rt_df_shifted <- Rt_df %>%
  left_join(peak_dates, by = "season_label") %>%
  mutate(days_since_peak = as.numeric(date - peak_date))  # x-axis variable


p_Rt_shifted_1 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean, color = pandemic_phase)) +
  geom_line(size = 1) +
  xlim(c(-50,50))+
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = pandemic_phase), alpha = 0.2, color = NA) +
  facet_wrap(~ season_label, ncol = 2, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Rt estimates aligned to peak pre-Christmas Rt",
       x = "Days since peak Rt", y = "Estimated Rt") +
 theme(legend.position = "bottom")+
  theme(legend.title = element_blank())

p_Rt_shifted_1

ggsave("rt_estimates_plots/pre-Christmas_rt_estimates_by_season_1.png",p_Rt_shifted_1,height =8 ,width = 10)
# Define current season
current_season <- "2025-2026"

Rt_df_shifted <- Rt_df_shifted %>%
  mutate(highlight = if_else(season_label == current_season, "Current", "Other"))

p_Rt_shifted_2 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean,
                                          color = highlight, fill = highlight)) +
  geom_line(size = 1) +
  xlim(c(-50, 50)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, color = NA) +
  facet_wrap(~ season_label, ncol = 2, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Current" = "#d95f02", "Other" = "grey70")) +
  scale_fill_manual(values = c("Current" = "#d95f02", "Other" = "grey85")) +
  labs(title = "Rt estimates aligned to peak pre-Christmas Rt",
       x = "Days since peak Rt", y = "Estimated Rt") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

p_Rt_shifted_2

ggsave("rt_estimates_plots/pre-Christmas_rt_estimates_by_season_2.png",p_Rt_shifted_2,height =8 ,width = 10)


# Define the current season (e.g., "2023-2024")
current_season <- "2025-2026"

p_Rt_shifted_stacked_1<-ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean,
                          color = season_label, fill = season_label)) +
  geom_line(size = 1) +
  xlim(c(-50, 50)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(title = "Rt estimates aligned to pre-Christmas peak Rt",
       x = "Days since peak Rt", y = "Estimated Rt") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_manual(values = ifelse(unique(Rt_df_shifted$season_label) == current_season,
                                     "#d95f02",  # Highlight color for current season
                                     "grey70")) +
  scale_fill_manual(values = ifelse(unique(Rt_df_shifted$season_label) == current_season,
                                    "#d95f02",  # Same highlight color for ribbon
                                    "grey80"))
p_Rt_shifted_stacked_1

ggsave("rt_estimates_plots/pre-Christmas_rt_estimates_by_season_staked_1.png",p_Rt_shifted_stacked_1,height =8 ,width = 10)

p_Rt_shifted_stacked_2 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean, color = season_label)) +
  geom_line(size = 1) +
  xlim(c(-50,50))+
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = season_label), alpha = 0.2, color = NA) +
#  facet_wrap(~ season_label, ncol = 2, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  #scale_color_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
#  scale_fill_manual(values = c("Pre-pandemic" = "#1b9e77", "Pandemic" = "#d95f02", "Post-pandemic" = "#7570b3")) +
  labs(title = "Rt estimates aligned to pre-Christmas peak Rt",
       x = "Days since peak Rt", y = "Estimated Rt") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())


p_Rt_shifted_stacked_2

ggsave("rt_estimates_plots/pre-Christmas_rt_estimates_by_season_staked_2.png",p_Rt_shifted_stacked_2,height =8 ,width = 10)

# Peak Rt per season
peak_table <- Rt_df %>%
  group_by(season_label) %>%
  slice_max(order_by = Rt_mean, n = 1) %>%
  ungroup() %>%
  select(season_label, Peak_Rt = Rt_mean, Peak_Date = date)

#  Peak Rt per season (restricted to pre-Christmas)
peak_pre_xmas_table <- Rt_df %>%
  filter(day_of_year <= 358) %>%  # Dec 24 or earlier
  group_by(season_label) %>%
  slice_max(order_by = Rt_mean, n = 1) %>%
  ungroup() %>%
  select(season_label, Peak_Rt_PreChristmas = Rt_mean, Peak_Date_PreChristmas = date)


# Save as CSV
write.csv(peak_table, "rt_estimates_plots/peak_Rt_per_season.csv", row.names = FALSE)
write.csv(peak_pre_xmas_table, "rt_estimates_plots/peak_Rt_preChristmas.csv", row.names = FALSE)

