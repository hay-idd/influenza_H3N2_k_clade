library(EpiEstim)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(patchwork)
library(tidyr)
library(ggthemes)
library(reshape2)

#set the working directory 
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
source("R/rt_estimation_help_functions.R")
source("R/funcs.R")

#read the data:fluenet data 
england_data<-read_csv("data/resp_datamart_influenza_cases_england.csv")

#covet the dates to a Date class in R 
#england_data$date=as.Date(england_data$date,format = "%d/%m/%Y")

#add the inluenza H3 and H1 types 
df<-data.frame("date"=england_data$date, "h3_cases"=england_data$H3_adj,
               "h1_cases"=england_data$H1_adj)


m_dat<-df %>% pivot_longer(-date)

m_dat<-m_dat %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    season_label = season_from_date(date,8),
    #season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )

#plot the data by strain: 
raw_data_plot<-ggplot(data=m_dat,aes(x=day_in_season,y=value,color=name))+
  geom_line()+
  facet_wrap(~season_label,ncol=4,scales="free_y")
raw_data_plot


#use the function to 
#convert the weekly data to daily data: 
#h3 cases data: 
h3_df<-data.frame("date"=df$date,"cases"=df$h3_cases)
h1_df<-data.frame("date"=df$date,"cases"=df$h1_cases)


#to compare the daily cases vs. weekly cases: 
#ggplot(data=daily_data,aes(x=date,y = rolling_mean))+
# geom_line()+
# geom_line(data=df,aes(x=date,y=cases), color="red")


#####estimate the Rt:
dat_h1<-h1_df %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    season_label = season_from_date(date,8),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )


dat_h3<-h3_df %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    season_label = season_from_date(date,8),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )

#remove tails;
weekly_data_h1<-remove_tails(dat_h1)
weekly_data_h1$subtype<-"A/H1N1pdm09"
weekly_data_h3<-remove_tails(dat_h3)
weekly_data_h3$subtype<-"A/H3N2"

ggplot(weekly_data_h1,aes(x=day_in_season,y=cases))+
  facet_wrap(~ season_label, scales = "free_y") +
  geom_line( linewidth = 1) 

#remove: 2011
yr_rm_h1<-c("2011/12","2016/17","2020/21","2021/22")

weekly_data_h1<-subset(weekly_data_h1 ,!(season_label %in% yr_rm_h1))

ggplot(weekly_data_h3,aes(x=day_in_season,y=cases))+
  facet_wrap(~ season_label, scales = "free_y") +
  geom_line( linewidth = 1) 

yr_rm_h3<-c("2008/09","2009/10","2010/11","2020/21","2021/22")
weekly_data_h3<-subset(weekly_data_h3 ,!(season_label %in% yr_rm_h3))



#turn to daily data:
#use the function to 
#convert the weekly data to daily data: 
#h3 cases data: 
h3_df<-data.frame("date"=weekly_data_h3$date,"cases"=weekly_data_h3$cases)
h1_df<-data.frame("date"=weekly_data_h1$date,"cases"=weekly_data_h1$cases)

daily_data_h1<-scale_weekly_to_daily(h1_df,7)
daily_data_h3<-scale_weekly_to_daily(h3_df,7)

##to check the estimation: 
# converted_weekly_h1<-rolling_mean_to_weekly(daily_data_h1,weekly_data_h1$date[1])
# converted_weekly_h3<-rolling_mean_to_weekly(daily_data_h3,weekly_data_h1$date[1])
# 
# converted_weekly_h1<-converted_weekly_h1 %>%
#   mutate(
#     season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
#     season_label = paste0(season_year, "-", season_year + 1),
#     day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
#   )
# 
# d_h1=weekly_data_h1
# d_h1$method<-"orginal"
# converted_weekly_h1$subtype<-"A/H1N1pdm09"
# converted_weekly_h1$method<-"back_calulated"
# 
# all_dd<-rbind(d_h1,converted_weekly_h1)
# 
# ggplot(all_dd,aes(x=date,y=cases,color=method))+
#   geom_line()

#smoothed data: 
dat_h1<- data.frame("date"=daily_data_h1$date,I = round(daily_data_h1$rolling_mean))
dat_h1 <- dat_h1[!is.na(dat_h1$I), ]

dat_h1<-dat_h1 %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    season_label = season_from_date(date,8), #paste0(season_year, "-", season_year + 1),
    #season_label = paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )

dat_h1$subtype="A/H1N1pdm09"

dat_h3<- data.frame("date"=daily_data_h3$date,I = round(daily_data_h3$rolling_mean))
dat_h3 <- dat_h3[!is.na(dat_h3$I), ]
dat_h3$subtype="A/H3N2"

dat_h3<-dat_h3 %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    season_label = season_from_date(date,8), #paste0(season_year, "-", season_year + 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01")))
  )

## Remove dates prior to 1% of total cases and after 99% of total cases
dat_h1_filter_date_start <- dat_h1 %>% group_by(season_label) %>% mutate(I_prop = I/sum(I)) %>% mutate(cum_sum=cumsum(I_prop))%>%group_by(season_label) %>% mutate(diff=abs(cum_sum - 0.01)) %>% filter(diff==min(diff)) %>%
  select(season_label,date) %>% rename(start_date=date)
dat_h1_filter_date_end <- dat_h1 %>% group_by(season_label) %>% mutate(I_prop = I/sum(I)) %>% mutate(cum_sum=cumsum(I_prop))%>%group_by(season_label) %>% mutate(diff=abs(cum_sum - 0.98)) %>% filter(diff==min(diff)) %>%
  select(season_label,date) %>% rename(end=date)
dat_h1 <- dat_h1 %>% left_join(dat_h1_filter_date_start,by="season_label") %>% 
  left_join(dat_h1_filter_date_end,by="season_label") %>%
  filter(date>=start_date & date<=end) %>% select(-start_date,-end)  

dat_h3_filter_date_start <- dat_h3 %>% group_by(season_label) %>% mutate(I_prop = I/sum(I)) %>% mutate(cum_sum=cumsum(I_prop))%>%group_by(season_label) %>% mutate(diff=abs(cum_sum - 0.01)) %>% filter(diff==min(diff)) %>%
  select(season_label,date) %>% rename(start_date=date)
dat_h3_filter_date_end <- dat_h3 %>% group_by(season_label) %>% mutate(I_prop = I/sum(I)) %>% mutate(cum_sum=cumsum(I_prop))%>%group_by(season_label) %>% mutate(diff=abs(cum_sum - 0.98)) %>% filter(diff==min(diff)) %>%
  select(season_label,date) %>% rename(end=date)
dat_h3 <- dat_h3 %>% 
  left_join(dat_h3_filter_date_start,by="season_label") %>% 
  left_join(dat_h3_filter_date_end,by="season_label") %>%
  mutate(season=season_from_date(date,8)) %>%
  filter((date>=start_date & date <= end) | (season == "2025/26" & date >= start_date))%>%
  select(-start_date,-end,-season)



daily_case_data<-rbind(dat_h1,dat_h3)

#plot the data;
p_incidence <- ggplot(daily_case_data, aes(x = date, y = I,color=subtype)) +
  geom_line( linewidth = 1) +
  #geom_line(data=m_dat,aes(x=day_in_season,y=value,color=variable)) +
  facet_wrap(~ season_label, scales = "free") +
  labs(title = "Daily Cases", x = "Date", y = "Cases") +
  theme_minimal()

#p_incidence


#incidence by dominant subtype:
rm_yrs_in_h1<-c("2011/12","2012/13","2014/15","2016/17","2017/18","2019/20","2022/23","2025/26")
keep_yrs_in_h1 <- c("2013/14","2015/16","2018/19","2023/24","2024/25")
dat_h1<-subset(dat_h1,(season_label %in% keep_yrs_in_h1))
#rm_yrs_in_h3<-c("2013/14","2015/16","2018/19","2024/25")
keep_yrs_in_h3<-c("2011/12","2012/13","2014/15","2016/17","2017/18","2022/23","2023/24","2025/26")
dat_h3<-subset(dat_h3,(season_label %in% keep_yrs_in_h3))


daily_case_data<-rbind(dat_h1,dat_h3)

#plot the data;

p_incidence <- ggplot(daily_case_data%>% 
                        mutate(season = season_from_date(date,8))%>% filter(!(season_label %in% c("2008/09","2009/10","2010/11",
                                                                "2019/20","2020/21","2021/22"))), aes(x = date, y = I, color = subtype)) +
  geom_line(linewidth = 1) +
  theme_bw()+    
  labs(#title = "Daily incidence by influenza season",
       x = "Date", y = "Daily incidence", color = "Dominant influenza A subtype") +
  facet_wrap(~ season_label, scales = "free",ncol=2) +
  theme(legend.position = "none")+
  scale_color_manual(values = c("#CC79A7","#0072B2"))


p_incidence


#define seriel interval, Mean SI = 3.6 days
#SD SI = 1.6 days
# Cowling et al., Epidemiology, 2009
mean_si <- 3.6 
std_si  <- 1.6 


window_size <- 14
dtout<-window_size 

Rt_df_h1<-estimate_Rt_by_season(dat_h1,window_size,dtout)
Rt_df_h1$subtype="A/H1N1pdm09"
Rt_df_h3<-estimate_Rt_by_season(dat_h3,window_size,dtout)
Rt_df_h3$subtype="A/H3N2"

Rt_df<-rbind(Rt_df_h1,Rt_df_h3)


Rt_df <- Rt_df %>%
  mutate(
    season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
    day_in_season = as.numeric(date - as.Date(paste0(season_year, "-09-01"))) # day 0 = September 1
  )


p_Rt <- ggplot(Rt_df, aes(x = day_in_season, y = Rt_mean)) +
  geom_line(color = "darkred", size = 1) +
  coord_cartesian(ylim=c(0,2))+
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.2, fill = "red") +
  facet_wrap(subtype~ season) +
  labs(title = "Estimated reproduction number (Rt) by influenza season",
       x = "Date", y = "Rt") +
  theme_minimal()

p_Rt

#Rt_df_2=Rt_df
#Rt_df_2$method="daily"

p <- 0.02  # proportion of cumulative cases (user-defined)

# cumulative incidence by season
cum_inc_df <- daily_case_data %>%
  group_by(season_label, subtype) %>%
  arrange(day_in_season) %>%
  mutate(
    cum_cases = cumsum(I),
    cum_prop  = cum_cases / max(cum_cases)
  ) %>%
  ungroup()

# join to Rt data and find max Rt after p% of cases
Rt_peak_df <- Rt_df %>%
  left_join(
    cum_inc_df %>% 
      select(season_label, subtype, day_in_season, cum_prop),
    by = c("season" = "season_label", "subtype", "day_in_season")
  ) %>%
  filter(cum_prop >= p) %>%
  group_by(season, subtype) %>%
  slice_max(Rt_mean, n = 1, with_ties = FALSE) %>%
  ungroup()



p_Rt <- ggplot(Rt_df %>% filter(!(season %in% c("2008/09","2009/10","2010/11",
                                                "2019/20","2021/22"))),
               aes(x = day_in_season, y = Rt_mean, color = subtype)) +
  geom_hline(yintercept = c(1), linetype = "dashed", color = "grey40",linewidth=0.5) +
  geom_line(size = 1) +
  #ylim(0,3)+
  geom_ribbon(
    aes(ymin = Rt_lower, ymax = Rt_upper, fill = subtype),
    alpha = 0.3,
    colour = NA
  ) +
  coord_cartesian(ylim=c(0,2)) +
  scale_y_continuous(breaks=seq(0,2,by=0.5)) +
  facet_wrap( ~ season, ncol=2) +
  scale_color_manual(name = "Dominant influenza A subtype", values = c("#CC79A7","#0072B2")) +
  scale_fill_manual(name = "Dominant influenza A subtype", values = c("#CC79A7","#0072B2")) +
  labs(
    #title = "Estimated Rt for dominant subtype by influenza season",
    x = "Date", y = "Rt for dominant influenza A subtype"
  ) +
  theme_bw()+                          # put minimal theme first
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  geom_vline(
    data = Rt_peak_df,
    aes(xintercept = day_in_season, color = subtype),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_text(
    data = Rt_peak_df %>% filter(!(season == "2023/24" & subtype == "A/H3N2")),
    aes(
      x = 150,
      y = 1.3,
      label = sprintf(
        "Rt = %.2f (%.2f–%.2f)\n%s",
        Rt_mean, Rt_lower, Rt_upper, date
      )
    ),
    hjust = -0.05,
    vjust = 0,
    size = 3,
    show.legend = FALSE
  )



pl<-(p_incidence | p_Rt) + plot_layout(guides="collect") & theme(legend.position="bottom")
pl

ggsave("figures/Fig3.png",p_Rt,height =12 ,width = 8)
ggsave("figures/Fig3.pdf",p_Rt,height =12 ,width = 8)
ggsave("figures/figS7.png",p_incidence + theme(legend.position="bottom"),height =10 ,width = 7)
ggsave("figures/figS7.pdf",p_incidence+ theme(legend.position="bottom"),height =10 ,width = 7)

#ggsave("figures/Fig3.png",pl,height =12 ,width = 12)
#ggsave("figures/Fig3.pdf",pl,height =12 ,width = 12)


write.csv(Rt_peak_df,"figures/peak_rt_estimates.csv",row.names = FALSE)




peak_single <- daily_case_data %>%
  filter(!is.na(I)) %>%
  group_by(subtype, season_label) %>%
  arrange(date, .by_group = TRUE) %>%
  slice_max(order_by = I, n = 1, with_ties = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    subtype,
    season_label,
    max_I = I,
    peak_date = date,
    day_in_season_peak = day_in_season,
    season_week = floor(day_in_season / 7) + 1
  )





peak_single <- peak_single %>% mutate(peak_date = as.Date(peak_date))
Rt_df       <- Rt_df       %>% mutate(date      = as.Date(date))


Rt_df <- Rt_df %>%
  rename(season_label= season)

# Join by subtype & season, then filter for rows at peak date
Rt_at_peak <- peak_single %>%
  select(subtype, season_label, peak_date) %>%
  distinct() %>%
  inner_join(Rt_df %>% select(subtype, season_label, date, Rt_mean,Rt_lower,Rt_upper),
             by = c("subtype", "season_label")) %>%
  filter(date == peak_date) %>%
  select(subtype, season_label, peak_date, Rt_mean,Rt_lower,Rt_upper)

Rt_at_peak



Rt_df_shifted <- Rt_df %>%
  left_join(Rt_at_peak,
            by = c("subtype", "season_label"),
            suffix = c("", "_peak")) %>%
  mutate(
    days_since_peak = as.numeric(date - peak_date)  # peak day = 0
  )

p_Rt_shifted_2 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean, color = subtype)) +
  geom_line(size = 1) +
  coord_cartesian(ylim=c(0,2.5))+
 xlim(c(-50,50))+
  # xlim(c(-100,100))+
  theme_bw()+
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkred") +
  # Add text labels for peak Rt  
  geom_text(data = Rt_at_peak, aes(x = -30, y = 1.5, 
                                  label = paste0(  "Peak on ", peak_date, " \n with Rt= ",
                                                 round(Rt_mean, 2)," (",round(Rt_lower, 2),",",round(Rt_upper, 2),")", "") ), 
          inherit.aes = FALSE,    color = "darkred",     vjust = -0.5,   size = 3  ) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = subtype), alpha = 0.4, color = NA) +
  facet_wrap(~ season_label, ncol = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(title = "Rt estimates aligned to peak incidence",
       x = "Days since peak Rt", y = expression("Estimated " * R[t])) +
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())

p_Rt_shifted_2

 ggsave("figures/rt_estimates_plot/pre-Christmas_rt_estimates_by_season_1_shifted_by_peak_incidence.png",p_Rt_shifted_2,height =8 ,width = 10)

 
 ###version 2 to avoid overlapping labels
 
 # install.packages("ggrepel")
 library(ggrepel)
 
 Rt_at_peak <- Rt_at_peak %>%
   mutate(
     lab = paste0("Peak on ", peak_date,
                  "\nRt = ", round(Rt_mean, 2), " (",
                  round(Rt_lower, 2), ", ", round(Rt_upper, 2), ")"),
     x_lab = -30 + as.numeric(factor(subtype)) * 1.5,
     y_lab = 1.5 + as.numeric(factor(subtype)) * 0.1
   )
 
 p_Rt_shifted_2 <- ggplot(Rt_df_shifted, aes(x = days_since_peak, y = Rt_mean, color = subtype)) +
   geom_line(size = 1) +
   geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper, fill = subtype), alpha = 0.4, color = NA) +
   coord_cartesian(ylim = c(0, 2.5)) +
   xlim(c(-50, 50)) +
   theme_bw() +
   geom_vline(xintercept = 0, linetype = "dotted", color = "darkred") +
   geom_text_repel(data = Rt_at_peak,
                   aes(x = x_lab, y = y_lab, label = lab, color = NULL),
                   inherit.aes = FALSE,
                   size = 3, box.padding = 0.3, point.padding = 0.1, max.overlaps = Inf,
                   min.segment.length = 0, color = "darkred") +
   facet_wrap(~ season_label, ncol = 3) +
   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
   scale_color_brewer(palette = "Dark2") +
   scale_fill_brewer(palette = "Dark2") +
   labs(title = "Rt estimates aligned to peak incidence",
        x = "Days since peak Rt", y = expression("Estimated " * R[t])) +
   theme(legend.position = "bottom", legend.title = element_blank())
 
 p_Rt_shifted_2
 
 
 ggsave("figures/rt_estimates_plot/pre-Christmas_rt_estimates_by_season_1_shifted_by_peak_incidence_2.png",p_Rt_shifted_2,height =8 ,width = 10)
 
