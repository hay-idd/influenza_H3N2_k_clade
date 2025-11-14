library(EpiStrainDynamics)
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)
library(data.table)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")

devtools::load_all("~/Documents/GitHub/EpiStrainDynamics/")
add_doubling_axis <- function(p,n_breaks=11){
  # Get the y-axis range from the built plot
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  
  # Calculate appropriate breaks
  #n_breaks <- 5
  breaks <- pretty(y_range, n = n_breaks)
  
  # Create labels
  times <- log(2) / abs(breaks)
  signs <- ifelse(breaks >= 0, "", "-")
  labels <- ifelse(abs(breaks) < 1e-6, Inf,
                   paste0(signs, round(times, 1)))
  
  # Add secondary axis with calculated breaks and labels
  p + scale_y_continuous(
    breaks=breaks,
    sec.axis = sec_axis(~ .,
                        breaks = breaks,
                        labels = labels,
                        name = "Doubling(+) / Halving(-) time")
  )
  
}

theme_use <- theme_bw() + theme(legend.text=element_text(size=6),
                               legend.title=element_text(size=6),
                               axis.text=element_text(size=8),
                               axis.title=element_text(size=8),
                               strip.text=element_text(size=8))

## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#############################################
## MAIN ILIplus dataset
#############################################
overall_ILIplus <- read_csv("data/final/ili_plus_datasets_by_age.csv")

gr_all <- NULL
index <- 1

overall_ILIplus <- overall_ILIplus %>% select(Year,Week,date,group,ILIplus,Season)
overall_ILIplus <- overall_ILIplus %>% bind_rows(overall_ILIplus %>% group_by(Year,Week,date,Season) %>% summarize(ILIplus = sum(ILIplus)) %>% mutate(group="All"))
for(agegroup in unique(overall_ILIplus$group)){
  all_dat_tmp <- overall_ILIplus %>% filter(group == agegroup)
  mod <- construct_model(
    #method = random_walk(), 
    method=p_spline(spline_degree = 3, days_per_knot = 3),
    pathogen_structure =single(
      case_timeseries = round((all_dat_tmp$ILIplus)),           # timeseries of case data
      time = as.numeric(as.factor(all_dat_tmp$date)),                       # date or time variable labels
      pathogen_name = 'Influenza A/H3'                # optional name of pathogen
    ),
    dow_effect = FALSE
  )
  
  fit <- fit_model(
    mod,
    iter = 2000,
    warmup = 1000,
    chains = 3
  )
  
  gr <- growth_rate(fit)
  gr_all[[index]] <- plot(gr) + xlab("Week") + ylab('Growth rate (per week)')+
    ggtitle(agegroup)
  index <- index + 1
}



tmp_dat_all<-NULL
for(i in 1:length(unique(overall_ILIplus$group))){
  tmp_dat <- gr_all[[i]]$data
  tmp_dat$agegroup <- unique(overall_ILIplus$group)[i]
  tmp_dat_all[[i]] <- tmp_dat
}
tmp_dat_all <- do.call("bind_rows",tmp_dat_all)
ggplot(tmp_dat_all) + geom_line(aes(x=time,y=y,col=agegroup))
ymin <- -Inf
ymax <- Inf




school_periods_oxford <- tribble(
  ~start,        ~end,         ~label,           ~acad_year,                ~type,
  
  # 2022/23 holidays
  "2022-10-24", "2022-10-28",  "Half-term Oct 2022",             "2022/23",   "holiday",
  "2022-12-19", "2023-01-02",  "Christmas Holiday 2022/23",      "2022/23",   "holiday",
  "2023-02-13", "2023-02-17",  "Half-term Feb 2023",             "2022/23",   "holiday",
  "2023-04-07", "2023-04-21",  "Easter Holiday 2023",            "2022/23",   "holiday",
  "2023-05-29", "2023-06-02",  "Late Spring Half-term 2023",     "2022/23",   "holiday",
  "2023-07-22", "2023-08-31",  "Summer Holiday 2023",            "2022/23",   "holiday",
  
  # 2023/24 holidays
  "2023-10-23", "2023-10-27",  "Half-term Oct 2023",             "2023/24",   "holiday",
  "2023-12-21", "2024-01-02",  "Christmas Holiday 2023/24",      "2023/24",   "holiday",
  "2024-02-12", "2024-02-16",  "Half-term Feb 2024",             "2023/24",   "holiday",
  "2024-03-29", "2024-04-12",  "Easter Holiday 2024",            "2023/24",   "holiday",
  "2024-05-27", "2024-05-31",  "Late Spring Half-term 2024",     "2023/24",   "holiday",
  "2024-07-20", "2024-08-31",  "Summer Holiday 2024",            "2023/24",   "holiday",
  
  # 2024/25 holidays
  "2024-10-25", "2024-10-29",  "Half-term Oct 2024",             "2024/25",   "holiday",
  "2024-12-20", "2025-01-03",  "Christmas Holiday 2024/25",      "2024/25",   "holiday",
  "2025-02-17", "2025-02-21",  "Half-term Feb 2025",             "2024/25",   "holiday",
  "2025-04-11", "2025-04-22",  "Easter Holiday 2025",            "2024/25",   "holiday",
  "2025-05-26", "2025-05-30",  "Late Spring Half-term 2025",     "2024/25",   "holiday",
  "2025-07-21", "2025-08-31",  "Summer Holiday 2025",            "2024/25",   "holiday",
  
  # 2025/26 holidays
  "2025-10-24", "2025-10-30",  "Half-term Oct 2025",             "2025/26",   "holiday",
  
  
) %>% mutate(start = as.Date(start), end = as.Date(end))
rects <- school_periods_oxford %>%
  mutate(xmin = start - 0, xmax = end + 1)

time_key <- overall_ILIplus %>% select(Year, Week, date,Season) %>%
  mutate(time = (Year-2023)*52 + Week) %>%
  mutate(time = as.numeric(as.factor(time))) %>%
  group_by(Season) %>%
  mutate(time_rel = min(Week) + (Week - min(Week))) %>%
  ungroup()

tmp_dat_all$agegroup <- factor(tmp_dat_all$agegroup, levels=c("All","1-4","5-14","15-44","45-64","65+"))


raw_data <- overall_ILIplus %>% group_by(group) %>% 
  mutate(gr = log(ILIplus/lag(ILIplus,1))) %>%
  mutate(gr = zoo::rollmean(gr,3,fill=NA,align='right')) %>%
  rename(agegroup=group)
raw_data$agegroup <- factor(raw_data$agegroup, levels=c("All","1-4","5-14","15-44","45-64","65+"))


tmp_dat_all <- tmp_dat_all %>% left_join(time_key) %>% distinct()
## Plot overall GR and by age
p_gr_by_age <- ggplot(tmp_dat_all) + 
  geom_rect(data = rects %>% filter(!(label %like% "Half-term Oct")), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="grey70",
            alpha = 0.4) +  
  geom_rect(data = rects %>% filter((label %like% "Half-term Oct")), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="red",
            alpha = 0.4) +  
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=date,y=y,col=agegroup),linewidth=0.75) +
  geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set2") +
  scale_color_brewer("Age group", palette = "Set2") +
  xlab("Week") + ylab('Growth rate (per week)') +
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1.2,1.2,by=0.2)) +
  coord_cartesian(ylim=c(-1.2,1.2)) +
  theme(legend.position="bottom",legend.direction="horizontal")
#p_gr_by_age <- add_doubling_axis(p_gr_by_age)


## Plot overall GR and by age, aligned by season

## Plot overall GR and by age
p_gr_by_age_shifted <- ggplot(tmp_dat_all%>% filter(Season != "2022 to 2023")) + 
  #geom_rect(data = rects, inherit.aes = FALSE,
    #        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="grey70",
   #         alpha = 0.4) +  
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=time_rel,ymin=lb_95,ymax=ub_95,fill=Season),alpha=0.5) +
  geom_line(aes(x=time_rel,y=y,col=Season),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set1") +
  scale_color_brewer("Age group", palette = "Set1") +
  xlab("Week (by season)") + ylab('Growth rate (per week)') +
  theme_use + facet_wrap(~agegroup,ncol=2) +
  #scale_y_continuous(breaks=seq(-1.2,1.2,by=0.2)) +
  #coord_cartesian(ylim=c(-1.2,1.2)) +
  theme(legend.position="bottom",legend.direction="horizontal")
#p_gr_by_age_shifted <- add_doubling_axis(p_gr_by_age_shifted)

## Also shift to line up by maximum growth rate
tmp_dat_shifted <-  tmp_dat_all %>% 
  filter(Season != "2022 to 2023") %>%
  left_join(
  tmp_dat_all %>% left_join(time_key) %>% distinct() %>% filter(Season != "2022 to 2023") %>%
  group_by(Season,agegroup) %>% filter(y == max(y)) %>% select(time_rel,Season,agegroup) %>% distinct() %>% rename(max_week = time_rel)
  ) %>%
  mutate(time_plot = time_rel - max_week)


p_gr_by_age_shifted_rel <- ggplot(tmp_dat_shifted) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=time_plot,ymin=lb_95,ymax=ub_95,fill=Season),alpha=0.5) +
  geom_line(aes(x=time_plot,y=y,col=Season),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set1") +
  scale_color_brewer("Age group", palette = "Set1") +
  xlab("Week (by season)") + ylab('Growth rate (per week)') +
  theme_use + facet_wrap(~agegroup,ncol=2) +
  #scale_y_continuous(breaks=seq(-1.2,1.2,by=0.2)) +
  #coord_cartesian(ylim=c(-1.2,1.2)) +
  theme(legend.position="bottom",legend.direction="horizontal")
#p_gr_by_age_shifted_rel <- add_doubling_axis(p_gr_by_age_shifted_rel)



#############################################
## SECONDARY FluNet dataset
#############################################
flunet <- read_csv("data/final/flunet_h3_cases_historic.csv")
flunet <- flunet %>% mutate(H3_sum = if_else(is.na(H3_sum),0,H3_sum))

mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    case_timeseries = round((flunet$H3_sum)),           # timeseries of case data
    time = as.numeric(as.factor(flunet$date)),                       # date or time variable labels
    pathogen_name = 'Influenza A/H3'                # optional name of pathogen
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  iter = 2000,
  warmup = 1000,
  chains = 3
)

gr <- growth_rate(fit)
p_flunet <- plot(gr)
gr_flunet_dat <- p_flunet$data

time_key_flunet <- flunet %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_flunet_dat <- gr_flunet_dat %>% left_join(time_key_flunet)

## For later years, just assign season manually
gr_flunet_dat <- gr_flunet_dat %>%
  mutate(season = case_when(
    date >= as.Date("2011-07-01") & date < as.Date("2012-07-01") ~ "2011/12",
    date >= as.Date("2012-07-01") & date < as.Date("2013-07-01") ~ "2012/13",
    date >= as.Date("2013-07-01") & date < as.Date("2014-07-01") ~ "2013/14",
    date >= as.Date("2014-07-01") & date < as.Date("2015-07-01") ~ "2014/15",
    date >= as.Date("2015-07-01") & date < as.Date("2016-07-01") ~ "2015/16",
    date >= as.Date("2016-07-01") & date < as.Date("2017-07-01") ~ "2016/17",
    date >= as.Date("2017-07-01") & date < as.Date("2018-07-01") ~ "2017/18",
    date >= as.Date("2018-07-01") & date < as.Date("2019-07-01") ~ "2018/19",
    date >= as.Date("2019-07-01") & date < as.Date("2020-07-01") ~ "2019/20",
    date >= as.Date("2020-07-01") & date < as.Date("2021-07-01") ~ "2020/21",
    date >= as.Date("2021-07-01") & date < as.Date("2022-07-01") ~ "2021/22",
    date >= as.Date("2022-07-01") & date < as.Date("2023-07-01") ~ "2022/23",
    date >= as.Date("2023-07-01") & date < as.Date("2024-07-01") ~ "2023/24",
    date >= as.Date("2024-07-01") & date < as.Date("2025-07-01") ~ "2024/25",
    date >= as.Date("2025-07-01") & date < as.Date("2026-07-01") ~ "2025/26",
    TRUE ~ NA
  ))



p_gr_flunet <- ggplot(data=gr_flunet_dat) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
  scale_y_continuous(limits=c(-1.2,1.2),breaks=seq(-1.2,1.2,by=0.2)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")


## Align by day of year
gr_flunet_dat$day_of_year <- yday(gr_flunet_dat$date)
gr_flunet_dat <- gr_flunet_dat %>% group_by(season) %>%
  mutate(day_of_year = if_else(day_of_year < 183, day_of_year + 365, day_of_year)) %>%
  ungroup()

## Create labels for plot, label seasons by pre-pandemic, during pandemic (2020 and 2021), post pandemic (2022-2024) and current season (2025/26)
gr_flunet_dat <- gr_flunet_dat %>%
  mutate(plot_label = case_when(
    season %in% c("2012/13","2013/14","2014/15","2015/16","2016/17","2017/18","2018/19","2019/20") ~ "Pre-pandemic",
    season %in% c("2020/21","2021/22") ~ "During pandemic",
    season %in% c("2022/23","2023/24","2024/25") ~ "Post-pandemic",
    season == "2025/26" ~ "Current season",
    TRUE ~ "Other"
  ))

p_gr_flunet_by_day<- ggplot(data=gr_flunet_dat) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
 # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label)) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_viridis_d("Time period") +
  scale_fill_viridis_d("Time period") +
  scale_y_continuous(limits=c(-1.2,1.2),breaks=seq(-1.2,1.2,by=0.2)) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")


## Align by peak data
max_gr <- gr_flunet_dat %>% group_by(season) %>% filter(y == max(y)) %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)
gr_flunet_dat <- gr_flunet_dat %>% left_join(max_gr)
gr_flunet_dat$day_shifted <- gr_flunet_dat$day_of_year - gr_flunet_dat$peak_time

p_gr_flunet_by_peak <- ggplot(data=gr_flunet_dat %>% filter(day_shifted <= 100, day_shifted >= -100)) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label)) +
  scale_y_continuous(limits=c(-1.2,1.2),breaks=seq(-1.2,1.2,by=0.2)) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day relative to peak") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_viridis_d("Time period") +
  scale_fill_viridis_d("Time period") +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_flunet <- add_doubling_axis(p_gr_flunet)
p_gr_flunet_by_day <- add_doubling_axis(p_gr_flunet_by_day)
p_gr_flunet_by_peak <- add_doubling_axis(p_gr_flunet_by_peak)

## Look at total number of H3 cases by peak time
gr_flunet_dat %>% left_join(flunet) %>% filter(day_shifted <= 0) %>% group_by(season) %>% summarize(total_cases = sum(H3_sum))

ggsave("figures/growth_rates_by_age.png",p_gr_by_age,width=8,height=8)
ggsave("figures/growth_rates_by_age_shifted.png",p_gr_by_age_shifted,width=8,height=8)
ggsave("figures/growth_rates_by_age_shifted_rel.png",p_gr_by_age_shifted_rel,width=8,height=8)

ggsave("figures/growth_rates_flunet.png",p_gr_flunet,width=8,height=4)
ggsave("figures/growth_rates_flunet_by_day.png",p_gr_flunet_by_day,width=6,height=4)
ggsave("figures/growth_rates_flunet_by_peak.png",p_gr_flunet_by_peak,width=6,height=4)
