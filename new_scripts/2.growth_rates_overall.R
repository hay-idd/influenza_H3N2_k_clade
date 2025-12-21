#library(EpiStrainDynamics)
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)
library(data.table)
library(mgcv)
library(ggrepel)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
source("R/funcs.R")

devtools::load_all("~/Documents/GitHub/EpiStrainDynamics/")

theme_use <- theme_bw() + theme(legend.text=element_text(size=6),
                               legend.title=element_text(size=6),
                               axis.text=element_text(size=8),
                               axis.title=element_text(size=8),
                               strip.text=element_text(size=8))

## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

##########################################################################################
## Plot all raw data, raw growth rates and GAM-smoothed growth rates for sanity check
##########################################################################################
## ILIplus
overall_ILIplus <- read_csv("new_data/ili_plus_datasets_by_age.csv")
overall_ILIplus <- overall_ILIplus %>% select(Year,Week,date,group,`Influenza A/H3N2`,`All influenza cases`,Season) %>%
  rename(ILIplus=`Influenza A/H3N2`,
         Influenza_cases=`All influenza cases`)
overall_ILIplus <- overall_ILIplus %>% bind_rows(overall_ILIplus %>% group_by(Year,Week,date,Season) %>% summarize(ILIplus = sum(ILIplus),Influenza_cases=sum(Influenza_cases)) %>% mutate(group="All"))
overall_ILIplus$group <- factor(overall_ILIplus$group, levels=c("All","0-4","5-18","19-64","65+"))
overall_ILIplus <- overall_ILIplus %>% 
  select(Season,Year,Week,date,group,ILIplus,Influenza_cases) %>% 
  mutate("Source"="RCGP") %>%
  rename(cases_H3=ILIplus,
         cases_total=Influenza_cases)
overall_ILIplus$Season <- season_from_date(overall_ILIplus$date)

## England case counts
#case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts <- read_csv("new_data/resp_datamart_influenza_cases_england.csv")
case_counts$date <- lubridate::ymd(case_counts$date)
case_counts <- case_counts %>% rename(cases_H3 = H3_adj) %>%
  rename(cases_total = total_cases) %>%
   mutate(Source = "UKHSA")

## WHO FluNet
flunet <- read_csv("new_data/WHO_flunet_cases.csv")
flunet <- flunet %>% mutate(H3_sum = if_else(is.na(H3_sum),0,H3_sum)) %>%
  rename(cases_H3=H3_sum, cases_total=flu_pos_sum) %>%
  mutate(Source="FluNet")
flunet$Season <- season_from_date(flunet$date)

## Plot H3 growth rates
case_counts <- calc_emp_growth(case_counts,"date","cases_H3",group2="Season",eps=0.5)
overall_ILIplus <- calc_emp_growth(overall_ILIplus,"date","cases_H3",group1="group",group2="Season",eps=0.5)
flunet <- calc_emp_growth(flunet,"date","cases_H3",group2="Season",eps=0.5)

all_grs_H3 <- bind_rows(overall_ILIplus%>% filter(group=="All")%>% select(Year,Week,Season,date,cases_H3,gr,gr_smooth,Source) ,
                         case_counts %>% select(Year,Week,Season,date,cases_H3,gr,gr_smooth,Source),
                         flunet %>% select(Year,Week,Season,date,cases_H3,gr,gr_smooth,Source))

fit_gams_H3 <- fit_gam_grs(all_grs_H3 %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_H3",time="date_int",group="Source")[[1]]

p_h3_gr_all <- ggplot(fit_gams_H3) + 
  geom_line(aes(x=date,y=gr,col=Source,group=interaction(Source,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=Source,group=interaction(Source,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=Source,group=interaction(Source,Season),alpha="GAM")) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  coord_cartesian(ylim=c(-1,1)) +
  facet_wrap(~Source,ncol=1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use

p_h3_gr_recent_compare <-  ggplot(fit_gams_H3 %>% filter(Year >= 2022)) + 
  geom_line(aes(x=date,y=gr_smooth,col=Source,group=interaction(Source,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=Source,group=interaction(Source,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=Source,group=interaction(Source,Season),alpha="GAM")) +
  coord_cartesian(ylim=c(-1,1)) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use
## Plot overall flu growth rates
case_counts <- calc_emp_growth(case_counts,"date","cases_total")
overall_ILIplus <- calc_emp_growth(overall_ILIplus,"date","cases_total",group1="group")
flunet <- calc_emp_growth(flunet,"date","cases_total")

all_grs_cases <- bind_rows(overall_ILIplus%>% filter(group=="All")%>% select(Year,Week,Season,date,
                                                                             cases_total,gr,gr_smooth,Source) ,
                     case_counts %>% select(Year,Week,Season,date,cases_total,gr,gr_smooth,Source),
                     flunet %>% select(Year,Week,Season,date,cases_total,gr,gr_smooth,Source))

fit_gams_all <- fit_gam_grs(all_grs_cases %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_total",time="date_int",group="Source")[[1]]

p_gr_all <- ggplot(fit_gams_all) + 
  geom_line(aes(x=date,y=gr,col=Source,group=interaction(Source,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=Source,group=interaction(Source,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=Source,group=interaction(Source,Season),alpha="GAM")) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  coord_cartesian(ylim=c(-1,1)) +
  facet_wrap(~Source,ncol=1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use



p_gr_recent_compare <-  ggplot(fit_gams_all %>% filter(Year >= 2022)) + 
  geom_line(aes(x=date,y=gr_smooth,col=Source,group=interaction(Source,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=Source,group=interaction(Source,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=Source,group=interaction(Source,Season),alpha="GAM")) +
  coord_cartesian(ylim=c(-1,1)) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use

## Plot by age group
rcgp_h3_gr <-  calc_emp_growth(overall_ILIplus,"date","cases_H3",group1="group")
rcgp_h3_gr <- fit_gam_grs(rcgp_h3_gr %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_H3",time="date_int",group="group")[[1]]


## Find difference between peak growth rates in each season for each age group and the 19-64 age group
peak_grs_h3_rcgp <- rcgp_h3_gr %>% filter(group == "19-64", Source=="RCGP") %>% group_by(Season) %>% filter(r_week == max(r_week)) %>%
  select(r_week,Season) %>% rename(peak_gr_adults=r_week)

rcgp_h3_gr %>% filter(Source=="RCGP") %>% group_by(Season, group) %>% filter(r_week==max(r_week)) %>% left_join(peak_grs_h3_rcgp) %>%
  mutate(diff = r_week - peak_gr_adults) %>% select(Season,group,diff)

p_gr_h3_age <- ggplot(rcgp_h3_gr) + 
  geom_line(aes(x=date,y=gr,col=group,group=interaction(group,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=group,group=interaction(group,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=group,group=interaction(group,Season),alpha="GAM")) +
  coord_cartesian(ylim=c(-1,1)) +
  facet_wrap(~group,ncol=2) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  scale_color_brewer("Age group",palette = "Set2") +
  scale_fill_brewer("Age group",palette = "Set2") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use

rcgp_all_gr <-  calc_emp_growth(overall_ILIplus,"date","cases_total",group1="group")
rcgp_all_gr <- fit_gam_grs(rcgp_all_gr %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_total",time="date_int",group="group")[[1]]

## Find difference between peak growth rates in each season for each age group and the 19-64 age group
peak_grs_rcgp <- rcgp_all_gr %>% filter(group == "19-64", Source=="RCGP") %>% group_by(Season) %>% filter(r_week == max(r_week)) %>%
  select(r_week,Season) %>% rename(peak_gr_adults=r_week)

diff_peak_gr <- rcgp_all_gr %>% filter(Source=="RCGP") %>% group_by(Season, group) %>% filter(r_week==max(r_week)) %>% left_join(peak_grs_rcgp) %>%
  mutate(diff = r_week - peak_gr_adults) %>% select(Season,group,diff) %>% pivot_wider(names_from=group,values_from=diff) %>%
  select(-All)
write.csv(diff_peak_gr,"new_figures/Table2.csv")


p_gr_all_age <- ggplot(rcgp_all_gr) + 
  geom_line(aes(x=date,y=gr,col=group,group=interaction(group,Season),alpha="Raw data"),linewidth=0.25) +
  geom_ribbon(aes(x=date,ymin=lower_r,ymax=upper_r,fill=group,group=interaction(group,Season)),alpha=0.5) +
  geom_line(aes(x=date,y=r_week,col=group,group=interaction(group,Season),alpha="GAM")) +
  coord_cartesian(ylim=c(-1,1)) +
  facet_wrap(~group,ncol=2) +
  scale_alpha_manual(values=c("Raw data"=0.25,"GAM"=1),name="Growth rate type")+
  scale_color_brewer("Age group",palette = "Set2") +
  scale_fill_brewer("Age group",palette = "Set2") +
  ylab("Weekly growth rate") +
  xlab("Date (start of Epi week)") +
  theme_use

## Look at the ratio of growth rates by age to overall
rcgp_all_gr %>% filter(group != "All") %>%
  left_join(rcgp_all_gr %>% filter(group == "All") %>% select(date,gr) %>% rename(gr_all=gr)) %>%
  mutate(gr_ratio = gr/gr_all) %>%
  ggplot() + 
  geom_line(aes(x=date,y=gr_ratio,col=group,group=interaction(group,Season))) +
  geom_hline(yintercept=1,linetype="dashed") +
  ylab("Ratio of age-specific growth rate to overall growth rate") +
  xlab("Date (start of Epi week)") +
  scale_color_brewer("Age group",palette = "Set2") +
  theme_use +
  facet_wrap(~group)

rcgp_h3_gr %>% filter(group != "19-64") %>%
  left_join(rcgp_h3_gr %>% filter(group == "19-64") %>% select(date,gr) %>% rename(gr_all=gr)) %>%
  mutate(gr_ratio = log(gr/gr_all)) %>%
  ggplot() + 
  geom_line(aes(x=date,y=gr_ratio,col=group,group=interaction(group,Season))) +
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Ratio of age-specific growth rate to overall growth rate") +
  xlab("Date (start of Epi week)") +
  scale_color_brewer("Age group",palette = "Set2") +
  theme_use +
  facet_wrap(~group)


## Incidence plots
p_inc_h3 <- plot_incidence(all_grs_H3,x="date",y="cases_H3",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed\n or suspected A/H3N2") + theme_use

p_inc_h3_recent <- plot_incidence(all_grs_H3 %>% filter(Year >= 2021),x="date",y="cases_H3",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use

p_inc_all <- plot_incidence(all_grs_cases,x="date",y="cases_total",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed or suspected influenza") + theme_use

p_inc_all_recent <- plot_incidence(all_grs_cases%>% filter(Year >= 2021),x="date",y="cases_total",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed or suspected influenza") + theme_use

p_inc_h3_age <- plot_incidence(overall_ILIplus,x="date",y="cases_H3",color="group",point=FALSE)+ 
  scale_color_brewer("Age group",palette = "Set2") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use
p_inc_all_age <- plot_incidence(overall_ILIplus,x="date",y="cases_total",color="group",point=FALSE)+ 
  scale_color_brewer("Age group",palette = "Set2") + xlab("Date (start of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use

## Save all plots
p_h3_gr_all
p_h3_gr_recent_compare
p_gr_all
p_gr_recent_compare
p_gr_h3_age
p_gr_all_age
p_inc_h3
p_inc_h3_recent
p_inc_all
p_inc_all_recent
p_inc_h3_age
p_inc_all_age

#write.csv(fit_gams_H3,file="results/fit_gams_H3.csv")
#write.csv(fit_gams_all,"results/fit_gams_all.csv")
#write.csv(rcgp_h3_gr,"results/rcgp_h3_gr.csv")
#write.csv(fit_gams_H3,"results/rcgp_all_gr.csv")
#ggsave("figures/raw_data/p_h3_gr_all.png",p_h3_gr_all,width=7,height=6)
#ggsave("figures/raw_data/p_h3_gr_recent.png",p_h3_gr_recent_compare,width=7,height=4)
#ggsave("figures/raw_data/p_gr_all.png",p_gr_all,width=7,height=6)
#ggsave("figures/raw_data/p_gr_recent_compare.png",p_gr_recent_compare,width=7,height=4)
#ggsave("figures/raw_data/p_gr_h3_age.png",p_gr_h3_age,width=8,height=7)
ggsave("new_figures/figS6.png",p_gr_all_age,width=8,height=6)
#ggsave("figures/raw_data/p_inc_h3.png",p_inc_h3,width=7,height=4)
#ggsave("figures/raw_data/p_inc_h3.pdf",p_inc_h3,width=7,height=4)
#ggsave("figures/raw_data/p_inc_h3.pdf",p_inc_h3 + theme(legend.position=c(0.2,0.7),
#axis.text = element_text(size=12),
#axis.title = element_text(size=14),
#strip.text = element_text(size=14),
#title = element_text(size=16)),width=5,height=3)
#ggsave("figures/raw_data/p_inc_h3_recent.png",p_inc_h3_recent,width=7,height=4)
#ggsave("figures/raw_data/p_inc_all.png",p_inc_all_recent,width=7,height=4)
#ggsave("figures/raw_data/p_inc_all_recent.png",p_inc_all,width=7,height=4)
#ggsave("figures/raw_data/p_inc_h3_age.png",p_inc_h3_age,width=7,height=4)
#ggsave("figures/raw_data/p_inc_all_age.png",p_inc_all_age,width=7,height=4)

#############################################
## MAIN ILIplus dataset
#############################################
if(FALSE){
overall_ILIplus <- read_csv("new_data/ili_plus_datasets_by_age.csv")
overall_ILIplus <- overall_ILIplus %>% select(Year,Week,date,group,`Influenza A/H3N2`,`All influenza cases`,Season)
overall_ILIplus <- overall_ILIplus %>% bind_rows(overall_ILIplus %>% group_by(Year,Week,date,Season) %>% summarize(ILIplus = sum(`Influenza A/H3N2`)) %>% mutate(group="All"))
overall_ILIplus <- overall_ILIplus %>% filter(group == "All")%>% filter(date >= "2022-09-01")

overall_ILIplus$cases <- round(overall_ILIplus$ILIplus)
overall_ILIplus$time <- as.numeric(as.factor(overall_ILIplus$date))
all_dat_tmp <- overall_ILIplus %>% select(cases, time) 
mod <- construct_model(
  method = random_walk(), 
  #method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    data=all_dat_tmp,
    case_timeseries = "cases",           # timeseries of case data
    time = "time"                      # date or time variable labels
  ), 
  smoothing_params = smoothing_structure(   # independent smoothing structure 
    tau_mean = c(0),           # parameter - one for each pathogen
    tau_sd = c(0.01)
  ),
  
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  n_iter = 2000,
  n_warmup = 1000,
  n_chain = 3
)

gr <- growth_rate(fit)

time_key <- overall_ILIplus %>% select(Year, Week, date,Season) %>%
  mutate(time = (Year-2023)*52 + Week) %>%
  mutate(time = as.numeric(as.factor(time))) %>%
  group_by(Season) %>%
  mutate(time_rel = time - min(time)) %>%
  ungroup()


raw_data <- overall_ILIplus %>% group_by(group) %>% 
  mutate(gr = log((ILIplus+0.1)/lag(ILIplus+0.1,1))) %>%
  #mutate(gr = zoo::rollmean(gr,3,fill=NA,align='right')) %>%
  rename(agegroup=group)
raw_data$agegroup <- factor(raw_data$agegroup, levels=c("All","0-4","5-14","15-44","45-64","65+"))


tmp_dat_all <- raw_data %>% left_join(time_key) %>% distinct()
## Plot overall GR and by age
p_gr_by_age <- ggplot(tmp_dat_all) + 
  geom_rect(data = rects %>% filter(!(label %like% "Half-term Oct" & acad_year == "2025/26")), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="grey70",
            alpha = 0.4) +  
  geom_rect(data = rects %>% filter((label %like% "Half-term Oct" & acad_year == "2025/26")), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="red",
            alpha = 0.4) +  
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=date,y=y,col=agegroup,alpha="Random walk"),linewidth=0.75) +
  geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup,alpha="Raw data")) +
  scale_fill_brewer("Age group", palette = "Set2") +
  scale_color_brewer("Age group", palette = "Set2") +
  xlab("Week") + ylab('Growth rate (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"Random walk"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1.2,1.2)) +
  theme(legend.position="bottom",legend.direction="horizontal")
}
#############################################
## SECONDARY FluNet dataset fit to H3
#############################################
flunet <- read_csv("new_data/WHO_flunet_cases.csv")
flunet <- flunet %>% mutate(H3_sum = if_else(is.na(H3_sum),0,H3_sum))


dat_tmp <- flunet %>% filter(Year >= 2011)
dat_tmp$H3_sum <- round(dat_tmp$H3_sum)
dat_tmp$time <- as.numeric(as.factor(dat_tmp$date))
mod <- construct_model(
  method = random_walk(),
  #method=p_spline(spline_degree = 3, days_per_knot = 3), ## I report 3 in the paper, but this was originally 4 if you see things change a lot
  pathogen_structure =single(
    data=dat_tmp,
    case_timeseries = "H3_sum",           # timeseries of case data
    time = "time"                       # date or time variable labels
  ),
  smoothing_params = smoothing_structure(   # independent smoothing structure 
    tau_mean = c(1),           # parameter - one for each pathogen
    tau_sd = c(0.05)
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  n_iter = 6000,
  n_warmup = 3000,
  n_chain = 3
)

gr <- growth_rate(fit)
p_flunet <- plot(gr)
gr_flunet_dat <- p_flunet$data

time_key_flunet <- flunet %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_flunet_dat <- gr_flunet_dat %>% left_join(time_key_flunet)

## Convert dates to flu season
gr_flunet_dat <- gr_flunet_dat %>%
  dplyr::mutate(season = flu_season(date))


p_gr_flunet_h3 <- ggplot(data=gr_flunet_dat) + 
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

max_t <- max(gr_flunet_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
min_t <- min(gr_flunet_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
breaks <- seq(150,max_t,by=50)
labels <- if_else(breaks > 365, breaks - 365, breaks)
p_gr_flunet_h3_by_day<- ggplot(data=gr_flunet_dat%>% filter(plot_label != "During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_vline(xintercept=365,linetype="dashed") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label),linewidth=0.75) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  coord_cartesian(ylim=c(-1.1,1.1)) +  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")


## Align by peak data
## Filter gr_flunet_dat to be between September and May each season
gr_flunet_dat1 <- gr_flunet_dat %>%
  filter(month(date) >= 9 | month(date) <= 4 | (month(date) == 5 & day(date) <= 1)) %>%
  group_by(season) %>%
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup()

max_gr <- gr_flunet_dat1 %>% group_by(season) %>% filter(y == max(y))# %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)

max_gr %>% select(season, Year, Week, date, y, lb_95, ub_95) %>% arrange(season) %>% tail(5) %>% print()


gr_flunet_dat <- gr_flunet_dat %>% left_join(max_gr %>% select(season, day_of_year) %>% rename(peak_time = day_of_year))
gr_flunet_dat$day_shifted <- gr_flunet_dat$day_of_year - gr_flunet_dat$peak_time

## Get max growth rate and week of max growth rate
## Merge in raw data
gr_flunet_dat <-gr_flunet_dat %>% left_join(flunet %>% mutate(gr = log((H3_sum+0.1)/lag(H3_sum+0.1,1)), gr_smooth = zoo::rollmean(gr, 4, fill=NA,align="right"))) 

p_gr_flunet_h3_by_peak <- ggplot(data=gr_flunet_dat %>% filter(day_shifted <= 100, day_shifted >= -100)) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label),linewidth=0.75) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1.2,1.2,by=0.2)) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day relative to peak") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_flunet_h3_by_day_faceted <- ggplot(data=gr_flunet_dat) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_line(aes(x=day_of_year,y=gr,col="Raw data")) +
  geom_line(aes(x=day_of_year,y=gr_smooth,group=season,col="Smoothed data")) +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,col="Estimate"),alpha=0.1) +
  geom_line(aes(x=day_of_year,y=y,group=season,col="Estimate")) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  #scale_y_continuous(limits=c(-1.2,1.2),breaks=seq(-1.2,1.2,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme_use + 
  facet_wrap(~season) +
  theme(legend.position="bottom",legend.direction="horizontal")


p_gr_flunet_h3 <- add_doubling_axis(p_gr_flunet_h3)
p_gr_flunet_h3_by_day <- add_doubling_axis(p_gr_flunet_h3_by_day)
p_gr_flunet_h3_by_peak <- add_doubling_axis(p_gr_flunet_h3_by_peak)

## Look at total number of H3 cases by peak time
gr_flunet_dat %>% left_join(flunet) %>% filter(day_shifted <= 0) %>% group_by(season) %>% summarize(total_cases = sum(H3_sum))


ggsave("figures/growth_rates/flunet_h3_growth_rate.png",p_gr_flunet_h3,width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_day.png",p_gr_flunet_h3_by_day,width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_day.pdf",p_gr_flunet_h3_by_day+
         theme(axis.text = element_text(size=12),
               axis.title = element_text(size=14),
               strip.text = element_text(size=14),
               legend.title=element_text(size=10),
               legend.text=element_text(size=10),
               title = element_text(size=16)),width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_peak.png",p_gr_flunet_h3_by_peak,width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_day_faceted_ps.png",p_gr_flunet_h3_by_day_faceted,width=7,height=7)

ggsave("figures/growth_rates/flunet_h3_growth_rate_comb.png",p_gr_flunet_h3_by_day/p_gr_flunet_h3_by_peak
,width=8,height=10)

write_csv(gr_flunet_dat,"results/flunet_h3_growth_rates.csv")


#############################################
## THIRD FluNet dataset fit to all cases
#############################################
flunet <- read_csv("new_data/WHO_flunet_cases.csv")
flunet <- flunet %>% mutate(flu_pos_sum = if_else(is.na(flu_pos_sum),0,flu_pos_sum))

dat_tmp <- flunet
dat_tmp$flu_pos_sum <- round(dat_tmp$flu_pos_sum)
dat_tmp$time <- as.numeric(as.factor(dat_tmp$date))

mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    data=dat_tmp,
    case_timeseries = "flu_pos_sum",           # timeseries of case data
    time = "time"                       # date or time variable labels
  ),
  smoothing_params = smoothing_structure(   # independent smoothing structure 
    tau_mean = c(1),           # parameter - one for each pathogen
    tau_sd = c(0.25)
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  n_iter = 2000,
  n_warmup = 1000,
  n_chain = 3
)

gr <- growth_rate(fit)
p_flunet <- plot(gr)
gr_flunet_dat <- p_flunet$data

time_key_flunet <- flunet %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_flunet_dat <- gr_flunet_dat %>% left_join(time_key_flunet)

gr_flunet_dat <- gr_flunet_dat %>%
  dplyr::mutate(season = flu_season(date))

p_gr_flunet_all <- ggplot(data=gr_flunet_dat) + 
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

max_t <- max(gr_flunet_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
min_t <- min(gr_flunet_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
breaks <- seq(150,max_t,by=50)
labels <- if_else(breaks > 365, breaks - 365, breaks)
p_gr_flunet_by_day<- ggplot(data=gr_flunet_dat%>% filter(plot_label != "During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_vline(xintercept=365,linetype="dashed") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label),linewidth=0.75) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  coord_cartesian(ylim=c(-1.1,1.1)) +  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")


## Align by peak data
## Filter gr_flunet_dat to be between September and May each season
gr_flunet_dat1 <- gr_flunet_dat %>%
  filter(month(date) >= 9 | month(date) <= 4 | (month(date) == 5 & day(date) <= 1)) %>%
  group_by(season) %>%
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup()

max_gr <- gr_flunet_dat1 %>% group_by(season) %>% filter(y == max(y))# %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)


gr_flunet_dat <- gr_flunet_dat %>% left_join(max_gr %>% select(season, day_of_year) %>% rename(peak_time = day_of_year))
gr_flunet_dat$day_shifted <- gr_flunet_dat$day_of_year - gr_flunet_dat$peak_time

max_gr %>% select(season, Year, Week, date, y, lb_95, ub_95) %>% arrange(season) %>% tail(5) %>% print()

## Merge in raw data
gr_flunet_dat <-gr_flunet_dat %>% left_join(flunet %>% mutate(gr = log((H3_sum+0.1)/lag(H3_sum+0.1,1)), gr_smooth = zoo::rollmean(gr, 4, fill=NA,align="right"))) 

p_gr_flunet_by_peak <- ggplot(data=gr_flunet_dat %>% filter(day_shifted <= 100, day_shifted >= -100)%>% filter(plot_label != "During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label),linewidth=0.75) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day relative to peak") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_flunet_by_day_faceted <- ggplot(data=gr_flunet_dat%>% filter(plot_label != "During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_line(aes(x=day_of_year,y=gr,col="Raw data")) +
  geom_line(aes(x=day_of_year,y=gr_smooth,group=season,col="Smoothed data")) +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,col="Estimate"),alpha=0.1) +
  geom_line(aes(x=day_of_year,y=y,group=season,col="Estimate")) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  #scale_y_continuous(limits=c(-1.2,1.2),breaks=seq(-1.2,1.2,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme_use + 
  facet_wrap(~season) +
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_flunet_all <- add_doubling_axis(p_gr_flunet_all)
p_gr_flunet_by_day <- add_doubling_axis(p_gr_flunet_by_day)
p_gr_flunet_by_peak <- add_doubling_axis(p_gr_flunet_by_peak)

ggsave("figures/growth_rates/flunet_all_growth_rate.png",p_gr_flunet_all,width=7,height=4)
ggsave("figures/growth_rates/flunet_all_growth_rate_by_day.png",p_gr_flunet_by_day,width=7,height=4)
ggsave("figures/growth_rates/flunet_all_growth_rate_by_peak.png",p_gr_flunet_by_peak,width=7,height=4)

ggsave("figures/growth_rates/flunet_all_growth_rate_by_day_faceted_ps.png",p_gr_flunet_by_day_faceted,width=7,height=7)

ggsave("figures/growth_rates/flunet_all_growth_rate_comb.png",p_gr_flunet_by_day/p_gr_flunet_by_peak
       ,width=8,height=10)

write_csv(gr_flunet_dat,"results/flunet_all_growth_rates.csv")

figSX <- (p_gr_flunet_by_day + labs(tag="A") + coord_cartesian(ylim=c(-1,1)) + theme(plot.tag=element_text(face="bold"))) / 
  (p_gr_flunet_h3_by_day + labs(tag="B")+ coord_cartesian(ylim=c(-1,1)) +theme(plot.tag=element_text(face="bold"))) + plot_layout(guides = 'collect') & theme(legend.position='bottom')
ggsave("new_figures/figS5.png",figSX,width=8,height=8)
ggsave("new_figures/figS5.pdf",figSX,width=8,height=8)


#############################################
## FOURTH fit to overall flu cases
#############################################
#case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts <- read_csv("data/final/england_h3_cases_historic_manual.csv")
case_counts$date <- lubridate::dmy(case_counts$date)
## Start from the 2011/12 season to miss the H1N1 pandemic
case_counts <- case_counts %>% filter(date >= as.Date("2011-07-01"))
## Drop last row
case_counts <- case_counts[-nrow(case_counts),]
case_counts <- case_counts %>% mutate(total = `Influenza A not subtyped` + `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`) 
case_counts <- case_counts %>% mutate(Week=`Week`,Year = year(date))


case_counts1 <- case_counts
case_counts1$total <- round(case_counts1$total)
case_counts1$time <- as.numeric(as.factor(case_counts1$date))

mod <- construct_model(
  method = random_walk(),
  #method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    data=case_counts1,
    case_timeseries = "total",           # timeseries of case data
    time ="time"                       # date or time variable labels
  ),
  smoothing_params = smoothing_structure(   # independent smoothing structure 
    tau_mean = c(0),           # parameter - one for each pathogen
    tau_sd = c(1)
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  n_iter = 2000,
  n_warmup = 1000,
  n_chain = 3
)


gr <- growth_rate(fit)
p_cases <- plot(gr)
gr_case_dat <- p_cases$data

time_key_cases <- case_counts %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_case_dat <- gr_case_dat %>% left_join(time_key_cases)

## Convert dates to flu season
gr_case_dat <- gr_case_dat %>%
  dplyr::mutate(season = flu_season(date))

gr_case_dat <- gr_case_dat %>% group_by(season) %>% mutate(days_from_start = date - min(date))
gr_case_dat <- gr_case_dat %>% drop_na()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(y) %>% hist()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(days_from_start) %>% as.numeric() %>% hist()


p_gr_cases_all <- ggplot(data=gr_case_dat %>% filter(!(season %in% c("2019/20","2020/21","2021/22")))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1.1,1.1)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

## Align by day of year
gr_case_dat$day_of_year <- yday(gr_case_dat$date)
gr_case_dat <- gr_case_dat %>% group_by(season) %>%
  mutate(day_of_year = if_else(day_of_year < 183, day_of_year + 365, day_of_year)) %>%
  ungroup()

## Create labels for plot, label seasons by pre-pandemic, during pandemic (2020 and 2021), post pandemic (2022-2024) and current season (2025/26)
gr_case_dat <- gr_case_dat %>%
  mutate(plot_label = case_when(
    season %in% c("2012/13","2013/14","2014/15","2015/16","2016/17","2017/18","2018/19") ~ "Pre-pandemic",
    season %in% c("2019/20","2020/21","2021/22") ~ "During pandemic",
    season %in% c("2022/23","2023/24","2024/25") ~ "Post-pandemic",
    season == "2025/26" ~ "Current season",
    TRUE ~ "Pre-pandemic"
  ))


max_t <- max(gr_case_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
min_t <- min(gr_case_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
breaks <- seq(150,max_t,by=50)
labels <- if_else(breaks > 365, breaks - 365, breaks)

p_gr_case_by_day<- ggplot(data=gr_case_dat %>% filter(plot_label != "During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_vline(xintercept=365,linetype="dashed") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label),size=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  
  coord_cartesian(ylim=c(-1.1,1.1)) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")



## Align by peak data
## Filter gr_flunet_dat to be between September and May each season
gr_case_dat1 <- gr_case_dat %>%
  filter(month(date) >= 9 | month(date) <= 4 | (month(date) == 5 & day(date) <= 1)) %>%
  group_by(season) %>%
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup()

max_gr <- gr_case_dat1 %>% group_by(season) %>% filter(y == max(y))# %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)

max_gr %>% select(season, Year, Week, date, y, lb_95, ub_95) %>% arrange(season) %>% tail(5) %>% print()


gr_case_dat <- gr_case_dat %>% left_join(max_gr %>% select(season, day_of_year) %>% rename(peak_time = day_of_year))
gr_case_dat$day_shifted <- gr_case_dat$day_of_year - gr_case_dat$peak_time

## Get max growth rate and week of max growth rate
## Merge in raw data
gr_case_dat <-gr_case_dat %>% left_join(case_counts %>% mutate(gr = log((total+0.1)/lag(total+0.1,1)), gr_smooth = zoo::rollmean(gr, 4, fill=NA,align="right"))) 

peak_2025_26 <- gr_case_dat %>%
  filter(season == "2025/26", day_shifted >= -100, day_shifted <= 100) %>%
  summarise(max_y = max(y, na.rm = TRUE)) %>%
  pull(max_y)

label_dat <- gr_case_dat %>%
  filter(day_shifted >= -100, day_shifted <= 100, plot_label != "During pandemic", season != "2019/20") %>%
  group_by(season) %>%
  filter(y == max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(y >= peak_2025_26)  # only keep peaks higher than 2025/26
label_dat_all <- label_dat
p_gr_case_by_peak1 <- ggplot(data=gr_case_dat %>% filter(day_shifted <= 100, day_shifted >= -100,plot_label!="During pandemic")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label),size=0.75) +
  scale_y_continuous(limits=c(-1.1,1.1),breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day relative to peak") + ylab('Growth rate (per week)') +
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00"))+
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00"))+
  #scale_color_viridis_d("Time period") +
  #scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")
p_gr_case_by_peak <- p_gr_case_by_peak1 +
  geom_label_repel(
    data = label_dat,
    aes(
      x = day_shifted,
      y = y,
      label = season
    ),
    size = 3,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0,
    show.legend = FALSE,
    nudge_x=10,
    nudge_y=0.08,
    direction="both"
  )

p_gr_case_by_peak

p_gr_cases_allb <- p_gr_cases_all
p_gr_cases_all <- add_doubling_axis(p_gr_cases_all)
p_gr_case_all_by_dayb <- p_gr_case_by_day
p_gr_case_all_by_day <- add_doubling_axis(p_gr_case_by_day)
p_gr_case_all_by_peakb <- p_gr_case_by_peak
p_gr_case_all_by_peak1b <- p_gr_case_by_peak1
p_gr_case_all_by_peak <- add_doubling_axis(p_gr_case_by_peak)

ggsave("figures/growth_rates/ukhsa_all_flu_growth_rate.png",p_gr_cases_all,width=7,height=4)
ggsave("figures/growth_rates/ukhsa_all_flu_growth_rate_by_day.png",p_gr_case_all_by_day,width=7,height=4)
ggsave("figures/growth_rates/ukhsa_all_flu_growth_rate_by_peak.png",p_gr_case_all_by_peak,width=7,height=4)
write_csv(gr_case_dat,"results/ukhsa_all_flu_growth_rates.csv")

## Pull out peak growth rates
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y)

gr_case_dat1 <- gr_case_dat %>% group_by(season) %>% mutate(total_rel = total/sum(total))

gr_peaks_all <- gr_case_dat1 %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y) %>%
  select(y,lb_95,ub_95,Year,Week,date,season,plot_label,total_rel) %>%
  mutate(Estimate = paste0(signif(y,3), " (", signif(lb_95,3),"-",signif(ub_95,3),")")) %>% 
  select(season,Week, Estimate,plot_label,total_rel) %>%
  rename(Season=season,`Time period`=plot_label, `Week of peak growth`=Week,`Peak growth rate of all influenza cases`=Estimate,
         `Cumulative incidence (% of all influenza cases)`=total_rel)

#############################################
## FIFTH H3 flu cases from England case data
#############################################
#case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts <- read_csv("data/final/england_h3_cases_historic_manual.csv")
case_counts$date <- lubridate::dmy(case_counts$date)
## Start from the 2011/12 season to miss the H1N1 pandemic
case_counts <- case_counts %>% filter(date >= as.Date("2011-07-01"))

## Drop most recent week
case_counts <- case_counts[-nrow(case_counts),]

case_counts <- case_counts %>%
  mutate(total_tested = `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`,
         inferred_H3 = `Influenza A not subtyped` * `Influenza A H3N2`/total_tested,
         inferred_H3 = if_else(is.na(inferred_H3),0,inferred_H3),
         inferred_H3 = inferred_H3 + `Influenza A H3N2`)

case_counts <- case_counts %>% mutate(Week=`Week`,Year = year(date))
case_counts1 <- case_counts
case_counts1$inferred_H3 <- round(case_counts1$inferred_H3)
case_counts1$date <- as.numeric(as.factor(case_counts1$date))
mod <- construct_model(
  method = random_walk(),
  #method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    data=case_counts1,
    case_timeseries="inferred_H3",
    time="date"
    #case_timeseries = round((case_counts$inferred_H3)),           # timeseries of case data
    #time = as.numeric(as.factor(case_counts$date))                       # date or time variable labels

  ),
  smoothing_params = smoothing_structure(   # independent smoothing structure 
    tau_mean = c(0),           # parameter - one for each pathogen
    tau_sd = c(1)
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  n_iter = 2000,
  n_warmup = 1000,
  n_chain = 3
)

gr <- growth_rate(fit)
p_cases <- plot(gr)
gr_case_dat <- p_cases$data

time_key_cases <- case_counts %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_case_dat <- gr_case_dat %>% left_join(time_key_cases)

## Convert dates to flu season
gr_case_dat <- gr_case_dat %>%
  dplyr::mutate(season = flu_season(date))
gr_case_dat <- gr_case_dat %>% group_by(season) %>% mutate(days_from_start = date - min(date))
gr_case_dat <- gr_case_dat %>% drop_na()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(y) %>% hist()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(days_from_start) %>% as.numeric() %>% hist()


p_gr_cases <- ggplot(data=gr_case_dat %>% filter(!(season %in% c("2019/20","2020/21","2021/22")))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
  coord_cartesian(ylim=c(-1.1,1.1)) +
 # scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")
## Align by day of year
gr_case_dat$day_of_year <- yday(gr_case_dat$date)
gr_case_dat <- gr_case_dat %>% group_by(season) %>%
  mutate(day_of_year = if_else(day_of_year < 183, day_of_year + 365, day_of_year)) %>%
  ungroup()

## Create labels for plot, label seasons by pre-pandemic, during pandemic (2020 and 2021), post pandemic (2022-2024) and current season (2025/26)
gr_case_dat <- gr_case_dat %>%
  mutate(plot_label = case_when(
    season %in% c("2012/13","2013/14","2014/15","2015/16","2016/17","2017/18","2018/19") ~ "Pre-pandemic",
    season %in% c("2019/20","2020/21","2021/22") ~ "During pandemic",
    season %in% c("2022/23","2023/24","2024/25") ~ "Post-pandemic",
    season == "2025/26" ~ "Current season",
    TRUE ~ "Pre-pandemic"
  ))

max_t <- max(gr_case_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))
min_t <- min(gr_case_dat %>% select(season, day_of_year) %>% distinct() %>% pull(day_of_year))

breaks <- seq(150,max_t,by=50)
labels <- if_else(breaks > 365, breaks - 365, breaks)

p_gr_case_by_day<- ggplot(data=gr_case_dat %>% filter(plot_label != "During pandemic")) + 
  geom_vline(xintercept=365,linetype="dashed") +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label),size=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day of year (start of Epi week)") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  coord_cartesian(ylim=c(-1.1,1.1)) +
  
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")



## Align by peak data
## Filter gr_flunet_dat to be between September and May each season
gr_case_dat1 <- gr_case_dat %>%
  filter(month(date) >= 9 | month(date) <= 4 | (month(date) == 5 & day(date) <= 1)) %>%
  group_by(season) %>%
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup()

max_gr <- gr_case_dat1 %>% group_by(season) %>% filter(y == max(y))# %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)

max_gr %>% select(season, Year, Week, date, y, lb_95, ub_95) %>% arrange(season) %>% tail(5) %>% print()


gr_case_dat <- gr_case_dat %>% left_join(max_gr %>% select(season, day_of_year) %>% rename(peak_time = day_of_year))
gr_case_dat$day_shifted <- gr_case_dat$day_of_year - gr_case_dat$peak_time

## Get max growth rate and week of max growth rate
## Merge in raw data
gr_case_dat <-gr_case_dat %>% left_join(case_counts %>% mutate(gr = log((inferred_H3+0.1)/lag(inferred_H3+0.1,1)), gr_smooth = zoo::rollmean(gr, 4, fill=NA,align="right"))) 

## Add labels for the seasons with peak growth rates over 0.6
peak_2025_26 <- gr_case_dat %>%
  filter(season == "2025/26", day_shifted >= -100, day_shifted <= 100) %>%
  summarise(max_y = max(y, na.rm = TRUE)) %>%
  pull(max_y)
peak_2025_26 <- 0.6
label_dat <- gr_case_dat %>%
  filter(day_shifted >= -100, day_shifted <= 100, plot_label != "During pandemic", season != "2019/20") %>%
  group_by(season) %>%
  filter(y == max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(y >= peak_2025_26)  # only keep peaks higher than 2025/26


p_gr_case_by_peak1 <- ggplot(data=gr_case_dat %>% filter(day_shifted <= 100, day_shifted >= -100,plot_label!="During pandemic" & season != "2019/20")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  # geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label),size=0.75) +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1.1,1.1)) +
  
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Day relative to peak") + ylab('Growth rate (per week)') +
  #scale_color_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  #scale_fill_manual(values=c("Pre-pandemic"="blue","During pandemic"="black","Post-pandemic"="orange","Current season"="red"))+
  scale_color_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  scale_fill_manual("Time period",values=c("Pre-pandemic"="grey","Post-pandemic"="#0072B2","Current season"="#D55E00")) +
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_case_by_peak <- p_gr_case_by_peak1 +
  geom_label_repel(
    data = label_dat,
    aes(
      x = day_shifted,
      y = y,
      label = season
    ),
    size = 3,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0,
    show.legend = FALSE,
    nudge_x=10,
    nudge_y=0.08,
    direction="both"
  )
p_gr_cases_h3b <- p_gr_cases_all
p_gr_cases_h3 <- add_doubling_axis(p_gr_cases_all)
p_gr_case_h3_by_dayb <- p_gr_case_by_day
p_gr_case_h3_by_day <- add_doubling_axis(p_gr_case_by_day)
p_gr_case_h3_by_peakb <- p_gr_case_by_peak
p_gr_case_h3_by_peak <- add_doubling_axis(p_gr_case_by_peak)


gr_h3_dat1 <- gr_case_dat %>% group_by(season) %>% mutate(total_rel = inferred_H3/sum(inferred_H3))

gr_peaks_h3 <- gr_h3_dat1 %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y) %>%
  select(y,lb_95,ub_95,Year,Week,date,season,plot_label,total_rel) %>%
  mutate(Estimate = paste0(signif(y,3), " (", signif(lb_95,3),"-",signif(ub_95,3),")")) %>% 
  select(season,Week, Estimate,plot_label,total_rel) %>%
  rename(Season=season,`Time period`=plot_label, `Week of peak A/H3N2 growth` =Week, `Peak growth rate of A/H3N2 cases`=Estimate,`Cumulative incidence (% of total A/H3N2 cases)`=total_rel)

peak_grs_comb <- left_join(gr_peaks_h3,gr_peaks_all) %>% filter(`Time period` != "During pandemic") %>% arrange(`Peak growth rate of all influenza cases`)


dominant_dat <- gr_case_dat %>% group_by(season) %>% 
  summarize(total_h3 = sum(`Influenza A H3N2`),total_h1 = sum(`Influenza A H1N1pdm09`),total_b = sum(`Influenza B`)) %>% pivot_longer(-season) %>% group_by(season) %>% filter(value == max(value))

dominant_key <- c("total_h3"="A/H3N2","total_h1"="A/H1N1pdm09","total_b"="Influenza B")
dominant_dat$Dominant <- dominant_key[dominant_dat$name]

peak_grs_comb <- peak_grs_comb %>% left_join(dominant_dat %>% select(season,Dominant) %>% rename(Season=season))

ggsave("figures/growth_rates/ukhsa_h3_flu_growth_rate.png",p_gr_cases_h3,width=7,height=4)
ggsave("figures/growth_rates/ukhsa_h3_flu_growth_rate_by_day.png",p_gr_case_h3_by_day,width=7,height=4)
ggsave("figures/growth_rates/ukhsa_h3_flu_growth_rate_by_peak.png",p_gr_case_h3_by_peak,width=7,height=4)

write_csv(gr_case_dat,"results/ukhsa_h3_flu_growth_rates.csv")
write_csv(peak_grs_comb,"results/ukhsa_peak_growth_rates.csv")
p1 <- p_gr_case_h3_by_peak + coord_cartesian(ylim=c(-1,1)) + theme(axis.title.x = element_blank(), axis.title.y=element_blank()) + ggtitle("A/H3N2 cases aligned by peak growth") + labs(tag="B")
p2 <- p_gr_case_all_by_peak + coord_cartesian(ylim=c(-1,1)) + theme(axis.title.y=element_blank(),axis.title.x=element_text(size=10))+ ggtitle("All influenza cases aligned by peak growth")+ labs(tag="D")
p3 <- p_gr_case_h3_by_day + coord_cartesian(ylim=c(-1,1))+ theme(axis.title.x = element_blank())+ theme(axis.title.y=element_blank())+ ggtitle("A/H3N2 cases aligned by day of year")+ labs(tag="A")
p4 <- p_gr_case_all_by_day + coord_cartesian(ylim=c(-1,1)) + theme(axis.title.y=element_blank(),axis.title.x=element_text(size=10))+ ggtitle("All influenza cases aligned by day of year")+ labs(tag="C")

comb_plot <- ((p3 / p4) | (p1 / p2))
comb_plot <- comb_plot+plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",legend.text = element_text(size=10),legend.title=element_text(size=10),plot.title=element_text(size=10),plot.tag=element_text(face="bold"))

comb_plot
ggsave("new_figures/main_figure.png",comb_plot,width=10,height=7)
ggsave("new_figures/main_figure.pdf",comb_plot,width=10,height=7)

figS4 <- (p_gr_case_all_by_peak + coord_cartesian(ylim=c(-1,1)) + labs(tag="A") + theme(plot.tag=element_text(face="bold"))) /
(p_gr_case_h3_by_peak + coord_cartesian(ylim=c(-1,1)) + labs(tag="B") + theme(plot.tag=element_text(face="bold"))) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",legend.text = element_text(size=10),legend.title=element_text(size=10),plot.title=element_text(size=10))

ggsave("new_figures/figS4.pdf",figS4,width=8,height=8)
ggsave("new_figures/figS4.png",figS4,width=8,height=8)

p_main_all <- p_gr_case_all_by_dayb + labs(tag="") + theme(axis.title=element_text(size=16),axis.text=element_text(size=12),legend.text=element_text(size=12),
                            legend.title=element_text(size=12),plot.title=element_text(size=16)) +
  ggtitle("Weekly growth rate of all influenza cases") +
  xlab("Day of year") +
  coord_cartesian(ylim=c(-1,1))+
  geom_label_repel(
    data = label_dat_all,
    aes(
      x = day_of_year,
      y = y,
      label = season
    ),
    size = 4,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0,
    show.legend = FALSE,
    nudge_x=15,
    nudge_y=0.08,
    direction="both"
  ) +
  geom_label(data=data.frame(x=366,y=-1.04,label="1st January"),aes(x=x,y=y,label=label))
ggsave("new_figures/Fig1.pdf",p_main_all,width=10,height=6)
ggsave("new_figures/Fig1.png",p_main_all,width=10,height=6)


p_main_all_shifted <- p_gr_case_all_by_peak1b + labs(tag="") + theme(axis.title=element_text(size=16),axis.text=element_text(size=12),legend.text=element_text(size=12),
                                                         legend.title=element_text(size=12),plot.title=element_text(size=16)) +
  ggtitle("Weekly growth rate of all influenza cases aligned to day of peak growth rate") +
  xlab("Days relative to peak growth rate") +
  coord_cartesian(ylim=c(-1,1))+
  geom_label_repel(
    data = label_dat_all,
    aes(
      x = day_shifted,
      y = y,
      label = season
    ),
    size = 4,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0,
    show.legend = FALSE,
    nudge_x=15,
    nudge_y=0.08,
    direction="both"
  ) 

ggsave("figures/growth_rates/all_cases_main_aligned.pdf",p_main_all_shifted,width=10,height=6)
ggsave("figures/growth_rates/all_cases_main_aligned.png",p_main_all_shifted,width=10,height=6)

gr_plot_data_comb <- bind_rows(p_gr_flunet_all$data %>% mutate(Scenario = "A/H3N2, WHO FluNet"),
p_gr_flunet_all$data %>% mutate(Scenario = "All influenza, WHO FluNet"),
p_gr_cases$data %>% mutate(Scenario = "A/H3N2, UKHSA"),
p_gr_cases_all$data %>% mutate(Scenario = "All influenza, UKHSA"))

gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "A/H3N2, WHO FluNet") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y) %>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "A/H3N2, WHO FluNet") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(y) %>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "All influenza, WHO FluNet") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y) %>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "All influenza, WHO FluNet") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(y) %>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "A/H3N2, UKHSA") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y)%>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "A/H3N2, UKHSA") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(y)%>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "All influenza, UKHSA") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(-y) %>% head()
gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2019/20","2021/22"))) %>% filter(Scenario == "All influenza, UKHSA") %>% group_by(season) %>% filter(y == max(y)) %>% arrange(y)%>% head()


p_gr_all_comb <- ggplot(data=gr_plot_data_comb %>% filter(!(season %in% c("2020/21","2021/22")))) + 
  ## Add rectange over pandemic
  geom_rect(xmin=as.Date("2020-01-01"),xmax=as.Date("2022-06-01"),ymin=-Inf,ymax=Inf,fill="grey",alpha=0.2) +
  ## Add label for pandemic
  geom_text(data=data.frame(x=as.Date("2021-03-01"),y=0.8,label="COVID-19\npandemic"),aes(x=x,y=y,label=label),size=3) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
  # scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_date(date_breaks = "5 years")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  facet_wrap(~Scenario)
ggsave("figures/growth_rates/all_growth_rates_combined.png",p_gr_all_comb,width=8,height=6)

