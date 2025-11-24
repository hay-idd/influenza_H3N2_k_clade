library(EpiStrainDynamics)
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)
library(data.table)
library(mgcv)

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
overall_ILIplus <- read_csv("data/final/ili_plus_datasets_by_age.csv")
overall_ILIplus <- overall_ILIplus %>% select(Year,Week,date,group,ILIplus,Influenza_cases,Season)
overall_ILIplus <- overall_ILIplus %>% bind_rows(overall_ILIplus %>% group_by(Year,Week,date,Season) %>% summarize(ILIplus = sum(ILIplus),Influenza_cases=sum(Influenza_cases)) %>% mutate(group="All"))
overall_ILIplus$group <- factor(overall_ILIplus$group, levels=c("All","1-4","5-14","15-44","45-64","65+"))
overall_ILIplus <- overall_ILIplus %>% 
  select(Season,Year,Week,date,group,ILIplus,Influenza_cases) %>% 
  mutate("Source"="RCGP") %>%
  rename(cases_H3=ILIplus,
         cases_total=Influenza_cases)
overall_ILIplus$Season <- season_from_date(overall_ILIplus$date)

## England case counts
case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts$date <- lubridate::dmy(case_counts$Date)
case_counts <- case_counts %>% mutate(total = `Influenza A not subtyped` + `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`) 
case_counts <- case_counts %>%
  mutate(total_tested = `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`,
         inferred_H3 = `Influenza A not subtyped` * `Influenza A H3N2`/total_tested,
         inferred_H3 = if_else(is.na(inferred_H3),0,inferred_H3),
         inferred_H3 = inferred_H3 + `Influenza A H3N2`) %>%
  rename(Week = `Week number`) %>%
  mutate(Year = lubridate::year(date)) %>%
  select(Week,Year,date,inferred_H3, total) %>%
  rename(cases_H3=inferred_H3,cases_total=total) %>% 
  mutate("Source"="UKHSA")
case_counts$Season <- season_from_date(case_counts$date)


## WHO FluNet
flunet <- read_csv("data/final/flunet_h3_cases_historic.csv")
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
  xlab("Date (end of Epi week)") +
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
  xlab("Date (end of Epi week)") +
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
  xlab("Date (end of Epi week)") +
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
  xlab("Date (end of Epi week)") +
  theme_use

## Plot by age group
rcgp_h3_gr <-  calc_emp_growth(overall_ILIplus,"date","cases_H3",group1="group")
rcgp_h3_gr <- fit_gam_grs(rcgp_h3_gr %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_H3",time="date_int",group="group")[[1]]

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
  xlab("Date (end of Epi week)") +
  theme_use

rcgp_all_gr <-  calc_emp_growth(overall_ILIplus,"date","cases_total",group1="group")
rcgp_all_gr <- fit_gam_grs(rcgp_all_gr %>% mutate(date_int = as.numeric(as.factor(date))),y="cases_total",time="date_int",group="group")[[1]]


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
  xlab("Date (end of Epi week)") +
  theme_use

## Look at the ratio of growth rates by age to overall
rcgp_all_gr %>% filter(group != "All") %>%
  left_join(rcgp_all_gr %>% filter(group == "All") %>% select(date,gr) %>% rename(gr_all=gr)) %>%
  mutate(gr_ratio = gr/gr_all) %>%
  ggplot() + 
  geom_line(aes(x=date,y=gr_ratio,col=group,group=interaction(group,Season))) +
  geom_hline(yintercept=1,linetype="dashed") +
  ylab("Ratio of age-specific growth rate to overall growth rate") +
  xlab("Date (end of Epi week)") +
  scale_color_brewer("Age group",palette = "Set2") +
  theme_use +
  facet_wrap(~group)

rcgp_h3_gr %>% filter(group != "15-44") %>%
  left_join(rcgp_h3_gr %>% filter(group == "15-44") %>% select(date,gr) %>% rename(gr_all=gr)) %>%
  mutate(gr_ratio = log(gr/gr_all)) %>%
  ggplot() + 
  geom_line(aes(x=date,y=gr_ratio,col=group,group=interaction(group,Season))) +
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Ratio of age-specific growth rate to overall growth rate") +
  xlab("Date (end of Epi week)") +
  scale_color_brewer("Age group",palette = "Set2") +
  theme_use +
  facet_wrap(~group)


## Incidence plots
p_inc_h3 <- plot_incidence(all_grs_H3,x="date",y="cases_H3",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (end of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use

p_inc_h3_recent <- plot_incidence(all_grs_H3 %>% filter(Year >= 2021),x="date",y="cases_H3",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (end of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use

p_inc_all <- plot_incidence(all_grs_cases,x="date",y="cases_total",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (end of Epi week)") +
  ylab("Weekly cases of confirmed or suspected influenza") + theme_use

p_inc_all_recent <- plot_incidence(all_grs_cases%>% filter(Year >= 2021),x="date",y="cases_total",color="Source",point=FALSE) + 
  scale_color_brewer("Source",palette = "Set1") + xlab("Date (end of Epi week)") +
  ylab("Weekly cases of confirmed or suspected influenza") + theme_use

p_inc_h3_age <- plot_incidence(overall_ILIplus,x="date",y="cases_H3",color="group",point=FALSE)+ 
  scale_color_brewer("Age group",palette = "Set2") + xlab("Date (end of Epi week)") +
  ylab("Weekly cases of confirmed or suspected A/H3N2") + theme_use
p_inc_all_age <- plot_incidence(overall_ILIplus,x="date",y="cases_total",color="group",point=FALSE)+ 
  scale_color_brewer("Age group",palette = "Set2") + xlab("Date (end of Epi week)") +
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

write.csv(fit_gams_H3,file="results/fit_gams_H3.csv")
write.csv(fit_gams_all,"results/fit_gams_all.csv")
write.csv(rcgp_h3_gr,"results/rcgp_h3_gr.csv")
write.csv(fit_gams_H3,"results/rcgp_all_gr.csv")


ggsave("figures/raw_data/p_h3_gr_all.png",p_h3_gr_all,width=7,height=6)
ggsave("figures/raw_data/p_h3_gr_recent.png",p_h3_gr_recent_compare,width=7,height=4)
ggsave("figures/raw_data/p_gr_all.png",p_gr_all,width=7,height=6)
ggsave("figures/raw_data/p_gr_recent_compare.png",p_gr_recent_compare,width=7,height=4)
ggsave("figures/raw_data/p_gr_h3_age.png",p_gr_h3_age,width=8,height=7)
ggsave("figures/raw_data/p_gr_all_age.png",p_gr_all_age,width=8,height=7)

ggsave("figures/raw_data/p_inc_h3.png",p_inc_h3,width=7,height=4)
ggsave("figures/raw_data/p_inc_h3_recent.png",p_inc_h3_recent,width=7,height=4)
ggsave("figures/raw_data/p_inc_all.png",p_inc_all_recent,width=7,height=4)
ggsave("figures/raw_data/p_inc_all_recent.png",p_inc_all,width=7,height=4)
ggsave("figures/raw_data/p_inc_h3_age.png",p_inc_h3_age,width=7,height=4)
ggsave("figures/raw_data/p_inc_all_age.png",p_inc_all_age,width=7,height=4)


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
  scale <- 1
  ## Upscale counts to get convergence
  if(agegroup != "All") scale <- 10
  mod <- construct_model(
    #method = random_walk(), 
    method=p_spline(spline_degree = 3, days_per_knot = 3),
    pathogen_structure =single(
      case_timeseries = round((all_dat_tmp$ILIplus*scale)),           # timeseries of case data
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


rects <- school_periods_oxford %>%
  mutate(xmin = start - 0, xmax = end + 1)

time_key <- overall_ILIplus %>% select(Year, Week, date,Season) %>%
  mutate(time = (Year-2023)*52 + Week) %>%
  mutate(time = as.numeric(as.factor(time))) %>%
  group_by(Season) %>%
  mutate(time_rel = time - min(time)) %>%
  ungroup()

tmp_dat_all$agegroup <- factor(tmp_dat_all$agegroup, levels=c("All","1-4","5-14","15-44","45-64","65+"))


raw_data <- overall_ILIplus %>% group_by(group) %>% 
  mutate(gr = log((ILIplus+0.1)/lag(ILIplus+0.1,1))) %>%
  #mutate(gr = zoo::rollmean(gr,3,fill=NA,align='right')) %>%
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
  geom_line(aes(x=date,y=y,col=agegroup,alpha="P-spline"),linewidth=0.75) +
  geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup,alpha="Raw data")) +
  scale_fill_brewer("Age group", palette = "Set2") +
  scale_color_brewer("Age group", palette = "Set2") +
  xlab("Week") + ylab('Growth rate (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1.2,1.2)) +
  theme(legend.position="bottom",legend.direction="horizontal")
#p_gr_by_age <- add_doubling_axis(p_gr_by_age)


## Plot overall GR and by age, aligned by season

## Plot overall GR and by age
p_gr_by_age_shifted <- ggplot(tmp_dat_all%>% filter(Season != "2022 to 2023") %>% group_by(Season) %>%
                                mutate(time = time - min(time))) + 
  #geom_rect(data = rects, inherit.aes = FALSE,
    #        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill="grey70",
   #         alpha = 0.4) +  
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=time_rel,ymin=lb_95,ymax=ub_95,fill=Season),alpha=0.5) +
  geom_line(aes(x=time_rel,y=y,col=Season),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set1") +
  scale_color_brewer("Age group", palette = "Set1") +
  xlab("Week (relative to start of season)") + ylab('Growth rate (per week)') +
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1.2,1.2)) +
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
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme(legend.position="bottom",legend.direction="horizontal")
#p_gr_by_age_shifted_rel <- add_doubling_axis(p_gr_by_age_shifted_rel)

## Compare growth rates by age to the 15-44 group
p_gr_by_age_diff <- ggplot(tmp_dat_all %>% filter(!(agegroup %in% c("15-44","All"))) %>%
  left_join(tmp_dat_all %>% filter(agegroup == "15-44") %>% select(time,Season,y) %>% rename(y_45_64 = y))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  #geom_ribbon(aes(x=date,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=date,y=y - y_45_64,col=agegroup),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set2") +
  scale_color_brewer("Age group", palette = "Set2") +
  xlab("Week") + ylab('Growth rate difference to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme(legend.position="bottom",legend.direction="horizontal")

## Compare growth rates by age to the 15-44 group, by date of season
p_gr_by_age_diff_season <- ggplot(tmp_dat_all %>% filter(!(agegroup %in% c("15-44","All")),Season != "2022 to 2023") %>%
                             left_join(tmp_dat_all %>% filter(agegroup == "15-44") %>% select(time_rel,Season,y) %>% rename(y_45_64 = y))%>% group_by(Season) %>%mutate(time = time - min(time))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  #geom_ribbon(aes(x=date,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=time_rel,y=y - y_45_64,col=Season,group=interaction(agegroup,Season)),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Season", palette = "Set1") +
  scale_color_brewer("Season", palette = "Set1") +
  xlab("Week") + ylab('Growth rate difference to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme(legend.position="bottom",legend.direction="horizontal")

## Align by peak growth rate
p_gr_by_age_diff_peak <- ggplot(tmp_dat_shifted %>% filter(!(agegroup %in% c("15-44","All")),Season != "2022 to 2023") %>%
                             left_join(tmp_dat_shifted %>% filter(agegroup == "15-44") %>% select(time,Season,y) %>% 
                                         rename(y_45_64 = y))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  #geom_ribbon(aes(x=time_plot,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=Season),alpha=0.5) +
  geom_line(aes(x=time_plot,y=y - y_45_64,col=Season),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Season", palette = "Set1") +
  scale_color_brewer("Season", palette = "Set1") +
  xlab("Week") + ylab('Growth rate difference to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-1,1)) +
  theme(legend.position="bottom",legend.direction="horizontal")

## Ratios
## Compare growth rates by age to the 15-44 group, ratios
p_gr_by_age_ratio <- ggplot(tmp_dat_all %>% filter(!(agegroup %in% c("15-44","All")),Season != "2022 to 2023") %>%
                             left_join(tmp_dat_all %>% filter(agegroup == "15-44") %>% select(time,Season,y) %>% rename(y_45_64 = y))) + 
  geom_hline(yintercept=1,linetype="dashed") +
  #geom_ribbon(aes(x=date,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=date,y=(y/y_45_64),col=agegroup),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Age group", palette = "Set2") +
  scale_color_brewer("Age group", palette = "Set2") +
  xlab("Week") + ylab('Growth rate ratio to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  #scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-5,5)) +
  theme(legend.position="bottom",legend.direction="horizontal")

## Compare growth rates by age to the 15-44 group, by date of season
p_gr_by_age_ratio_season <- ggplot(tmp_dat_all %>% filter(!(agegroup %in% c("15-44","All")),Season != "2022 to 2023") %>%
                                    left_join(tmp_dat_all %>% filter(agegroup == "15-44") %>% select(time_rel,Season,y) %>% rename(y_45_64 = y))) + 
  geom_hline(yintercept=1,linetype="dashed") +
  #geom_ribbon(aes(x=date,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=agegroup),alpha=0.5) +
  geom_line(aes(x=time_rel,y=y/y_45_64,col=Season,group=interaction(agegroup,Season)),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Season", palette = "Set1") +
  scale_color_brewer("Season", palette = "Set1") +
  xlab("Week") + ylab('Growth rate ratio to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  #scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-5,5)) +
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_by_age_ratio_peak <- ggplot(tmp_dat_shifted %>% filter(!(agegroup %in% c("15-44","All")),Season != "2022 to 2023") %>%
                                  left_join(tmp_dat_shifted %>% filter(agegroup == "15-44") %>% select(time,Season,y) %>% 
                                              rename(y_45_64 = y))) + 
  geom_hline(yintercept=1,linetype="dashed") +
  #geom_ribbon(aes(x=time_plot,ymin=lb_95 - y_45_64,ymax=ub_95 - y_45_64,fill=Season),alpha=0.5) +
  geom_line(aes(x=time_plot,y=y /y_45_64,col=Season),linewidth=0.75) +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  scale_fill_brewer("Season", palette = "Set1") +
  scale_color_brewer("Season", palette = "Set1") +
  xlab("Week") + ylab('Growth rate ratio to 15-44 (per week)') +
  scale_alpha_manual(values=c("Raw data"=0.25,"P-spline"=1),name="Growth rate type")+
  theme_use + facet_wrap(~agegroup,ncol=2) +
  #scale_y_continuous(breaks=seq(-1,1,by=0.2)) +
  coord_cartesian(ylim=c(-5,5)) +
  theme(legend.position="bottom",legend.direction="horizontal")



ggsave("figures/growth_rates/RCGP_p_gr_by_age_diff.png",p_gr_by_age_diff,width=8,height=6)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_diff_season.png",p_gr_by_age_diff_season,width=8,height=6)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_diff_peak.png",p_gr_by_age_diff_peak,width=8,height=6)


ggsave("figures/growth_rates/RCGP_p_gr_by_age_ratio.png",p_gr_by_age_ratio,width=8,height=6)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_ratio_season.png",p_gr_by_age_ratio_season,width=8,height=6)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_ratio_peak.png",p_gr_by_age_ratio_peak,width=8,height=6)


ggsave("figures/growth_rates/RCGP_p_gr_by_age.png",p_gr_by_age,width=8,height=7)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_shifted.png",p_gr_by_age_shifted,width=8,height=7)
ggsave("figures/growth_rates/RCGP_p_gr_by_age_shifted_rel.png",p_gr_by_age_shifted_rel,width=8,height=7)

write_csv(tmp_dat_all,"results/RCGP_growth_rates_by_age.csv")
write_csv(tmp_dat_shifted,"results/RCGP_growth_rates_by_age_shifted.csv")

#############################################
## SECONDARY FluNet dataset fit to H3
#############################################
flunet <- read_csv("data/final/flunet_h3_cases_historic.csv")
flunet <- flunet %>% mutate(H3_sum = if_else(is.na(H3_sum),0,H3_sum))

mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3), ## I report 3 in the paper, but this was originally 4 if you see things change a lot
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


ggsave("figures/growth_rates/flunet_h3_growth_rate.png",p_gr_flunet,width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_day.png",p_gr_flunet_by_day,width=7,height=4)
ggsave("figures/growth_rates/flunet_h3_growth_rate_by_peak.png",p_gr_flunet_by_peak,width=7,height=4)
write_csv(gr_flunet_dat,"results/flunet_h3_growth_rates.csv")


#############################################
## THIRD FluNet dataset fit to all cases
#############################################
flunet <- read_csv("data/final/flunet_h3_cases_historic.csv")
flunet <- flunet %>% mutate(flu_pos_sum = if_else(is.na(flu_pos_sum),0,flu_pos_sum))

mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    case_timeseries = round((flunet$flu_pos_sum)),           # timeseries of case data
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

p_gr_flunet_all <- add_doubling_axis(p_gr_flunet_all)
p_gr_flunet_by_day <- add_doubling_axis(p_gr_flunet_by_day)
p_gr_flunet_by_peak <- add_doubling_axis(p_gr_flunet_by_peak)

ggsave("figures/growth_rates/flunet_all_growth_rate.png",p_gr_flunet,width=7,height=4)
ggsave("figures/growth_rates/flunet_all_growth_rate_by_day.png",p_gr_flunet_by_day,width=7,height=4)
ggsave("figures/growth_rates/flunet_all_growth_rate_by_peak.png",p_gr_flunet_by_peak,width=7,height=4)
write_csv(gr_flunet_dat,"results/flunet_all_growth_rates.csv")


#############################################
## FOURTH fit to overall flu cases
#############################################
case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts$date <- lubridate::dmy(case_counts$Date)
case_counts <- case_counts %>% mutate(total = `Influenza A not subtyped` + `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`) 
case_counts <- case_counts %>% mutate(Week=`Week number`,Year = year(date))
mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    case_timeseries = round((case_counts$total)),           # timeseries of case data
    time = as.numeric(as.factor(case_counts$date)),                       # date or time variable labels
    pathogen_name = 'Influenza'                # optional name of pathogen
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
p_cases <- plot(gr)
gr_case_dat <- p_cases$data

time_key_cases <- case_counts %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_case_dat <- gr_case_dat %>% left_join(time_key_cases)

## For later years, just assign season manually
gr_case_dat <- gr_case_dat %>%
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
gr_case_dat <- gr_case_dat %>% group_by(season) %>% mutate(days_from_start = date - min(date))
gr_case_dat <- gr_case_dat %>% drop_na()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(y) %>% hist()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(days_from_start) %>% as.numeric() %>% hist()


p_gr_cases_all <- ggplot(data=gr_case_dat %>% filter(!(season %in% c("2020/21","2021/22")))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

p_gr_cases_all <- add_doubling_axis(p_gr_cases_all)
ggsave("figures/growth_rates/ukhsa_all_flu_growth_rate.png",p_gr_cases_all,width=7,height=4)
write_csv(gr_case_dat,"results/ukhsa_all_flu_growth_rates.csv")


#############################################
## FIFTH H3 flu cases from England case data
#############################################
case_counts <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
case_counts$date <- lubridate::dmy(case_counts$Date)
case_counts <- case_counts %>%
  mutate(total_tested = `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza B`,
         inferred_H3 = `Influenza A not subtyped` * `Influenza A H3N2`/total_tested,
         inferred_H3 = if_else(is.na(inferred_H3),0,inferred_H3),
         inferred_H3 = inferred_H3 + `Influenza A H3N2`)

case_counts <- case_counts %>% mutate(Week=`Week number`,Year = year(date))
mod <- construct_model(
  #method = random_walk(),
  method=p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    case_timeseries = round((case_counts$inferred_H3)),           # timeseries of case data
    time = as.numeric(as.factor(case_counts$date)),                       # date or time variable labels
    pathogen_name = 'Influenza H3'                # optional name of pathogen
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
p_cases <- plot(gr)
gr_case_dat <- p_cases$data

time_key_cases <- case_counts %>% select(Year, Week, date) %>%
  mutate(time = as.numeric(as.factor(date)))
gr_case_dat <- gr_case_dat %>% left_join(time_key_cases)

## For later years, just assign season manually
gr_case_dat <- gr_case_dat %>%
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
gr_case_dat <- gr_case_dat %>% group_by(season) %>% mutate(days_from_start = date - min(date))
gr_case_dat <- gr_case_dat %>% drop_na()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(y) %>% hist()
gr_case_dat %>% group_by(season) %>% filter(y == max(y)) %>% pull(days_from_start) %>% as.numeric() %>% hist()


p_gr_cases <- ggplot(data=gr_case_dat %>% filter(!(season %in% c("2020/21","2021/22")))) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_ribbon(aes(x=date,ymin=lb_50,ymax=ub_50,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=date,ymin=lb_95,ymax=ub_95,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=date,y=y,group=season),col="black") +
  #geom_line(data=raw_data,aes(x=date,y=gr,col=agegroup),alpha=0.4) +
  xlab("Date (by week)") + ylab('Growth rate (per week)') +
 # scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,by=0.2)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year")+
  theme_use + 
  theme(legend.position="bottom",legend.direction="horizontal")

#p_gr_cases <- add_doubling_axis(p_gr_cases)
ggsave("figures/growth_rates/ukhsa_h3_flu_growth_rate.png",p_gr_cases,width=7,height=4)
write_csv(gr_case_dat,"results/ukhsa_h3_flu_growth_rates.csv")

gr_plot_data_comb <- bind_rows(p_gr_flunet$data %>% mutate(Scenario = "A/H3N2, WHO FluNet"),
p_gr_flunet_all$data %>% mutate(Scenario = "All influenza, WHO FluNet"),
p_gr_cases$data %>% mutate(Scenario = "A/H3N2, UKHSA"),
p_gr_cases_all$data %>% mutate(Scenario = "All influenza, UKHSA"))

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
