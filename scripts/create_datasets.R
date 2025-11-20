## Combines the various data streams that Punya and Alex collated into ILI+ indicators

library(tidyverse)
library(purrr)
library(patchwork)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
source("R/funcs.R")

save_wd <- "figures/raw_data/"

theme_use <- theme_bw()

desired_age_groups <- c("1-4","5-14","15-44","45-64","65+")

## Demography data
age_cap <- 90L
age_groups <- read_csv("data/ons_ages_2024.csv")
age_groups <- age_groups %>% select(Age_group, `ENGLAND`)
age_groups <- age_groups[2:nrow(age_groups),]
colnames(age_groups) <- c("age_group","N")
age_groups <- age_groups %>% mutate(prop = N/sum(N))
N_england <- sum(age_groups$N)
age_groups <- age_groups %>% 
  mutate(age_group = if_else(age_group=="90+", "90", age_group)) %>% 
  mutate(age_group=as.numeric(age_group))
## Test function
get_age_group_n(age_groups, 0, 4)

######################################################
## First dataset -- ILI overall and by age
######################################################
## Read in ILI per 100,000 by age group
ili_cases_recent <- read_csv("data/rcgp_by_age/data_25_24_23.csv") %>% fill(Year,.direction="down") %>% arrange(Year,Week) %>% drop_na()
ili_cases_older <- read_csv("data/rcgp_by_age/data_23_22_21.csv")%>% 
  fill(Year,.direction="down")%>% 
  arrange(Year,Week)%>% drop_na() %>%
  ## Remove last week as it's also in recent data
  filter(!(Year==2023 & Week==26))

## Get overall ILI cases
ili_cases_all <- bind_rows(ili_cases_recent %>% select(Year, Week, `All ages`), ili_cases_older %>% select(Year, Week, `All ages`)) %>% arrange(Year, Week)
ili_cases_all <- ili_cases_all %>% mutate(
  date = epiweek_to_date_formula(Year,Week,FALSE)
)
## Get actual start and end dates of the Epi weeks
p_ili_all <- ggplot(ili_cases_all, aes(x=date, y=`All ages`*N_england/100000)) + geom_line() + theme_use + xlab("Date (start of Epi week)") + ylab("ILI cases (RCGP, England, all ages)") +
  scale_x_date(date_labels="%b %Y",date_breaks="1 year")
ggsave(filename=paste0(save_wd,"ili_all_rcgp.png"),p_ili_all,width=8,height=4,units="in",dpi=300)

## Create age-structured data. Need to combine different age groups and re-weight
## Recent data
ili_cases_recent <- ili_cases_recent %>% mutate(
  date = epiweek_to_date_formula(Year,Week,FALSE)
)
ili_cases_recent <- ili_cases_recent %>% select(-`All ages`) %>% pivot_longer(-c(Year,Week,date))
ili_cases_recent <- extract_age_range(ili_cases_recent)
ili_cases_recent <- ili_cases_recent %>% rename(age_group=name,ILI_per_100k=value)
ili_cases_recent_expanded <- expand_age_group_ili(ili_cases_recent)
ili_cases_recent_expanded <- ili_cases_recent_expanded %>% left_join(age_groups %>% rename(age=age_group)) %>% select(-c(age_start,age_end,age_group))

## Older data
ili_cases_older <- ili_cases_older %>% mutate(
  date = epiweek_to_date_formula(Year,Week,FALSE)
)
ili_cases_older <- ili_cases_older %>% select(-`All ages`) %>% pivot_longer(-c(Year,Week,date))
ili_cases_older <- extract_age_range(ili_cases_older)
ili_cases_older <- ili_cases_older %>% rename(age_group=name,ILI_per_100k=value)
ili_cases_older_expanded <- expand_age_group_ili(ili_cases_older)
ili_cases_older_expanded <- ili_cases_older_expanded %>% left_join(age_groups %>% rename(age=age_group)) %>% select(-c(age_start,age_end,age_group))

ili_cases_comb_expanded <- bind_rows(ili_cases_older_expanded, ili_cases_recent_expanded) %>% mutate(ILI = ILI_per_100k*N/100000) %>% select(-prop)

## Be very careful, as older data does not have any data on 0 year olds, so need to request 1-14 only
ili_cases_comb_expanded_grouped <- combine_age_groups_ILI(ili_cases_comb_expanded,desired_age_groups)

p_ili_by_age <- ggplot(ili_cases_comb_expanded_grouped) + 
  geom_line(aes(x=date,y=ILI*100000/N,colour=group)) + 
  theme_use + 
  xlab("Date (start of Epi week)") + 
  ylab("ILI cases (RCGP, England, by age group)") +
  scale_x_date(date_labels="%b %Y",date_breaks="1 year") +
  scale_colour_brewer(palette="Set1",name="Age group")

######################################################
## Second dataset -- % of tests which are positive for influenza
######################################################
flu_pos_recent <- read_csv("data/ukhsa/datamart_weekly_positivity.csv")
flu_pos_recent <- flu_pos_recent %>% pivot_longer(-c(Date,`Week number`)) %>% rename(age_group=name,positivity=value) %>%
  mutate(positivity=positivity/100)
flu_pos_recent <- parse_age_groups_pos(flu_pos_recent)
flu_pos_recent$Date <- lubridate::dmy(flu_pos_recent$Date)
flu_pos_recent <- flu_pos_recent %>% rename(date=Date,Week=`Week number`) %>% select(-age_group) %>% mutate(Year=lubridate::year(date))
flu_pos_recent <- expand_age_group_pos(flu_pos_recent) %>% left_join(age_groups %>% rename(age=age_group))
flu_pos_recent <- combine_age_groups_pos(flu_pos_recent,desired_age_groups)

## Need to shift days +6 as this is start of reporting period
flu_pos_recent$date <- flu_pos_recent$date + 6

flu_pos_old <- read_csv("data/ukhsa/datamart_weekly_positivity_2.csv") %>% 
  rename(Year=year,Week=week) %>%
  fill(Year,.direction="down") %>% 
  arrange(Year,Week) %>% drop_na()
flu_pos_old <- flu_pos_old %>% mutate(
  date = epiweek_to_date_formula(Year,Week,FALSE)
) 

flu_pos_old <- flu_pos_old %>% pivot_longer(-c(date,Week,Year)) %>% rename(age_group=name,positivity=value) %>%
  mutate(positivity=positivity/100)
flu_pos_old <- parse_age_groups_pos(flu_pos_old) %>% select(-age_group)
flu_pos_old <- expand_age_group_pos(flu_pos_old) %>% left_join(age_groups %>% rename(age=age_group))
flu_pos_old <- combine_age_groups_pos(flu_pos_old,desired_age_groups)

flu_pos_comb <- bind_rows(flu_pos_old,flu_pos_recent)

p_influenza_cases_age <- left_join(ili_cases_comb_expanded_grouped,flu_pos_comb) %>% drop_na() %>% 
  mutate(ILIplus = ILI*positivity) %>%
  ggplot() + geom_line(aes(x=date,y=ILIplus,col=group))+ 
  theme_use + 
  xlab("Date (start of Epi week)") + 
  ylab("Estimated influenza cases (RCGP and Resp DataMart,\n England, by age group)") +
  scale_x_date(date_labels="%b %Y",date_breaks="1 year") +
  scale_colour_brewer(palette="Set1",name="Age group")


######################################################
## Third dataset -- % influenza tests which are H3
######################################################
flu_H3 <- read_csv("data/ukhsa/overall_subtype_percentages.csv")
flu_H3$Date <- lubridate::dmy(flu_H3$Date)
flu_H3 <- flu_H3 %>% rename(date=Date,Week=`Week number`) %>% mutate(Year=lubridate::year(date))
flu_H3 <- flu_H3 %>% 
  select(-`Influenza A(not subtype)`) %>%
  pivot_longer(-c(date,Week,Season,Year)) %>% 
  rename(subtype=name,percentage=value) %>% 
  group_by(date,Week,Season,Year) %>%
  mutate(percentage=percentage/sum(percentage)) %>% 
  filter(subtype == "Influenza A(H3N2)") %>%
  select(-subtype) %>%
  rename(percentage_h3=percentage)

## Shift to end of reporting period
flu_H3$date <- flu_H3$date + 6

final_dataset <- left_join(ili_cases_comb_expanded_grouped,flu_pos_comb %>% select(-N)) %>% drop_na() %>%
  left_join(flu_H3) %>% drop_na() %>%
  mutate(
    Influenza_cases = ILI*positivity,
    ILIplus=ILI*positivity*percentage_h3)

final_dataset$group <- factor(final_dataset$group, levels=desired_age_groups)


ggplot(final_dataset) + 
  geom_line(aes(x=date,y=ILI,col="ILI")) +
  geom_line(aes(x=date,y=ILI*positivity,col="Influenza cases")) +
  geom_line(aes(x=date,y=ILIplus,col="A/H3N2 cases")) +
  facet_wrap(~group)

p_final <- final_dataset %>% 
  select(date,group,ILI,Influenza_cases,ILIplus) %>% 
  pivot_longer(-c(group,date)) %>% 
  mutate(name=factor(name,levels=c("ILI","Influenza_cases","ILIplus"))) %>%
  ggplot() + 
  geom_line(aes(x=date,y=value,col=group)) + 
  facet_wrap(~name,scales="free",ncol=1) + 
  theme_use +
  scale_colour_brewer(palette="Set1",name="Age group") 

write_csv(final_dataset,"data/final/ili_plus_datasets_by_age.csv")


## Compare to UKHSA influenza case counts England 2009-2025
influenza_cases_eng <- read_csv("data/ukhsa/case_counts_england_2009_2025.csv")
influenza_cases_eng$Date <- lubridate::dmy(influenza_cases_eng$Date)
influenza_cases_eng <- influenza_cases_eng %>% rename(date=Date,Week=`Week number`) %>% mutate(Year=lubridate::year(date))
write_csv(influenza_cases_eng,"data/final/england_h3_cases_historic.csv")

## Compare to WHO FluNet data
flunet_data <- read_csv("data/WHO_FluNet/England_All Sites_02Jan2012_27Oct2025.csv") %>% select(-1)
colnames(flunet_data) <- c("country","surv_type","year_week","week_start","N","flu_pos","flu_neg","H1N1pdm","H3","not_subtyped")
## Create integer time index from the weeks
t_start <- flunet_data$week_start[1]
time_key <- data.frame(week_start=seq(from=as.Date(t_start),
                                      to=as.Date("2026-01-01"),
                                      by="week"))
time_key$index <- 1:nrow(time_key)
flunet_data <- flunet_data %>% left_join(time_key)

## Redistribute unsubtyped flu counts proportionally to subtyped counts
flunet_data <- flunet_data %>%
  mutate(H3 = if_else(is.na(H3),0,H3),
         H1n1pdm = if_else(is.na(H1N1pdm),0,H1N1pdm)) %>%
  mutate(H3_prop = H3/(H1N1pdm+H3)) %>%
  mutate(H3_prop = if_else(is.na(H3_prop),0,H3_prop)) %>%
  mutate(H3_adj = H3 + round(not_subtyped * H3_prop))

flunet_data <- flunet_data %>% group_by(week_start) %>% summarize(H3_sum = sum(H3_adj),flu_pos_sum=sum(flu_pos))
flunet_data <- flunet_data %>% rename(date=week_start) %>% mutate(Year=lubridate::year(date),
                                                        Week=lubridate::isoweek(date))
flunet_data$date <- flunet_data$date + 6
write_csv(flunet_data,"data/final/flunet_h3_cases_historic.csv")

p1 <- ggplot(flunet_data) + geom_line(aes(x=date,y=H3_sum,col="WHO FluNet"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng,aes(x=date,y=`Influenza A H3N2`,col="UKHSA"),linewidth=0.65) +
  geom_line(data=final_dataset %>% group_by(date) %>% summarize(ILIplus=sum(ILIplus)),aes(x=date,y=ILIplus,col="RCGP ILI+"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week") +
  ylab("Estimated number of H3\n cases or samples") +
  scale_colour_brewer(palette="Set1",name="Data source")
p2 <- ggplot(flunet_data %>% filter(Year >= 2023)) + geom_line(aes(x=date,y=H3_sum,col="WHO FluNet"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng%>% filter(Year >= 2023),aes(x=date,y=`Influenza A H3N2`,col="UKHSA"),linewidth=0.65) +
  geom_line(data=final_dataset %>% filter(Year >= 2023) %>% group_by(date) %>% summarize(ILIplus=sum(ILIplus)),aes(x=date,y=ILIplus,col="RCGP ILI+"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week") +
  ylab("Estimated number of H3\n cases or samples") +
  scale_colour_brewer(palette="Set1",name="Data source")
p_flunet <- p1/p2



## Plots
p_ili_all
p_ili_by_age
p_influenza_cases_age
p_final

p_flunet

ggsave("figures/raw_data/all_indicators.png",p_final,width=8,height=8,units="in",dpi=300)
ggsave("figures/raw_data/p_flunet.png",p_flunet,width=8,height=6,units="in",dpi=300)
