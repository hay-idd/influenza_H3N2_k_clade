## Combines the various data streams that Punya and Alex collated into ILI+ indicators

library(tidyverse)
library(purrr)
library(patchwork)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
source("R/funcs.R")

save_wd <- "figures/raw_data/"

theme_use <- theme_bw()

date_max <- as.Date("2025-12-08")

desired_age_groups <- c("0-4","5-18","19-64","65+")

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
ili_cases_recent <- read_csv("data/rcgp_ili_2023_2025.csv") %>% fill(Year,.direction="down") %>% arrange(Year,Week) %>% drop_na()
ili_cases_older <- read_csv("data/rcgp_ili_2021_2023.csv")%>% 
  fill(Year,.direction="down")%>% 
  arrange(Year,Week)

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

ili_cases_comb_expanded_grouped <- ili_cases_comb_expanded_grouped %>%
  mutate(Season = flu_season(date))
write_csv(ili_cases_comb_expanded_grouped,"data/rcgp_ili_by_age_detailed.csv")

######################################################
## Second dataset -- % of tests which are positive for influenza and by subtype from RCGP
######################################################
## Read in flu positivity data by age overall
## 2022-23 
rcgp_flu <- read_csv("data/rcgp_influenza_pos_2022_2023.csv")
rcgp_flu <- rcgp_flu %>% mutate(date=lubridate::dmy(WeekCommencing))
rcgp_flu_22_23 <- rcgp_flu %>% mutate(total_tests = `Influenza` + `COVID-19` + `RSV` + `hMPV` + `Entero/Rhinovirus` + `Adenovirus` + `Other Coronavirus` + `Awaiting results` + `None detected`) %>%
  select(`Influenza`,total_tests,date,`age group`) %>% 
  rename(number_of_influenza_samples = `Influenza`,total_samples=total_tests,age=`age group`)
age_key1 <- c("less than 5"="<5", "19-64"="19-64",   "05to 18"="5-18", "65+"="65+")
rcgp_flu_22_23$age <- age_key1[rcgp_flu_22_23$age]
## 2024-25 
rcgp_flu <- read_csv("data/rcgp_influenza_pos_2024_2025.csv")
rcgp_flu <- rcgp_flu %>% mutate(date=lubridate::dmy(date_start))
rcgp_flu <- rcgp_flu %>% mutate(age=if_else(age == "19-65","19-64",age))
rcgp_flu_24_25 <- rcgp_flu %>% 
  select(number_of_influenza_samples,total_samples, date, age)
rcgp_flu <- bind_rows(rcgp_flu_22_23, rcgp_flu_24_25)

rcgp_flu <- rcgp_flu %>% 
  mutate(flu_smooth = zoo::rollmean(number_of_influenza_samples,k=3,fill=NA,align="right")) %>%
  mutate(total_smooth = zoo::rollmean(total_samples,k=3,fill=NA,align="right")) %>%
  mutate(flu_prop  = number_of_influenza_samples/total_samples) %>% 
  mutate(flu_prop_smooth = flu_smooth/total_smooth) %>%
  select(flu_prop,flu_prop_smooth,date,age) %>%
  mutate(flu_prop = if_else(is.na(flu_prop),0,flu_prop))

## Get subtype proportions
## 2022/23 data
rcgp_subtype_22_23 <- read_csv("data/rcgp_subtype_pos_2022_2023.csv")
rcgp_subtype_22_23$date <- lubridate::dmy(rcgp_subtype_22_23$Week_commencing)
rcgp_subtype_22_23 <- rcgp_subtype_22_23 %>% rename(age=age_group)
## Get proportion of flu samples within each subtype
rcgp_subtype_22_23 <- rcgp_subtype_22_23 %>% rowwise() %>% 
  mutate(N = sum(H1N1, H3, A, B,na.rm=TRUE)) %>% 
  mutate(A = if_else(is.na(A),0,A))%>%
  mutate(N = if_else(is.na(N),0,N))%>%
  mutate(H3 = if_else(is.na(H3),0,H3))%>%
  mutate(H1N1 = if_else(is.na(H1N1),0,H1N1))%>%
  mutate(B = if_else(is.na(B),0,B)) %>%
  ## Convert to percentages
  mutate(percentage_h3 = if_else(N==0,0,H3/N)) %>%
  mutate(percentage_h1n1 = if_else(N==0,0,H1N1/N)) %>%
  mutate(percentage_a = if_else(N==0,0,A/N)) %>%
  mutate(percentage_b= if_else(N==0,0,B/N)) %>%
  mutate(extra_h3 = percentage_a * percentage_h3 / (percentage_h3 + percentage_h1n1)) %>%
  mutate(extra_h3 = if_else(is.na(extra_h3),0,extra_h3)) %>%
  mutate(extra_h1n1 = percentage_a * percentage_h1n1 / (percentage_h3 + percentage_h1n1)) %>%
  mutate(extra_h1n1 = if_else(is.na(extra_h1n1),0,extra_h1n1)) %>%
  mutate(percentage_h3_new = percentage_h3 + extra_h3) %>%
  mutate(percentage_h1n1_new = percentage_h1n1 + extra_h1n1) %>%
  select(date,age,percentage_h3,percentage_h3_new,percentage_h1n1,percentage_h1n1_new,percentage_b)

## Get 2024/25 data
rcgp_subtype_24_25 <- read_csv("data/rcgp_subtype_pos_2024_2025.csv")
rcgp_subtype_24_25 <- rcgp_subtype_24_25 %>% mutate(date=lubridate::ymd(date_start))
rcgp_subtype_24_25 <- rcgp_subtype_24_25 %>% 
  pivot_wider(names_from=virus,values_from=percent) %>%
  rowwise() %>% 
  mutate(N = sum(H1, H3, A, fluB,na.rm=TRUE)) %>% 
  mutate(A = if_else(is.na(A),0,A))%>%
  mutate(N = if_else(is.na(N),0,N))%>%
  mutate(H3 = if_else(is.na(H3),0,H3))%>%
  mutate(H1 = if_else(is.na(H1),0,H1))%>%
  mutate(fluB = if_else(is.na(fluB),0,fluB)) %>%
  ## Convert to percentages
  mutate(percentage_h3 = if_else(N==0,0,H3/N)) %>%
  mutate(percentage_h1n1 = if_else(N==0,0,H1/N)) %>%
  mutate(percentage_a = if_else(N==0,0,A/N)) %>%
  mutate(percentage_b= if_else(N==0,0,fluB/N)) %>%
  mutate(extra_h3 = percentage_a * percentage_h3 / (percentage_h3 + percentage_h1n1)) %>%
  mutate(extra_h3 = if_else(is.na(extra_h3),0,extra_h3)) %>%
  mutate(extra_h1n1 = percentage_a * percentage_h1n1 / (percentage_h3 + percentage_h1n1)) %>%
  mutate(extra_h1n1 = if_else(is.na(extra_h1n1),0,extra_h1n1)) %>%
  mutate(percentage_h3_new = percentage_h3 + extra_h3) %>%
  mutate(percentage_h1n1_new = percentage_h1n1 + extra_h1n1) %>%
  select(date,age,percentage_h3,percentage_h3_new,percentage_h1n1,percentage_h1n1_new,percentage_b)
age_key1 <- c("<5"="0-4","5"="5-18","18"="19-64","65+"="65+")
rcgp_subtype_24_25$age <- age_key1[rcgp_subtype_24_25$age]

rcgp_subtype <- bind_rows(rcgp_subtype_22_23, rcgp_subtype_24_25)%>%
  mutate(age = if_else(age == "<5","0-4",age)) %>%
  ungroup() %>%
  mutate(percentage_h3_new = zoo::rollmean(percentage_h3_new,k=3,fill=NA,align="right"))

final_dataset <- left_join(rcgp_subtype, rcgp_flu%>%
                              mutate(age = if_else(age == "<5","0-4",age)) ) %>% rename(group=age) %>% left_join(ili_cases_comb_expanded_grouped %>% mutate(date = date + 1)) %>%
  mutate(Influenza = ILI * flu_prop_smooth) %>%
  mutate(ILI_plus = ILI * flu_prop_smooth * percentage_h3_new) 
final_dataset_app <- final_dataset %>% 
  filter(date >= "2022-01-01") %>%
  filter(date <= "2023-12-31") %>%
  select(
  Year, Week, date, group, ILI, ILI_per_100k, N, Influenza, flu_prop, flu_prop_smooth,N, ILI_plus
) %>%
  rename(age_group=group,flu_samples=Influenza,prop_flu=flu_prop,smooth_prop=flu_prop_smooth,total_samples=N,ILI_flu=ILI_plus)
write_csv(final_dataset_app,"data/rcgp_ili_flu_by_age_for_app.csv")

## Recode Influenza and ILI_plus to All influenza cases and A/H3N2 cases
final_dataset <- final_dataset %>% select(Year,Week, date,group,ILI_per_100k, N, ILI, flu_prop,Influenza,,ILI_plus) %>%
  rename(`Proportion influenza positive`=flu_prop,
         `All influenza cases`=Influenza,
         `Influenza A/H3N2`=ILI_plus)

final_dataset <- final_dataset %>% mutate(Season = flu_season(date))

p_final <- final_dataset %>% 
  filter(date < date_max) %>%
  select(date,group,ILI,`All influenza cases`,`Influenza A/H3N2`) %>% 
  pivot_longer(-c(group,date)) %>% 
  mutate(name=factor(name,levels=c("ILI","All influenza cases","Influenza A/H3N2"))) %>%
  ggplot() + 
  geom_line(aes(x=date,y=value,col=group)) + 
  facet_wrap(~name,scales="free",ncol=1) + 
  xlab("Date") +
  ylab("Imputed cases") +
  theme_use +
  scale_colour_brewer(palette="Set1",name="Age group") 


write_csv(final_dataset,"data/ili_plus_datasets_by_age.csv")

######################################################
## Third dataset -- Respiratory DataMart cases
######################################################
## Compare to UKHSA influenza case counts England 2009-2025
influenza_cases_eng <- read_csv("data/raw/resp_datamart_all_flu.csv")
influenza_cases_eng$date <- lubridate::dmy(influenza_cases_eng$date)
influenza_cases_eng <- influenza_cases_eng %>% filter(date < date_max) %>%
  mutate(total_cases = `Influenza A H1N1pdm09` + `Influenza A H3N2` + `Influenza A not subtyped` + `Influenza B`)

## Distribute not subtyped influenza A cases proportionally to H1N1pdm09 and H3N2
influenza_cases_eng <- influenza_cases_eng %>%
  mutate(H1_prop = `Influenza A H1N1pdm09` / (`Influenza A H1N1pdm09` + `Influenza A H3N2`)) %>%
  mutate(H1_prop = if_else(is.na(H1_prop),0,H1_prop)) %>%
  mutate(H3_prop = `Influenza A H3N2` / (`Influenza A H1N1pdm09` + `Influenza A H3N2`)) %>%
  mutate(H3_prop = if_else(is.na(H3_prop),0,H3_prop)) %>%
  mutate(H1_adj = `Influenza A H1N1pdm09` + round(`Influenza A not subtyped` * H1_prop)) %>%
  mutate(H3_adj = `Influenza A H3N2` + round(`Influenza A not subtyped` * H3_prop))

influenza_cases_eng <- influenza_cases_eng %>% mutate(Season = flu_season(date)) 

write_csv(influenza_cases_eng,"data/resp_datamart_influenza_cases_england.csv")

######################################################
## Fourth dataset -- Respiratory DataMart cases
######################################################
## Compare to WHO FluNet data
flunet_data <- read_csv("data/raw/FlunetData_United Kingdom, England_All Sites_for_03 January 2011 to 15 December 2025.csv") %>% select(-1) %>% select(-c(11))

colnames(flunet_data) <- c("country","surv_type","year_week","week_start","N","flu_pos","flu_neg","H1N1pdm","H3","not_subtyped")
## Create integer time index from the weeks
t_start <- flunet_data$week_start[1]
time_key <- data.frame(week_start=seq(from=as.Date(t_start),
                                      to=as.Date("2026-01-01"),
                                      by="week"))
time_key$index <- 1:nrow(time_key)
flunet_data <- flunet_data %>% left_join(time_key)

## Create flu seasons from June to June
flunet_data <- flunet_data %>% mutate(Year = lubridate::year(week_start),
                                      Month = lubridate::month(week_start)) %>%
  mutate(Season = if_else(Month >=6, paste0(Year,"-",Year+1), paste0(Year-1,"-",Year))) %>%
  select(-c(Year,Month))

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
flunet_data <- flunet_data %>% filter(date < date_max)
flunet_data <- flunet_data %>% mutate(Season = flu_season(date))

write_csv(flunet_data,"data/WHO_flunet_cases.csv")

p1 <- ggplot(flunet_data) + geom_line(aes(x=date,y=H3_sum,col="WHO FluNet A/H3N2"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng,aes(x=date,y=`Influenza A H3N2`,col="UKHSA A/H3N2"),linewidth=0.65) +
  geom_line(data=final_dataset %>% group_by(date) %>% summarize(ILIplus=sum(`Influenza A/H3N2`)),aes(x=date,y=ILIplus,col="RCGP ILI+ A/H3N2"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week)") +
  ylab("Estimated number of A/H3N2 cases") +
  scale_colour_brewer(palette="Set1",name="Data source")

p2 <- ggplot(flunet_data %>% filter(Year >= 2025)) + geom_line(aes(x=date,y=H3_sum,col="WHO FluNet A/H3N2"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng%>% filter(Year >= 2025),aes(x=date,y=`Influenza A H3N2`,col="UKHSA A/H3N2"),linewidth=0.65) +
  geom_line(data=final_dataset %>% filter(Year >= 2025) %>% group_by(date) %>% summarize(ILIplus=sum(`Influenza A/H3N2`)),aes(x=date,y=ILIplus,col="RCGP ILI+ A/H3N2"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week)") +
  ylab("Estimated number of A/H3N2 cases") +
  scale_colour_brewer(palette="Set1",name="Data source")


p3 <- ggplot(flunet_data) + geom_line(aes(x=date,y=flu_pos_sum,col="WHO FluNet total"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng,aes(x=date,y=total_cases,col="UKHSA total"),linewidth=0.65) +
  geom_line(data=final_dataset %>% group_by(date) %>% summarize(ILIplus=sum(`All influenza cases`)),aes(x=date,y=ILIplus,col="RCGP ILI+"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week)") +
  ylab("Estimated number of influenza cases") +
  scale_colour_brewer(palette="Set1",name="Data source")

p4 <- ggplot(flunet_data %>% filter(Year >= 2025)) + geom_line(aes(x=date,y=H3_sum,col="WHO FluNet total"),linewidth=0.65) +
  geom_line(data=influenza_cases_eng%>% filter(Year >= 2025),aes(x=date,y=total_cases,col="UKHSA total"),linewidth=0.65) +
  geom_line(data=final_dataset %>% filter(Year >= 2025) %>% group_by(date) %>% summarize(ILIplus=sum(`All influenza cases`)),aes(x=date,y=ILIplus,col="RCGP ILI+ total"),linewidth=0.65) +
  theme_use +
  xlab("Date (start of Epi week)") +
  ylab("Estimated number of influenza cases") +
  scale_colour_brewer(palette="Set1",name="Data source")

figS1 <- p1/p2
figS2 <- p3/p4
figS3 <- p_final

## Plots
ggsave("figures/figS1.png",figS1,width=8,height=6,units="in",dpi=300)
ggsave("figures/figS2.png",figS2,width=8,height=6,units="in",dpi=300)
ggsave("figures/figS3.png",figS3,width=8,height=8,units="in",dpi=300)
