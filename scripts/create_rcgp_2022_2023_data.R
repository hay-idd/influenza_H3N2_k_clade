library(tidyverse)
library(purrr)
library(patchwork)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
source("R/funcs.R")
theme_use <- theme_bw()
## Read in RCGP flu sample data
dat_flu <- read_csv("data/rcgp_digitised/raw_2022:23/virology_wide_all_extended_2022_23.csv")
dat_flu$date <- lubridate::dmy(dat_flu$WeekCommencing)
## Get total samples collected
dat_flu <- dat_flu %>% mutate(total_samples = Influenza + `COVID-19` + RSV + hMPV + `Entero/Rhinovirus` + 
                                Adenovirus + `Other Coronavirus` + `Awaiting results` + `None detected`)
## Proportion of these which were positive for flu
dat_flu <- dat_flu %>% mutate(prop_flu = Influenza/total_samples)

## Relabel age groups
age_key1 <- c("less than 5"="0-4","05to 18"="5-18","19-64"="19-64","65+"="65+")
dat_flu$age_group <- age_key1[dat_flu$`age group`]

ggplot(dat_flu) + geom_line(aes(x=date,y=prop_flu,col=age_group)) +
  labs(x="Date",y="Proportion of samples positive for Influenza",
       title="Proportion of RCGP samples positive for Influenza by age group, 2022-23 season",
       caption="Data source: RCGP RSC via UKHSA") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",name="Age group") +
  theme(text = element_text(size=16))

## Pull out key columns and do some smoothing, also report by end of week rather than start of week
dat_flu <- dat_flu %>% select(date, Influenza, total_samples, prop_flu, age_group) %>% rename(flu_samples=Influenza) %>% 
  mutate(smooth_total_samples = zoo::rollmean(total_samples,k=3,fill=NA,align="center"), 
         smooth_flu_samples = zoo::rollmean(flu_samples,k=3,fill=NA,align="center"),
         smooth_prop = smooth_flu_samples/smooth_total_samples)  %>% 
  mutate(date = date + 7)

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

desired_age_groups <- c("0-4","5-18","19-64","65+")
ili_cases_comb_expanded_grouped <- combine_age_groups_ILI(ili_cases_comb_expanded,desired_age_groups)

final_2022_2023_dat <- ili_cases_comb_expanded_grouped %>% select(Year,Week,date,group,ILI,ILI_per_100k) %>% rename(age_group=group) %>% #filter(date >= "2025-06-01", date <= "2026-06-01")  %>%
  left_join(dat_flu %>% select(date, age_group,flu_samples,prop_flu,smooth_total_samples,smooth_flu_samples,smooth_prop) %>% mutate(date = date - 1)) %>% drop_na() %>% mutate(ILI_flu = ILI * smooth_prop)

ggplot(final_2022_2023_dat) + geom_line(aes(x=date,y=ILI_flu,colour=age_group)) + theme_use + 
  xlab("Date") + ylab("Influenza-attributable ILI cases (RCGP, England)") +
  scale_x_date(date_labels="%b %Y",date_breaks="1 year") +
  scale_colour_brewer(palette="Set1",name="Age group")
write_csv(final_2022_2023_dat,"data/final/flu_2022_2023_new.csv")
## Old data
old_data <- read_csv("data/final/flu_2022_2023.csv")
head(old_data)

comb <- bind_rows(
  final_2022_2023_dat %>% select(date,age_group,ILI_flu) %>% mutate(ver="New"),
  old_data %>% select(date,age_group,ILI_flu) %>% mutate(ver="Old")
)

ggplot(comb) + geom_line(aes(x=date,y=ILI_flu,colour=ver)) + theme_use + 
  xlab("Date") + ylab("Influenza-attributable ILI cases (RCGP, England)") +
  scale_x_date(date_labels="%b %Y",date_breaks="1 year")  + facet_wrap(~age_group, scales="free_y")
