library(tidyverse)
library(ggplot2)

flu_age <- read_csv("data/respmart_flu_pos_by_age.csv") %>% rename(date=Date)
head(flu_age)


flu_age <- flu_age %>% mutate(date = dmy(date))

flu_age <- flu_age %>% mutate(Season = if_else(date >= "2025-07-01","2025 to 2026",
                                               if_else(date >= "2024-07-01","2024 to 2025",
                                                       if_else(date >= "2023-07-01","2023 to 2024",
                                                               "2022 to 2023")))) %>% filter(Season != "2023 to 2024")

flu_age <- flu_age %>% pivot_longer(-c(Season,date,`Week number`)) %>% rename(Age=name,Positivity=value)

flu_age$day <- yday(flu_age$date)
flu_age <- flu_age %>% group_by(Season) %>% mutate(day_from_start = date - min(date))
flu_age$Positivity <- pmax(flu_age$Positivity, 0.01)  ## Avoid zeros for log growth rate calculation
flu_age <-  flu_age %>% group_by(Season,Age) %>% arrange(date) %>%
  mutate(growth_rate = log(Positivity/lag(Positivity))) %>%
  mutate(growth_rate_smooth1 = log(zoo::rollmean(Positivity, k=5, fill=NA, align="center")/
                                     lag(zoo::rollmean(Positivity, k=5, fill=NA, align="center")))) %>%
  mutate(growth_rate_smooth2 = zoo::rollmean(growth_rate_smooth1, k=5, fill=NA, align="center"))


flu_age <- flu_age %>% left_join(day_match) %>% mutate(day_rel = date - date_match)
flu_age <- flu_age %>% group_by(Season) %>% mutate(day_from_start = date - min(date))
flu_age$Age <- factor(flu_age$Age,levels=c("Up to 5 years", "5 to 14 years", "15 to 44 years", "45 to 64 years", 
                                           "65 to 79 years", "80 years and above"))

flu_age <- flu_age %>% filter(Season != "2024 to 2025", Age %in% c("Up to 5 years", "5 to 14 years", "15 to 44 years", "45 to 64 years"))

p_by_age_from_start <- ggplot(flu_age) + 
  geom_line(aes(x=day_from_start,y=Positivity,col=Season)) + 
  xlab("Days from start of season (1st July)") +
  ylab("Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)

save_plots(p_by_age_from_start,"figures/","flu_positivity_by_age_from_start",width=8,height=6)


p_by_age_rel <- ggplot(flu_age %>% filter(day_rel > -100)) + geom_line(aes(x=day_rel,y=Positivity,col=Season)) +  
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)
save_plots(p_by_age_rel,"figures/","flu_positivity_by_age_relative",width=8,height=6)


p_gr_by_age_rel <- ggplot(flu_age) + 
  geom_line(aes(x=day_from_start,y=growth_rate_smooth1,col=Season)) +  
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)
save_plots(p_gr_by_age_rel,"figures/","flu_growth_rate_by_age_relative",width=8,height=6)
