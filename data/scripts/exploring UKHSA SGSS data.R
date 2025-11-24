library(tidyverse)
library(ggplot2)

flu_pos <- read_csv("data/influenza_positivity_england.csv")
head(flu_pos)

flu_pos$date <- dmy(flu_pos$Date)
flu_pos$day <- yday(flu_pos$date)
flu_pos <- flu_pos %>% group_by(Season) %>% mutate(day_from_start = date - min(date))
flu_pos$smooth_pos <- zoo::rollmean(flu_pos$`Positivity (%)`, k=5, fill=NA, align="center")
flu_pos$growth_rate <- log(flu_pos$`Positivity (%)`/lag(flu_pos$`Positivity (%)`))
flu_pos$growth_rate_smooth1 <- log(flu_pos$smooth_pos/lag(flu_pos$smooth_pos))
flu_pos <- flu_pos %>% mutate(growth_rate_smooth2 = zoo::rollmean(growth_rate_smooth1, k=5, fill=NA, align="center"))

## Look at % positivity over time by same day of year
p_pos_all <- ggplot(flu_pos %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=`Positivity (%)`,x=day_from_start,col=Season,linewidth=is_2025)) +
  xlab("Days from start of season (30th June - 4th July)") +
  ylab("Influenza Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2023 to 2024"="#3382CD","2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical")

## Try to align when % positivity is most similar
day_match <- flu_pos %>% group_by(Season) %>% mutate(diff_pos = abs(`Positivity (%)` - 11.8)) %>% filter(day > 300 )%>% filter(diff_pos == min(diff_pos), na.rm=TRUE) %>% select(Season, day,date) %>% rename(day_match = day, date_match=date)

flu_pos <- flu_pos %>% left_join(day_match) %>% mutate(day_rel = day - day_match)
p_pos_rel <- ggplot(flu_pos %>% filter(day_rel > -100) %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=`Positivity (%)`,x=day_rel,col=Season,linewidth=is_2025)) +
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Influenza Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                             "2023 to 2024"="#3382CD","2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical")

## Look at growth rates
p_gr_rel <- ggplot(flu_pos %>% filter(day_rel > -100) %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=growth_rate_smooth1,x=day_rel,col=Season,linewidth=is_2025))+
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Smoothed growth rate") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2023 to 2024"="#3382CD","2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical")


p_gr_all <- ggplot(flu_pos %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=growth_rate_smooth1,x=day_from_start,col=Season,linewidth=is_2025))+
  xlab("Days from start of season (30th June - 4th July)") +
  ylab("Smoothed growth rate") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2023 to 2024"="#3382CD","2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical")

save_plots(p_pos_all,"figures/","flu_positivity_all_seasons",width=8,height=4)
save_plots(p_pos_rel,"figures/","flu_positivity_aligned",width=8,height=4)
save_plots(p_gr_all,"figures/","flu_growth_rate_all_seasons",width=8,height=4)
save_plots(p_gr_rel,"figures/","flu_growth_rate_aligned",width=8,height=4)

## By age SGSS
flu_age <- read_csv("data/influenza_by_age_england.csv")
flu_age <- flu_age %>% mutate(Date = dmy(Date)) %>% pivot_longer(-c(Date)) %>% rename(Age = name, Positivity = value)

head(flu_age)
## Label season
flu_age <- flu_age %>% mutate(Season = if_else(Date >= "2025-07-01","2025 to 2026",
                                        if_else(Date >= "2024-07-01","2024 to 2025",
                                        if_else(Date >= "2023-07-01","2023 to 2024",
                                        "2022 to 2023")))) %>% filter(Season != "2023 to 2024")
flu_age <- flu_age %>% left_join(day_match) %>% mutate(day_rel = Date - date_match)
flu_age <- flu_age %>% group_by(Season) %>% mutate(day_from_start = Date - min(Date))
flu_age$Age <- factor(flu_age$Age,levels=c("0 to 4 years", "5 to 14 years", "15 to 24 years", "25 to 44 years", 
                                              "45 to 54 years", "55 to 64 years", "65 to 74 years", "75 to 84 years", 
                                              "85 years and above"))

flu_age %>% ggplot() + geom_line(aes(y=Positivity,x=day_from_start,col=Season)) + facet_wrap(~Age)
p_age_pos_rel <- ggplot(flu_age %>% filter(day_rel > -100) %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=Positivity,x=day_rel,col=Season,linewidth=is_2025)) +
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Influenza Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)
save_plots(p_age_pos_rel,"figures/","flu_positivity_by_age_aligned",width=8,height=6)

p_age_pos_all <- ggplot(flu_age %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=Positivity,x=day_from_start,col=Season,linewidth=is_2025)) +
  xlab("Days from start of season (1st July)") +
  ylab("Influenza Positivity (%)") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)
save_plots(p_age_pos_all,"figures/","flu_positivity_by_age_all_seasons",width=8,height=6)
## Find growth rates by age
flu_age <- flu_age %>% group_by(Season,Age) %>% arrange(Date) %>%
  mutate(growth_rate = log(Positivity/lag(Positivity))) %>%
  mutate(growth_rate_smooth1 = log(zoo::rollmean(Positivity, k=5, fill=NA, align="center")/
                                     lag(zoo::rollmean(Positivity, k=5, fill=NA, align="center")))) %>%
  mutate(growth_rate_smooth2 = zoo::rollmean(growth_rate_smooth1, k=5, fill=NA, align="center"))
p_age_gr_rel <- ggplot(flu_age %>% filter(day_rel > -100) %>% mutate(is_2025 = Season == "2025 to 2026")) +
  geom_line(aes(y=growth_rate_smooth1,x=day_rel,col=Season))+
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Smoothed growth rate") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)
p_age_gr_rel_zoomed <- ggplot(flu_age %>% filter(day_rel > -50, day_rel <= 0) %>% mutate(is_2025 = Season == "2025 to 2026")) +
  geom_line(aes(y=growth_rate_smooth1,x=day_rel,col=Season))+
  xlab("Day relative to day with Positivity % closest to 11.8%") +
  ylab("Smoothed growth rate") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)

save_plots(p_age_gr_rel,"figures/","flu_growth_rate_by_age_aligned",width=8,height=6)
save_plots(p_age_gr_rel_zoomed,"figures/","flu_growth_rate_by_age_aligned_zoomed",width=8,height=6)
p_age_gr_all <- ggplot(flu_age %>% mutate(is_2025 = Season == "2025 to 2026")) +
  geom_line(aes(y=growth_rate_smooth1,x=day_from_start,col=Season))+
  xlab("Days from start of season (1st July)") +
  ylab("Smoothed growth rate") +
  theme_psi() + 
  scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
                                        "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Age)

save_plots(p_age_gr_all,"figures/","flu_growth_rate_by_age_all_seasons",width=8,height=6)


p_age_pos_alt <- ggplot(flu_age %>% mutate(is_2025 = Season == "2025 to 2026")) + 
  geom_line(aes(y=Positivity,x=day_from_start,col=Age)) +
  xlab("Days from start of season (1st July)") +
  ylab("Influenza Positivity (%)") +
  theme_psi() + 
  scale_color_viridis_d() +
  #scale_color_manual("Season", values=c("2025 to 2026"="#EB5F17","2024 to 2025"="#A32CA5",
  #                                      "2022 to 2023"="#009786")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.5)) +
  theme(panel.grid.major.x = element_line(colour='grey90'),legend.position="bottom",legend.box="vertical") +
  facet_wrap(~Season,ncol=1)
save_plots(p_age_pos_alt,"figures/","flu_positivity_by_age_all_seasons_alt",width=8,height=6)

