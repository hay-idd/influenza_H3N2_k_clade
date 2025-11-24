## Check if peak, timing or second derivative of growth rate correlates with cumulative hospitalisations for each season

library(tidyverse)
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")

hosp_dat <- read_csv('data/hospitalisations/weekly_hospital_admissions_flu_overall.csv')

hosp_dat$date <- lubridate::dmy(hosp_dat$Date)
ggplot(hosp_dat) + geom_line(aes(x=date,y=Rate,group=Season))
hosp_dat_summ <- hosp_dat %>% group_by(Season) %>% summarize(total_hosp = sum(Rate))

## H3 seasons
Season <- c("2016/17", "2017/18", "2018/19", "2022/23", "2023/24", "2024/25", "2025/26")
Season1 <- c("2016 to 2017", "2017 to 2018", "2018 to 2019", "2022 to 2023", "2023 to 2024", "2024 to 2025", "2025 to 2026")
Dominant <- c("H3", "B", "H1", "H3", "Mixed", "H1", "H3")
Season_key <- data.frame(Season=Season1,season=Season,Dominant=Dominant)

## Compare to flunet growth rates
flunet_grs <- read_csv("results/flunet_all_growth_rates.csv")
get_peak_grs <- flunet_grs %>% group_by(season) %>% filter(y == max(y)) %>% select(y,season,day_of_year)

## Does peak growth rate predict total hospitalisations?
fit_dat <- left_join(hosp_dat_summ, Season_key) %>% select(-Season) %>% left_join(get_peak_grs) %>%
  filter(Season != "2025/26")
p_hosp_regression_flunet_peak <-   ggplot(fit_dat) + geom_point(aes(x=y,y=log(total_hosp)))+ geom_smooth(aes(x=y,y=log(total_hosp)),method="lm") +
    theme_bw() +
    xlab("Peak growth rate") +
    ylab("Log cumulative hospitalisation rate (per 100,000)")

summary(lm(log(total_hosp) ~ y, data=fit_dat))

## Does timing of peak growth rate predict total hospitalisations?
p_hosp_regression_flunet_peak_time <- left_join(hosp_dat_summ, Season_key) %>% select(-Season) %>% left_join(get_peak_grs) %>%
  filter(Season != "2025/26") %>%
  ggplot() + geom_point(aes(x=day_of_year,y=log(total_hosp)))+ geom_smooth(aes(x=day_of_year,y=log(total_hosp)),method="lm") +
  theme_bw() +
  xlab("Peak growth rate") +
  ylab("Log cumulative hospitalisation rate (per 100,000)")

summary(lm(log(total_hosp) ~ day_of_year, data=fit_dat))

ggsave(p_hosp_regression_flunet_peak, filename="figures/hospitalisations_vs_flunet_peak_gr.png", width=5, height=3)
ggsave(p_hosp_regression_flunet_peak_time, filename="figures/hospitalisations_vs_flunet_peak_gr_time.png", width=5, height=3)



## Compare to UKHSA growth rates
flunet_grs <- read_csv("results/ukhsa_all_flu_growth_rates.csv")
get_peak_grs <- flunet_grs %>% group_by(season) %>% filter(y == max(y)) %>% mutate(day_of_year=yday(date)) %>% select(y,season,day_of_year)

## Does peak growth rate predict total hospitalisations?
fit_dat <- left_join(hosp_dat_summ, Season_key) %>% select(-Season) %>% left_join(get_peak_grs) %>%
  filter(Season != "2025/26")
p_hosp_regression_ukhsa_peak <-   ggplot(fit_dat) + geom_point(aes(x=y,y=log(total_hosp)))+ geom_smooth(aes(x=y,y=log(total_hosp)),method="lm") +
  theme_bw() +
  xlab("Peak growth rate") +
  ylab("Log cumulative hospitalisation rate (per 100,000)")

summary(lm(log(total_hosp) ~ y, data=fit_dat))

## Does timing of peak growth rate predict total hospitalisations?
p_hosp_regression_ukhsa_peak_time <- left_join(hosp_dat_summ, Season_key) %>% select(-Season) %>% left_join(get_peak_grs) %>%
  filter(Season != "2025/26") %>%
  ggplot() + geom_point(aes(x=day_of_year,y=log(total_hosp)))+ geom_smooth(aes(x=day_of_year,y=log(total_hosp)),method="lm") +
  theme_bw() +
  xlab("Peak growth rate") +
  ylab("Log cumulative hospitalisation rate (per 100,000)")

summary(lm(log(total_hosp) ~ day_of_year, data=fit_dat))

ggsave(p_hosp_regression_ukhsa_peak, filename="figures/hospitalisations_vs_ukhsa_peak_gr.png", width=5, height=3)
ggsave(p_hosp_regression_ukhsa_peak_time, filename="figures/hospitalisations_vs_ukhsa_peak_gr_time.png", width=5, height=3)