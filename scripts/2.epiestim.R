library(EpiStrainDynamics)
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)
library(EpiEstim)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")

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
overall_ILIplus <- overall_ILIplus %>% select(Year,Week,date,group,ILIplus,Season)
overall_ILIplus <- overall_ILIplus %>% bind_rows(overall_ILIplus %>% group_by(Year,Week,date,Season) %>% summarize(ILIplus = sum(ILIplus)) %>% mutate(group="All"))
overall_ILIplus <- overall_ILIplus %>% filter(Season != "2022 to 2023")


mean_si <- 3
std_si <- 1.5
method <- "parametric_si"
config <- make_config(list(mean_si = mean_si,
                           std_si = std_si))
R_est <- NULL
index <- 1
for(tmp_season in unique(overall_ILIplus$Season)){
  print(tmp_season)
  inc <- overall_ILIplus %>% filter(Season==tmp_season & group=="All") %>%
    arrange(date) %>% pull(ILIplus) %>% round()

  output <- EpiEstim::estimate_R(incid = inc,
                                 dt = 7L,
                                 recon_opt = "match",
                                 method = method,
                                 iter=20L,
                                 tol=1e-6,
                                 grid = list(precision = 0.001, min = -1, max = 1),
                                 config = config)
  p <- plot(output)
  p[[2]]$data %>% filter(start > 20, start < 50) %>% pull(meanR) %>% mean() %>% print()
  R_est[[index]] <- output$R %>%  mutate(Season=tmp_season)
  index <- index + 1
}
R_est_comb <- do.call('bind_rows',R_est)

p_R_epiestim <- ggplot(R_est_comb) + 
  geom_hline(yintercept=1,linetype="dashed") +
  geom_hline(yintercept=1.2,linetype="dashed") +
  geom_hline(yintercept=1.4,linetype="dashed") +
  geom_ribbon(aes(x=t_start,ymin=`Quantile.0.025(R)`,ymax=`Quantile.0.975(R)`,fill=Season),alpha=0.2) +
  geom_line(aes(x=t_start,y=`Mean(R)`,col=Season)) +
 
  #scale_x_date(breaks="month") +
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  ylab("Effective reproductive number") + xlab("Day (starting from 1st August)") +
  theme(legend.position="bottom")

## Try FluNet data instead

flunet <- read_csv("data/final/flunet_h3_cases_historic.csv")
flunet <- flunet %>% mutate(H3_sum = if_else(is.na(H3_sum),0,H3_sum))

flunet <- flunet %>%
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


R_est <- NULL
index <- 1
for(tmp_season in unique(flunet$season)){
  print(tmp_season)
  inc <- flunet %>% filter(season==tmp_season) %>%
    arrange(date) %>% pull(H3_sum) %>% round()
  
  output <- EpiEstim::estimate_R(incid = inc,
                                 dt = 7L,
                                 dt_out=7L,
                                 recon_opt = "naive",
                                 method = method,
                                 iter=20L,
                                 tol=1e-6,
                                 grid = list(precision = 0.001, min = -1, max = 1),
                                 config = config)
  
  R_est[[index]] <- output$R %>%  mutate(Season=tmp_season)
  index <- index + 1
}
R_est_comb <- do.call('bind_rows',R_est)
