## Fit growth rates and Rt to WHO FluNet counts back to 2012

library(tidyverse)
library(ggplot2)
devtools::load_all("~/Documents/GitHub/EpiStrainDynamics/")
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")

## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

flunet <- read_csv("data/WHO_FluNet/England_All Sites_02Jan2012_27Oct2025.csv") %>% select(-1)

head(flunet)

colnames(flunet) <- c("country","surv_type","year_week","week_start","N","flu_pos","flu_neg","H1N1pdm","H3","not_subtyped")

## Create integer time index from the weeks
t_start <- flunet$week_start[1]
time_key <- data.frame(week_start=seq(from=as.Date(t_start),
                            to=as.Date("2026-01-01"),
                            by="week"))
time_key$index <- 1:nrow(time_key)
flunet <- flunet %>% left_join(time_key)
ggplot(flunet %>% filter(surv_type == "Non-sentinel") %>% group_by(index) %>% summarize(N=sum(flu_pos))) + geom_line(aes(x=index,y=N)) 

## Redistribute unsubtyped flu counts proportionally to subtyped counts
flunet <- flunet %>%
  mutate(H3 = if_else(is.na(H3),0,H3),
         H1n1pdm = if_else(is.na(H1N1pdm),0,H1N1pdm)) %>%
  mutate(H3_prop = H3/(H1N1pdm+H3)) %>%
  mutate(H3_adj = H3 + round(not_subtyped * H3_prop))

flunet_tmp <- flunet %>% filter(surv_type == "Non-sentinel") %>% filter(week_start > "2012-01-01")
mod <- construct_model(
  #method = p_spline(spline_degree = 3, days_per_knot = 3),
  method = random_walk(),
  pathogen_structure =single(
    case_timeseries = round((flunet_tmp$H3)),           # timeseries of case data
    time = flunet_tmp$index,                       # date or time variable labels
    pathogen_name = 'flu'                # optional name of pathogen
  ),
  dow_effect = FALSE
)

fit <- fit_model(
  mod,
  iter = 2000,
  warmup = 1000,
  chains = 3
)



## Start time is 1, but we want to go back to 2013-01-01 as 0, so subtract accordingly
start_time <- as.numeric((as.Date("2013-09-30") - as.Date("2013-01-01"))/365)
date_key <- date_key %>% mutate(i = time - start_time)

t_start <- as.numeric((as.Date("2014-01-01") - min(flunet_tmp$week_start))/7)

rt <- Rt(fit, tau_max = 7, gi_dist = function(x) 4*x*exp(-2*x))

gr <- growth_rate(fit)
gr_dat <- gr$measure
gr_dat <- gr_dat %>% left_join(flunet_tmp %>% select(index,week_start),by=c("time"="index")) %>% mutate(year=year(week_start))
ggplot(gr_dat) + 
  geom_ribbon(aes(x=week_start,ymin=lb_95,ymax=ub_95,group=year),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=week_start,ymin=lb_50,ymax=ub_50,group=year),alpha=0.5,fill="blue") +
  geom_line(aes(x=week_start,y=y,group=year)) +
  scale_x_date(breaks="1 year")

library(dplyr)
library(lubridate)

# parameters
gap_threshold_weeks <- 4   # change this if you want a different gap size to break a season
season_cutoff_month <- 7   # months >= this are considered "start year" (Jul 1 cutoff -> season spans Jul-Jun)

label_seasons <- function(df, date_col = "week_start", pathogen_col = "pathogen",
                          gap_weeks = gap_threshold_weeks, cutoff_month = season_cutoff_month) {
  df %>%
    arrange(.data[[pathogen_col]], .data[[date_col]]) %>%
    group_by(.data[[pathogen_col]]) %>%
    mutate(
      # difference in weeks to previous row
      diff_weeks = as.numeric(difftime(.data[[date_col]], lag(.data[[date_col]]), units = "days")) / 7,
      # mark start of a new contiguous run (first row OR gap >= threshold)
      new_run = if_else(is.na(diff_weeks) | diff_weeks >= gap_weeks, 1L, 0L),
      # integer run id per pathogen
      run_id = cumsum(new_run)
    ) %>%
    group_by(.data[[pathogen_col]], run_id) %>%
    mutate(
      season_start = min(.data[[date_col]]),
      # derive start year: if month >= cutoff_month, use that year, else use year-1
      season_start_year = if_else(month(season_start) >= cutoff_month,
                                  year(season_start),
                                  year(season_start) - 1L),
      # label like "2013/14"
      season = sprintf("%d/%02d", season_start_year, (season_start_year + 1) %% 100)
    ) %>%
    ungroup() %>%
    select(-diff_weeks, -new_run, -run_id, -season_start_year)
}

# Apply only to flu (or to entire df if you want seasons for every pathogen)
gr_dat_labeled <- gr_dat %>%
  filter(pathogen == "flu") %>%
  label_seasons()

head(gr_dat_labeled, 10)

## For later years, just assign season manually
gr_dat_labeled <- gr_dat_labeled %>%
  mutate(season = case_when(
    week_start >= as.Date("2021-07-01") & week_start < as.Date("2022-07-01") ~ "2021/22",
    week_start >= as.Date("2022-07-01") & week_start < as.Date("2023-07-01") ~ "2022/23",
    week_start >= as.Date("2023-07-01") & week_start < as.Date("2024-07-01") ~ "2023/24",
    week_start >= as.Date("2024-07-01") & week_start < as.Date("2025-07-01") ~ "2024/25",
    week_start >= as.Date("2025-07-01") & week_start < as.Date("2026-07-01") ~ "2025/26",
    TRUE ~ season
  ))


ggplot(gr_dat_labeled) + 
  geom_ribbon(aes(x=week_start,ymin=lb_95,ymax=ub_95,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=week_start,ymin=lb_50,ymax=ub_50,group=season),alpha=0.5,fill="blue") +
  geom_line(aes(x=week_start,y=y,group=season)) +
  scale_x_date(breaks="1 year")


## Align by day of year
gr_dat_labeled$day_of_year <- yday(gr_dat_labeled$week_start)
gr_dat_labeled <- gr_dat_labeled %>% group_by(season) %>%
  mutate(day_of_year = if_else(day_of_year < 183, day_of_year + 365, day_of_year)) %>%
  ungroup()

## Create labels for plot, label seasons by pre-pandemic, during pandemic (2020 and 2021), post pandemic (2022-2024) and current season (2025/26)
gr_dat_labeled <- gr_dat_labeled %>%
  mutate(plot_label = case_when(
    season %in% c("2012/13","2013/14","2014/15","2015/16","2016/17","2017/18","2018/19","2019/20") ~ "Pre-pandemic",
    season %in% c("2020/21","2021/22") ~ "During pandemic",
    season %in% c("2022/23","2023/24") ~ "Post-pandemic",
    season == "2024/25" ~ "Recent season",
    season == "2025/26" ~ "Current season",
    TRUE ~ "Other"
  ))

ggplot(gr_dat_labeled) + 
  #geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.5) +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label))


## Align by peak data
max_gr <- gr_dat_labeled %>% group_by(season) %>% filter(y == max(y)) %>% select(season, day_of_year) %>% rename(peak_time = day_of_year)
gr_dat_labeled <- gr_dat_labeled %>% left_join(max_gr)
gr_dat_labeled$day_shifted <- gr_dat_labeled$day_of_year - gr_dat_labeled$peak_time


## 
p2 <- ggplot(gr_dat_labeled) + 
  geom_rect(aes(xmin=300,xmax=300-7*4,ymin=-Inf,ymax=Inf),fill="lightgray",alpha=0.5) +
  #geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=day_of_year,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_line(aes(x=day_of_year,y=y,group=season,col=plot_label)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw()+
  ylab("Growth rate") +
  xlab("Day of the year")

p1 <- ggplot(gr_dat_labeled) + 
  #geom_ribbon(aes(x=day_of_year,ymin=lb_95,ymax=ub_95,group=season),alpha=0.25,fill="blue") +
  geom_ribbon(aes(x=day_shifted,ymin=lb_50,ymax=ub_50,group=season,fill=plot_label),alpha=0.1) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_line(aes(x=day_shifted,y=y,group=season,col=plot_label)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +
  ylab("Growth rate") +
  xlab("Day of the year")

add_doubling_axis <- function(p){
  # Get the y-axis range from the built plot
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  
  # Calculate appropriate breaks
  n_breaks <- 5
  breaks <- pretty(y_range, n = n_breaks)
  
  # Create labels
  times <- log(2) / abs(breaks)
  signs <- ifelse(breaks >= 0, "", "-")
  labels <- ifelse(abs(breaks) < 1e-6, Inf,
                   paste0(signs, round(times, 1)))
  
  # Add secondary axis with calculated breaks and labels
  p + scale_y_continuous(
    sec.axis = sec_axis(~ .,
                        breaks = breaks,
                        labels = labels,
                        name = "Doubling(+) / Halving(-) time")
  )
  
}
add_doubling_axis(p2) 
add_doubling_axis(p1)
