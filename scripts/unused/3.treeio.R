library(ape)
library(phangorn)
library(ips)
library(treeio)
library(tidyverse)
tr <- treeio::read.beast("~/Downloads/nextstrain_seasonal-flu_h3n2_ha_2y_timetree (1).nexus")
tr1 <- as_tibble(tr)
unique(tr1$subclade)

tr1$date <- as.numeric(tr1$num_date)

## split epiweek variable into year (first 4 digits) and epi week (last 2 digits)
tr1$epiweek <- as.numeric(tr1$epiweek)
tr1$epiweek_year <- floor(tr1$epiweek/100)


tr1 <- tr1 %>% mutate(year = as.integer(substr(epiweek, 1, 4)), week = as.integer(substr(epiweek, 5, 6)))

tr1$raw_date <- tr1$year + (tr1$week - 1)/52.1775

tr1_subclade <- tr1 %>% 
  filter(year >= 2025) %>%
  group_by(subclade, week) %>%
  tally()

expand_grid(subclade = unique(tr1_subclade$subclade), week = 1:52) %>%
  left_join(tr1_subclade) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(week) %>% mutate(n_tot=sum(n)) %>%
  mutate(prop = n/n_tot) %>%
  ggplot() + geom_line(aes(x=week,y=prop,col=subclade))

expand_grid(subclade = unique(tr1_subclade$subclade), week = 1:52) %>%
  left_join(tr1_subclade) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(week) %>% mutate(n_tot=sum(n)) %>%
  mutate(prop = n/n_tot) %>%
  ggplot() + geom_bar(aes(x=week,y=prop,fill=subclade),position="fill",stat="identity")

tr1 %>% 
  filter(year >= 2025, week > 20) %>%
  group_by(subclade, week) %>%
  tally() %>%  
  group_by(week) %>% mutate(n_tot=sum(n)) %>%
  mutate(prop = n/n_tot) %>%
  ggplot() + geom_line(aes(x=week,y=prop,col=subclade))

h3_dat <- tr1 %>% 
  filter(year >= 2024, week > 20,region=="Europe") %>%
  group_by(subclade, week) %>%
  tally() %>%  
  group_by(week) %>% mutate(n_tot=sum(n)) %>%
  mutate(prop = n/n_tot)

