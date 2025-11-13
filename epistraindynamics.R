library(EpiStrainDynamics)
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)

setwd("~/Documents/GitHub/EpiStrainDynamics/")

devtools::load_all("~/Documents/GitHub/EpiStrainDynamics/")


## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

gr_all <- NULL
index <- 1
for(agegroup in unique(all_dat_comb$age_new)){
all_dat_tmp <- all_dat_comb %>% filter(age_new == agegroup)
mod <- construct_model(
  method = p_spline(spline_degree = 3, days_per_knot = 3),
  pathogen_structure =single(
    case_timeseries = round((all_dat_tmp$ILI_plus)*10000),           # timeseries of case data
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
  gr_all[[index]] <- plot(gr) + scale_x_continuous(limits=c(0,35),breaks=seq(0,35,by=5)) + xlab("Week") + ylab('Growth rate (per week)')+
    ggtitle(agegroup)
  index <- index + 1
}

gr_all[[1]] + scale_y_continuous(limits=c(-1.2,1.2)) +
  gr_all[[2]] + scale_y_continuous(limits=c(-1.2,1.2))+
  gr_all[[3]] + scale_y_continuous(limits=c(-1.2,1.2))+
  gr_all[[4]] + scale_y_continuous(limits=c(-1.2,1.2)) +
  plot_layout(ncol=2)

tmp_dat_all<-NULL
for(i in 1:4){
  tmp_dat <- gr_all[[i]]$data
  tmp_dat$agegroup <- unique(all_dat_comb$age_new)[i]
  tmp_dat_all[[i]] <- tmp_dat
}
tmp_dat_all <- do.call("bind_rows",tmp_dat_all)
ggplot(tmp_dat_all) + geom_line(aes(x=time,y=y,col=agegroup))
ggplot(tmp_dat_all) + geom_ribbon(aes(x=time,ymin=lb_95,ymax=ub_95,fill=agegroup),alpha=0.3) +
  geom_line(aes(x=time,y=y,col=agegroup)) +
  xlab("Week") + ylab('Growth rate (per week)') +
  theme_minimal()
