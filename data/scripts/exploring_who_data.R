library(tidyverse)
library(ggplot2)
source("~/Documents/GitHub/cop_framework/R/psi_theme.R")
save_plots <- function (p, wd, save_name, width, height) 
{
  message("Saving to ", paste0(wd, save_name, ".pdf/png"))
  
  showtext_auto(enable = TRUE)
  ggsave(p, filename = paste0(wd, save_name, ".pdf"), width = width, 
         height = height)
  
  showtext_auto(enable = FALSE)
  ggsave(p, filename = paste0(wd, save_name, ".png"), width = width, 
         height = height, units = "in", dpi = 300)
  showtext_auto()
  message("Success")
  NULL
}


dat <- read_csv("data/FlunetData_United Kingdom, England_All Sites_for_01 August 2012 to 20 October 2025 (1).csv") %>% filter(`Surveillance site type` == "Non-sentinel")
  # %>% filter(`Week start date (ISO 8601 calendar)` >"2023-05-01")
head(dat)

## Calculate growth rates and smoothed growth rates
dat <- dat %>% mutate(date=`Week start date (ISO 8601 calendar)`,y=`Influenza positive`,gr=log(y/lag(y,1))) %>% mutate(smooth_gr = zoo::rollmean(gr, k=5, fill=NA, align="center"))


p_inc <- ggplot(dat) + geom_ribbon(aes(x=date,ymax=y,ymin=0),fill=psi_core["oxford_blue"]) + 
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,15000,by=1000)) +
  ylab("Number of specimens") + xlab("Week start date") +
  scale_x_date(breaks="year") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 
p_inc

p_gr <- ggplot(dat) + 
  geom_hline(yintercept=0,linetype="dashed") +
  geom_path(aes(x=date,y=gr,group=year,col="Raw data"),linewidth=0.5)+
  geom_path(aes(x=date,y=smooth_gr,group=year,col="Smoothed (5 weeks)"),linewidth=0.75)+ 
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  scale_color_manual("", values=c("Raw data"="grey70","Smoothed (5 weeks)"="#4A4AB2"))+
  scale_y_continuous(expand=c(0,0),breaks=seq(-2,2,by=0.25)) +
  coord_cartesian(ylim=c(-2,2))+
  ylab("Log growth rate of specimens (per week)") + xlab("Week start date") +
  scale_x_date(breaks="year")+
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom") 

save_plots(p_inc,"figures/","incidence",width=8,height=4)
save_plots(p_gr,"figures/","gr",width=8,height=4)
## Convert to flu seasons
## Get influenza season from September to April
dat <- dat %>% mutate(year = ifelse(month(date) >= 8, year(date), year(date)-1))
dat <- dat %>% filter(month(date) >= 8 | month(date) <=4)
dat <- dat %>% group_by(year) %>% mutate(t = 1:n()) %>% mutate(t = t * 7) %>% mutate(t_date = as.Date("2025-08-01") + t)
dat$year_comb <- paste0(dat$year,"/",dat$year+1)

## Convert date to day of the year
dat <- dat %>% mutate(date_rel = yday(date),year_rel=year(date))
dat <- dat %>% mutate(
  year_label = case_when(
    year == 2025 ~ "2025",
    year == 2024 ~ "Post-pandemic (2022-24)",
    year == 2023 ~ "Post-pandemic (2022-24)",
    year == 2022 ~ "Post-pandemic (2022-24)",
    year == 2021 ~ "Pandemic (2020-21)",
    year == 2020 ~ "Pandemic (2020-21)",
    .default = "Pre-pandemic"
  )
)
p_growth_rates_aligned <- ggplot(dat %>% mutate(is_2025=year == 2025)) + 
  geom_hline(yintercept=0,linetype="dashed") +
  #geom_line(aes(x=date_rel,y=gr),col='red') +
  geom_line(aes(x=t_date,y=smooth_gr,col=as.factor(year_label),
                linewidth=is_2025,
                group=interaction(year_label,year)))  +
  scale_color_manual("Time period", values=c("2025"="#EB5F17","Post-pandemic (2022-24)"="#A32CA5",
                              "Pandemic (2020-21)"="#3382CD","Pre-pandemic"="grey")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.75)) +
  coord_cartesian(ylim=c(-1,1)) +
  scale_x_date(breaks="month") +
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  ylab("Smoothed log growth rate of specimens\n (per week)") + xlab("Week start date") +
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom")

save_plots(p_growth_rates_aligned,"figures/","p_growth_rates_aligned",width=8,height=4)


p_inc_alt <- ggplot(dat) + 
  #geom_line(aes(x=date_rel,y=gr),col='red') +
  geom_line(aes(x=t_date,y=y,col=as.factor(year_label),group=interaction(year_label,year)))  +
  scale_color_manual(values=c("2025"="red","Post-pandemic (2022-24)"="blue","Pandemic (2020-21)"="orange","Pre-pandemic"="grey")) +
  scale_color_manual("Time period", values=c("2025"="#EB5F17","Post-pandemic (2022-24)"="#A32CA5",
                                             "Pandemic (2020-21)"="#3382CD","Pre-pandemic"="grey")) +
  #scale_linewidth_manual(values=c("2025"=1.2,"Post-pandemic (2022-24)"=0.8,"Pandemic (2020-21)"=0.5,"Pre-pandemic"=0.5)) +
  scale_x_date(breaks="month") +
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  ylab("Number of specimens") + xlab("Week start date") +
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom")

save_plots(p_inc_alt,"figures/","p_inc_alt",width=8,height=6)


## Is the growth rate unusually high for the incidence?
p_gr_inc_compare <- dat %>% ggplot() + geom_point(aes(x=smooth_gr,y=log10(y),col=as.factor(year_label))) +
  scale_color_manual("Time period", values=c("2025"="red","Post-pandemic (2022-24)"="blue","Pandemic (2020-21)"="orange","Pre-pandemic"="grey")) +
  theme_psi() + xlab("Smoothed growth rate") + ylab("log10(Number of specimens)") + theme(legend.position="bottom")
save_plots(p_gr_inc_compare,"figures/","p_gr_inc_compare",width=6,height=4)


install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev', 'https://cloud.r-project.org'))
library("EpiEstim")

mean_si <- 2.5
std_si <- 1.5
method <- "parametric_si"
config <- make_config(list(mean_si = mean_si,
                           std_si = std_si))
R_est <- NULL
index <- 1
for(tmp_year in unique(dat$year)){

  inc <- dat %>% filter(year == tmp_year) %>% pull(y)

  output <- EpiEstim::estimate_R(incid = inc,
                                 dt = 7L,
                                 recon_opt = "match",
                                 method = method,
                                 iter=20L,
                                 config = config)

  R_est[[index]] <- output$R %>%  mutate(Year=tmp_year)
  index <- index + 1
}
R_est_comb <- do.call('bind_rows',R_est) %>% rename(year=Year)

R_est_comb <- R_est_comb %>% mutate(
  year_label = case_when(
    year == 2025 ~ "2025",
    year == 2024 ~ "Post-pandemic (2022-24)",
    year == 2023 ~ "Post-pandemic (2022-24)",
    year == 2022 ~ "Post-pandemic (2022-24)",
    year == 2021 ~ "Pandemic (2020-21)",
    year == 2020 ~ "Pandemic (2020-21)",
    .default = "Pre-pandemic"
  )
)

p_R_epiestim <- ggplot(R_est_comb %>% mutate(is_2025=year == 2025)) + 
  geom_hline(yintercept=1,linetype="dashed") +
  geom_hline(yintercept=1.2,linetype="dashed") +
  geom_hline(yintercept=1.4,linetype="dashed") +
  geom_line(aes(x=t_start,y=`Mean(R)`,col=year_label,group=interaction(year_label,year),linewidth=is_2025)) +
  scale_color_manual("Time period", values=c("2025"="#EB5F17","Post-pandemic (2022-24)"="#A32CA5",
                                             "Pandemic (2020-21)"="#3382CD","Pre-pandemic"="grey")) +
  scale_linewidth_manual(values=c("TRUE"=1.5,"FALSE"=0.75)) +
  #scale_x_date(breaks="month") +
  theme_psi() +
  theme(panel.grid.major.x = element_line(colour='grey90')) +
  ylab("Effective reproductive number") + xlab("Day (starting from 1st August)") +
  theme(legend.position="bottom")
save_plots(p_R_epiestim,"figures/","p_R_epiestim",width=8,height=4)

                          