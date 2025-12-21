# run_model_grid.R — minimal driver to run model over scenarios
library(tidyverse)
library(socialmixr)
library(data.table)
library(pracma)
library(deSolve)
library(patchwork)
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/ModelFluUk-H3N2//")
# ---- SOURCE helper files (must exist in same project) ----
source("auxiliary_funcs.R")   # contains setup_period_tibble, setup_holiday_tibble, etc.
source("sir_functions.R")     # contains get_beta, setup_C, epi_ode_size
source("helper_funcs.R")
source("vars.R")
# ---- CONSTANTS (copied from app) ----
half_term_start <- as.Date("2022-10-21")
half_term_end   <- as.Date("2022-10-30")
shopping_period_start <- as.Date("2022-12-01")
shopping_period_end   <- as.Date("2022-12-21")
winter_holiday_start  <- as.Date("2022-12-21")
winter_holiday_end    <- as.Date("2023-01-04")
N_tot <- 60000000
age_breaks <- c(0,5,18,65)

data("polymod")
contacts_all <- polymod$contacts
polymod_base <- polymod


# ---- core runner: run_model_once (takes a named list of parameters) ----
run_model_once <- function(params) {
  # params expected keys:
  # sim_start_date, sim_end_date, seed_size, R0, gamma, reporting_rate,
  # overall_immune_escape, prop_immune_younger, prop_immune_younger2, prop_immune_older, prop_immune_oldest,
  # prop_home_contacts_in_hols, prop_work_contacts_in_hols, prop_rest_contacts_in_hols,
  # prop_all_contacts_christmas, prop_rest_contacts_in_christmas, prop_home_contacts_in_christmas,
  # symp_fracs = c(symp_1,..)
  #
  # keep defaults for missing entries
  defaults <- list(
    sim_start_date = as.Date("2022-09-19"),
    sim_end_date = as.Date("2023-03-01"),
    seed_size = 1000,
    R0 = 1.9, gamma = 4.5,
    reporting_rate = 0.008,
    overall_immune_escape = 1,
    prop_immune_younger = 0.35, prop_immune_younger2 = 0.6,
    prop_immune_older = 0.7, prop_immune_oldest = 0.7,
    prop_home_contacts_in_hols = 1, prop_work_contacts_in_hols = 0.75, prop_rest_contacts_in_hols = 1.1,
    prop_all_contacts_christmas = 1.25, prop_rest_contacts_in_christmas = 0.7, prop_home_contacts_in_christmas = 3,
    symp_fracs = c(0.05, 0.2, 0.3, 0.5)
  )
  p <- modifyList(defaults, params)
  
  # dates & time grid (daily)
  start_date <- as.Date(p$sim_start_date)
  end_date <- as.Date(p$sim_end_date)
  dates <- seq(start_date, end_date, by = "1 week")
  start_day <- as.numeric(start_date)
  end_day_num <- as.numeric(end_date)
  ts <- seq(start_day, end_day_num, by = 1) - start_day + 1
  
  # contact matrices
  contact_params <- list(
    prop_home_contacts_in_hols = p$prop_home_contacts_in_hols,
    prop_work_contacts_in_hols = p$prop_work_contacts_in_hols,
    prop_rest_contacts_in_hols = p$prop_rest_contacts_in_hols,
    prop_all_contacts_christmas = p$prop_all_contacts_christmas,
    prop_rest_contacts_in_christmas = p$prop_rest_contacts_in_christmas,
    prop_home_contacts_in_christmas = p$prop_home_contacts_in_christmas
  )
  cm <- build_contact_matrices(contact_params)
  C_term <- cm$term
  C_holidays <- cm$term_break
  C_christmas_break <- cm$christmas_period
  C_xmas_holidays <- cm$christmas_holiday
  polymod_c_term <- cm$polymod_term
  
  # weights
  weights_res <- make_weights_df(start_date, end_date, smooth_time = 7)
  weights_df <- weights_res$weights
  
  # immunity vector (apply overall multiplier)
  raw_prop_immune <- c(p$prop_immune_younger, p$prop_immune_younger2, p$prop_immune_older, p$prop_immune_oldest)
  prop_immune <- pmin(1, p$overall_immune_escape * raw_prop_immune)
  
  # population N (two immunity compartments)
  N_props <- polymod_c_term$demography$proportion
  N_age_classes <- length(N_props)
  N_props_long <- c(N_props * (1 - prop_immune), N_props * prop_immune)
  alphas <- c(1, 0) # keep immune class susceptibility code consistent with app (alpha2 input not used here)
  alphas <- alphas/mean(alphas)
  N <- matrix(N_props_long * N_tot, ncol = length(alphas), nrow = N_age_classes)
  
  # beta & C matrices adjusted for N
  beta_par <- get_beta(C_term, polymod_c_term$participants$proportion, p$gamma, p$R0)
  beta_scales <- rep(1, N_age_classes)
  C_term_use <- setup_C(C_term, N, beta_scales)
  C_holidays_use <- setup_C(C_holidays, N, beta_scales)
  C_christmas_break_use <- setup_C(C_christmas_break, N, beta_scales)
  C_xmas_holidays_use <- setup_C(C_xmas_holidays, N, beta_scales)
  
  C_list <- vector("list", nrow(weights_df))
  for (i in seq_len(nrow(weights_df))) {
    w_term <- weights_df$w_term[i]
    h1 <- weights_df$h_termBreak[i]
    h2 <- weights_df$h_ChristmasPeriod[i]
    h3 <- weights_df$h_ChristmasBreak[i]
    C_list[[i]] <- C_term_use * w_term +
      C_holidays_use * h1 +
      C_christmas_break_use * h2 +
      C_xmas_holidays_use * h3
  }
  
  # run model
  y_base <- epi_ode_size(C_list, beta_par, p$gamma, N, ts = ts,
                         alphas = alphas,
                         age_seed = 3, immunity_seed = 1,
                         seed_size = p$seed_size, return_compartments = TRUE)
  
  # extract incidence columns and convert to weekly reported (same as app)
  use_cols <- which(colnames(y_base) %like% "inc")
  ret <- as.matrix(y_base[, use_cols])
  ret_daily <- apply(ret, 2, function(x) c(0, diff(x))) # daily incidence
  # combine immunity columns like app: sum pairs to get per-age daily incidence
  ret_age_daily <- ret_daily[, seq(1, ncol(ret_daily), 2)] + ret_daily[, seq(2, ncol(ret_daily), 2)]
  # apply symptomatic fractions and reporting
  symp_frac <- p$symp_fracs
  ret_age_daily_reported <- t(p$reporting_rate * symp_frac * t(ret_age_daily))
  
  # aggregate to weekly (7-day) sums
  row_groups <- gl(nrow(ret_age_daily_reported) %/% 7 + (nrow(ret_age_daily_reported) %% 7 > 0), 7, nrow(ret_age_daily_reported))
  ret_weekly_reported <- aggregate(ret_age_daily_reported, by = list(row_groups), FUN = sum)[, -1]  # matrix weeks x ages
  
  date_key <- data.frame(t = 1:nrow(ret_weekly_reported), date = dates)
  
  # compute weekly totals (reported)
  weekly_totals_reported <- rowSums(as.matrix(ret_weekly_reported), na.rm = TRUE)
  peak_idx <- if (length(weekly_totals_reported)) which.max(weekly_totals_reported) else NA
  peak_date <- if (!is.na(peak_idx) && length(peak_idx)) date_key$date[peak_idx] else NA
  reported_peak <- if (!is.na(peak_idx)) weekly_totals_reported[peak_idx] else NA
  # estimate underlying symptomatic at peak (divide by reporting rate)
  peak_symptomatic_est <- if (!is.na(reported_peak) && p$reporting_rate > 0) reported_peak / p$reporting_rate else NA
  # cumulative incidence in model state: sum 'inc' compartment at final time
  y_base_long <- y_base %>% pivot_longer(-time)
  y_base_long <- y_base_long %>%
    mutate(compartment = stringr::str_split(name, "_", simplify = TRUE)[,1],
           age = as.integer(stringr::str_split(name, "_", simplify = TRUE)[,2]),
           immunity = as.integer(stringr::str_split(name, "_", simplify = TRUE)[,3]))
  
  ## Calculate overall weekly growth rate
  weekly_growth_rate_total <- y_base_long %>% 
    filter(compartment == "inc") %>% 
    mutate(time=floor(time/7)) %>% 
    group_by(time) %>% summarize(inc=sum(value)) %>% 
    mutate(gr=log(inc/lag(inc,1))) 
  
  ## Get peak after the first 4 weeks
  peak_gr <- weekly_growth_rate_total %>% filter(time >= 5) %>% filter(gr == max(gr))
  
  ## By age
  weekly_growth_rate_age <- y_base_long %>% 
    filter(compartment == "inc") %>%
    mutate(time=floor(time/7)) %>% 
    group_by(time,age) %>% 
    summarize(inc=sum(value)) %>%
    group_by(age) %>% 
    mutate(gr=log(inc/lag(inc,1))) 
  
  peak_gr_age <- weekly_growth_rate_age %>% filter(time >= 5) %>% group_by(age) %>% filter(gr == max(gr))
  
  
  # total infections by age at final time
  cumulative_incidence <- y_base_long %>% filter(time == max(time)) %>%
    filter(compartment == "inc") %>%
    group_by(age) %>%
    summarize(total_inf = sum(value), .groups = "drop")
  
  # match age group labels
  age_groups_all <- c("[0,5)", "[5,18)", "[18,65)", "65+")
  cumulative_incidence <- cumulative_incidence %>% mutate(age_group = age_groups_all[age])
  
  total_infections_final <- sum(cumulative_incidence$total_inf)/sum(N)
  total_infections_65plus <- cumulative_incidence %>% filter(age_group == "65+") %>% pull(total_inf)
  total_infections_65plus <- ifelse(length(total_infections_65plus)==1, total_infections_65plus, NA)
  total_infections_65plus <- total_infections_65plus/sum(N[4,])
  # also compute total symptomatic final (apply symp_frac)
  symptomatic_by_age <- cumulative_incidence %>% mutate(symp = total_inf * symp_frac)
  total_symptomatic_final <- sum(symptomatic_by_age$symp, na.rm = TRUE)
  
  # prepare outputs
  out <- list(
    res = list(inc = as.data.frame(ret_weekly_reported), date_key = date_key, cumulative_incidence = cumulative_incidence, weekly_growth_rate_total=weekly_growth_rate_total, peak_gr=peak_gr),
    metrics = tibble(
      sim_start_date = start_date,
      sim_end_date = end_date,
      seed_size = p$seed_size,
      R0 = p$R0,
      gamma = p$gamma,
      overall_immune_escape = p$overall_immune_escape,
      prop_immune_younger = p$prop_immune_younger,
      prop_immune_younger2 = p$prop_immune_younger2,
      prop_immune_older = p$prop_immune_older,
      prop_immune_oldest = p$prop_immune_oldest,
      final_size_infections = total_infections_final,
      final_size_symptomatic = total_symptomatic_final,
      date_of_peak = peak_date,
      reported_peak_weekly = reported_peak,
      estimated_peak_symptomatic = peak_symptomatic_est,
      final_size_65plus = total_infections_65plus
    )
  )
  out
}

# ---- Example parameter grid (customise as needed) ----
# ---- Run model over different scenarios ----

## Loop over R0 values
R0_vals <- seq(1.05,3,by=0.05)
results_list <- vector("list", length(R0_vals))
for (i in seq_along(R0_vals)) {
  params <- list(R0=R0_vals[i])
  
  cat(sprintf("[%d/%d] running: R0=%.2f\n",
              i, length(R0_vals),params$R0))
  
  res <- run_model_once(params)
  tmp_gr <- res$res$weekly_growth_rate_total %>% rename(t=time)
  ## Pull out peak growth rate after 1% of cases have elapsed
  tmp_gr <- tmp_gr %>% mutate(cumu_inc=cumsum(inc)/sum(inc)) %>% filter(cumu_inc > 0.01) %>% filter(gr == max(gr)) %>%
    left_join(res$res$date_key) %>%
    select(-c(t,inc)) %>%
    rename(gr_early_nov=gr,nov_date=date)
  
  results_list[[i]] <- bind_cols(res$metrics,tmp_gr)
  ## Pull out growth rates on 1st Oct, 1st Nov, and 1st Dec
}
results_tbl_R0 <- bind_rows(results_list)

## Loop over Seed values
seed_vals <- seq(as.Date("2022-08-01"),as.Date("2022-11-01"),by="7 days")
results_list <- vector("list", length(seed_vals))
for (i in seq_along(seed_vals)) {
  params <- list(sim_start_date=seed_vals[i])
  
  cat(sprintf("[%d/%d] running: Seed time=%.2f\n",
              i, length(seed_vals),params$sim_start_date))
  
  res <- run_model_once(params)
  tmp_gr <- res$res$weekly_growth_rate_total %>% rename(t=time)
  tmp_gr <- tmp_gr %>% mutate(cumu_inc=cumsum(inc)/sum(inc)) %>% filter(cumu_inc > 0.01) %>% filter(gr == max(gr)) %>%
    left_join(res$res$date_key) %>%
    select(-c(t,inc)) %>%
    rename(gr_early_nov=gr,nov_date=date)
  
  results_list[[i]] <- bind_cols(res$metrics,tmp_gr)
}
results_tbl_seed <- bind_rows(results_list)

## Loop over immune escape vals
overall_escape_vals <- seq(0,1,by=0.05)
results_list <- vector("list", length(overall_escape_vals))
for (i in seq_along(overall_escape_vals)) {
  params <- list(overall_immune_escape=overall_escape_vals[i])
  
  cat(sprintf("[%d/%d] running: overall_immune_escape=%.2f\n",
              i, length(overall_escape_vals),params$overall_immune_escape))
  
  res <- run_model_once(params)
  tmp_gr <- res$res$weekly_growth_rate_total %>% rename(t=time)
  tmp_gr <- tmp_gr %>% mutate(cumu_inc=cumsum(inc)/sum(inc)) %>% filter(cumu_inc > 0.01) %>% filter(gr == max(gr)) %>%
    left_join(res$res$date_key) %>%
    select(-c(t,inc)) %>%
    rename(gr_early_nov=gr,nov_date=date)
  results_list[[i]] <- bind_cols(res$metrics,tmp_gr)
}
results_tbl_escape <- bind_rows(results_list)

## Loop over immune escape in younger ages
overall_escape_vals <- seq(0,1,by=0.05)
results_list <- vector("list", length(overall_escape_vals))
for (i in seq_along(overall_escape_vals)) {
  prop_immune_younger <- 0.3
  prop_immune_younger2 <- 0.6
  params <- list(prop_immune_younger=overall_escape_vals[i]*prop_immune_younger,
                 prop_immune_younger2=overall_escape_vals[i]*prop_immune_younger2)
  
  cat(sprintf("[%d/%d] running: in children overall_immune_escape=%.2f\n",
              i, length(overall_escape_vals),overall_escape_vals[i]))
  
  res <- run_model_once(params)
  tmp_gr <- res$res$weekly_growth_rate_total %>% rename(t=time)
  tmp_gr <- tmp_gr %>% mutate(cumu_inc=cumsum(inc)/sum(inc)) %>% filter(cumu_inc > 0.01) %>% filter(gr == max(gr)) %>%
    left_join(res$res$date_key) %>%
    select(-c(t,inc)) %>%
    rename(gr_early_nov=gr,nov_date=date)
  results_list[[i]] <- bind_cols(res$metrics,tmp_gr)
  results_list[[i]]$index <- overall_escape_vals[i]
}
results_tbl_young_escape <- bind_rows(results_list)

## Create data for varying R0
results_tbl_R0_long <- results_tbl_R0 %>% mutate(Scenario = "R0") %>% select(R0, estimated_peak_symptomatic, final_size_infections, final_size_65plus,gr_early_nov,
 Scenario) %>% rename(variable=R0) %>% 
  mutate(estimated_peak_symptomatic = 100000*estimated_peak_symptomatic/60000000) %>%
  rename(`Final size`=final_size_infections,
         `Final size 65+`=final_size_65plus,
         `Peak infection incidence \nper 100,000`=estimated_peak_symptomatic,
         `Peak growth rate`=gr_early_nov) %>%
  pivot_longer(-c(variable,Scenario))

## Create data for varying Seed time
results_tbl_seed_long <- results_tbl_seed %>% mutate(Scenario = "Seed time") %>% select(sim_start_date, estimated_peak_symptomatic, final_size_infections, final_size_65plus, gr_early_nov, Scenario) %>% rename(variable=sim_start_date) %>% 
  mutate(estimated_peak_symptomatic = 100000*estimated_peak_symptomatic/60000000) %>%
  rename(`Final size`=final_size_infections,
         `Final size 65+`=final_size_65plus,
         `Peak infection incidence \nper 100,000`=estimated_peak_symptomatic,
         `Peak growth rate`=gr_early_nov) %>%
  pivot_longer(-c(variable,Scenario))


## Create data for varying Overall immune escape
results_tbl_escape_long <- results_tbl_escape %>% mutate(Scenario = "Immune escape") %>% select(overall_immune_escape, estimated_peak_symptomatic, final_size_infections, final_size_65plus, gr_early_nov,  Scenario) %>% rename(variable=overall_immune_escape) %>% 
  mutate(estimated_peak_symptomatic = 100000*estimated_peak_symptomatic/60000000) %>%
  rename(`Final size`=final_size_infections,
         `Final size 65+`=final_size_65plus,
         `Peak infection incidence \nper 100,000`=estimated_peak_symptomatic,
         `Peak growth rate`=gr_early_nov) %>%
  pivot_longer(-c(variable,Scenario))

## Create data for varying Immune escape in younger ages
results_tbl_young_escape_long <- results_tbl_young_escape %>% mutate(Scenario = "Immune escape in children") %>% select(index, estimated_peak_symptomatic, final_size_infections, final_size_65plus,gr_early_nov,  Scenario) %>% rename(variable=index) %>% 
  mutate(estimated_peak_symptomatic = 100000*estimated_peak_symptomatic/60000000) %>%
  rename(`Final size`=final_size_infections,
         `Final size 65+`=final_size_65plus,
         `Peak infection incidence \nper 100,000`=estimated_peak_symptomatic,
         `Peak growth rate`=gr_early_nov) %>%
  pivot_longer(-c(variable,Scenario))


## Read in estimated growth rates from RCGP and WHO FluNet
gr_flunet <- read_csv("../results/flunet_h3_growth_rates.csv") %>% filter(season=='2025/26') %>%
  mutate(month = month(date)) %>% filter(month %in% c(11)) %>% group_by(month) %>% filter(date == min(date)) %>%
  select(y,lb_95,ub_95) %>% mutate(Dataset="WHO FluNet") %>%
  
  mutate(Label="WHO FluNet") %>%
  mutate(name = "Peak growth rate")

rcgp_dat <- read_csv("../results/RCGP_growth_rates_by_age.csv") %>% filter(agegroup=="All",Season=='2025 to 2026') %>%
  mutate(month = month(date)) %>% filter(month %in% c(11)) %>% group_by(month) %>% filter(date == min(date)) %>%
  select(y,lb_95,ub_95) %>% mutate(Dataset="RCGP") %>% 
  mutate(Label="RCGP") %>%
  mutate(name = "Peak growth rate")


gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)
gr_est_dat$x <- c(3,3.5)


ylims <- tibble(
  name = c(
    "Final size",
    "Final size 65+",
    "Peak growth rate",
    "Peak infection incidence \nper 100,000"
  ),
  ymin = c(0, 0, 0, 0),
  ymax = c(0.35, 0.35, 1.5, 5000)
)
p_r0 <- ggplot(results_tbl_R0_long) + 
  geom_blank(
    data = ylims,
    aes(y = ymin)
  ) +
  geom_blank(
    data = ylims,
    aes(y = ymax)
  ) +
  geom_vline(xintercept=1.9, linetype="dashed", color="grey") +
  geom_line(aes(x=variable, y=value, color=name),linewidth=0.75) +
  #geom_hline(data=gr_est_dat,aes(yintercept=y,linetype=Dataset)) +
  #geom_pointrange(data=gr_est_dat,aes(y=y,ymin=lb_95,ymax=ub_95,group=Dataset,x=x,col=Label),linewidth=0.5,size=0.25) +
  #geom_text(data=gr_est_dat,aes(y=y + 0.04,label=Label,x=x2),size=3,hjust=0) +
  facet_wrap(~name,scales="free_y",nrow=1) +
  xlab("R0") +
  ylab("Value") +
  scale_color_brewer("Metric", palette="Set2") +
  #scale_linetype_manual(values=c("WHO FluNet"="dashed","RCGP"="dotted")) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position="bottom")
p_r0

gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)
gr_est_dat$x <- c(0.7,0.8)

#gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)%>% mutate(x1=if_else(Dataset=="WHO FluNet",0.7,0.8)) %>% mutate(x2 = if_else(Dataset=="WHO FluNet",.05,.05))

ylims <- tibble(
  name = c(
    "Final size",
    "Final size 65+",
    "Peak growth rate",
    "Peak infection incidence \nper 100,000"
  ),
  ymin = c(0, 0, 0, 0),
  ymax = c(1, 1, 1.5, 20000)
)

p_escape <- ggplot(results_tbl_escape_long) +
  geom_vline(xintercept=1, linetype="dashed", color="grey") +
  geom_line(aes(x=variable, y=value, color=name),linewidth=0.75) + 
  
  geom_blank(
    data = ylims,
    aes(y = ymin)
  ) +
  geom_blank(
    data = ylims,
    aes(y = ymax)
  ) +
  #geom_hline(data=gr_est_dat,aes(yintercept=y,linetype=Dataset)) +
  #geom_pointrange(data=gr_est_dat,aes(y=y,ymin=lb_95,ymax=ub_95,group=Dataset,x=x,col=Label),linewidth=0.5,size=0.25,alpha=1) +
  #geom_text(data=gr_est_dat,aes(y=y + 0.04,label=Dataset,x=x2),size=3,hjust=0) +
  
  facet_wrap(~name,scales="free_y",nrow=1) +
  xlab("Proportion immune relative\n to baseline") +
  scale_linetype_manual(values=c("WHO FluNet"="dashed","RCGP"="dotted")) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Value") +
  scale_color_brewer("Metric", palette="Set2") +
  theme_bw() +
  theme(legend.position="bottom")
p_escape


gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)
gr_est_dat$x <- c(0.7,0.8)
#gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)%>% mutate(x1=if_else(Dataset=="WHO FluNet",0.7,0.8)) %>% mutate(x2 = if_else(Dataset=="WHO FluNet",.05,.05))
ylims <- tibble(
  name = c(
    "Final size",
    "Final size 65+",
    "Peak growth rate",
    "Peak infection incidence \nper 100,000"
  ),
  ymin = c(0, 0, 0, 0),
  ymax = c(0.35, 0.35, 1.5, 5000)
)

p_escape_young <- ggplot(results_tbl_young_escape_long) +
  geom_vline(xintercept=1, linetype="dashed", color="grey") +
  
  geom_blank(
    data = ylims,
    aes(y = ymin)
  ) +
  geom_blank(
    data = ylims,
    aes(y = ymax)
  ) +
  geom_line(aes(x=variable, y=value, color=name),linewidth=0.75) +
  #geom_hline(data=gr_est_dat,aes(yintercept=y,linetype=Dataset)) +
 # geom_pointrange(data=gr_est_dat,aes(y=y,ymin=lb_95,ymax=ub_95,group=Dataset,x=x,col=Label),linewidth=0.5,size=0.25,alpha=1) +
  
#  geom_text(data=gr_est_dat,aes(y=y + 0.04,label=Dataset,x=x2),size=3,hjust=0) +
  
  facet_wrap(~name,scales="free_y",nrow=1) +
  xlab("Proportion immune relative\n to baseline in <18") +
  ylab("Value") +
  scale_color_brewer("Metric", palette="Set2") +
  scale_linetype_manual(values=c("WHO FluNet"="dashed","RCGP"="dotted")) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position="bottom")
p_escape_young


gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)
gr_est_dat$x <- c(as.Date("2022-10-01"),as.Date("2022-10-15"))
#gr_est_dat <- bind_rows(gr_flunet, rcgp_dat)%>% 
#  mutate(x1=if_else(Dataset=="WHO FluNet",as.Date("2022-10-01"),as.Date("2022-10-15"))) %>%
#  mutate(x2 = if_else(Dataset=="WHO FluNet",as.Date("2022-08-01"),as.Date("2022-08-01")))

ylims <- tibble(
  name = c(
    "Final size",
    "Final size 65+",
    "Peak growth rate",
    "Peak infection incidence \nper 100,000"
  ),
  ymin = c(0, 0, 0, 0),
  ymax = c(0.2, 0.2, 1.5, 1000)
)

p_seed <- ggplot(results_tbl_seed_long) +
  geom_vline(xintercept=as.Date("2022-09-19"), linetype="dashed", color="grey") +
  geom_blank(
    data = ylims,
    aes(y = ymin)
  ) +
  geom_blank(
    data = ylims,
    aes(y = ymax)
  ) +
  geom_line(aes(x=variable, y=value, color=name),linewidth=0.75) +
  #geom_hline(data=gr_est_dat,aes(yintercept=y,linetype=Dataset)) +
  #geom_pointrange(data=gr_est_dat,aes(y=y,ymin=lb_95,ymax=ub_95,group=Dataset,x=x,col=Label),linewidth=0.5,size=0.25,alpha=1) +
  #geom_text(data=gr_est_dat,aes(y=y + 0.04,label=Dataset,x=x2),size=3,hjust=0) +
  
  facet_wrap(~name,scales="free_y",nrow=1) +
  xlab("Seed date (1000 infections)") +
  ylab("Value") +
  scale_color_brewer("Metric", palette="Set2") +
  scale_linetype_manual(values=c("WHO FluNet"="dashed","RCGP"="dotted")) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position="bottom")
p_seed
p_main <- p_r0/
  p_escape/
  p_escape_young/
  p_seed + plot_layout(guides="collect")&
  theme(legend.position='none',axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        strip.text = element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        title = element_text(size=16))

ggsave("~/Documents/GitHub/influenza_H3N2_k_clade/new_figures/fig_model_grid_results.png", p_main, width=10, height=10)
ggsave("~/Documents/GitHub/influenza_H3N2_k_clade/new_figures/fig_model_grid_results.pdf", p_main, width=10, height=10)
