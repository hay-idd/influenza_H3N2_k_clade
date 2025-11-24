# Variables and data to load

half_term_start <- as.Date("2022-10-21")
half_term_end   <- as.Date("2022-10-30")
shopping_period_start <- as.Date("2022-12-01")
shopping_period_end   <- as.Date("2022-12-15")
winter_holiday_start  <- as.Date("2022-12-21")
winter_holiday_end    <- as.Date("2023-01-04")
age_breaks <- c(0,5,18,65)

data("polymod")

defaults <- list(
  
  ## Core parameter (main tab)
  R0 = 2,
  gamma = 4,
  overall_immune_escape = 1,
  seed_size = 1000,
  sim_start_date = as.Date("2022-09-10"),
  reporting_rate = 0.007,
  
  ## Immunity and disease tab
  prop_immune_younger = 0.3,
  prop_immune_younger2 = 0.6,
  prop_immune_older = 0.7,
  prop_immune_oldest = 0.75,
  
  symp_1 = 0.05,
  symp_2 = 0.2,
  symp_3 = 0.35,
  symp_4 = 0.45,
  
  ## Contact patterns tab
  N_tot = 60e6,
  prop_home_contacts_in_hols = 1,
  prop_work_contacts_in_hols = 0.75,
  prop_rest_contacts_in_hols = 1.1,
  prop_all_contacts_christmas = 1.3,
  prop_rest_contacts_in_christmas = 0.67,
  prop_home_contacts_in_christmas = 3,
  
  ## Plot options tab
  y_lim_max = 4000,
  plot_start_date = as.Date("2022-10-01"),
  sim_end_date = as.Date("2023-03-01"),
  
  ## Growth plot scale
  growth_y_scale = 1.5,
  
  ## Remove/hide
  initial_immune_frac = 0,
  alpha2 = 0
)