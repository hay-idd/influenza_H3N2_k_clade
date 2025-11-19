## Demography data
## Population of England
N_tot <- 60000000

# read input parameters
prop_work_contacts_in_hols <- input$prop_work_contacts_in_hols
prop_home_contacts_in_hols <- input$prop_home_contacts_in_hols

prop_immune_younger <- input$prop_immune_younger
prop_immune_younger2 <- input$prop_immune_younger2
prop_immune_older <- input$prop_immune_older
prop_immune_oldest <- input$prop_immune_oldest

alphas <- c(1, input$alpha2)
R0 <- input$R0
gamma <- input$gamma
seed_size <- input$seed_size
reporting_rate <- input$reporting_rate

R0 <- input$R0
symp_frac <- c(input$symp_1,input$symp_2,input$symp_3,input$symp_4)

start_day <- as.numeric(start_date)
end_day <- as.numeric(end_date)
date_seq <- seq(start_date,end_date,by="1 day")
ts <- seq(start_day,end_day,by=1) - start_day + 1

########################################################
## CREATE CONTACT MATRICES
########################################################
## POLYMOD contact matrix
data("polymod")
contacts <- polymod$contacts
## Create new contact survey with no school contacts, decreased work contacts, and increased home contacts
contacts_no_schools <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home==0)# & cnt_leisure == 0)
## Resample to give frac of work contacts
contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(prop_work_contacts_in_hols)
## Resample with replacement to give more home contacts
contacts_home <- contacts %>% filter(cnt_home == 1)
contacts_home <- contacts_home[sample(1:nrow(contacts_home), size = nrow(contacts_home)*prop_home_contacts_in_hols, replace=TRUE), ]
contacts_holidays <- bind_rows(contacts_no_schools,contacts_work,contacts_home)

## Create holiday contact matrix
age_breaks <- c(0,5,18,65)
polymod1 <- polymod
polymod1$contacts <- contacts_holidays
polymod_c_holidays <- contact_matrix(polymod1,
                                     countries="United Kingdom",
                                     #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                     age.limits=age_breaks,
                                     symmetric=TRUE,
                                     missing.contact.age = "sample",
                                     missing.participant.age = "remove")
C_holidays <- polymod_c_holidays$matrix
row.names(C_holidays) <- colnames(C_holidays)

## Create term time contact matrix
polymod_c_term <- contact_matrix(polymod,
                                 countries="United Kingdom",
                                 #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                 age.limits=age_breaks,
                                 symmetric=TRUE,
                                 missing.contact.age = "sample",
                                 missing.participant.age = "remove")
C_term <- polymod_c_term$matrix
row.names(C_term) <- colnames(C_term)


########################################################
## SETUP MODEL FOR 2023/24 SEASON
########################################################

## From https://www.kennedycollegeoxford.org.uk/term-dates
## Holiday dates
## Setup contact matrix list
school_days_weighted <- setup_holiday_tibble(start_date,end_date, half_term_start,half_term_end,
                                             winter_holiday_start,winter_holiday_end, smooth_time = 7)
prop_immune <- c(prop_immune_younger,prop_immune_younger2,prop_immune_older,prop_immune_oldest)
## Setup population proportions
N_props <- polymod_c_term$demography$proportion
N_age_classes <- length(N_props)
N_immunity_classes <- 2

## Number of people in each age group and age class
N <- matrix(N_props*N_tot,ncol=1,nrow=N_age_classes)
beta_scales <- rep(1, N_age_classes)
N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
beta_par <- get_beta(C_term,polymod_c_term$participants$proportion,gamma, R0)

beta_scales <- rep(1,N_age_classes)
C_use_holiday <- setup_C(C_holidays, N, beta_scales)
C_use_term <- setup_C(C_term, N, beta_scales)

C_list <- NULL
## Create list of interpolated contact matrices
for(i in 1:nrow(school_days_weighted)){
  C_list[[i]] <- C_use_term*school_days_weighted$weight[i] + C_use_holiday*(1 - school_days_weighted$weight[i])
}



## Get proporiton immune vector
prop_immune <- c(prop_immune_youngest,prop_immune_younger, prop_immune_older,prop_immune_oldest)

## Housekeeping for setting up contact matrix
N_age_classes <- length(N_props)
N_immunity_classes <- 2

N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
#N_props_long <- N_props
N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)

## Get beta parameter
beta_par <- get_beta(C_term,N_props,gamma, R0)
beta_scales <- rep(1,N_age_classes)
y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                       alphas=alphas, age_seed=3,immunity_seed=1,seed_size=seed_size,return_compartments=TRUE)
use_cols <- which(colnames(y_base) %like% "inc")


ret <- as.matrix(y_base[,use_cols])

## Get daily incidence
ret <- apply(ret, 2, function(x) c(0, diff(x)))

## Get expected symptomatic cases scaled by reporting rate
# sum every 2 columns for the two immune states
ret2 <- ret[, seq(1, ncol(ret), 2)] + ret[, seq(2, ncol(ret), 2)]
ret2 <- t(reporting_rate*symp_frac*t(ret2))

# then sum every 7 rows to get weekly data
row_groups <- gl(nrow(ret2) %/% 7 + (nrow(ret2) %% 7 > 0), 7, nrow(ret2))
ret_sum <- aggregate(ret2, by = list(row_groups), FUN = sum)[, -1]
ret_sum
