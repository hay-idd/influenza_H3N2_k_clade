# app.R
library(shiny)
library(socialmixr)
library(pracma)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(data.table)
library(tidyr)
library(DT)

#library(shinycssloaders) # optional, shows spinner while plotting
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/shiny_app/")
# UI ----
ui <- fluidPage(
  titlePanel("Influenza H3N2 model"),
  sidebarLayout(
    sidebarPanel(
      # A small note / run button
      #actionButton("run_model", "Run model"),
      downloadButton("download_plot", "Download Plot (PNG)"),
      downloadButton("download_table", "Download Table (CSV)"),
      hr(),
      
      h4("Plot settings"),
      numericInput("y_lim_max", "Maximum y-axis value", value = 4000, step = 100),
      
      hr(),
      
      h4("Model parameters"),
      # Other epidemiological parameters
      numericInput("alpha2", "Relative susceptibility susceptible class", value = 0.65, step = 0.05),
      numericInput("R0", "Basic reproductive number R0", value = 2, step = 0.01),
      numericInput("gamma", "Infectious period (days) gamma", value = 4.6, step = 0.1),
      numericInput("seed_size", "Seed size (number of initial infections)", value = 150, step = 10),
      # New: simulation start date slider (user requested a slider)
      sliderInput("sim_start_date", "Seed date",
                  min = as.Date("2022-07-01"),
                  max = as.Date("2022-12-31"),
                  value = as.Date("2022-09-01"),
                  timeFormat = "%Y-%m-%d"),
      sliderInput("plot_start_date", "Plot start date",
                  min = as.Date("2022-07-01"),
                  max = as.Date("2022-12-31"),
                  value = as.Date("2022-10-01"),
                  timeFormat = "%Y-%m-%d"),
      
      sliderInput("sim_end_date", "Simulation end date",
                  min = as.Date("2023-01-01"),
                  max = as.Date("2023-07-01"),
                  value = as.Date("2023-03-01"),
                  timeFormat = "%Y-%m-%d"),
      hr(),
      # Reporting / symptomatic fraction
      sliderInput("reporting_rate", "Reporting rate (fraction of symptomatic reported)",
                  min = 0, max = 0.01, value = 0.0012, step = 0.0001),
      hr(),
      h4("School holiday contacts"),
      # Contacts / holiday behavior
      ## Break sliders
      sliderInput("prop_home_contacts_in_hols", "Multiplier for home contacts in breaks",
                  min = 0, max = 3, value = 1, step = 0.01),
      sliderInput("prop_work_contacts_in_hols", "Proportion of work contacts kept in school holidays",
                  min = 0, max = 1, value = 0.75, step = 0.01),
      sliderInput("prop_rest_contacts_in_hols", "Multiplier for non-school or work contacts in breaks",
                  min = 0.5, max = 3, value = 1.1, step = 0.01),
      hr(),
      h4("Christmas period contacts (before school holiday)"),
            ## Christmas period sliders
      sliderInput("prop_all_contacts_christmas", "Multiplier for all non-school contacts in Christmas period",
                  min = 0.5, max = 3, value = 1.3, step = 0.01),
      hr(),
      ## Christmas holiday sliders
      h4("Christmas holiday contacts"),
      #sliderInput("prop_work_contacts_in_christmas", "Proportion of work contacts kept in Christmas holidays",
      #            min = 0, max = 1, value = 1, step = 0.01),
      sliderInput("prop_rest_contacts_in_christmas", "Proportion of all contacts kept over Christmas",
                  min = 0, max = 1, value = 0.67, step = 0.01),
      sliderInput("prop_home_contacts_in_christmas", "Multiplier for home contacts over Christmas holiday",
                  min = 1, max = 5, value = 3, step = 0.01),
      
      
      hr(),
      # Immunity proportions by age-groups
      sliderInput("prop_immune_younger", "Proportion immune (0-4 yrs)",
                  min = 0, max = 1, value = 0.3, step = 0.01),
      sliderInput("prop_immune_younger2", "Proportion immune (5-18 yrs)",
                  min = 0, max = 1, value = 0.4, step = 0.01),
      sliderInput("prop_immune_older", "Proportion immune (19-64 yrs)",
                  min = 0, max = 1, value = 0.8, step = 0.01),
      sliderInput("prop_immune_oldest", "Proportion immune (65+ yrs)",
                  min = 0, max = 1, value = 0.8, step = 0.01),
      hr(),
      # Symptomatic fraction by age-groups
      sliderInput("symp_1", "Symptomatic frac incorporating reporting rate (0-4 yrs)",
                  min = 0, max = 1, value = 0.25, step = 0.01),
      sliderInput("symp_2", "Symptomatic frac incorporating reporting rate (5-18 yrs)",
                  min = 0, max = 1, value = 0.4, step = 0.01),
      sliderInput("symp_3", "Symptomatic frac incorporating reporting rate (19-64 yrs)",
                  min = 0, max = 1, value = 0.5, step = 0.01),
      sliderInput("symp_4", "Symptomatic frac incorporating reporting rate (65+ yrs)",
                  min = 0, max = 1, value = 0.5, step = 0.01),
      width = 3
    ),
    mainPanel(
      # Top area can show a quick summary
    
      # Plot at bottom right in its own panel
      fluidRow(
        column(12,
               wellPanel(
                 h4("Reported symptomatic cases (weekly)"),
                 plotOutput("inc_plot", height = "600px"),
                 hr(),
                 h4("Cumulative cases"),
                 # scrollable container so big tables don't blow up layout
                 div(style = "max-height:300px; overflow-y:auto; overflow-x:auto;",
                     DT::dataTableOutput("cum_table")
                 ),
                 h4("Contact matrices (term / term-break / Christmas period / Christmas break)"),
                 plotOutput("contact_matrices_plot", height = "900px")
               ))
      ),
      width = 9
    )
  )
)

# Server ----
server <- function(input, output, session) {
  
  # Reactive values for fixed dates and baseline constants copied from your script:
  #end_date <- as.Date("2023-06-01")
  half_term_start <- as.Date("2022-10-21")
  half_term_end <- as.Date("2022-10-30")
  shopping_period_start <- as.Date("2022-12-01")
  shopping_period_end <- as.Date("2022-12-15")
  winter_holiday_start <- as.Date("2022-12-21")
  winter_holiday_end <- as.Date("2023-01-04")
  N_tot <- 60000000
  
  age_breaks <- c(0,5,18,65)
  
  # Preload polymod contact matrices once (these calls can be slow)
  # We'll construct term and holiday contact matrices inside the reactive model run
  data("polymod")
  contacts_all <- polymod$contacts
  polymod_base <- polymod
  
  ## Comparison data
  flu_dat <- read_csv("data/final/flu_2022_2023.csv")

  # Reactive wrapper that runs when the button is clicked or inputs change:
  #run_model <- eventReactive(input$run_model, {
  run_model <- reactive({
    library(shiny)
    library(socialmixr)
    library(pracma)
    library(deSolve)
    library(ggplot2)
    library(patchwork)
    library(tidyverse)
    library(ggpubr)
    library(data.table)
    library(lazymcmc)
    library(tidyr)
    
    age_groups_all <- c("[0,5)","[5,18)","[18,65)","65+")
    
    
    # Wrap in tryCatch so app doesn't crash if user environment lacks functions
    tryCatch({
      # source helper files (assume they're in the working directory or project folder)
      # If not found, error will be shown in app.
      source("auxiliary_funcs.R")
      source("sir_functions.R")
      
      start_date <- input$sim_start_date
      end_date <- input$sim_end_date
      dates <- seq(start_date,end_date,by="1 week")
      # read input parameters

      prop_immune_younger <- input$prop_immune_younger
      prop_immune_youngest <- input$prop_immune_younger2
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
      contacts <- contacts_all
      ## Create new contact survey with no school contacts, decreased work contacts, and increased home contacts
      ## Resample to give frac of work contacts
      ## Resample with replacement to give more home contacts
      ## Term time contact matrix
      ## Create term time contact matrix
      polymod_c_term <- contact_matrix(polymod_base,
                                       countries="United Kingdom",
                                       #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                       age.limits=age_breaks,
                                       symmetric=TRUE,
                                       missing.contact.age = "sample",
                                       missing.participant.age = "remove")
      C_term <- polymod_c_term$matrix
      row.names(C_term) <- colnames(C_term)
      
      ## Create school holiday contact matrix
      prop_home_contacts_in_hols <- input$prop_home_contacts_in_hols
      prop_work_contacts_in_hols <- input$prop_work_contacts_in_hols
      prop_rest_contacts_in_hols <- input$prop_rest_contacts_in_hols
      
      contacts <- contacts_all

      
      #contacts_schools <- contacts %>% filter(cnt_school == 1)  %>% sample_frac(0)
      ## Get all non-school, non-work and non-home contacts
      contacts_others <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home == 0) %>% sample_frac(prop_rest_contacts_in_hols,replace=TRUE)
      ## Subsample work contacts
      contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(prop_work_contacts_in_hols)
      ## Up-sample at home contacts
      contacts_home <- contacts %>% filter(cnt_home == 1) %>% sample_frac(prop_home_contacts_in_hols,replace=TRUE)
      ## Combine other, home and work contacts. No school contacts
      contacts_overall <- bind_rows(contacts_others, #contacts_schools, 
                                    contacts_work, contacts_home)
    
      ## Up-sample all contacts by some proportion
      #contacts_TermBreak <- contacts_overall[sample(1:nrow(contacts_overall), size = nrow(contacts_overall) * prop_rest_contacts_in_hols, replace = TRUE), ]
      
      polymodTermBreak <- polymod_base
      polymodTermBreak$contacts <- contacts_overall
      
      polymod_c_holidays <- contact_matrix(polymodTermBreak,
                                           countries="United Kingdom",
                                           #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                           age.limits=age_breaks,
                                           symmetric=TRUE,
                                           missing.contact.age = "sample",
                                           missing.participant.age = "remove")
      C_holidays <- polymod_c_holidays$matrix
      row.names(C_holidays) <- colnames(C_holidays)
      
      ## Create Christmas break (pre school holidays)
      prop_all_contacts_christmas <- input$prop_all_contacts_christmas
      contacts <- contacts_all
      contacts_schools <- contacts %>% filter(cnt_school == 1)#  %>% sample_frac(prop_school_contacts_in_Christmasperiod)
      #contacts_others <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home == 0)
      #contacts_work <- contacts %>% filter(cnt_work == 1) #%>% sample_frac(prop_work_contacts_in_Christmasperiod)
      #contacts_home <- contacts %>% filter(cnt_home == 1)
      contacts_non_school <- contacts %>% filter(cnt_school == 0) %>% sample_frac(prop_all_contacts_christmas, replace = TRUE)
      #contacts_home <- contacts_home[sample(1:nrow(contacts_home), size = nrow(contacts_home) * prop_home_contacts_in_Christmasperiod, replace = TRUE), ]
      ## Not school contacts!
      #contacts_overall <- bind_rows(contacts_others, contacts_work, contacts_home) %>% sample_frac(prop_all_contacts_christmas, replace = TRUE)
      contacts_overall <- contacts_non_school %>% bind_rows(contacts_schools)
      #contacts_ChristmasPeriod <- contacts_overall[sample(1:nrow(contacts_overall), size = nrow(contacts_overall) * prop_all_contacts_christmas, replace = TRUE), ]
      
      polymodChristmasPeriod <- polymod_base
      polymodChristmasPeriod$contacts <- contacts_overall
  
      polymod_c_christmas_break <- contact_matrix(polymodChristmasPeriod,
                                           countries="United Kingdom",
                                           #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                           age.limits=age_breaks,
                                           symmetric=TRUE,
                                           missing.contact.age = "sample",
                                           missing.participant.age = "remove")
      C_christmas_break <- polymod_c_christmas_break$matrix
      row.names(C_christmas_break) <- colnames(C_christmas_break)
      
      
      ## Create Christmas school holiday contact matrix
      prop_work_contacts_in_christmas <- input$prop_work_contacts_in_christmas
      prop_rest_contacts_in_christmas <- input$prop_rest_contacts_in_christmas
      prop_home_contacts_in_christmas <- input$prop_home_contacts_in_christmas
      
      contacts <- contacts_all
      
      ## No school contacts
      # contacts_schools <- contacts %>% filter(cnt_school == 1)  %>% sample_frac(0)
      #contacts_others <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home == 0)
      #contacts_work <- contacts %>% filter(cnt_work == 1)# %>% sample_frac(prop_work_contacts_in_christmas)
      contacts_non_school <- contacts %>% filter(cnt_school == 0) %>% sample_frac(prop_rest_contacts_in_christmas)
      contacts_home <- contacts_non_school %>% filter(cnt_home == 1) %>% sample_frac(prop_home_contacts_in_christmas, replace = TRUE)
      contacts_overall <- bind_rows(contacts_non_school, contacts_home)
      
      #contacts_overall <- bind_rows(contacts_others, contacts_schools, contacts_work, contacts_home) %>% sample_frac(prop_rest_contacts_in_christmas,replace=TRUE)
      #contacts_ChristmasBreak <- contacts_overall[sample(1:nrow(contacts_overall), size = nrow(contacts_overall) * prop_rest_contacts_in_christmas, replace = TRUE), ]
      
      polymodChristmasHoliday <- polymod_base
      polymodChristmasHoliday$contacts <- contacts_overall
      
      ## Create term time contact matrix
      polymod_c_xmas_holidays <- contact_matrix(polymodChristmasHoliday,
                                       countries="United Kingdom",
                                       #age.limits = c(0,5,15,25,35,45,55,65,75,85),
                                       age.limits=age_breaks,
                                       symmetric=TRUE,
                                       missing.contact.age = "sample",
                                       missing.participant.age = "remove")
      C_xmas_holidays <- polymod_c_xmas_holidays$matrix
      row.names(C_xmas_holidays) <- colnames(C_xmas_holidays)
      
   
      ########################################################
      ## SETUP MODEL FOR 2023/24 SEASON
      ########################################################
      
      ## From https://www.kennedycollegeoxford.org.uk/term-dates
      ## Holiday dates
      ## Setup contact matrix list
      school_days_weighted <- setup_holiday_tibble(start_date,end_date, half_term_start,half_term_end,
                                                   shopping_period_start,winter_holiday_end, smooth_time = 7)
      
      
      half_term_weighting <- setup_period_tibble(start_date, end_date,
                                                 half_term_start, half_term_end,
                                                 smooth_time = 7)
      
      christmas_holiday_weighting <- setup_period_tibble(start_date, end_date,
                                                         winter_holiday_start, winter_holiday_end,
                                                         smooth_time = 7)
      
      christmas_shopping_period <- setup_period_tibble(start_date, end_date,
                                                       shopping_period_start, shopping_period_end,
                                                       smooth_time = 7)
      ref_dates <- school_days_weighted$date
      
      # join weights by date to ensure alignment (will error if missing dates)
      weights_df <- tibble::tibble(date = ref_dates) %>%
        dplyr::left_join(half_term_weighting %>% dplyr::select(date, weight) %>% dplyr::rename(w_termBreak = weight), by = "date") %>%
        dplyr::left_join(christmas_shopping_period %>% dplyr::select(date, weight) %>% dplyr::rename(w_ChristmasPeriod = weight), by = "date") %>%
        dplyr::left_join(christmas_holiday_weighting %>% dplyr::select(date, weight) %>% dplyr::rename(w_ChristmasBreak = weight), by = "date")
      
      #  convert to holidayness, compute sum, normalise if > 1, compute final term weight
      weights_df <- weights_df %>%
        dplyr::mutate(
          h_termBreak       = 1 - w_termBreak,
          h_ChristmasPeriod = 1 - w_ChristmasPeriod,
          h_ChristmasBreak  = 1 - w_ChristmasBreak,
          sum_h = h_termBreak + h_ChristmasPeriod + h_ChristmasBreak,
          norm_factor = ifelse(sum_h > 1, 1 / sum_h, 1),
          # scaled holiday contributions (share when overlapping)
          h_termBreak = h_termBreak * norm_factor,
          h_ChristmasPeriod = h_ChristmasPeriod * norm_factor,
          h_ChristmasBreak = h_ChristmasBreak * norm_factor,
          # final term weight ensures weights sum to 1
          w_term = 1 - (h_termBreak + h_ChristmasPeriod + h_ChristmasBreak)
        ) %>%
        dplyr::select(date, w_term, h_termBreak, h_ChristmasPeriod, h_ChristmasBreak)
      
      prop_immune <- c(prop_immune_younger,prop_immune_youngest,prop_immune_older,prop_immune_oldest)
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
      C_term_use <- setup_C(C_term, N, beta_scales)
      C_holidays_use <- setup_C(C_holidays, N, beta_scales)
      C_christmas_break_use <- setup_C(C_christmas_break, N, beta_scales)
      C_xmas_holidays_use <- setup_C(C_xmas_holidays, N, beta_scales)
      
      C_list <- vector("list", nrow(weights_df))
      ## Create list of interpolated contact matrices
      for(i in 1:nrow(school_days_weighted)){
        w_term <- weights_df$w_term[i]
        h1     <- weights_df$h_termBreak[i]
        h2     <- weights_df$h_ChristmasPeriod[i]
        h3     <- weights_df$h_ChristmasBreak[i]
        
        C_list[[i]] <- C_term_use * w_term +
          C_holidays_use * h1 +
          C_christmas_break_use * h2 +
          C_xmas_holidays_use * h3
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
      date_key <- data.frame(t=1:nrow(ret_sum),date=dates)
      
      ## Pull out cumulative incidence numbers

      y_base_long <- y_base %>% pivot_longer(-time)
      y_base_long <- y_base_long %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
                                            age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
                                            immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
      total_inf <- y_base_long %>% filter(time == max(time),compartment == "inc")  %>% pull(value) %>% sum()

      N_from_sim <- y_base_long %>% filter(time == 1) %>% 
        filter(compartment != "inc") %>% 
        group_by(age) %>% 
        summarize(total_by_age=sum(value))

      cumulative_incidence <- y_base_long %>% filter(time == max(time)) %>% 
        filter(compartment == "inc")  %>% 
        group_by(age) %>% 
        summarize(total_inf=sum(value))

            cumulative_incidence <- cumulative_incidence %>% 
        left_join(N_from_sim) %>% 
        mutate(prop=total_inf/total_by_age) %>%
        rename("Age group"=age,"Cumulative influenza infections"=total_inf,"N"=total_by_age,"Proportion infected"=prop)
      cumulative_incidence$`Symptomatic` <- cumulative_incidence$`Cumulative influenza infections`*symp_frac
      cumulative_incidence$`Symptomatic and reported` <- cumulative_incidence$`Symptomatic`*reporting_rate
      
      cumulative_incidence$`Age group` <- age_groups_all
      
      cumulative_incidence_all <- tibble("Age group"="Total",
                                         "Cumulative influenza infections"=sum(cumulative_incidence$`Cumulative influenza infections`),
                                         "N"=sum(cumulative_incidence$N)) %>% 
        mutate(`Proportion infected`=`Cumulative influenza infections`/N) %>%
        mutate(`Symptomatic` = sum(cumulative_incidence$`Symptomatic`)) %>%
        mutate("Symptomatic and reported"=sum(cumulative_incidence$`Symptomatic`*reporting_rate))

      
      cumulative_incidence <- bind_rows(cumulative_incidence,cumulative_incidence_all)
      inc <- y_base_long %>% filter(compartment == "inc")
      inc <- inc %>% group_by(age, immunity) %>% mutate(value = value - lag(value, 1))
    
      # return list of objects needed for plotting / summary
      list(inc = ret_sum,
           date_key=date_key,
           cumulative_incidence=cumulative_incidence,
           contact_matrices=list(term=C_term,half_term=C_holidays, christmas_shopping=C_christmas_break,christmas_holidays=C_xmas_holidays),
           
           meta = list(start_date = start_date, end_date = end_date,plot_start_date=input$plot_start_date,
                       half_term_start = half_term_start, half_term_end = half_term_end,
                       reporting_rate=input$reporting_rate,
                       shopping_period_start = shopping_period_start,
                       winter_holiday_start = winter_holiday_start, winter_holiday_end = winter_holiday_end,
                       y_lim_max=input$y_lim_max))
    }, error = function(e) {
      # return the error to display in UI
      list(error = TRUE, message = paste("Error running model:", e$message))
    })
  })#, ignoreNULL = FALSE)
  
  # Summary text
  output$summary_text <- renderText({
    res <- run_model()
    if (isTRUE(res$error)) {
      return(res$message)
    } else {
      paste0("Dates: ", as.character(res$meta$start_date), " to ", as.character(res$meta$end_date),
             " | Estimated total infections (fraction of population): ",
             format(res$total_infections_fraction, digits = 4))
    }
  })
  
  # Plot
  output$inc_plot <- renderPlot({
    res <- run_model()
    if (isTRUE(res$error)) {
      plot.new()
      text(0.5, 0.5, res$message)
      return()
    }
    inc <- as.data.frame(res$inc)
    inc$t <- 1:nrow(inc)
    inc <- inc %>% pivot_longer(-t)
    
    colnames(inc) <- c("t","age_group","incidence")
    age_group_key <- c("inc_1_1"="[0,5)","inc_2_1"="[5,18)","inc_3_1"="[18,65)","inc_4_1"="65+")
    inc$age_group <- age_group_key[inc$age_group]
    inc$age_group <- factor(inc$age_group,levels=age_group_key)
    
    date_key <- res$date_key
    
    age_group_key1 <- c("0-4"="[0,5)","5-18"="[5,18)","19-64"="[18,65)","65+"="65+")
    
    
    flu_dat1 <- flu_dat %>% filter(date >= res$meta$start_date,date <= res$meta$end_date)
    flu_dat1$age_group <- age_group_key1[flu_dat1$age_group]
    flu_dat1$age_group <- factor(flu_dat1$age_group,levels = age_group_key1)
    
    
    rects <- data.frame(
      xmin = as.Date(c(half_term_start, shopping_period_start, winter_holiday_start)),
      xmax = as.Date(c(half_term_end,   shopping_period_end,   winter_holiday_end)),
      label = c("Half term", "Christmas period\n (pre holidays)", "Christmas\nholidays"),
      fill  = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
      stringsAsFactors = FALSE
    )
    # midpoints for label placement and y position at 95% of plot max
    rects$mid <- as.Date( (as.numeric(rects$xmin) + as.numeric(rects$xmax)) / 2, origin = "1970-01-01" )
    y_pos <- res$meta$y_lim_max * 0.95
    rects$y <- y_pos
    
    
    p <- ggplot(inc %>% left_join(date_key)) +
      geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                inherit.aes = FALSE, alpha = 0.25, show.legend = TRUE) +
      geom_text(data = rects, aes(x = mid, y = y, label = label), inherit.aes = FALSE,
                size = 3, fontface = "bold", vjust = 1) +
      geom_line(data=flu_dat1,aes(x=date,y=ILI_flu,col=age_group,linetype="2023/23 season data"),alpha=0.5,linewidth=1) +
      geom_line(aes(x = date, y = incidence, col = age_group,linetype="Model"),linewidth=1) +
      scale_color_brewer("Age group",palette="Set1") +
      scale_fill_brewer("Holiday period",palette="Set2") +
      scale_linetype_manual("Data source",values=c("2023/23 season data"="dashed","Model"="solid")) +
      #geom_vline(xintercept = c(half_term_start, half_term_end,shopping_period_start, winter_holiday_start, winter_holiday_end),
        #         linetype = "dashed", color = "red") +
      scale_x_date(breaks = "1 month",expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      coord_cartesian(ylim=c(0,res$meta$y_lim_max),xlim=c(res$meta$plot_start_date,res$meta$end_date)) +
      theme_bw() +
      theme(#panel.background=element_rect(fill="grey10"),
            #panel.grid.major=element_line(color="white"),
            axis.text=element_text(size=14),
            axis.title=element_text(size=16),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14)) +
      xlab("Date (end of reporting week)") +
      ylab("Reported influenza cases (weekly)") +
      labs(color = "Age group")
    
    
    
    # inc here is ret_sum (weekly reported), as a data.frame or matrix with rows = weeks
    weekly_totals <- rowSums(as.matrix(res$inc), na.rm = TRUE)   # weekly reported (all ages)
    peak_val <- if(length(weekly_totals)) max(weekly_totals, na.rm = TRUE) else NA
    peak_idx <- if(length(weekly_totals)) which.max(weekly_totals) else NA
    peak_date <- if(!is.na(peak_idx) && length(peak_idx)) res$date_key$date[peak_idx] else NA
    reported_peak <- peak_val
    peak_symptomatic <- reported_peak / res$meta$reporting_rate
    
    # cumulative symptomatic total (model symptomatic, not necessarily reported)
    cum_symp_total <- NA
    if(!is.null(res$cumulative_incidence)){
      if("Age group" %in% colnames(res$cumulative_incidence)){
        row_total <- res$cumulative_incidence %>% filter(`Age group` == "Total")
        if(nrow(row_total) == 1 && "Symptomatic" %in% colnames(row_total)){
          cum_symp_total <- row_total$Symptomatic
        } else {
          # fallback: sum per-age Symptomatic if Total row not present
          if("Symptomatic" %in% colnames(res$cumulative_incidence)){
            cum_symp_total <- sum(res$cumulative_incidence$Symptomatic, na.rm = TRUE)
          }
        }
      }
    }
    
    # formatted label text
    fmt_num <- function(x) if(is.na(x)) "n/a" else format(round(x), big.mark = ",", scientific = FALSE)
    label_text <- paste0(
      "Peak (weekly): ", fmt_num(peak_symptomatic), "\n",
      "Peak (weekly reported): ", fmt_num(peak_val), "\n",
      "Date of peak: ", if(is.na(peak_date)) "n/a" else format(peak_date, "%Y-%m-%d"), "\n",
      "Cumulative symptomatic: ", fmt_num(cum_symp_total)
    )
    
    # choose placement: top-right inside plot area
    x_pos <- as.Date(res$meta$end_date) - 3       # 7 days left of end; tweak if needed
    y_pos <- 0.98 * res$meta$y_lim_max            # near top
    
    # add box label to plot
    p <- p + annotate("label",
                      x = x_pos, y = y_pos,
                      label = label_text,
                      hjust = 1, vjust = 1,
                      size = 3.5, fontface = "bold",
                      fill = "white", alpha = 0.8)
    
    
    print(p)
  })
  
  output$cum_table <- DT::renderDataTable({
    res <- run_model()
    if (isTRUE(res$error)) {
      return(data.frame(Message = res$message))
    }
    
    df <- res$cumulative_incidence
    
    # Identify numeric columns except "Proportion infected"
    numcols <- sapply(df, is.numeric)
    df[numcols] <- lapply(df[numcols], function(x) signif(x, 3))
    
    # Build the DT table
    datatable(df,
              rownames = FALSE,
              options = list(
                pageLength = 10,
                dom = 't',
                ordering = FALSE
              )
    ) %>%
      # Highlight 65+ and Total rows
      DT::formatStyle(
        'Age group',
        target = 'row',
        backgroundColor = DT::styleEqual(
          c("65+", "Total"),
          c("#ffe0e0", "#e0ffe0")   # light red for 65+, light green for Total
        )
      )
  })
  
  output$contact_matrices_plot <- renderPlot({
  res <- run_model()
  if (isTRUE(res$error)) {
    plot.new(); text(0.5, 0.5, res$message); return()
  }
  mats <- res$contact_matrices
  if (is.null(mats) || length(mats) == 0) {
    plot.new(); text(0.5, 0.5, "No contact matrices returned by model"); return()
  }

  # Convert each matrix to long format and tag with name
  df_list <- lapply(names(mats), function(nm) {
    m <- mats[[nm]]
    # ensure matrix has rownames / colnames
    rns <- rownames(m); cns <- colnames(m)
    if (is.null(rns)) rns <- paste0("r", seq_len(nrow(m)))
    if (is.null(cns)) cns <- paste0("c", seq_len(ncol(m)))
    # as.table -> data.frame gives Var1,Var2,Freq
    d <- as.data.frame(as.table(m))
    colnames(d) <- c("row", "col", "value")
    d$matrix_name <- nm
    # preserve ordering so heatmaps look consistent
    d$row <- factor(as.character(d$row), levels = unique(rns))
    d$col <- factor(as.character(d$col), levels = unique(cns)
                    ) 
    d
  })

  plot_df <- do.call(rbind, df_list)

  # Draw 2x2 facetted heatmap (matrix rows on y-axis, cols on x-axis)
  p <- ggplot(plot_df, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    facet_wrap(~ matrix_name, ncol = 2) +
    labs(x = NULL, y = NULL, fill = "Contact\nrate") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(size = 14),
      axis.title = element_text(size=14),
      axis.text = element_text(size = 12)
    ) +
    coord_fixed()
  print(p)
})
  
  # ---- Download plot ----
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("model_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      res <- run_model()
      if (isTRUE(res$error)) stop("Model error: ", res$message)
      
      start_date <- res$meta$start_date
      end_date   <- res$meta$end_date
      half_term_start <- res$meta$half_term_start
      half_term_end   <- res$meta$half_term_end
      winter_holiday_start <- res$meta$winter_holiday_start
      winter_holiday_end   <- res$meta$winter_holiday_end
      shopping_period_start <- res$meta$shopping_period_start
      
      inc <- as.data.frame(res$inc)
      inc$t <- 1:nrow(inc)
      inc <- inc %>% pivot_longer(-t)
      colnames(inc) <- c("t","age_group","incidence")
      age_group_key <- c("inc_1_1"="[0,5)","inc_2_1"="[5,18)","inc_3_1"="[18,65)","inc_4_1"="65+")
      inc$age_group <- age_group_key[inc$age_group]
      inc$age_group <- factor(inc$age_group, levels=age_group_key)
      date_key <- res$date_key
      
      flu_dat1 <- flu_dat %>% filter(date >= start_date, date <= end_date)
      age_group_key1 <- c("0-4"="[0,5)","5-18"="[5,18)","19-64"="[18,65)","65+"="65+")
      flu_dat1$age_group <- age_group_key1[flu_dat1$age_group]
      flu_dat1$age_group <- factor(flu_dat1$age_group, levels = age_group_key1)
      
      
      rects <- data.frame(
        xmin = as.Date(c(half_term_start, shopping_period_start, winter_holiday_start)),
        xmax = as.Date(c(half_term_end,   shopping_period_end,   winter_holiday_end)),
        label = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
        fill  = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
        stringsAsFactors = FALSE
      )
      rects$mid <- as.Date( (as.numeric(rects$xmin) + as.numeric(rects$xmax)) / 2, origin = "1970-01-01" )
      y_pos <- input$y_lim_max * 0.95
      rects$y <- y_pos
      
      
      p <- ggplot(inc %>% left_join(date_key)) +
        geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                  inherit.aes = FALSE, alpha = 0.25, show.legend = TRUE) +
        geom_text(data = rects, aes(x = mid, y = y, label = label), inherit.aes = FALSE,
                  size = 3, fontface = "bold", vjust = 1) +
        geom_line(data=flu_dat1,aes(x=date,y=ILI_flu,col=age_group,linetype="2023/23 season data"),alpha=0.5,linewidth=1) +
        geom_line(aes(x = date, y = incidence, col = age_group,linetype="Model"),linewidth=1) +
        scale_color_brewer("Age group",palette="Set1") +
        scale_fill_brewer("Holiday period",palette="Set2") +
        scale_linetype_manual("Data source",values=c("2023/23 season data"="dashed","Model"="solid")) +
        #geom_vline(xintercept = c(half_term_start, half_term_end,shopping_period_start, winter_holiday_start, winter_holiday_end),
        #         linetype = "dashed", color = "red") +
        scale_x_date(breaks = "1 month",expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        coord_cartesian(ylim=c(0,res$meta$y_lim_max),xlim=c(res$meta$plot_start_date,res$meta$end_date)) +
        theme_bw() +
        theme(#panel.background=element_rect(fill="grey10"),
          #panel.grid.major=element_line(color="white"),
          axis.text=element_text(size=8),
          axis.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.title=element_text(size=8)) +
        xlab("Date (end of reporting week)") +
        ylab("Reported influenza cases (weekly)") +
        labs(color = "Age group")
      
      
      # inc here is ret_sum (weekly reported), as a data.frame or matrix with rows = weeks
      weekly_totals <- rowSums(as.matrix(res$inc), na.rm = TRUE)   # weekly reported (all ages)
      peak_val <- if(length(weekly_totals)) max(weekly_totals, na.rm = TRUE) else NA
      peak_idx <- if(length(weekly_totals)) which.max(weekly_totals) else NA
      peak_date <- if(!is.na(peak_idx) && length(peak_idx)) res$date_key$date[peak_idx] else NA
      
      reported_peak <- peak_val
      peak_symptomatic <- reported_peak / input$reporting_rate
      
      # cumulative symptomatic total (model symptomatic, not necessarily reported)
      cum_symp_total <- NA
      if(!is.null(res$cumulative_incidence)){
        if("Age group" %in% colnames(res$cumulative_incidence)){
          row_total <- res$cumulative_incidence %>% filter(`Age group` == "Total")
          if(nrow(row_total) == 1 && "Symptomatic" %in% colnames(row_total)){
            cum_symp_total <- row_total$Symptomatic
          } else {
            # fallback: sum per-age Symptomatic if Total row not present
            if("Symptomatic" %in% colnames(res$cumulative_incidence)){
              cum_symp_total <- sum(res$cumulative_incidence$Symptomatic, na.rm = TRUE)
            }
          }
        }
      }
      
      # formatted label text
      fmt_num <- function(x) if(is.na(x)) "n/a" else format(round(x), big.mark = ",", scientific = FALSE)
      label_text <- paste0(
        "Peak (weekly): ", fmt_num(peak_symptomatic), "\n",
        "Peak (weekly reported): ", fmt_num(peak_val), "\n",
        "Date of peak: ", if(is.na(peak_date)) "n/a" else format(peak_date, "%Y-%m-%d"), "\n",
        "Cumulative symptomatic: ", fmt_num(cum_symp_total)
      )
      
      # choose placement: top-right inside plot area
      x_pos <- as.Date(res$meta$end_date) - 3       # 7 days left of end; tweak if needed
      y_pos <- 0.98 * res$meta$y_lim_max            # near top
      
      # add box label to plot
      p <- p + annotate("label",
                        x = x_pos, y = y_pos,
                        label = label_text,
                        hjust = 1, vjust = 1,
                        size = 3.5, fontface = "bold",
                        fill = "white", alpha = 0.8)
      
      ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
    }
  )
  
  # ---- Download table ----
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("cumulative_incidence_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- run_model()
      if (isTRUE(res$error)) stop("Model error: ", res$message)
      write.csv(res$cumulative_incidence, file, row.names = FALSE)
    }
  )
}




# Run app
shinyApp(ui = ui, server = server)
