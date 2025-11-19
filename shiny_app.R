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
library(lazymcmc)
library(tidyr)
#library(shinycssloaders) # optional, shows spinner while plotting

# UI ----
ui <- fluidPage(
  titlePanel("Influenza H3N2 k-clade model (interactive)"),
  sidebarLayout(
    sidebarPanel(
      h4("Model parameters"),
      # Contacts / holiday behavior
      sliderInput("prop_work_contacts_in_hols", "Proportion of work contacts kept in holidays",
                  min = 0, max = 1, value = 0.4, step = 0.05),
      sliderInput("prop_home_contacts_in_hols", "Multiplier for home contacts in holidays",
                  min = 1, max = 5, value = 2, step = 0.1),
      hr(),
      # Immunity proportions by age-groups
      sliderInput("prop_immune_younger", "Proportion immune (0-4 yrs)",
                  min = 0, max = 1, value = 0.5, step = 0.01),
      sliderInput("prop_immune_younger2", "Proportion immune (5-14 yrs)",
                  min = 0, max = 1, value = 0.5, step = 0.01),
      sliderInput("prop_immune_older", "Proportion immune (15-64 yrs)",
                  min = 0, max = 1, value = 0.7, step = 0.01),
      sliderInput("prop_immune_oldest", "Proportion immune (65+ yrs)",
                  min = 0, max = 1, value = 0.7, step = 0.01),
      hr(),
      # Other epidemiological parameters
      numericInput("alpha1", "Relative susceptibility alpha1 (immune class 1)", value = 1.2, step = 0.05),
      numericInput("alpha2", "Relative susceptibility alpha2 (immune class 2)", value = 1, step = 0.05),
      numericInput("R0", "Basic reproductive number R0", value = 1.2, step = 0.01),
      numericInput("gamma", "Infectious period (days) gamma", value = 3, step = 0.1),
      numericInput("seed_size", "Seed size (number of initial infections)", value = 200, step = 10),
      hr(),
      # Reporting / symptomatic fraction
      sliderInput("reporting_rate", "Reporting rate (fraction of symptomatic reported)",
                  min = 0, max = 1, value = 0.01, step = 0.001),
      # A small note / run button
      actionButton("run_model", "Run model"),
      width = 3
    ),
    mainPanel(
      # Top area can show a quick summary
      fluidRow(
        column(12,
               wellPanel(
                 h4("Run summary"),
                 textOutput("summary_text")
               ))
      ),
      # Plot at bottom right in its own panel
      fluidRow(
        column(12,
               wellPanel(
                 h4("Model output: reported symptomatic daily cases"),
                 plotOutput("inc_plot", height = "600px")
               ))
      ),
      width = 9
    )
  )
)

# Server ----
server <- function(input, output, session) {

  # Reactive values for fixed dates and baseline constants copied from your script:
  start_date <- as.Date("2023-09-01")
  end_date <- as.Date("2024-05-01")
  half_term_start <- as.Date("2023-10-21")
  half_term_end <- as.Date("2023-10-30")
  winter_holiday_start <- as.Date("2023-12-15")
  winter_holiday_end <- as.Date("2024-01-04")
  N_tot <- 60000000
  # Symptomatic fraction per age group (kept fixed as in your script)
  symp_frac <- c(0.75,0.5,0.4,0.4,0.4,0.5,0.5,0.75,0.75)
  
  # Preload polymod contact matrices once (these calls can be slow)
  # We'll construct term and holiday contact matrices inside the reactive model run
  data("polymod")
  contacts_all <- polymod$contacts
  polymod_base <- polymod
  
  # read input parameters
  prop_work_contacts_in_hols <-0.4
  prop_home_contacts_in_hols <- 2
  
  prop_immune_younger <- 0.5
  prop_immune_younger2 <- 0.5
  prop_immune_older <- 0.5
  prop_immune_oldest <- 0.5
  
  alphas <- c(1,0.5)
  R0 <- 2
  gamma <-3
  seed_size <- 200
  reporting_rate <- 0.001
  
  
  # Reactive wrapper that runs when the button is clicked or inputs change:
  run_model <- eventReactive(input$run_model, {
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
    # Wrap in tryCatch so app doesn't crash if user environment lacks functions
    tryCatch({
      # source helper files (assume they're in the working directory or project folder)
      # If not found, error will be shown in app.
      source("auxiliary_funcs.R")
      source("sir_functions.R")
      
      # read input parameters
      prop_work_contacts_in_hols <- input$prop_work_contacts_in_hols
      prop_home_contacts_in_hols <- input$prop_home_contacts_in_hols
      
      prop_immune_younger <- input$prop_immune_younger
      prop_immune_younger2 <- input$prop_immune_younger2
      prop_immune_older <- input$prop_immune_older
      prop_immune_oldest <- input$prop_immune_oldest
      
      alphas <- c(input$alpha1, input$alpha2)
      R0 <- input$R0
      gamma <- input$gamma
      seed_size <- input$seed_size
      reporting_rate <- input$reporting_rate
      
      # Build holiday-polymod contacts as in your script
      contacts <- contacts_all
      
      contacts_no_schools <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home == 0)
      contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(prop_work_contacts_in_hols)
      contacts_home <- contacts %>% filter(cnt_home == 1)
      contacts_home <- contacts_home[sample(1:nrow(contacts_home), size = nrow(contacts_home) * prop_home_contacts_in_hols, replace = TRUE), ]
      contacts_holidays <- bind_rows(contacts_no_schools, contacts_work, contacts_home)
      
      polymod1 <- polymod_base
      polymod1$contacts <- contacts_holidays
      
      # Age limits used in your script
      age.limits <- c(0,5,15,25,35,45,55,65,75,85)
      
      polymod_c_holidays <- contact_matrix(polymod1,
                                           countries = "United Kingdom",
                                           age.limits = age.limits,
                                           symmetric = TRUE,
                                           missing.contact.age = "sample",
                                           missing.participant.age = "remove")
      C_holidays <- polymod_c_holidays$matrix
      row.names(C_holidays) <- colnames(C_holidays)
      
      polymod_c_term <- contact_matrix(polymod_base,
                                       countries = "United Kingdom",
                                       age.limits = age.limits,
                                       symmetric = TRUE,
                                       missing.contact.age = "sample",
                                       missing.participant.age = "remove")
      C_term <- polymod_c_term$matrix
      row.names(C_term) <- colnames(C_term)
      
      # setup holiday weighting
      school_days_weighted <- setup_holiday_tibble(start_date, end_date,
                                                   half_term_start, half_term_end,
                                                   winter_holiday_start, winter_holiday_end,
                                                   smooth_time = 7)
      
      prop_immune <- c(prop_immune_younger,
                       prop_immune_younger2,
                       rep(prop_immune_older, 5),
                       prop_immune_oldest,
                       prop_immune_oldest)
      
      # population props and N
      N_props <- polymod_c_term$participants$proportion
      N_age_classes <- length(N_props)
      N_immunity_classes <- length(alphas)
      
      N_props_long <- c(N_props * (1 - prop_immune), N_props * prop_immune)
      N <- matrix(N_props_long * N_tot, ncol = N_immunity_classes, nrow = N_age_classes)
      
      beta_scales <- rep(1, N_age_classes)
      beta_par <- get_beta(C_term, polymod_c_term$participants$proportion, gamma, R0)
      
      C_use_holiday <- setup_C(C_holidays, N, beta_scales)
      C_use_term <- setup_C(C_term, N, beta_scales)
      
      C_list <- list()
      for (i in 1:nrow(school_days_weighted)) {
        C_list[[i]] <- C_use_term * school_days_weighted$weight[i] + C_use_holiday * (1 - school_days_weighted$weight[i])
      }
      
      # time indexing
      start_day <- as.numeric(start_date)
      end_day <- as.numeric(end_date)
      date_seq <- seq(start_date, end_date, by = "1 day")
      ts <- seq(start_day, end_day, by = 1) - start_day + 1
      
      # Run your ODE model (epi_ode_size) — you used return_compartments = TRUE
      y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts = ts,
                             alphas = alphas, age_seed = 4, immunity_seed = 1,
                             seed_size = seed_size, return_compartments = TRUE)
      
      # Transform to long and compute daily reported symptomatic cases per age group
      y_base <- y_base %>% pivot_longer(-time)
      y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify = TRUE)[, 1],
                                  age = as.integer(str_split(name, "_", simplify = TRUE)[, 2]),
                                  immunity = as.integer(str_split(name, "_", simplify = TRUE)[, 3]))
      
      inc <- y_base %>% filter(compartment == "inc")
      inc <- inc %>% group_by(age, immunity) %>% mutate(value = value - lag(value, 1))
      inc <- inc %>% left_join(polymod_c_term$participants %>% mutate(age = 1:nrow(polymod_c_term$participants)))
      inc$age.group <- factor(inc$age.group, levels = c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", "[55,65)", "[65,75)", "75+"))
      
      # sum across immunity classes to get total infections by age group
      inc <- inc %>% group_by(age.group, time) %>% summarize(value = sum(value, na.rm = TRUE)) %>% ungroup()
      
      symp_frac_dat <- data.frame(age.group = c("[0,5)", "[5,15)", "[15,25)", "[25,35)", "[35,45)", "[45,55)", "[55,65)", "[65,75)", "75+"),
                                  symp_frac = symp_frac)
      
      inc <- inc %>% left_join(symp_frac_dat, by = "age.group")
      
      date_key <- data.frame(time = ts, date = as.Date("2023-09-01") + ts - 1)
      
      inc <- inc %>% left_join(date_key, by = "time")
      
      # compute reported symptomatic cases
      inc <- inc %>% mutate(reported = value * symp_frac * reporting_rate)
      
      # return list of objects needed for plotting / summary
      list(inc = inc,
           total_infections_fraction = sum(y_base %>% filter(compartment == "inc") %>% pull(value), na.rm = TRUE) / sum(N),
           meta = list(start_date = start_date, end_date = end_date,
                       half_term_start = half_term_start, half_term_end = half_term_end,
                       winter_holiday_start = winter_holiday_start, winter_holiday_end = winter_holiday_end))
    }, error = function(e) {
      # return the error to display in UI
      list(error = TRUE, message = paste("Error running model:", e$message))
    })
  }, ignoreNULL = FALSE)
  
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
    
    inc <- res$inc
    
    p <- ggplot(inc %>% ungroup()) +
      geom_line(aes(x = date, y = reported, col = age.group)) +
      scale_color_viridis_d() +
      geom_vline(xintercept = c(half_term_start, half_term_end, winter_holiday_start, winter_holiday_end),
                 linetype = "dashed", color = "red") +
      scale_x_date(breaks = "1 month") +
      theme_bw() +
      xlab("Date") +
      ylab("Reported symptomatic cases (daily)") +
      labs(color = "Age group")
    
    print(p)
  })
}

# Run app
shinyApp(ui = ui, server = server)
