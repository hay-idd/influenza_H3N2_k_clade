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
library(shinycssloaders)
library(bslib)

source("vars.R")
source("helper_funcs.R")

ui <- fluidPage(
  theme = bs_theme(version = 5,
                   # bootswatch = "litera"),
                   # bootswatch = "flatly"),
                   bootswatch = "cosmo"),
  
  # put this near the top of fluidPage
  tags$style(HTML("
    /* Accordion titles + bodies */
    .accordion .accordion-button {
      font-size: 14px !important;
      padding: 0.35rem 0.6rem !important;
    }
    .accordion .accordion-body {
      font-size: 12px;
      padding: 0.4rem 0.6rem !important;
    }

    /* Make labels + inputs inside the sidebar smaller */
    .compact-sidebar .form-label,
    .compact-sidebar .control-label,
    .compact-sidebar .shiny-input-container {
      font-size: 12px;
    }

    .compact-sidebar .form-control,
    .compact-sidebar .form-select {
      font-size: 12px;
      padding: 0.1rem 0.35rem;
      height: calc(1.4em + 0.4rem + 2px);
    }

    /* Sliders */
    .compact-sidebar .form-range {
      height: 0.25rem;
    }
  ")),
  
  
  
  # ---- UI ----
  titlePanel("Influenza H3N2 model"),
  fluidRow(
    column(9, tags$div(style = "margin-bottom:10px;")),
    column(3, align = "right",
           actionButton("reset_defaults", "Reset defaults", class = "btn-primary")
    )
  ),
  tabsetPanel(
    tabPanel("Model",
             sidebarLayout(
               sidebarPanel(
                 div(class = "compact-sidebar",
                     
                     actionButton("run_model", "Run Model"),
                     hr(),
                     
                     actionButton("set_ref", 
                                  "Set current scenario as new reference"),
                     hr(),
                     
                     fileInput("file", 
                               "[Optional] Upload custom data as a csv file. See About tab for format.", 
                               accept = ".csv"),
                     hr(),
                     
                     # single combined download
                     downloadButton("download_all", "Download results (zip)"),
                     hr(),
                     
                     # ---- Collapsible sections ----
                     accordion(
                       open = c("Model parameters (core)"),
                       id = "sidebar_accordion",
                       
                       # ---- Core parameters ----
                       accordion_panel(
                         title = "Model parameters (core)",
                         numericInput("R0", "Basic reproductive number R0", value = defaults$R0, step = 0.01),
                         numericInput("gamma", "Infectious period (days) gamma", value = defaults$gamma, step = 0.1),
                         br(),
                         
                         sliderInput("overall_immune_escape", "Overall immune escape (multiplier on immunity - lower value means more antigenic drift )",
                                     min = 0, max = 2, value = defaults$overall_immune_escape, step = 0.01),
                         
                         
                         numericInput("seed_size", "Seed size (number of initial infections)", 
                                      value = defaults$seed_size, step = 10),
                         br(),
                         
                         sliderInput("sim_start_date", "Seed date",
                                     min = as.Date("2022-07-01"),
                                     max = as.Date("2022-12-31"),
                                     value = defaults$sim_start_date,
                                     timeFormat = "%b-%d"),
                         sliderInput("sim_end_date", "Simulation end date",
                                     min = as.Date("2023-01-01"),
                                     max = as.Date("2023-07-01"),
                                     value = defaults$sim_end_date,
                                     timeFormat = "%b-%d"),
                         br(),
                         
                         sliderInput("reporting_rate", "Reporting rate (fraction of symptomatic reported)",
                                     min = 0, max = 0.01, value = defaults$reporting_rate, step = 0.0001)
                         
                         # numericInput("initial_immune_frac", 
                         #              "Initial proportion of the population fully immune", 
                         #              value = defaults$initial_immune_frac, step = 0.01),
                         # numericInput("alpha2", 
                         #              "Relative susceptibility of immune class", 
                         #              value = defaults$alpha2, step = 0.05)
                       ),
                       
                       hr(),
                       
                       # ---- Plot options (matches defaults$Plot options tab) ----
                       accordion_panel(
                         title = "Plot options",
                         numericInput("y_lim_max", "Maximum y-axis value", value = defaults$y_lim_max, step = 100),
                         numericInput("growth_y_scale", "Growth plot y-scale (symmetric)", value = defaults$growth_y_scale, step = 0.1),
                         sliderInput("plot_start_date", "Plot start date",
                                     min = as.Date("2022-07-01"),
                                     max = as.Date("2022-12-31"),
                                     value = defaults$plot_start_date,
                                     timeFormat = "%b-%d")
                       ),
                       hr(),
                       
                       
                       
                       # ---- Contact patterns ----
                       accordion_panel(
                         title="Demographics & Contact patterns",
                         
                         sliderInput("N_tot", "Population Size",
                                     min = 10e6,
                                     max = 150e6,
                                     value = defaults$N_tot),
                         
                         h5("School holiday contacts"),
                         sliderInput("prop_home_contacts_in_hols", "Multiplier for home contacts in breaks",
                                     min = 0, max = 3, value = defaults$prop_home_contacts_in_hols, step = 0.01),
                         sliderInput("prop_work_contacts_in_hols", "Proportion of work contacts kept in school holidays",
                                     min = 0, max = 1, value = defaults$prop_work_contacts_in_hols, step = 0.01),
                         sliderInput("prop_rest_contacts_in_hols", "Multiplier for non-school or work contacts in breaks",
                                     min = 0.5, max = 3, value = defaults$prop_rest_contacts_in_hols, step = 0.01),
                         hr(),
                         h5("Christmas period contacts (before school holiday)"),
                         sliderInput("prop_all_contacts_christmas", "Multiplier for all non-school contacts in Christmas period",
                                     min = 0.5, max = 3, value = defaults$prop_all_contacts_christmas, step = 0.01),
                         hr(),
                         h5("Christmas holiday contacts"),
                         sliderInput("prop_rest_contacts_in_christmas", "Proportion of all contacts kept over Christmas",
                                     min = 0, max = 1, value = defaults$prop_rest_contacts_in_christmas, step = 0.01),
                         sliderInput("prop_home_contacts_in_christmas", "Multiplier for home contacts over Christmas holiday",
                                     min = 1, max = 5, value = defaults$prop_home_contacts_in_christmas, step = 0.01)
                       ),
                       
                       hr(),
                       
                       # ---- Immunity and disease ----
                       accordion_panel(
                         title="Immunity and disease",
                         sliderInput("prop_immune_younger", "Proportion immune (0-4 yrs)",
                                     min = 0, max = 1, value = defaults$prop_immune_younger, step = 0.01),
                         sliderInput("prop_immune_younger2", "Proportion immune (5-18 yrs)",
                                     min = 0, max = 1, value = defaults$prop_immune_younger2, step = 0.01),
                         sliderInput("prop_immune_older", "Proportion immune (19-64 yrs)",
                                     min = 0, max = 1, value = defaults$prop_immune_older, step = 0.01),
                         sliderInput("prop_immune_oldest", "Proportion immune (65+ yrs)",
                                     min = 0, max = 1, value = defaults$prop_immune_oldest, step = 0.01),
                         
                         hr(),
                         h4("Symptomatic fraction by age-groups"),
                         sliderInput("symp_1", "Symptomatic frac (0-4 yrs)", min = 0, max = 1, value = defaults$symp_1, step = 0.01),
                         sliderInput("symp_2", "Symptomatic frac (5-18 yrs)", min = 0, max = 1, value = defaults$symp_2, step = 0.01),
                         sliderInput("symp_3", "Symptomatic frac (19-64 yrs)", min = 0, max = 1, value = defaults$symp_3, step = 0.01),
                         sliderInput("symp_4", "Symptomatic frac (65+ yrs)", min = 0, max = 1, value = defaults$symp_4, step = 0.01)
                       )
                     )
                 ),
                 
                 width = 3
               ),
               
               mainPanel(
                 fluidRow(
                   column(12,
                          wellPanel(
                            h4("Status"),
                            textOutput("status"),
                            h4(),
                            
                            h4("Reported symptomatic cases (weekly)"),
                            plotOutput("combined_plot", height = "760px"), # combined plot includes growth plot
                            hr(),
                            
                            h4("Cumulative cases"),
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
    ),
    tabPanel("About",
             fluidRow(
               column(9,
                      h3("About this model"),
                      p("Placeholder: description of the model and data goes here. (You can edit this later.)"),
                      br(),
                      h3("Format for custom CSV file."),
                      p("Please ensure your file conforms to the following format."),
                      DT::dataTableOutput("about_csv")
               ),
               column(3,
                      h4("Contact / Disclaimer"),
                      p(HTML("<strong>Contact:</strong> James Hay &lt;james.hay@ndm.ox.ac.uk&gt;")),
                      p(HTML("<strong>Disclaimer:</strong> This tool is a work in progress. Use results for exploration only; not for operational decision-making.")),
                      br(),
                      p("Last updated: ", Sys.Date())
               )
             )
    )
  )
)


server <- function(input, output, session) {
  
  # ensure helper functions from your files are available
  source("auxiliary_funcs.R")
  source("sir_functions.R")
  
  ref_path <- reactiveVal(NULL)
  
  
  output$about_csv <- DT::renderDataTable({
    df <- read.csv("data/final/flu_2022_2023.csv", 
                   stringsAsFactors = FALSE) %>%
      head(10)
    
    DT::datatable(df, options = list(scrollX = TRUE)) |>
      DT::formatRound(
        columns = 5:12,
        digits = 2
      )
  }, 
  options = list(pageLength = 10, scrollX = TRUE))
  
  
  flu_dat <- reactive({
    if (is.null(input$file)) {
      read_csv("data/final/flu_2022_2023.csv")
    } else {
      read_csv(input$file$datapath)
    }
  })
  
  model_running <- reactiveVal(FALSE)
  
  # helper that builds & runs the model; made reactive so UI changes update automatically
  run_model <- eventReactive(input$run_model,{
    model_running(TRUE)
    on.exit(model_running(FALSE), add = TRUE)
    
    
    # local inputs
    start_date <- input$sim_start_date
    end_date   <- input$sim_end_date
    dates <- seq(start_date, end_date, by = "1 week")
    
    prop_immune_younger <- input$prop_immune_younger
    prop_immune_youngest <- input$prop_immune_younger2
    prop_immune_older <- input$prop_immune_older
    prop_immune_oldest <- input$prop_immune_oldest
    
    alphas <- c(1, defaults$alpha2)
    alphas <- alphas/mean(alphas)
    R0 <- input$R0
    gamma <- input$gamma
    seed_size <- input$seed_size
    reporting_rate <- input$reporting_rate
    symp_frac <- c(input$symp_1, input$symp_2, input$symp_3, input$symp_4)
    
    start_day <- as.numeric(start_date)
    end_day_num <- as.numeric(end_date)
    ts <- seq(start_day, end_day_num, by = 1) - start_day + 1
    
    # contact matrices
    cm <- build_contact_matrices(input)
    C_term <- cm$term
    C_holidays <- cm$term_break
    C_christmas_break <- cm$christmas_period
    C_xmas_holidays <- cm$christmas_holiday
    polymod_c_term <- cm$polymod_term
    
    # weights
    weights_res <- make_weights_df(start_date, end_date, smooth_time = 7)
    weights_df <- weights_res$weights
    
    # apply overall immune escape multiplier to age-group immunity sliders
    overall_mult <- input$overall_immune_escape
    raw_prop_immune <- c(prop_immune_younger, prop_immune_youngest, prop_immune_older, prop_immune_oldest)
    prop_immune <- pmin(1, overall_mult * raw_prop_immune)
    
    # population & N
    N_props <- polymod_c_term$demography$proportion
    N_age_classes <- length(N_props)
    N_props_long <- c(N_props * (1 - prop_immune), N_props * (prop_immune))
    N <- matrix(N_props_long * defaults$N_tot, ncol = length(alphas), nrow = N_age_classes)
    
    beta_par <- get_beta(C_term, polymod_c_term$participants$proportion, gamma, R0)
    beta_scales <- rep(1, N_age_classes)
    C_term_use <- setup_C(C_term, N, beta_scales)
    C_holidays_use <- setup_C(C_holidays, N, beta_scales)
    C_christmas_break_use <- setup_C(C_christmas_break, N, beta_scales)
    C_xmas_holidays_use <- setup_C(C_xmas_holidays, N, beta_scales)
    
    # create list of contact matrices for each day
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
    y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts = ts,
                           alphas = alphas,
                           initial_immune_frac = defaults$initial_immune_frac,
                           age_seed = 3, immunity_seed = 1,
                           seed_size = seed_size, return_compartments = TRUE)
    
    use_cols <- which(colnames(y_base) %like% "inc")
    ret <- as.matrix(y_base[, use_cols])
    ret <- apply(ret, 2, function(x) c(0, diff(x))) # daily incidence
    ret2 <- ret[, seq(1, ncol(ret), 2)] + ret[, seq(2, ncol(ret), 2)]
    ret2 <- t(reporting_rate * symp_frac * t(ret2)) # apply symptomatic * reporting
    
    row_groups <- gl(nrow(ret2) %/% 7 + (nrow(ret2) %% 7 > 0), 7, nrow(ret2))
    ret_sum <- aggregate(ret2, by = list(row_groups), FUN = sum)[, -1]
    date_key <- data.frame(t = 1:nrow(ret_sum), date = dates)
    
    # cumulative incidence & formatting
    y_base_long <- y_base %>% pivot_longer(-time)
    y_base_long <- y_base_long %>% mutate(compartment = str_split(name, "_", simplify = TRUE)[,1],
                                          age = as.integer(str_split(name, "_", simplify = TRUE)[,2]),
                                          immunity = as.integer(str_split(name, "_", simplify = TRUE)[,3]))
    N_from_sim <- y_base_long %>% filter(time == 1) %>%
      filter(compartment != "inc") %>%
      group_by(age) %>%
      summarize(total_by_age = sum(value), .groups = "drop")
    cumulative_incidence <- y_base_long %>% filter(time == max(time)) %>%
      filter(compartment == "inc") %>%
      group_by(age) %>%
      summarize(total_inf = sum(value), .groups = "drop")
    cumulative_incidence <- cumulative_incidence %>%
      left_join(N_from_sim) %>%
      mutate(prop = total_inf / total_by_age) %>%
      rename("Age group" = age, "Cumulative influenza infections" = total_inf, "N" = total_by_age, "Proportion infected" = prop)
    cumulative_incidence$Symptomatic <- cumulative_incidence$`Cumulative influenza infections` * symp_frac
    cumulative_incidence$`Symptomatic and reported` <- cumulative_incidence$Symptomatic * reporting_rate
    age_groups_all <- c("[0,5)", "[5,18)", "[18,65)", "65+")
    cumulative_incidence$`Age group` <- age_groups_all
    cumulative_incidence_all <- tibble("Age group" = "Total",
                                       "Cumulative influenza infections" = sum(cumulative_incidence$`Cumulative influenza infections`),
                                       "N" = sum(cumulative_incidence$N)) %>%
      mutate(`Proportion infected` = `Cumulative influenza infections` / N) %>%
      mutate(`Symptomatic` = sum(cumulative_incidence$Symptomatic)) %>%
      mutate("Symptomatic and reported" = sum(cumulative_incidence$Symptomatic * reporting_rate))
    cumulative_incidence <- bind_rows(cumulative_incidence, cumulative_incidence_all)
    
    contact_matrices <- list(term = C_term, half_term = C_holidays,
                             christmas_shopping = C_christmas_break, christmas_holidays = C_xmas_holidays)
    
    meta <- list(start_date = start_date,
                 end_date = end_date,
                 plot_start_date = input$plot_start_date,
                 half_term_start = half_term_start,
                 half_term_end = half_term_end,
                 shopping_period_start = shopping_period_start,
                 winter_holiday_start = winter_holiday_start,
                 winter_holiday_end = winter_holiday_end,
                 reporting_rate = reporting_rate,
                 y_lim_max = input$y_lim_max,
                 prop_immune = prop_immune) # include immunity info for annotation
    
    list(inc = ret_sum,
         date_key = date_key,
         cumulative_incidence = cumulative_incidence,
         contact_matrices = contact_matrices,
         meta = meta)
  }) # end run_model
  
  scenario_snapshot <- reactiveVal(NULL)
  
  # Update ONLY when button is clicked
  observeEvent(input$set_ref, {
    scenario_snapshot(run_model()) 
  })
  
  
  output$status <- renderText({
    if (input$run_model == 0) {
      "Click 'Run Model' to begin."
    } else if (model_running()) {
      "Model is running..."
    } else {
      "Model finished."
    }
  })
  
  # Reset defaults
  observeEvent(input$reset_defaults, {
    updateNumericInput(session, "y_lim_max", value = defaults$y_lim_max)
    updateNumericInput(session, "alpha2", value = defaults$alpha2)
    updateNumericInput(session, "R0", value = defaults$R0)
    updateNumericInput(session, "gamma", value = defaults$gamma)
    updateNumericInput(session, "seed_size", value = defaults$seed_size)
    updateSliderInput(session, "sim_start_date", value = defaults$sim_start_date)
    updateSliderInput(session, "plot_start_date", value = defaults$plot_start_date)
    updateSliderInput(session, "sim_end_date", value = defaults$sim_end_date)
    updateSliderInput(session, "reporting_rate", value = defaults$reporting_rate)
    updateSliderInput(session, "prop_home_contacts_in_hols", value = defaults$prop_home_contacts_in_hols)
    updateSliderInput(session, "prop_work_contacts_in_hols", value = defaults$prop_work_contacts_in_hols)
    updateSliderInput(session, "prop_rest_contacts_in_hols", value = defaults$prop_rest_contacts_in_hols)
    updateSliderInput(session, "prop_all_contacts_christmas", value = defaults$prop_all_contacts_christmas)
    updateSliderInput(session, "prop_rest_contacts_in_christmas", value = defaults$prop_rest_contacts_in_christmas)
    updateSliderInput(session, "prop_home_contacts_in_christmas", value = defaults$prop_home_contacts_in_christmas)
    updateSliderInput(session, "overall_immune_escape", value = defaults$overall_immune_escape)
    updateNumericInput(session, "growth_y_scale", value = defaults$growth_y_scale)
    updateSliderInput(session, "prop_immune_younger", value = defaults$prop_immune_younger)
    updateSliderInput(session, "prop_immune_younger2", value = defaults$prop_immune_younger2)
    updateSliderInput(session, "prop_immune_older", value = defaults$prop_immune_older)
    updateSliderInput(session, "prop_immune_oldest", value = defaults$prop_immune_oldest)
    updateSliderInput(session, "symp_1", value = defaults$symp_1)
    updateSliderInput(session, "symp_2", value = defaults$symp_2)
    updateSliderInput(session, "symp_3", value = defaults$symp_3)
    updateSliderInput(session, "symp_4", value = defaults$symp_4)
  })
  
  # Combined plot render
  output$combined_plot <- renderPlot({
    res <- run_model()
    if (isTRUE(res$error)) {
      plot.new(); text(0.5, 0.5, res$message); return()
    }
    
    last_scenario <- scenario_snapshot()
    combined <- build_combined_plot(res, 
                                    input, 
                                    flu_dat(),
                                    last_scenario)
    print(combined)
  })
  
  # Contact matrices
  output$contact_matrices_plot <- renderPlot({
    res <- run_model()
    if (isTRUE(res$error)) {
      plot.new(); text(0.5, 0.5, res$message); return()
    }
    p <- build_contact_matrices_plot(res$contact_matrices)
    print(p)
  })
  
  # Table (DT)
  output$cum_table <- DT::renderDataTable({
    res <- run_model()
    if (isTRUE(res$error)) {
      return(datatable(data.frame(Message = res$message), rownames = FALSE))
    }
    df <- res$cumulative_incidence
    numcols <- sapply(df, is.numeric)
    df[numcols] <- lapply(df[numcols], function(x) signif(x, 3))
    dat <- datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = 't', ordering = FALSE)) %>%
      formatStyle('Age group', target = 'row',
                  backgroundColor = styleEqual(c("65+", "Total"), c("#ffe0e0", "#e0ffe0")))
    dat
  })
  
  
  
  # SERVER: single combined download handler
  output$download_all <- downloadHandler(
    filename = function() {
      paste0("model_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # run model
      res <- run_model()
      if (isTRUE(res$error)) stop("Model error: ", res$message)
      
      # build combined plot (dynamics + growth)
      # assumes you have build_combined_plot defined (uses build_incidence_plot and build_growth_plot)
      p_combined <- build_combined_plot(res, input, flu_dat())
      
      # prepare data tables to include
      # cumulative incidence table (from model)
      cum_df <- res$cumulative_incidence
      
      # incidence long (weekly) used in plots
      inc_df <- as.data.frame(res$inc)
      inc_df$t <- seq_len(nrow(inc_df))
      inc_long <- inc_df %>%
        pivot_longer(-t, names_to = "age_col", values_to = "weekly_reported") %>%
        mutate(date = res$date_key$date[t]) %>%
        select(date, t, age_col, weekly_reported)
      
      # compute growth rates (same logic as build_growth_plot)
      age_group_key <- c("inc_1_1"="[0,5)","inc_2_1"="[5,18)","inc_3_1"="[18,65)","inc_4_1"="65+")
      growth_df <- inc_long %>%
        mutate(age_group = age_group_key[age_col]) %>%
        group_by(age_group) %>%
        arrange(date) %>%
        mutate(lag_weekly = lag(weekly_reported),
               log_growth = log((weekly_reported + 1) / (coalesce(lag_weekly, 0) + 1))) %>%
        ungroup() %>%
        select(date, age_col, age_group, t, weekly_reported, lag_weekly, log_growth)
      
      # filenames to appear at top level of zip
      png_name   <- "model_combined_plot.png"
      rdata_name <- "model_results.RData"
      cum_name   <- "cumulative_incidence.csv"
      inc_name   <- "weekly_incidence_long.csv"
      growth_name<- "weekly_log_growth_by_age.csv"
      
      # write into a short-named tempdir so zip has top-level files
      tmpdir <- tempfile("model_results_zip")
      dir.create(tmpdir)
      oldwd <- setwd(tmpdir)
      on.exit({
        setwd(oldwd)
        unlink(tmpdir, recursive = TRUE, force = TRUE)
      }, add = TRUE)
      
      # save PNG of combined plot
      # pick size so both panels render clearly
      ggsave(filename = png_name, plot = p_combined, width = 12, height = 9, dpi = 300)
      
      # save RData with objects useful for later inspection
      # save combined plot, component plots (if available), the model result 'res' and tables
      # if you want p1/p2 separately, rebuild them quickly:
      p1 <- build_incidence_plot(res, 
                                 input,
                                 flu_dat())
      p2 <- build_growth_plot(res, input, age_palette = "Set1")
      save(p_combined, p1, p2, res, cum_df, inc_long, growth_df, file = rdata_name)
      
      # write CSVs
      write.csv(cum_df, cum_name, row.names = FALSE)
      write.csv(inc_long, inc_name, row.names = FALSE)
      write.csv(growth_df, growth_name, row.names = FALSE)
      
      utils::zip(zipfile = file, files = c(png_name, rdata_name, cum_name, inc_name, growth_name))
      
      # tmpdir and files removed by on.exit
    }
  )
  
} # end server

# --------------------------
# Run app
# --------------------------
shinyApp(ui = ui, server = server)
