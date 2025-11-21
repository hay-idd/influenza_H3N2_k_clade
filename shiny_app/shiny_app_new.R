# app.R (refactored)
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

# working dir (adjust if needed)
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/shiny_app/")

# --------------------------
# Global constants & data
# --------------------------
half_term_start <- as.Date("2022-10-21")
half_term_end   <- as.Date("2022-10-30")
shopping_period_start <- as.Date("2022-12-01")
shopping_period_end   <- as.Date("2022-12-15")
winter_holiday_start  <- as.Date("2022-12-21")
winter_holiday_end    <- as.Date("2023-01-04")
N_tot <- 60000000
age_breaks <- c(0,5,18,65)

data("polymod")
contacts_all <- polymod$contacts
polymod_base <- polymod

flu_dat <- read_csv("data/final/flu_2022_2023.csv")

# --------------------------
# Defaults (single source)
# --------------------------
defaults <- list(
  y_lim_max = 4000,
  alpha2 = 0.65,
  R0 = 2,
  gamma = 4.6,
  seed_size = 150,
  initial_immune_frac = 0,
  sim_start_date = as.Date("2022-09-01"),
  plot_start_date = as.Date("2022-10-01"),
  sim_end_date = as.Date("2023-03-01"),
  reporting_rate = 0.0012,
  prop_home_contacts_in_hols = 1,
  prop_work_contacts_in_hols = 0.75,
  prop_rest_contacts_in_hols = 1.1,
  prop_all_contacts_christmas = 1.3,
  prop_rest_contacts_in_christmas = 0.67,
  prop_home_contacts_in_christmas = 3,
  prop_immune_younger = 0.3,
  prop_immune_younger2 = 0.4,
  prop_immune_older = 0.8,
  prop_immune_oldest = 0.8,
  symp_1 = 0.25,
  symp_2 = 0.4,
  symp_3 = 0.5,
  symp_4 = 0.5
)

# --------------------------
# Helper functions
# --------------------------

# Build contact matrices (returns list of four raw matrices and the polymod term object)
build_contact_matrices <- function(input) {
  # term
  polymod_c_term <- contact_matrix(polymod_base,
                                   countries = "United Kingdom",
                                   age.limits = age_breaks,
                                   symmetric = TRUE,
                                   missing.contact.age = "sample",
                                   missing.participant.age = "remove")
  C_term <- polymod_c_term$matrix
  row.names(C_term) <- colnames(C_term)
  
  # term-break (no school contacts; adjust home/work/other)
  prop_home_contacts_in_hols <- input$prop_home_contacts_in_hols
  prop_work_contacts_in_hols <- input$prop_work_contacts_in_hols
  prop_rest_contacts_in_hols <- input$prop_rest_contacts_in_hols
  
  contacts <- contacts_all
  contacts_others <- contacts %>% filter(cnt_school == 0, cnt_work == 0, cnt_home == 0) %>% sample_frac(prop_rest_contacts_in_hols, replace = TRUE)
  contacts_work <- contacts %>% filter(cnt_work == 1) %>% sample_frac(prop_work_contacts_in_hols)
  contacts_home <- contacts %>% filter(cnt_home == 1) %>% sample_frac(prop_home_contacts_in_hols, replace = TRUE)
  contacts_overall <- bind_rows(contacts_others, contacts_work, contacts_home)
  polymodTermBreak <- polymod_base
  polymodTermBreak$contacts <- contacts_overall
  polymod_c_holidays <- contact_matrix(polymodTermBreak,
                                       countries = "United Kingdom",
                                       age.limits = age_breaks,
                                       symmetric = TRUE,
                                       missing.contact.age = "sample",
                                       missing.participant.age = "remove")
  C_holidays <- polymod_c_holidays$matrix
  row.names(C_holidays) <- colnames(C_holidays)
  
  # christmas period (pre-holidays): amplify non-school contacts
  prop_all_contacts_christmas <- input$prop_all_contacts_christmas
  contacts <- contacts_all
  contacts_schools <- contacts %>% filter(cnt_school == 1)
  contacts_non_school <- contacts %>% filter(cnt_school == 0) %>% sample_frac(prop_all_contacts_christmas, replace = TRUE)
  contacts_overall2 <- bind_rows(contacts_non_school, contacts_schools)
  polymodChristmasPeriod <- polymod_base
  polymodChristmasPeriod$contacts <- contacts_overall2
  polymod_c_christmas_break <- contact_matrix(polymodChristmasPeriod,
                                              countries = "United Kingdom",
                                              age.limits = age_breaks,
                                              symmetric = TRUE,
                                              missing.contact.age = "sample",
                                              missing.participant.age = "remove")
  C_christmas_break <- polymod_c_christmas_break$matrix
  row.names(C_christmas_break) <- colnames(C_christmas_break)
  
  # christmas holidays: no school, adjust rest/home
  prop_rest_contacts_in_christmas <- input$prop_rest_contacts_in_christmas
  prop_home_contacts_in_christmas <- input$prop_home_contacts_in_christmas
  contacts <- contacts_all
  contacts_non_school2 <- contacts %>% filter(cnt_school == 0) %>% sample_frac(prop_rest_contacts_in_christmas, replace = TRUE)
  contacts_home_ch <- contacts_non_school2 %>% filter(cnt_home == 1) %>% sample_frac(prop_home_contacts_in_christmas, replace = TRUE)
  contacts_overall_ch <- bind_rows(contacts_non_school2, contacts_home_ch)
  polymodChristmasHoliday <- polymod_base
  polymodChristmasHoliday$contacts <- contacts_overall_ch
  polymod_c_xmas_holidays <- contact_matrix(polymodChristmasHoliday,
                                            countries = "United Kingdom",
                                            age.limits = age_breaks,
                                            symmetric = TRUE,
                                            missing.contact.age = "sample",
                                            missing.participant.age = "remove")
  C_xmas_holidays <- polymod_c_xmas_holidays$matrix
  row.names(C_xmas_holidays) <- colnames(C_xmas_holidays)
  
  list(term = C_term,
       term_break = C_holidays,
       christmas_period = C_christmas_break,
       christmas_holiday = C_xmas_holidays,
       polymod_term = polymod_c_term)
}

# Create normalized weights dataframe (combines three period weightings into term + holidayness)
make_weights_df <- function(start_date, end_date, smooth_time = 7) {
  school_days_weighted <- setup_holiday_tibble(start_date, end_date,
                                               half_term_start, half_term_end,
                                               shopping_period_start, winter_holiday_end,
                                               smooth_time = smooth_time)
  half_term_weighting <- setup_period_tibble(start_date, end_date,
                                             half_term_start, half_term_end,
                                             smooth_time = smooth_time)
  christmas_holiday_weighting <- setup_period_tibble(start_date, end_date,
                                                     winter_holiday_start, winter_holiday_end,
                                                     smooth_time = smooth_time)
  christmas_shopping_period <- setup_period_tibble(start_date, end_date,
                                                   shopping_period_start, shopping_period_end,
                                                   smooth_time = smooth_time)
  ref_dates <- school_days_weighted$date
  weights_df <- tibble::tibble(date = ref_dates) %>%
    dplyr::left_join(half_term_weighting %>% dplyr::select(date, weight) %>% dplyr::rename(w_termBreak = weight), by = "date") %>%
    dplyr::left_join(christmas_shopping_period %>% dplyr::select(date, weight) %>% dplyr::rename(w_ChristmasPeriod = weight), by = "date") %>%
    dplyr::left_join(christmas_holiday_weighting %>% dplyr::select(date, weight) %>% dplyr::rename(w_ChristmasBreak = weight), by = "date") %>%
    dplyr::mutate(
      h_termBreak       = 1 - w_termBreak,
      h_ChristmasPeriod = 1 - w_ChristmasPeriod,
      h_ChristmasBreak  = 1 - w_ChristmasBreak,
      sum_h = h_termBreak + h_ChristmasPeriod + h_ChristmasBreak,
      norm_factor = ifelse(sum_h > 1, 1 / sum_h, 1),
      h_termBreak = h_termBreak * norm_factor,
      h_ChristmasPeriod = h_ChristmasPeriod * norm_factor,
      h_ChristmasBreak = h_ChristmasBreak * norm_factor,
      w_term = 1 - (h_termBreak + h_ChristmasPeriod + h_ChristmasBreak)
    ) %>%
    dplyr::select(date, w_term, h_termBreak, h_ChristmasPeriod, h_ChristmasBreak)
  list(schedule = school_days_weighted, weights = weights_df)
}

# Build ggplot object for incidence (used by UI and download)
build_incidence_plot <- function(res, input, flu_dat) {
  inc <- as.data.frame(res$inc); inc$t <- 1:nrow(inc)
  inc <- inc %>% pivot_longer(-t)
  colnames(inc) <- c("t","age_group","incidence")
  age_group_key <- c("inc_1_1"="[0,5)","inc_2_1"="[5,18)","inc_3_1"="[18,65)","inc_4_1"="65+")
  inc$age_group <- age_group_key[inc$age_group]
  inc$age_group <- factor(inc$age_group, levels = age_group_key)
  date_key <- res$date_key
  
  age_group_key1 <- c("0-4"="[0,5)","5-18"="[5,18)","19-64"="[18,65)","65+"="65+")
  flu_dat1 <- flu_dat %>% filter(date >= res$meta$start_date, date <= res$meta$end_date)
  flu_dat1$age_group <- age_group_key1[flu_dat1$age_group]
  flu_dat1$age_group <- factor(flu_dat1$age_group, levels = age_group_key1)
  
  rects <- data.frame(
    xmin = as.Date(c(half_term_start, shopping_period_start, winter_holiday_start)),
    xmax = as.Date(c(half_term_end,   shopping_period_end,   winter_holiday_end)),
    label = c("Half term", "Christmas period\n (pre holidays)", "Christmas\nholidays"),
    fill  = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
    stringsAsFactors = FALSE
  )
  rects$mid <- as.Date((as.numeric(rects$xmin) + as.numeric(rects$xmax)) / 2, origin = "1970-01-01")
  y_pos <- res$meta$y_lim_max * 0.95
  rects$y <- y_pos
  
  p <- ggplot(inc %>% left_join(date_key)) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.25, show.legend = TRUE) +
    geom_text(data = rects, aes(x = mid, y = y, label = label), inherit.aes = FALSE,
              size = 3, fontface = "bold", vjust = 1) +
    geom_line(data = flu_dat1, aes(x = date, y = ILI_flu, col = age_group, linetype = "2023/23 season data"),
              alpha = 0.5, linewidth = 1) +
    geom_line(aes(x = date, y = incidence, col = age_group, linetype = "Model"), linewidth = 1) +
    scale_color_brewer("Age group", palette = "Set1") +
    scale_fill_brewer("Holiday period", palette = "Set2") +
    scale_linetype_manual("Data source", values = c("2023/23 season data" = "dashed", "Model" = "solid")) +
    scale_x_date(breaks = "1 month", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, res$meta$y_lim_max), xlim = c(res$meta$plot_start_date, res$meta$end_date)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) +
    xlab("Date (end of reporting week)") +
    ylab("Reported influenza cases (weekly)") +
    labs(color = "Age group")
  
  # Compute peak & cumulative labels
  weekly_totals <- rowSums(as.matrix(res$inc), na.rm = TRUE)
  peak_val <- if (length(weekly_totals)) max(weekly_totals, na.rm = TRUE) else NA
  peak_idx <- if (length(weekly_totals)) which.max(weekly_totals) else NA
  peak_date <- if (!is.na(peak_idx) && length(peak_idx)) res$date_key$date[peak_idx] else NA
  reported_peak <- peak_val
  peak_symptomatic <- if (!is.na(reported_peak) && res$meta$reporting_rate > 0) reported_peak / res$meta$reporting_rate else NA
  
  cum_symp_total <- NA
  if (!is.null(res$cumulative_incidence)) {
    if ("Age group" %in% colnames(res$cumulative_incidence)) {
      row_total <- res$cumulative_incidence %>% filter(`Age group` == "Total")
      if (nrow(row_total) == 1 && "Symptomatic" %in% colnames(row_total)) {
        cum_symp_total <- row_total$Symptomatic
      } else if ("Symptomatic" %in% colnames(res$cumulative_incidence)) {
        cum_symp_total <- sum(res$cumulative_incidence$Symptomatic, na.rm = TRUE)
      }
    }
  }
  fmt_num <- function(x) if (is.na(x)) "n/a" else format(round(x), big.mark = ",", scientific = FALSE)
  label_text <- paste0(
    "Peak (weekly symptomatic): ", fmt_num(peak_symptomatic), "\n",
    "Peak (weekly reported): ", fmt_num(peak_val), "\n",
    "Date of peak: ", if (is.na(peak_date)) "n/a" else format(peak_date, "%Y-%m-%d"), "\n",
    "Cumulative symptomatic: ", fmt_num(cum_symp_total)
  )
  x_pos <- as.Date(res$meta$end_date) - 3
  y_pos <- 0.98 * res$meta$y_lim_max
  p <- p + annotate("label",
                    x = x_pos, y = y_pos,
                    label = label_text,
                    hjust = 1, vjust = 1,
                    size = 3.5, fontface = "bold",
                    fill = "white", alpha = 0.8)
  p
}

# Contact matrices plot builder
build_contact_matrices_plot <- function(contact_matrices) {
  df_list <- lapply(names(contact_matrices), function(nm) {
    m <- contact_matrices[[nm]]
    rns <- rownames(m); cns <- colnames(m)
    if (is.null(rns)) rns <- paste0("r", seq_len(nrow(m)))
    if (is.null(cns)) cns <- paste0("c", seq_len(ncol(m)))
    d <- as.data.frame(as.table(m))
    colnames(d) <- c("row", "col", "value")
    d$matrix_name <- nm
    d$row <- factor(as.character(d$row), levels = unique(rns))
    d$col <- factor(as.character(d$col), levels = unique(cns))
    d
  })
  plot_df <- do.call(rbind, df_list)
  p <- ggplot(plot_df, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    facet_wrap(~ matrix_name, ncol = 2) +
    labs(x = NULL, y = NULL, fill = "Contact\nrate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    coord_fixed()
  p
}

# Save plot + data into a zip (png + RData at top level)
save_plot_and_data_zip <- function(plot_obj, plot_data, flu_dat1, res, zip_path) {
  tmpdir <- tempfile("model_plot_zip")
  dir.create(tmpdir)
  oldwd <- setwd(tmpdir)
  on.exit({
    setwd(oldwd)
    unlink(tmpdir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  png_name <- "model_plot.png"
  rdata_name <- "model_plot_data.RData"
  ggsave(filename = png_name, plot = plot_obj, width = 12, height = 6, dpi = 300)
  save(plot = plot_obj, plot_data, flu_dat1, res, file = rdata_name)
  utils::zip(zipfile = zip_path, files = c(png_name, rdata_name))
}

# --------------------------
# UI
# --------------------------
ui <- fluidPage(
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
                 downloadButton("download_plot", "Download Plot (zip)"),
                 downloadButton("download_table", "Download Table (csv)"),
                 hr(),
                 h4("Plot settings"),
                 numericInput("y_lim_max", "Maximum y-axis value", value = defaults$y_lim_max, step = 100),
                 sliderInput("plot_start_date", "Plot start date",
                             min = as.Date("2022-07-01"),
                             max = as.Date("2022-12-31"),
                             value = defaults$plot_start_date,
                             timeFormat = "%Y-%m-%d"),
                 hr(),
                 h4("Model parameters"),
                 numericInput("alpha2", "Relative susceptibility of immune class", value = defaults$alpha2, step = 0.05),
                 numericInput("R0", "Basic reproductive number R0", value = defaults$R0, step = 0.01),
                 numericInput("gamma", "Infectious period (days) gamma", value = defaults$gamma, step = 0.1),
                 numericInput("seed_size", "Seed size (number of initial infections)", value = defaults$seed_size, step = 10),
                 numericInput("initial_immune_frac", "Initial proportion of the population fully immune", value = defaults$initial_immune_frac, step = 0.01),
                 sliderInput("sim_start_date", "Seed date",
                             min = as.Date("2022-07-01"),
                             max = as.Date("2022-12-31"),
                             value = defaults$sim_start_date,
                             timeFormat = "%Y-%m-%d"),
                 sliderInput("sim_end_date", "Simulation end date",
                             min = as.Date("2023-01-01"),
                             max = as.Date("2023-07-01"),
                             value = defaults$sim_end_date,
                             timeFormat = "%Y-%m-%d"),
                 hr(),
                 sliderInput("reporting_rate", "Reporting rate (fraction of symptomatic reported)",
                             min = 0, max = 0.01, value = defaults$reporting_rate, step = 0.0001),
                 hr(),
                 h4("School holiday contacts"),
                 sliderInput("prop_home_contacts_in_hols", "Multiplier for home contacts in breaks",
                             min = 0, max = 3, value = defaults$prop_home_contacts_in_hols, step = 0.01),
                 sliderInput("prop_work_contacts_in_hols", "Proportion of work contacts kept in school holidays",
                             min = 0, max = 1, value = defaults$prop_work_contacts_in_hols, step = 0.01),
                 sliderInput("prop_rest_contacts_in_hols", "Multiplier for non-school or work contacts in breaks",
                             min = 0.5, max = 3, value = defaults$prop_rest_contacts_in_hols, step = 0.01),
                 hr(),
                 h4("Christmas period contacts (before school holiday)"),
                 sliderInput("prop_all_contacts_christmas", "Multiplier for all non-school contacts in Christmas period",
                             min = 0.5, max = 3, value = defaults$prop_all_contacts_christmas, step = 0.01),
                 hr(),
                 h4("Christmas holiday contacts"),
                 sliderInput("prop_rest_contacts_in_christmas", "Proportion of all contacts kept over Christmas",
                             min = 0, max = 1, value = defaults$prop_rest_contacts_in_christmas, step = 0.01),
                 sliderInput("prop_home_contacts_in_christmas", "Multiplier for home contacts over Christmas holiday",
                             min = 1, max = 5, value = defaults$prop_home_contacts_in_christmas, step = 0.01),
                 hr(),
                 h4("Immunity proportions by age-groups"),
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
                 sliderInput("symp_4", "Symptomatic frac (65+ yrs)", min = 0, max = 1, value = defaults$symp_4, step = 0.01),
                 width = 3
               ),
               mainPanel(
                 fluidRow(
                   column(12,
                          wellPanel(
                            h4("Reported symptomatic cases (weekly)"),
                            plotOutput("inc_plot", height = "600px"),
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
               column(8,
                      h3("About this model"),
                      p("Placeholder: description of the model and data goes here. (You can edit this later.)")
               ),
               column(4,
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

# --------------------------
# Server
# --------------------------
server <- function(input, output, session) {
  # reactive model (raw) then debounced
  run_model_raw <- reactive({
    # ensure helper functions from your files are available
    source("auxiliary_funcs.R")
    source("sir_functions.R")
    
    # local inputs
    start_date <- input$sim_start_date
    end_date   <- input$sim_end_date
    dates <- seq(start_date, end_date, by = "1 week")
    
    prop_immune_younger <- input$prop_immune_younger
    prop_immune_youngest <- input$prop_immune_younger2
    prop_immune_older <- input$prop_immune_older
    prop_immune_oldest <- input$prop_immune_oldest
    
    alphas <- c(1, input$alpha2)
    ## Set alphas to average 1
    alphas <- alphas/mean(alphas)
    print(alphas)
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
    
    # population & N
    prop_immune <- c(prop_immune_younger, prop_immune_youngest, prop_immune_older, prop_immune_oldest)
    N_props <- polymod_c_term$demography$proportion
    N_age_classes <- length(N_props)
    N_props_long <- c(N_props * (1 - prop_immune), N_props * (prop_immune))
    N <- matrix(N_props_long * N_tot, ncol = length(alphas), nrow = N_age_classes)
    
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
                           initial_immune_frac = input$initial_immune_frac,
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
                 y_lim_max = input$y_lim_max)
    
    list(inc = ret_sum,
         date_key = date_key,
         cumulative_incidence = cumulative_incidence,
         contact_matrices = contact_matrices,
         meta = meta)
  }) # end run_model_raw
  
  # Debounced reactive
  run_model <- shiny::debounce(run_model_raw, millis = 700)
  
  # Reset defaults
  observeEvent(input$reset_defaults, {
    # numeric inputs
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
    updateSliderInput(session, "prop_immune_younger", value = defaults$prop_immune_younger)
    updateSliderInput(session, "prop_immune_younger2", value = defaults$prop_immune_younger2)
    updateSliderInput(session, "prop_immune_older", value = defaults$prop_immune_older)
    updateSliderInput(session, "prop_immune_oldest", value = defaults$prop_immune_oldest)
    updateSliderInput(session, "symp_1", value = defaults$symp_1)
    updateSliderInput(session, "symp_2", value = defaults$symp_2)
    updateSliderInput(session, "symp_3", value = defaults$symp_3)
    updateSliderInput(session, "symp_4", value = defaults$symp_4)
  })
  
  # --------------------------
  # Incidence plot render
  # --------------------------
  output$inc_plot <- renderPlot({
    res <- run_model()
    if (isTRUE(res$error)) {
      plot.new(); text(0.5, 0.5, res$message); return()
    }
    p <- build_incidence_plot(res, input, flu_dat)
    print(p)
  })
  
  # --------------------------
  # Table (DT) render
  # --------------------------
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
  
  # --------------------------
  # Contact matrices render
  # --------------------------
  output$contact_matrices_plot <- renderPlot({
    res <- run_model()
    if (isTRUE(res$error)) {
      plot.new(); text(0.5, 0.5, res$message); return()
    }
    p <- build_contact_matrices_plot(res$contact_matrices)
    print(p)
  })
  
  # --------------------------
  # Download (plot + data) -> zip
  # --------------------------
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("model_plot_and_data_", Sys.Date(), ".zip")
    },
    content = function(file) {
      res <- run_model()
      if (isTRUE(res$error)) stop("Model error: ", res$message)
      p <- build_incidence_plot(res, input, flu_dat)
      plot_data <- (as.data.frame(res$inc) %>% mutate(t = 1:nrow(res$inc)) %>% pivot_longer(-t))
      # small wrapper that filters/renames to match earlier saved objects
      flu_dat1 <- flu_dat %>% filter(date >= res$meta$start_date, date <= res$meta$end_date)
      save_plot_and_data_zip(p, plot_data, flu_dat1, res, zip_path = file)
    }
  )
  
  # --------------------------
  # Download table CSV
  # --------------------------
  output$download_table <- downloadHandler(
    filename = function() paste0("cumulative_incidence_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- run_model()
      if (isTRUE(res$error)) stop("Model error: ", res$message)
      write.csv(res$cumulative_incidence, file, row.names = FALSE)
    }
  )
  
} # end server

# --------------------------
# Run app
# --------------------------
shinyApp(ui = ui, server = server)
