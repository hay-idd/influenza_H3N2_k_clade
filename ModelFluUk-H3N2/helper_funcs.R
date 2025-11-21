# Helper functions

contacts_all <- polymod$contacts
polymod_base <- polymod

build_contact_matrices <- function(input) {
  polymod_c_term <- contact_matrix(polymod_base,
                                   countries = "United Kingdom",
                                   age.limits = age_breaks,
                                   symmetric = TRUE,
                                   missing.contact.age = "sample",
                                   missing.participant.age = "remove")
  C_term <- polymod_c_term$matrix
  row.names(C_term) <- colnames(C_term)
  
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

build_incidence_plot <- function(res, input, flu_dat, last_scenario) {
  inc <- as.data.frame(res$inc)
  inc$t <- 1:nrow(inc)
  inc <- inc %>% pivot_longer(-t)
  colnames(inc) <- c("t","age_group","incidence")
  age_group_key <- c("inc_1_1"="[0,5)","inc_2_1"="[5,18)","inc_3_1"="[18,65)","inc_4_1"="65+")
  inc$age_group <- age_group_key[inc$age_group]
  inc$age_group <- factor(inc$age_group, levels = age_group_key)
  date_key <- res$date_key
  
  if (!is.null(last_scenario)){
    
    flu_subset_dat <- as.data.frame(last_scenario$inc)
    flu_subset_dat$t <- 1:nrow(flu_subset_dat)
    flu_subset_dat <- flu_subset_dat %>% pivot_longer(-t)
    colnames(flu_subset_dat) <- c("date","age_group","ILI_flu")
    flu_subset_dat$age_group <- age_group_key[flu_subset_dat$age_group]
    flu_subset_dat$age_group <- factor(flu_subset_dat$age_group, levels = age_group_key)
    
  } else{
    
    age_group_key1 <- c("0-4"="[0,5)","5-18"="[5,18)","19-64"="[18,65)","65+"="65+")
    flu_subset_dat <- flu_dat %>% filter(date >= res$meta$start_date, date <= res$meta$end_date)
    flu_subset_dat$age_group <- age_group_key1[flu_subset_dat$age_group]
    flu_subset_dat$age_group <- factor(flu_subset_dat$age_group, levels = age_group_key1)

  }
  
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
    
    geom_line(data = flu_subset_dat, 
              aes(x = date, y = ILI_flu, col = age_group, linetype = "2022/23 season data"),
              alpha = 0.5, linewidth = 1) +
    
    geom_line(aes(x = date, y = incidence, col = age_group, linetype = "Model"), linewidth = 1) +
    scale_color_brewer("Age group", palette = "Set1") +
    scale_fill_brewer("Holiday period", palette = "Set2") +
    scale_linetype_manual("Data source", values = c("2022/23 season data" = "dashed", "Model" = "solid")) +
    scale_x_date(breaks = "1 month", date_labels = "%b-%d", expand = c(0, 0)) +
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
  
  # Compute peak & cumulative labels (unchanged)
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
    "Date of peak: ", if (is.na(peak_date)) "n/a" else format(peak_date, "%b-%d"), "\n",
    "Cumulative symptomatic: ", fmt_num(cum_symp_total)
  )
  x_pos_peak <- as.Date(res$meta$end_date) - 3
  y_pos_peak <- 0.98 * res$meta$y_lim_max
  
  p <- p + annotate("label",
                    x = x_pos_peak, y = y_pos_peak,
                    label = label_text,
                    hjust = 1, vjust = 1,
                    size = 3.5, fontface = "bold",
                    fill = "white", alpha = 0.8)
  
  # --- New: proportion immune annotation (stacked, left-aligned, same style as peak label) ---
  # compute the effective immune proportions (same logic as in server)
  overall_mult <- if (!is.null(input$overall_immune_escape)) input$overall_immune_escape else 1
  raw_prop_immune <- c(
    if (!is.null(input$prop_immune_younger)) input$prop_immune_younger else NA,
    if (!is.null(input$prop_immune_younger2)) input$prop_immune_younger2 else NA,
    if (!is.null(input$prop_immune_older)) input$prop_immune_older else NA,
    if (!is.null(input$prop_immune_oldest)) input$prop_immune_oldest else NA
  )
  prop_immune <- pmin(1, overall_mult * raw_prop_immune)
  # friendly percent formatting
  pct <- function(x) if (is.na(x)) "n/a" else paste0(format(round(100 * x, 1), trim = TRUE), "%")
  immune_lines <- paste0(
    "Proportion immune:\n",
    "0-4: ", pct(prop_immune[1]), "\n",
    "5-18: ", pct(prop_immune[2]), "\n",
    "19-64: ", pct(prop_immune[3]), "\n",
    "65+: ", pct(prop_immune[4])
  )
  # position immune label on the left, vertically near the top, and use same font/size
  x_pos_immune <- as.Date(res$meta$plot_start_date) + 3
  y_pos_immune <- 0.98 * res$meta$y_lim_max
  
  p <- p + annotate("label",
                    x = x_pos_immune, y = y_pos_immune,
                    label = immune_lines,
                    hjust = 0, vjust = 1,
                    size = 3.5, fontface = "bold",
                    fill = "white", alpha = 0.8, lineheight = 0.95)
  
  p
}


build_growth_plot <- function(res, input, age_palette = "Set1") {
  inc_df <- as.data.frame(res$inc)
  inc_df$t <- 1:nrow(inc_df)
  inc_long <- inc_df %>% pivot_longer(-t, names_to = "age_col", values_to = "incidence")
  age_group_key <- c("inc_1_1"="[0,5)", "inc_2_1"="[5,18)", "inc_3_1"="[18,65)", "inc_4_1"="65+")
  inc_long$age_group <- age_group_key[inc_long$age_col]
  inc_long$age_group <- factor(inc_long$age_group, levels = age_group_key)
  date_key <- res$date_key
  inc_long <- inc_long %>% left_join(date_key, by = "t")
  inc_long <- inc_long %>%
    arrange(age_group, t) %>%
    group_by(age_group) %>%
    mutate(lag_inc = lag(incidence),
           log_growth = log((incidence + 1) / (coalesce(lag_inc, 0) + 1))) %>%
    ungroup()
  
  rects <- data.frame(
    xmin = as.Date(c(half_term_start, shopping_period_start, winter_holiday_start)),
    xmax = as.Date(c(half_term_end,   shopping_period_end,   winter_holiday_end)),
    label = c("Half term", "Christmas period\n (pre holidays)", "Christmas\nholidays"),
    fill  = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
    stringsAsFactors = FALSE
  )
  rects$mid <- as.Date((as.numeric(rects$xmin) + as.numeric(rects$xmax)) / 2, origin = "1970-01-01")
  
  # set y-limits from user input growth_y_scale (symmetric about zero)
  gscale <- if (!is.null(input$growth_y_scale)) input$growth_y_scale else defaults$growth_y_scale
  ylim_val <- gscale
  
  p <- ggplot(inc_long, aes(x = date, y = log_growth, colour = age_group, group = age_group)) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.12, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(size = 1) +
    scale_color_brewer("Age group", palette = age_palette) +
    scale_fill_brewer("Holiday period", palette = "Set2") +
    scale_x_date(breaks = "1 month", date_labels = "%b-%d", expand = c(0, 0)) +
    coord_cartesian(ylim = c(-ylim_val, ylim_val), xlim = c(res$meta$plot_start_date, res$meta$end_date)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),          # match incidence text sizes
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(6, 10, 6, 10)) +
    xlab(NULL) +
    ylab("Weekly growth rate")
  
  # return both the plot and the prepared growth data (useful for download)
  attr(p, "growth_data") <- inc_long
  p
}

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

# Combined plot builder (returns a single patchwork object)
build_combined_plot <- function(res, input, flu_dat, last_scenario) {
  p1 <- build_incidence_plot(res, input, flu_dat, last_scenario)
  p2 <- build_growth_plot(res, input, age_palette = "Set1")
  
  # remove x-axis from top plot and shrink its bottom margin; nudge top margin of bottom plot
  p1_clean <- p1 + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(6, 10, 2, 10)
  )
  p2_clean <- p2 + theme(plot.margin = margin(2, 10, 6, 10))
  
  combined <- (p1_clean / p2_clean) +
    plot_layout(ncol = 1, heights = c(3, 1.9), guides = "collect") +
    plot_annotation(
      title = "Reported symptomatic cases (weekly) and log weekly growth",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  # apply shared legend placement and stacking + use incidence plot's text settings
  combined <- combined & theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
  
  combined
}

save_plot_and_data_zip <- function(plot_obj, plot_data, growth_data, flu_subset_dat, res, zip_path) {
  tmpdir <- tempfile("model_plot_zip")
  dir.create(tmpdir)
  oldwd <- setwd(tmpdir)
  on.exit({
    setwd(oldwd)
    unlink(tmpdir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  png_name <- "model_combined_plot.png"
  rdata_name <- "model_plot_and_data.RData"
  
  # Save higher-resolution PNG to accommodate combined plot
  ggsave(filename = png_name, plot = plot_obj, width = 14, height = 9, dpi = 300)
  
  # Save R objects: include combined plot object, plot_data (incidence), growth_data, flu_subset_dat and res
  save(plot = plot_obj, plot_data, growth_data, flu_subset_dat, res, file = rdata_name)
  
  utils::zip(zipfile = zip_path, files = c(png_name, rdata_name))
}

