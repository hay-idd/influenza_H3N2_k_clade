epiweek_to_date_formula <- function(year, week,start=TRUE) {
  jan4 <- as.Date(paste0(year, "-01-04"))
  week1_monday <- jan4 - (lubridate::wday(jan4, week_start = 1) - 1)
  start_date <- week1_monday + (week - 1) * 7
  end_date <- start_date + 6
  if(start) return(start_date)
  else return(end_date)
}

extract_age_range <- function(df){
 df %>%
    mutate(
      age_start = case_when(
        str_detect(name, "^<\\d+") ~ 0,
        TRUE ~ as.numeric(str_extract(name, "^\\d+"))
      ),
      age_end = case_when(
        str_detect(name, "\\+") ~ Inf,
        str_detect(name, "^<\\d+") ~ as.numeric(str_extract(name, "\\d+")) - 1,
        TRUE ~ as.numeric(str_extract(name, "(?<=-)\\d+"))
      )
    )
}

get_age_group_n <- function(age_groups, age_min, age_max){
  if(age_max==Inf){
    n <- age_groups %>% filter(age_group >= age_min) %>% summarise(N=sum(N)) %>% pull(N)
  } else {
    n <- age_groups %>% filter(age_group >= age_min & age_group <= age_max) %>% summarise(N=sum(N)) %>% pull(N)
  }
  return(n)
}

expand_age_group_ili <- function(df,age_cap=90L){
  df %>%
    # convert Inf upper bounds to the cap we want (e.g. 90)
    mutate(age_end2 = ifelse(is.infinite(age_end), age_cap, age_end)) %>%
    # create a list-column of ages covered by each age-group row
    mutate(ages = map2(age_start, age_end2, ~ seq(from = .x, to = .y, by = 1))) %>%
    # expand into one row per age
    unnest(ages) %>%
    rename(age = ages) %>%
    # keep only ages within 0:age_cap (defensive)
    filter(between(age, 0, age_cap)) %>%
    # optional: drop helper column
    select(-age_end2)
}


combine_age_groups_ILI <- function(data, groups) {
  # groups: e.g., c("0-14", "15-64", "65+")
  map_dfr(groups, function(g) {
    start <- as.numeric(str_extract(g, "^\\d+"))
    end   <- ifelse(str_detect(g, "\\+"), Inf, as.numeric(str_extract(g, "(?<=-)\\d+")))
    
    data %>%
      filter(age >= start, age <= end) %>%
      summarise(
        group = g,
        ILI = sum(ILI, na.rm = TRUE),
        N = sum(N, na.rm = TRUE),
        ILI_per_100k = 1e5 * ILI / N,
        .by = c(Year, Week, date)
      )
  })
}

combine_age_groups_pos <- function(data, groups) {
  # groups: e.g., c("0-14", "15-64", "65+")
  map_dfr(groups, function(g) {
    start <- as.numeric(str_extract(g, "^\\d+"))
    end   <- ifelse(str_detect(g, "\\+"), Inf, as.numeric(str_extract(g, "(?<=-)\\d+")))
    
    data %>%
      filter(age >= start, age <= end) %>%
      summarise(
        group = g,
        # weighted mean positivity
        positivity = weighted.mean(positivity, N, na.rm = TRUE),
        N = sum(N, na.rm = TRUE),
        .by = c(Year, Week, date)
      )
  })
}


parse_age_groups_pos <- function(groups) {
  groups %>%
    mutate(
      age_start = case_when(
        str_detect(age_group, regex("up to", ignore_case = TRUE)) ~ 0,
        str_detect(age_group, regex("and above", ignore_case = TRUE)) ~ as.numeric(str_extract(age_group, "\\d+")),
        TRUE ~ as.numeric(str_extract(age_group, "^\\d+"))
      ),
      age_end = case_when(
        str_detect(age_group, regex("up to", ignore_case = TRUE)) ~ as.numeric(str_extract(age_group, "\\d+")),
        str_detect(age_group, regex("and above", ignore_case = TRUE)) ~ Inf,
        str_detect(age_group, regex("and over", ignore_case = TRUE)) ~ Inf,
        TRUE ~ as.numeric(str_extract(age_group, "(?<=to )\\d+"))
      )
    )
}

expand_age_group_pos <- function(df,age_cap=90L){
  df %>%
    # convert Inf upper bounds to the cap we want (e.g. 90)
    mutate(age_end2 = ifelse(is.infinite(age_end), age_cap, age_end)) %>%
    # create a list-column of ages covered by each age-group row
    mutate(ages = map2(age_start, age_end2, ~ seq(from = .x, to = .y, by = 1))) %>%
    # expand into one row per age
    unnest(ages) %>%
    rename(age = ages) %>%
    # keep only ages within 0:age_cap (defensive)
    filter(between(age, 0, age_cap)) %>%
    # optional: drop helper column
    select(-age_end2)
}



#' Fit a smoothed GAM and estimate weekly growth rates
#'
#' Fits a penalized spline GAM to log-transformed case counts over time,
#' computes fitted values with 95% confidence intervals, and derives
#' weekly growth rates and their approximate bounds.
#'
#' @param dat A data frame containing at least:
#'   \describe{
#'     \item{y}{Log-transformed counts (numeric).}
#'     \item{time_w}{Time in weeks (numeric).}
#'     \item{date}{Calendar date (Date).}
#'   }
#' @param epsilon A small constant to avoid log(0); default is 0.1.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{1}{The input data with fitted values, CIs, and growth rates.}
#'     \item{2}{A ggplot of the fitted log-counts with 95\% CI ribbon.}
#'     \item{3}{A ggplot of estimated weekly growth rates with 95\% CI ribbon.}
#'   }
#'
#' @details
#' The function fits \code{y ~ s(time_w, bs = "ps", k = nrow(dat)/8)} using
#' REML smoothing with term selection (\code{select = TRUE}). Growth rates are
#' computed as finite differences of fitted log-counts per week.
#'
#' @importFrom mgcv gam
#' @importFrom ggplot2 ggplot geom_line geom_ribbon aes
#' @importFrom dplyr mutate
#'
#' @examples
#' \dontrun{
#' fit <- fit_gam_grs(mydata)
#' fit[[2]]  # fitted values plot
#' fit[[3]]  # growth rate plot
#' }
#' @export
fit_gam_grs <- function(dat, y = "y", time = "time_w", group = NULL, epsilon = 0.1) {
  # prepare: add log outcome
  dat2 <- dat %>% mutate(log_y = log(.data[[y]] + epsilon))
  
  fit_one <- function(df) {
    k <- max(5, round(nrow(df) / 8))
    fm <- as.formula(paste0("log_y ~ s(", time, ", bs = 'tp', k = ", k, ")"))
    fit <- gam(fm, data = df, method = "REML", select = TRUE)
    pr <- predict(fit, newdata = df, se.fit = TRUE, unconditional = TRUE)
    df$fitted_log <- pr$fit
    df$fitted_se  <- pr$se.fit
    df$lower <- df$fitted_log - 1.96 * df$fitted_se
    df$upper <- df$fitted_log + 1.96 * df$fitted_se
    # derivative per unit of 'time' (assumes time is numeric)
    tvec <- df[[time]]
    df$r_week <- c(NA, diff(df$fitted_log) / diff(tvec))
    df$lower_r <- c(NA, diff(df$lower) / diff(tvec))
    df$upper_r <- c(NA, diff(df$upper) / diff(tvec))
    df
  }
  
  if (!is.null(group) && length(group) > 0) {
    dat_out <- dat2 %>%
      group_by(across(all_of(group))) %>%
      group_modify(~ fit_one(.x)) %>%
      ungroup()
  } else {
    dat_out <- fit_one(dat2)
  }
  
  # plots (use aes_string to keep string interface)
  p_fit <- ggplot(dat_out) +
    geom_ribbon(aes_string(x = "date", ymin = "lower", ymax = "upper"), alpha = 0.25) +
    geom_line(aes_string(x = "date", y = "fitted_log")) +
    labs(y = "fitted log") +
    theme_minimal()
  
  p_gr <- ggplot(dat_out) +
    geom_ribbon(aes_string(x = "date", ymin = "lower_r", ymax = "upper_r"), alpha = 0.25) +
    geom_line(aes_string(x = "date", y = "r_week")) +
    labs(y = "growth (d log / unit time)") +
    theme_minimal()
  
  list(dat_out, p_fit, p_gr)
}
add_doubling_axis <- function(p,n_breaks=11){
  # Get the y-axis range from the built plot
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  
  # Calculate appropriate breaks
  #n_breaks <- 5
  breaks <- pretty(y_range, n = n_breaks)
  
  # Create labels
  times <- log(2) / abs(breaks)
  signs <- ifelse(breaks >= 0, "", "-")
  labels <- ifelse(abs(breaks) < 1e-6, Inf,
                   paste0(signs, round(times, 1)))
  
  # Add secondary axis with calculated breaks and labels
  p + scale_y_continuous(
    breaks=breaks,
    sec.axis = sec_axis(~ .,
                        breaks = breaks,
                        labels = labels,
                        name = "Doubling(+) / Halving(-) time")
  )
  
}

plot_incidence <- function(data, x = "date", y = "ILIplus",
                           color = NULL, facet = NULL, 
                           point = TRUE, line = TRUE) {

  # build aes (use aes_string for concision)
  aes_args <- if (is.null(color)) 
    aes_string(x = x, y = y) 
  else 
    aes_string(x = x, y = y, color = color, group = color)
  
  p <- ggplot(data, aes_args)
  if (line)  p <- p + geom_line()
  if (point) p <- p + geom_point(size = 0.6)
  if (!is.null(facet)) p <- p + facet_wrap(as.formula(paste("~", facet)))
  p + labs(x = x, y = y) 
}

align_x_by_facet <- function(data, x = "date", facet = "Season") {
  data %>%
    group_by(.data[[facet]]) %>%
    mutate(day_of_year = lubridate::yday(.data[[x]])) %>%
    mutate(min_date = first(day_of_year),
           min_int = first(.data[[x]])) %>%
    mutate(diff_int = as.numeric(.data[[x]] - min_int)) %>%
    mutate(day_of_year_shifted = min_date+ diff_int) %>%
    ungroup() 
}

calc_emp_growth <- function(df, date = "date", count = "count",
                            group1 = NULL, group2 = NULL, eps = 0.5,
                            additional_smooth=3) {

  if(!is.null(group1)){
    df <- df %>% group_by(.data[[group1]])
  }
  if(!is.null(group1) & !is.null(group2)){
    df <- df %>% group_by(.data[[group1]],.data[[group2]])
  }
  df <- df %>% arrange(.data[[date]])
  df <- df %>%
    mutate(log_count = log(.data[[count]] + eps)) %>%
    mutate(dt = as.numeric(.data[[date]] - lag(.data[[date]]))) %>%
    mutate(gr = log_count - lag(log_count)) %>%
    mutate(gr = 7*gr/dt)
  
  if(!is.null(additional_smooth)){
    df <- df %>% mutate(gr_smooth = zoo::rollmean(gr, k=additional_smooth,fill=NA,align="right"))
  }
  df %>% ungroup()
}

library(tidyr)
library(dplyr)
school_periods_oxford <- tribble(
  ~start,        ~end,         ~label,           ~acad_year,                ~type,
  
  # 2022/23 holidays
  "2022-10-24", "2022-10-28",  "Half-term Oct 2022",             "2022/23",   "holiday",
  "2022-12-19", "2023-01-02",  "Christmas Holiday 2022/23",      "2022/23",   "holiday",
  "2023-02-13", "2023-02-17",  "Half-term Feb 2023",             "2022/23",   "holiday",
  "2023-04-07", "2023-04-21",  "Easter Holiday 2023",            "2022/23",   "holiday",
  "2023-05-29", "2023-06-02",  "Late Spring Half-term 2023",     "2022/23",   "holiday",
  "2023-07-22", "2023-08-31",  "Summer Holiday 2023",            "2022/23",   "holiday",
  
  # 2023/24 holidays
  "2023-10-23", "2023-10-27",  "Half-term Oct 2023",             "2023/24",   "holiday",
  "2023-12-21", "2024-01-02",  "Christmas Holiday 2023/24",      "2023/24",   "holiday",
  "2024-02-12", "2024-02-16",  "Half-term Feb 2024",             "2023/24",   "holiday",
  "2024-03-29", "2024-04-12",  "Easter Holiday 2024",            "2023/24",   "holiday",
  "2024-05-27", "2024-05-31",  "Late Spring Half-term 2024",     "2023/24",   "holiday",
  "2024-07-20", "2024-08-31",  "Summer Holiday 2024",            "2023/24",   "holiday",
  
  # 2024/25 holidays
  "2024-10-25", "2024-10-29",  "Half-term Oct 2024",             "2024/25",   "holiday",
  "2024-12-20", "2025-01-03",  "Christmas Holiday 2024/25",      "2024/25",   "holiday",
  "2025-02-17", "2025-02-21",  "Half-term Feb 2025",             "2024/25",   "holiday",
  "2025-04-11", "2025-04-22",  "Easter Holiday 2025",            "2024/25",   "holiday",
  "2025-05-26", "2025-05-30",  "Late Spring Half-term 2025",     "2024/25",   "holiday",
  "2025-07-21", "2025-08-31",  "Summer Holiday 2025",            "2024/25",   "holiday",
  
  # 2025/26 holidays
  "2025-10-24", "2025-10-30",  "Half-term Oct 2025",             "2025/26",   "holiday",
  
  
) %>% mutate(start = as.Date(start), end = as.Date(end))


season_from_date <- function(dates, start_month = 7) {
  dates <- as.Date(dates)
  na_idx <- is.na(dates)
  y <- as.integer(format(dates, "%Y"))
  m <- as.integer(format(dates, "%m"))
  start_y <- ifelse(m >= start_month, y, y - 1)
  end_y   <- start_y + 1
  res <- sprintf("%d/%02d", start_y, end_y %% 100)
  res[na_idx] <- NA_character_
  res
}

flu_season <- function(date) {
  date <- as.Date(date)
  
  start_year <- ifelse(
    format(date, "%m") >= "07",
    as.integer(format(date, "%Y")),
    as.integer(format(date, "%Y")) - 1
  )
  
  paste0(
    start_year,
    "/",
    substr(start_year + 1, 3, 4)
  )
}
