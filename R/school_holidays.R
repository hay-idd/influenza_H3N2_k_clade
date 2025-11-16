library(lubridate)
library(dplyr)

# helper: Monday on/after a date
next_monday_on_or_after <- function(date) {
  date <- as.Date(date)
  w <- wday(date, week_start = 1)
  date + ((1 - w) %% 7)
}
# helper: Monday on/before a date
prev_monday_on_or_before <- function(date) {
  date <- as.Date(date)
  w <- wday(date, week_start = 1)
  date - ((w - 1) %% 7)
}
# compute Easter (Meeus algorithm) for given year (Gregorian)
easter_date <- function(y) {
  a <- y %% 19
  b <- y %/% 100
  c <- y %% 100
  d <- b %/% 4
  e <- b %% 4
  f <- (b + 8) %/% 25
  g <- (b - f + 1) %/% 3
  h <- (19*a + b - d - g + 15) %% 30
  i <- c %/% 4
  k <- c %% 4
  l <- (32 + 2*e + 2*i - h - k) %% 7
  m <- (a + 11*h + 22*l) %/% 451
  month <- (h + l - 7*m + 114) %/% 31
  day   <- ((h + l - 7*m + 114) %% 31) + 1
  as.Date(sprintf("%04d-%02d-%02d", y, month, day))
}

years <- 2009:2025  # academic-year starting year e.g. 2009 -> 2009/10
rows <- lapply(years, function(y) {
  ay <- sprintf("%d/%02d", y, (y+1) %% 100)
  # Oct half-term (week containing Oct 24 of start year)
  oct_anchor <- as.Date(sprintf("%04d-10-24", y))
  oct_start <- prev_monday_on_or_before(oct_anchor)
  oct_end   <- oct_start + 4
  
  # Christmas: Dec 20 (start) -> Jan 2 (next year)
  xmas_start <- as.Date(sprintf("%04d-12-20", y))
  xmas_end   <- as.Date(sprintf("%04d-01-02", y+1))
  
  # Feb half-term: in next calendar year (y+1), week containing Feb 14
  feb_anchor <- as.Date(sprintf("%04d-02-14", y+1))
  feb_start  <- prev_monday_on_or_before(feb_anchor)
  feb_end    <- feb_start + 4
  
  # Easter: compute Easter for calendar year y+1; start = Monday one week before Easter; two weeks
  eas <- easter_date(y+1)
  eas_start <- prev_monday_on_or_before(eas - 7)
  eas_end   <- eas_start + 13
  
  # Late May half-term: last Monday of May (y+1)
  last_may_day <- as.Date(sprintf("%04d-05-31", y+1))
  may_start <- prev_monday_on_or_before(last_may_day)
  may_end   <- may_start + 4
  
  # Summer: Monday on/after 20 July (y+1) -> 31 Aug (y+1)
  july_anchor <- as.Date(sprintf("%04d-07-20", y+1))
  summer_start <- next_monday_on_or_after(july_anchor)
  summer_end   <- as.Date(sprintf("%04d-08-31", y+1))
  
  tibble::tribble(
    ~start, ~end, ~label, ~acad_year, ~type,
    as.character(oct_start), as.character(oct_end), paste("Half-term Oct", y), ay, "holiday",
    as.character(xmas_start), as.character(xmas_end), paste("Christmas", ay), ay, "holiday",
    as.character(feb_start), as.character(feb_end), paste("Half-term Feb", y+1), ay, "holiday",
    as.character(eas_start), as.character(eas_end), paste("Easter", y+1), ay, "holiday",
    as.character(may_start), as.character(may_end), paste("Late May Half-term", y+1), ay, "holiday",
    as.character(summer_start), as.character(summer_end), paste("Summer Holiday", y+1), ay, "holiday"
  )
})

school_periods_oxford_2009_2025 <- dplyr::bind_rows(rows) %>%
  mutate(start = as.Date(start), end = as.Date(end))

make_daily_holiday_calendar <- function(hol, buffer_days = 7) {
  hol <- hol %>% mutate(start = as.Date(start), end = as.Date(end))
  # all calendar days
  all_dates <- tibble(date = seq(min(hol$start), max(hol$end), by = "day"))
  # explicit holiday days with label
  holiday_days <- hol %>%
    rowwise() %>%
    mutate(date = list(seq(start, end, by = "day"))) %>%
    unnest(cols = c(date)) %>%
    select(date, label, acad_year)
  # buffered days (± buffer_days) excluding true holiday days
  buffer_days_df <- hol %>%
    rowwise() %>%
    mutate(date = list(seq(start - buffer_days, end + buffer_days, by = "day"))) %>%
    unnest(cols = c(date)) %>%
    select(date, label, acad_year) %>%
    anti_join(holiday_days, by = "date")
  # join and assign priority: in_holiday > within_1week > term_time
  all_dates %>%
    left_join(holiday_days, by = "date") %>%
    left_join(buffer_days_df %>% rename(buf_label = label, buf_acad_year = acad_year),
              by = "date") %>%
    mutate(
      type = case_when(
        !is.na(label) ~ "in_holiday",
        !is.na(buf_label) ~ "within_1week",
        TRUE ~ "term_time"
      ),
      label = coalesce(label, buf_label),
      acad_year = coalesce(acad_year, buf_acad_year)
    ) %>%
    select(date, type, label, acad_year)
}

# run it
school_days_full <- make_daily_holiday_calendar(school_periods_oxford_2009_2025, buffer_days = 7)


make_daily_holiday_calendar2 <- function(hol, buffer_days = 7) {
  hol <- hol %>% mutate(start = as.Date(start), end = as.Date(end))
  all_dates <- tibble(date = seq(min(hol$start) - buffer_days, max(hol$end), by = "day"))
  # holiday days with original start (needed to compute time_since)
  holiday_days <- hol %>%
    rowwise() %>%
    mutate(date = list(seq(start, end, by = "day"))) %>%
    unnest(cols = c(date)) %>%
    select(date, label, acad_year, start_hol = start)
  # buffer days = week before holiday only (start - buffer_days ... start-1)
  buffer_days_df <- hol %>%
    rowwise() %>%
    mutate(date = list(seq(start - buffer_days, start - 1, by = "day"))) %>%
    unnest(cols = c(date)) %>%
    select(date, label, acad_year, start_hol = start)
  # join and assign priority: in_holiday > within_1week > term_time
  df <- all_dates %>%
    left_join(holiday_days, by = "date") %>%
    left_join(buffer_days_df %>% rename(buf_label = label, buf_acad_year = acad_year, buf_start = start_hol),
              by = "date") %>%
    mutate(
      type = case_when(
        !is.na(label)         ~ "in_holiday",
        !is.na(buf_label)     ~ "within_1week",
        TRUE                  ~ "term_time"
      ),
      label = coalesce(label, buf_label),
      acad_year = coalesce(acad_year, buf_acad_year),
      # compute time since (start - buffer_days) for both buffer and holiday days
      time_since_preweek = case_when(
        !is.na(start_hol) ~ as.integer(date - (start_hol - buffer_days)),   # holiday days
        !is.na(buf_start) ~ as.integer(date - (buf_start - buffer_days)),   # buffer days
        TRUE              ~ 0L
      )
    ) %>%
    select(date, type, label, acad_year, time_since_preweek)
  df
}

# run it
school_days_full2 <- make_daily_holiday_calendar2(school_periods_oxford_2009_2025, buffer_days = 7)

# quick checks
school_days_full2 %>% filter(type != "term_time") %>% slice(1:12)
school_days_full2 %>% filter(type == "term_time") %>% slice(1:6)


