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
