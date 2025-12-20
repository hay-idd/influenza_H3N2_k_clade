# combine_sizes_from_subfolders.R
library(tidyverse)

# === user config ===
csv_dir   <- "~/Documents/GitHub/influenza_H3N2_k_clade/new_scenarios/"   # top-level folder containing many subfolders
out_table <- "~/Documents/GitHub/influenza_H3N2_k_clade/new_figures/scenario_table_65_and_total.csv"

# === read subfolders' cumulative_incidence.csv ===
subdirs <- list.dirs(csv_dir, recursive = FALSE, full.names = TRUE)
all_tbls <- map_dfr(subdirs, function(d) {
  f <- file.path(d, "cumulative_incidence.csv")
  if (!file.exists(f)) return(tibble())          # skip missing files
  read_csv(f, show_col_types = FALSE) %>%
    mutate(Scenario = basename(d)) %>%
    relocate(Scenario)
})

if (nrow(all_tbls) == 0) stop("No tables read from subfolders")

# try to identify the cumulative column (first reasonable match)
cum_col_candidates <- names(all_tbls)[str_detect(names(all_tbls),
                                                 regex("Cumulative|total_inf|cumulative influenza|final_size|total", ignore_case = TRUE))]
if (length(cum_col_candidates) == 0) stop("No cumulative-like column found")
cum_col <- cum_col_candidates[1]

# make sure it's numeric (coerce safely)
all_tbls <- all_tbls %>% mutate(!!cum_col := as.numeric(!!sym(cum_col)))

# --- compute final_size_total per scenario ---
# prefer explicit 'Total' row, otherwise sum numeric rows by scenario
df_total_explicit <- all_tbls %>%
  filter(str_to_lower(`Age group`) %in% c("total")) %>%
  select(Scenario, final_size_total = !!sym(cum_col)) %>%
  distinct()

if (nrow(df_total_explicit) > 0) {
  df_total <- df_total_explicit
} else {
  df_total <- all_tbls %>%
    filter(!is.na(!!sym(cum_col))) %>%
    group_by(Scenario) %>%
    summarise(final_size_total = sum(!!sym(cum_col), na.rm = TRUE), .groups = "drop")
}

# --- compute final_size_65 per scenario (tolerant matching) ---
df_65 <- all_tbls %>%
  filter(str_detect(`Age group`, regex("^(65\\+|65|65\\s*years|65\\s*yrs)$", ignore_case = TRUE))) %>%
  group_by(Scenario) %>%
  summarise(final_size_65 = sum(!!sym(cum_col), na.rm = TRUE), .groups = "drop")

# join
sizes_df <- df_total %>%
  left_join(df_65, by = "Scenario")

# --- find base values (use raw scenario names) ---
find_base_value <- function(df, value_col) {
  # search for common base names
  idx <- which(str_detect(tolower(df$Scenario), "\\bbase_counts\\b|\\bbase\\b"))
  if (length(idx) == 0) idx <- which(str_detect(tolower(df$Scenario), "^base"))
  if (length(idx) == 0) idx <- 1                # fallback to first
  val <- df[[value_col]][idx[1]]
  if (is.na(val) || val == 0) warning("Base value is NA or zero; check scenarios")
  val
}
base_total_value <- find_base_value(sizes_df, "final_size_total")
base_65_value    <- if ("final_size_65" %in% names(sizes_df)) find_base_value(sizes_df, "final_size_65") else NA_real_

# --- compute ratios safely ---
result_df <- sizes_df %>%
  mutate(
    final_size_total = as.numeric(final_size_total),
    final_size_65    = as.numeric(final_size_65),
    final_size_total_rel = ifelse(!is.na(base_total_value) & base_total_value > 0,
                                  final_size_total / base_total_value, NA_real_),
    final_size_65_rel = ifelse(!is.na(base_65_value) & base_65_value > 0,
                               final_size_65 / base_65_value, NA_real_)
  ) %>%
  arrange(desc(final_size_total))

# optional: pretty labels mapping (do AFTER computing ratios)
scenario_labels <- c(
  "base" = "A. Base case (loosely based on 2022/23)",
  "base_counts" = "A. Base case (loosely based on 2022/23)",
  "moderate_immune_escape" = "B. Moderate immune escape (5% overall)",
  "high_immune_escape" = "C. High immune escape (10% overall)",
  "highest_immune_escape" = "D. Highest immune escape (20% overall)",
  "fifty_percent_immune_kids" = "E. Reduce immune fraction in <18 by 50%",
  "youngest_20_older_50_immune_escape" = "F. Reduce immune fraction to 20% in\n <5 and 50% in 5-18",
  "higher_r0" = "G. Higher transmissibility (R0=2.2)",
  "highest_r0" = "H. Highest transmissibility (R0=2.4)",
  "two_weeks_early" = "I. Two weeks earlier seeding",
  "month_early" = "J. One month earlier seeding",
  "two_weeks_moderate_escape" = "K. Two weeks earlier seeding and 5% immune escape"
)

result_df <- result_df %>%
  mutate(Scenario_pretty = if_else(Scenario %in% names(scenario_labels),
                                   scenario_labels[Scenario], Scenario)) %>%
  select(Scenario = Scenario_pretty, everything(), -Scenario)

# write out

result_df <- result_df %>% arrange(Scenario)
result_df$final_size_total_rel <- result_df$final_size_total/result_df$final_size_total[result_df$Scenario == "A. Base case (loosely based on 2022/23)"]
result_df$final_size_65_rel <- result_df$final_size_65/result_df$final_size_65[result_df$Scenario == "A. Base case (loosely based on 2022/23)"]

write.csv(result_df, out_table)
message("Wrote final-size table to: ", out_table)

# return in interactive sessions
invisible(result_df)

