# combine_files_and_plots.R
library(tidyverse)
library(patchwork)   # for composing ggplots
library(ggplot2)

# === user config ===
csv_dir   <- "~/Documents/GitHub/influenza_H3N2_k_clade/scenarios/"       # folder containing many .csv files
plots_dir <- "~/Documents/GitHub/influenza_H3N2_k_clade/scenarios/"     # folder containing subfolders, each has one .RData with object `plot_obj`
out_csv   <- "combined_all_csvs.csv"
out_png   <- "combined_plots.png"
ncol_out  <- 2                          # number of columns in final multi-panel

# === combine CSVs ===
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
combined_tbl <- csv_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE), .id = "source_file") %>%
  mutate(Scenario = tools::file_path_sans_ext(basename(source_file))) %>%
  select(-source_file) %>%
  relocate(Scenario)

# write combined csv
write_csv(combined_tbl, out_csv)
message("Wrote combined CSV to: ", out_csv)


# find the cumulative column (first column name that contains "Cumulative")
cum_col <- names(combined_tbl)[str_detect(names(combined_tbl), regex("Cumulative", ignore_case = TRUE))][1]
if (is.na(cum_col)) stop("No column with 'Cumulative' found in combined_tbl")

# extract 65+ rows and keep only Scenario + cumulative value
df65 <- combined_tbl %>%
  filter(`Age group` %in% c("65+", "65+ years", "65 yrs")) %>%   # tolerant matching
  select(Scenario, cum_65 = !!sym(cum_col)) %>%
  mutate(cum_65 = as.numeric(cum_65))

# try to find the base case row (several fallbacks)
base_row <- df65 %>%
  filter(str_detect(tolower(Scenario), "\\bbase_counts\\b|\\bbase\\b")) %>%
  slice_head(n = 1)

if (nrow(base_row) == 0) {
  # fallback: look for scenario starting with "base" or use the first row
  base_row <- df65 %>%
    filter(str_detect(tolower(Scenario), "^base")) %>%
    slice_head(n = 1)
}
if (nrow(base_row) == 0) {
  base_row <- df65 %>% slice_head(n = 1)
  warning("Could not find explicit 'base' scenario; using first Scenario as base: ", base_row$Scenario)
}

base_value <- base_row$cum_65[1]

# produce result table with ratio
result_df <- df65 %>%
  mutate(ratio_to_base = cum_65 / base_value) %>%
  arrange(desc(cum_65))

# optional: nicely formatted numeric columns if desired
result_df_formatted <- result_df %>%
  mutate(
    cum_65 = signif(cum_65, 6),
    ratio_to_base = round(ratio_to_base, 3)
  )

# outputs
result_df      # raw numeric table (mergeable)
result_df_formatted  # human-friendly formatted table

scenario_labels <- c(
  "base_counts"                               = "A. Base case (loosely based on the 2022/23 season)",
  "month_early"                         = "B. Earlier seeding",
  "moderate_immune_escape"              = "C. Moderate immune escape",
  "severe_immune_escape"                = "D. Severe immune escape",
  "r0_10_percent"                       = "E. Moderately higher transmissibility",
  "r0_20_percent"                       = "F. Severely higher transmissibility",
  "moderate_immune_escape_and_earlier"  = "G. Moderate immune escape and earlier seeding",
  "immune_escape_by_age"                = "H. Increased immune escape in younger ages",
  "more_xmas_mixing"                    = "I. Increased mixing in the Christmas period",
  "moderate_immune_escape_earlier_r0_10pct" =
    "J. Moderate immune escape + earlier seeding\n + higher transmissibility"
)

result_df_formatted$Scenario <- scenario_labels[result_df_formatted$Scenario]
result_df_formatted <- result_df_formatted %>% arrange(Scenario)
write.csv(result_df_formatted, "~/Documents/GitHub/influenza_H3N2_k_clade/figures/scenario_table_65.csv")
