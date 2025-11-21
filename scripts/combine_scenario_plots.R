# combine_plots_from_RData_facets.R
library(tidyverse)
library(ggplot2)

# ---- CONFIG ----
plots_dir <- "~/Documents/GitHub/influenza_H3N2_k_clade/scenarios/" # each subfolder must contain one .RData with 'res' (and optionally 'plot_obj'/'plot_data')
out_png   <- "combined_facet_plot.png"
ncol_out  <- 2
age_group_key_plot <- c("inc_1_1" = "[0,5)", "inc_2_1" = "[5,18)", "inc_3_1" = "[18,65)", "inc_4_1" = "65+")

# ---- helpers ----
safe_load_res <- function(rdata_path) {
  e <- new.env()
  tryCatch({
    load(rdata_path, envir = e)
    if (exists("res", envir = e)) return(get("res", envir = e))
    NULL
  }, error = function(err) {
    warning("Failed loading ", rdata_path, ": ", conditionMessage(err))
    NULL
  })
}

# ---- iterate folders and extract canonical dataframes from res ----
folders <- list.dirs(plots_dir, recursive = FALSE, full.names = TRUE)
if (length(folders) == 0) stop("No subfolders found in ", plots_dir)

all_inc <- list()
all_rects <- list()
all_annots <- list()
scenarios <- character()

for (fld in folders) {
  rfiles <- list.files(fld, pattern = "\\.RData$", full.names = TRUE)
  if (length(rfiles) == 0) next
  res <- safe_load_res(rfiles[[1]])
  if (is.null(res)) next
  
  scen <- basename(fld)
  scenarios <- c(scenarios, scen)
  
  # build weekly incidence dataframe (same logic as app)
  inc_mat <- as.matrix(res$inc)
  inc_df <- as.data.frame(inc_mat)
  inc_df$t <- seq_len(nrow(inc_df))
  inc_long <- inc_df %>%
    pivot_longer(-t, names_to = "age_group_raw", values_to = "incidence") %>%
    mutate(age_group = age_group_key_plot[age_group_raw],
           date = res$date_key$date[t],
           Scenario = scen)
  
  all_inc[[scen]] <- inc_long
  
  # build rectangle shading for that scenario using res$meta dates
  meta <- res$meta
  meta$shopping_period_end <- "2022-12-15"
  rects <- tibble(
    xmin = as.Date(c(meta$half_term_start, meta$shopping_period_start, meta$winter_holiday_start)),
    xmax = as.Date(c(meta$half_term_end,   meta$shopping_period_end,   meta$winter_holiday_end)),
    fill  = c("Half term", "Christmas period (pre holidays)", "Christmas holidays"),
    Scenario = scen
  )
  all_rects[[scen]] <- rects
  
  # compute annotation values (peak reported, peak date, cumulative symptomatic if available)
  weekly_totals <- rowSums(as.matrix(res$inc), na.rm = TRUE) # reported weekly sum across age groups
  peak_reported <- if (length(weekly_totals)) max(weekly_totals, na.rm = TRUE) else NA_real_
  peak_idx <- if (length(weekly_totals)) which.max(weekly_totals) else NA_integer_
  peak_date <- if (!is.na(peak_idx)) res$date_key$date[peak_idx] else as.Date(NA)
  
  # cumulative symptomatic: prefer res$cumulative_incidence$Symptomatic Total row, else sum age Symptomatic
  cum_symp <- NA_real_
  if (!is.null(res$cumulative_incidence)) {
    if ("Age group" %in% colnames(res$cumulative_incidence)) {
      total_row <- res$cumulative_incidence %>% filter(`Age group` == "Total")
      if (nrow(total_row) == 1 && "Symptomatic" %in% colnames(total_row)) {
        cum_symp <- total_row$Symptomatic
      } else if ("Symptomatic" %in% colnames(res$cumulative_incidence)) {
        cum_symp <- sum(res$cumulative_incidence$Symptomatic, na.rm = TRUE)
      }
    }
  }
  reporting_rate <- if (!is.null(meta$reporting_rate)) meta$reporting_rate else NA_real_
  peak_symp_est <- if (!is.na(reporting_rate) && reporting_rate > 0) peak_reported / reporting_rate else NA_real_
  
  label_text <- paste0(
    "Peak (reported): ", ifelse(is.na(peak_reported), "n/a", format(round(peak_reported), big.mark = ",")), "\n",
    "Date of peak: ", ifelse(is.na(peak_date), "n/a", format(peak_date, "%Y-%m-%d")), "\n",
    "Cumulative symptomatic (model): ", ifelse(is.na(cum_symp), "n/a", format(round(cum_symp), big.mark = ",")), "\n",
    ifelse(!is.na(peak_symp_est), paste0("Peak symptomatic (est): ", format(round(peak_symp_est), big.mark = ",")), "")
  )
  all_annots[[scen]] <- tibble(Scenario = scen, label = label_text)
}

if (length(all_inc) == 0) stop("No 'res' data extracted from any .RData files.")

combined_inc <- bind_rows(all_inc)
combined_inc$Scenario <- factor(combined_inc$Scenario, levels = unique(scenarios))

# combine rects and annots
rects_df <- bind_rows(all_rects)
annots_df <- bind_rows(all_annots)

# compute ymax per scenario for annotation placement
ymax_df <- combined_inc %>% group_by(Scenario) %>% summarize(ymax = max(incidence, na.rm = TRUE), .groups = "drop")
annots_df <- left_join(annots_df, ymax_df, by = "Scenario") %>%
  mutate(x = max(combined_inc$date, na.rm = TRUE) - 3, y = ymax * 0.98)

# ensure age_group factor ordering (shared)
age_levels <- c("[0,5)","[5,18)","[18,65)","65+")
combined_inc$age_group <- factor(combined_inc$age_group, levels = age_levels)

scenario_labels <- c(
  "base"                               = "A. Base case (loosely based on the 2022/23 season)",
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
keep_levels <- names(scenario_labels)
keep_levels_labels <- scenario_labels
combined_inc$Scenario <- factor(combined_inc$Scenario, levels=keep_levels)

# apply mapping
combined_inc$Scenario <- factor(
  scenario_labels[as.character(combined_inc$Scenario)],
  levels = scenario_labels[levels(combined_inc$Scenario)]
)

# Apply same mapping to rectangles and annotation data frames:
rects_df$Scenario  <- scenario_labels[rects_df$Scenario]
annots_df$Scenario <- scenario_labels[annots_df$Scenario]

# ensure ordering is preserved:
combined_inc$Scenario <- factor(combined_inc$Scenario,levels = scenario_labels)
rects_df$Scenario  <- factor(rects_df$Scenario,  levels = levels(combined_inc$Scenario))
annots_df$Scenario <- factor(annots_df$Scenario, levels = levels(combined_inc$Scenario))

# ---- build shared-facet plot ----
p <- ggplot() +
  geom_rect(data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = combined_inc, aes(x = date, y = incidence, color = age_group, group = age_group), size = 0.9) +
  facet_wrap(~ Scenario, ncol = ncol_out) +
  scale_color_brewer("Age group", palette = "Set1") +
  scale_fill_brewer("Holiday period", palette = "Set2") +
  #geom_label(data = annots_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
  #           hjust = 1, vjust = 1, size = 3.2, fill = "white", alpha = 0.85) +
  theme_bw() +
  labs(x = "Date", y = "Reported symptomatic influenza (weekly)") +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size=10),
    legend.position = "bottom",
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.box = "vertical",              # <-- stack legends vertically
    legend.spacing.y = unit(4, "pt"),     # <-- spacing between legend boxes
    panel.spacing.x = unit(1, "cm")       # (your earlier spacing adjustment)
  ) +
  coord_cartesian(ylim=c(0,5000), expand = FALSE,xlim=as.Date(c("2022-09-01","2023-02-01")))

# Save
ggsave(out_png, p, width = 9, height = 11, dpi = 300)
message("Saved combined plot: ", out_png)
break
# --- Clean all quotes from Scenario & label ---
annots_clean <- annots_df %>%
  mutate(
    Scenario = str_replace_all(Scenario, "[\"'`’‘“”]", ""),
    label    = str_replace_all(label,    "[\"'`’‘“”]", "")
  )

# --- Regex pattern for your exact label structure ---
pat <- paste0(
  "Peak \\(reported\\):\\s*([0-9,]+)\\s*\n",
  "Date of peak:\\s*([0-9]{4}-[0-9]{2}-[0-9]{2})\\s*\n",
  "Cumulative symptomatic \\(model\\):\\s*([0-9,]+)\\s*\n",
  "Peak symptomatic \\(est\\):\\s*([0-9,]+)"
)

# --- Parse into tidy explicit columns ---
parsed_annots <- annots_clean %>%
  mutate(matches = str_match(label, pat)[, -1]) %>%
  mutate(
    peak_reported           = matches[,1] %>% str_remove_all(",") %>% as.numeric(),
    peak_date               = matches[,2] %>% ymd(),
    cumulative_symptomatic  = matches[,3] %>% str_remove_all(",") %>% as.numeric(),
    peak_symptomatic_est    = matches[,4] %>% str_remove_all(",") %>% as.numeric()
  ) %>%
  select(Scenario, peak_reported, peak_date,
         cumulative_symptomatic, peak_symptomatic_est)

# --- Identify Base scenario safely ---
base_idx <- which(str_detect(tolower(parsed_annots$Scenario), "\\bbase case\\b") |
                    str_detect(tolower(parsed_annots$Scenario), "\\bbase\\b"))

if (length(base_idx) == 0) {
  warning("Could not find Base scenario. Using first row instead.")
  base_idx <- 1
}

base_vals <- parsed_annots[base_idx,
                           c("peak_reported","cumulative_symptomatic","peak_symptomatic_est")] %>%
  unlist()

# --- Create ratio table (retain peak_date) ---
ratio_table <- parsed_annots %>%
  mutate(
    ratio_peak_reported          = peak_reported / base_vals["peak_reported"],
    ratio_cumulative_symptomatic = cumulative_symptomatic / base_vals["cumulative_symptomatic"],
    ratio_peak_symptomatic_est   = peak_symptomatic_est / base_vals["peak_symptomatic_est"]
  ) %>%
  select(Scenario, peak_date,
         ratio_peak_reported,
         ratio_cumulative_symptomatic,
         ratio_peak_symptomatic_est)

# --- Print output ---
parsed_annots
ratio_table

write.csv(ratio_table %>% arrange(Scenario),"~/Documents/GitHub/influenza_H3N2_k_clade/figures/scenario_table.csv")
