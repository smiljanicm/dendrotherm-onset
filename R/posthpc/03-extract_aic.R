library(tidyverse)
library(duckdb)
library(DBI)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) stop("Usage: compile_boots file_suffix")
  suff <- args[1]
} else {
  suff <- 'full'
}

thresh_file <- paste0('./output/tables/thresholds_', suff, '.csv')

clean_df <- read_csv(thresh_file)

con <- dbConnect(duckdb())

diag_file <- paste0('./hpc/output/diag_', suff, '.csv')

aic_raw_all <- dbGetQuery(con, sprintf("
  SELECT species, model_name, boot_id, datafile, value AS aic
  FROM read_csv_auto('%s')
  WHERE metric = 'AIC'
    AND (model_name LIKE '%%_base_10'
         OR model_name LIKE '%%_doy_10'
         OR model_name LIKE '%%gddtb%%_10')
", diag_file))

dbDisconnect(con, shutdown = TRUE)

aic_wide_all <- aic_raw_all %>%
  mutate(
    datarange = case_when(
      grepl('gs_', datafile)    ~ '60-172',
      grepl('wide_', datafile)  ~ '32-172',
      grepl('narrow_', datafile) ~ '91-172'
    ),
    metric_fam = case_when(
      grepl('omin', model_name) ~ 'omin',
      grepl('max',  model_name) ~ 'max',
      TRUE                      ~ 'mean'
    ),
    model_type = case_when(
      grepl('_base_', model_name) ~ 'base',
      grepl('_doy_', model_name)  ~ 'doy',
      grepl('gddtb', model_name)  ~ paste0("gdd_tb", str_extract(model_name, "(?<=gddtb)\\d")),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(model_type)) %>%
  select(species, metric_fam, datarange, boot_id, model_type, aic) %>%
  pivot_wider(names_from = model_type, values_from = aic)

# One row per species × metric_fam × datarange × Tb, delta vs base
aic_summary_all <- aic_wide_all %>%
  pivot_longer(
    cols = starts_with("gdd_tb"),
    names_to = "tb", values_to = "aic_gdd", names_prefix = "gdd_tb"
  ) %>%
  mutate(
    delta_aic_gdd = aic_gdd - base,
    delta_aic_doy = doy - base
  ) %>%
  group_by(species, metric_fam, datarange, tb) %>%
  summarize(
    n_boot                = n(),
    median_aic_base        = median(base, na.rm = TRUE),
    median_aic_doy          = median(doy, na.rm = TRUE),
    median_delta_aic_doy    = median(delta_aic_doy, na.rm = TRUE),
    frac_doy_better         = mean(delta_aic_doy < -2, na.rm = TRUE),
    median_aic_gdd          = median(aic_gdd, na.rm = TRUE),
    median_delta_aic_gdd    = median(delta_aic_gdd, na.rm = TRUE),
    frac_gdd_better         = mean(delta_aic_gdd < -2, na.rm = TRUE),
    .groups = "drop"
  )

aic_file <- paste0('./output/tables/aic_comparison_all_sweep_', suff, '.csv')
aic_fileDE <- paste0('./output/tables/german/aic_comparison_all_sweep_', suff, '.csv')
write_csv(aic_summary_all, aic_file)
write_csv2(aic_summary_all, aic_fileDE)
