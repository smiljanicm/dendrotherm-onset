library(tidyverse)
library(duckdb)
library(DBI)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) stop("Usage: compile_boots file_suffix")
  suff <- args[1]
} else {
  suff <- 'full_gdd'
}

thresh_file <- paste0('./output/thresholds_', suff, '.csv')

clean_df <- read_csv(thresh_file)

con <- dbConnect(duckdb())

diag_file <- paste0('./postproc/output/diag_', suff, '.csv')

conc_raw_all <- dbGetQuery(con, sprintf("
  SELECT species, model_name, boot_id, datafile, term, metric, value
  FROM read_csv_auto('%s')
  WHERE (model_name LIKE '%%gddtb%%_10' OR model_name LIKE '%%_doy_10')
    AND term IN ('s(Temperature)', 's(TempOverallMin)', 's(TempMax)')
    AND (
      metric LIKE 'worst:s(gdd%%'    OR metric LIKE 'observed:s(gdd%%'    OR
      metric LIKE 'worst:s(doy)'     OR metric LIKE 'observed:s(doy)'
    )
", diag_file))

dbDisconnect(con, shutdown = TRUE)

conc_summary_all <- conc_raw_all %>%
  mutate(
    datarange = case_when(
      grepl('gs_', datafile)     ~ '60-172',
      grepl('wide_', datafile)   ~ '32-172',
      grepl('narrow_', datafile) ~ '91-172'
    ),
    metric_fam = case_when(
      grepl('omin', model_name) ~ 'omin',
      grepl('max',  model_name) ~ 'max',
      TRUE                      ~ 'mean'
    ),
    confound_type = if_else(grepl('_doy_', model_name), 'doy', 'gdd'),
    tb        = if_else(confound_type == 'gdd', str_extract(model_name, "(?<=gddtb)\\d"), NA_character_),
    conc_type = str_extract(metric, "^[a-z]+")
  ) %>%
  group_by(species, metric_fam, datarange, confound_type, tb, conc_type) %>%
  summarize(n_boot = n(),
            concurvity_q05 = quantile(value, 0.05, na.rm=T),
            concurvity_med = median(value, na.rm = TRUE),
            concurvity_q95 = quantile(value, 0.95, na.rm=T),
            .groups = "drop") %>%
  pivot_wider(
    names_from  = conc_type,
    values_from = c(concurvity_q05, concurvity_med, concurvity_q95, n_boot),
    names_glue  = "{conc_type}_{.value}"
  )

conc_file <- paste0('./output/concurvity_all_sweep_', suff, '.csv')
conc_fileDE <- paste0('./output/german/concurvity_all_sweep_', suff, '.csv')
write_csv(conc_summary_all, conc_file)
write_csv2(conc_summary_all, conc_fileDE)
