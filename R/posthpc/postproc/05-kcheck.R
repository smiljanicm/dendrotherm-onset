library(duckdb)
library(DBI)
library(tidyverse)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Usage: compile_boots file_suffix")
  }
  suff <- args[1]
} else {
  suff <- 'full_gdd'
}

diag_file <- paste0('./postproc/output/diag_', suff, '.csv')

con <- dbConnect(duckdb())

kcheck_raw <- dbGetQuery(con, sprintf("
  SELECT species, model_name, boot_id, datafile, term, metric, value
  FROM read_csv_auto('%s')
  WHERE model_name LIKE '%%_10'
    AND metric IN ('p-value', 'k-index', 'edf', 'k''')
", diag_file))

dbDisconnect(con, shutdown = TRUE)

kcheck_summary <- kcheck_raw %>%
  mutate(
    datarange = case_when(
      grepl('gs_', datafile)     ~ '60-172',
      grepl('wide_', datafile)   ~ '32-172',
      grepl('narrow_', datafile) ~ '91-172',
      TRUE ~ NA_character_
    ),
    metric_fam = case_when(
      grepl('omin', model_name) ~ 'omin',
      grepl('max',  model_name) ~ 'max',
      TRUE                      ~ 'mean'
    ),
    model_type = case_when(
      grepl('_base_', model_name) ~ 'base',
      grepl('_doy_',  model_name) ~ 'doy',
      grepl('gddtb',  model_name) ~ paste0('gdd_tb', str_extract(model_name, "(?<=gddtb)\\d")),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(model_type), !is.na(datarange)) %>%
  pivot_wider(
    id_cols = c(species, metric_fam, datarange, model_type, term, boot_id),
    names_from = metric,
    values_from = value
  ) %>%
  group_by(species, metric_fam, datarange, model_type, term) %>%
  summarize(
    n_boot        = n(),
    median_pvalue = median(`p-value`, na.rm = TRUE),
    frac_fail     = mean(`p-value` < 0.05, na.rm = TRUE),
    median_edf    = median(edf, na.rm = TRUE),
    median_kprime = median(`k'`, na.rm = TRUE),
    edf_ratio     = median_edf / median_kprime,
    .groups = "drop"
  )

kcheck_file <- paste0('./output/kcheck_summary_', suff, '.csv')
kcheck_fileDE <- paste0('./output/german/kcheck_summary_', suff, '.csv')
write_csv(kcheck_summary, kcheck_file)
write_csv2(kcheck_summary, kcheck_fileDE)

