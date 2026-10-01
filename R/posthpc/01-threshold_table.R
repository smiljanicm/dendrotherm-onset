library(tidyverse)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check that at least one argument was provided
  if (length(args) != 1) {
    stop("Usage: compile_boots file_suffix")
  }
  
  # Use the first argument
  suff <- args[1]
} else {
  suff <- 'full'
}

diag_file <- paste0('./hpc/output/diag_', suff, '.csv')
seg_file <- paste0('./hpc/output/seg_', suff, '.csv')
thresh_file <- paste0('./output/tables/thresholds_', suff, '.csv')
threshDE_file <- paste0('./output/tables/german/thresholds_', suff, '.csv')

df <- read_csv(seg_file)
clean_df <- df

thresh <- 
  clean_df %>%
  filter(grepl('_10', model_name)) %>%
  mutate(datarange = case_when(grepl('gs_', datafile) ~ '60-172',
                               grepl('wide_', datafile) ~ '32-172',
                               grepl('narrow_', datafile) ~ '91-172')) %>%
  select(-datafile) %>%
  group_by(species, model_name, datarange) %>%
  pivot_longer(c(-boot_id, -temp_min, -temp_max, -species, -model_name,-t_q05, -t_q10, -t_q90, -t_q95, -obs_max_prob, -datarange)) %>% 
  mutate(value = case_when(value == temp_max ~ NA_real_,
                           T ~ value),
         value = case_when(value == temp_min ~ NA_real_,
                           T ~ value)) %>% 
  pivot_wider() %>% 
  summarize(
    n          = n(),
    t_min      = min(temp_min),
    t_q05      = median(t_q05),
    t_q10      = median(t_q10),
    t_q90      = median(t_q90),
    t_q95      = median(t_q95),
    t_max      = max(temp_max),
    obs_max_prob = median(obs_max_prob),
    # Absolute thresholds (Rossi-comparable)
    onset_p50_low  = as.numeric(quantile(onset_p50, 0.05, na.rm = TRUE)),
    onset_p50_med      = as.numeric(quantile(onset_p50, 0.50, na.rm = TRUE)),
    onset_p50_high = as.numeric(quantile(onset_p50, 0.95, na.rm = TRUE)),
    onset_p25_low  = as.numeric(quantile(onset_p25, 0.05, na.rm = TRUE)),
    onset_p25_med      = as.numeric(quantile(onset_p25, 0.50, na.rm = TRUE)),
    onset_p25_high = as.numeric(quantile(onset_p25, 0.95, na.rm = TRUE)),
    onset_p10_low  = as.numeric(quantile(onset_p10, 0.05, na.rm = TRUE)),
    onset_p10_med      = as.numeric(quantile(onset_p10, 0.50, na.rm = TRUE)),
    onset_p10_high = as.numeric(quantile(onset_p10, 0.95, na.rm = TRUE)),
    onset_p05_low  = as.numeric(quantile(onset_p05, 0.05, na.rm = TRUE)),
    onset_p05_med      = as.numeric(quantile(onset_p05, 0.50, na.rm = TRUE)),
    onset_p05_high = as.numeric(quantile(onset_p05, 0.95, na.rm = TRUE)),
    # Relative thresholds
    onset_rel50_low  = as.numeric(quantile(onset_rel50, 0.05, na.rm = TRUE)),
    onset_rel50_med      = as.numeric(quantile(onset_rel50, 0.50, na.rm = TRUE)),
    onset_rel50_high = as.numeric(quantile(onset_rel50, 0.95, na.rm = TRUE)),
    onset_rel25_low  = as.numeric(quantile(onset_rel25, 0.05, na.rm = TRUE)),
    onset_rel25_med      = as.numeric(quantile(onset_rel25, 0.50, na.rm = TRUE)),
    onset_rel25_high = as.numeric(quantile(onset_rel25, 0.95, na.rm = TRUE)),
    onset_rel10_low  = as.numeric(quantile(onset_rel10, 0.05, na.rm = TRUE)),
    onset_rel10_med      = as.numeric(quantile(onset_rel10, 0.50, na.rm = TRUE)),
    onset_rel10_high = as.numeric(quantile(onset_rel10, 0.95, na.rm = TRUE)),
    onset_rel05_low  = as.numeric(quantile(onset_rel05, 0.05, na.rm = TRUE)),
    onset_rel05_med      = as.numeric(quantile(onset_rel05, 0.50, na.rm = TRUE)),
    onset_rel05_high = as.numeric(quantile(onset_rel05, 0.95, na.rm = TRUE)),
    .groups = "drop"
  )

thresh %>% print(n=Inf)

thresh %>%
  write_csv(thresh_file)

thresh %>%
  write_csv2(threshDE_file)