library(tidyverse)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) stop("Usage: compile_boots file_suffix")
  suff <- args[1]
} else {
  suff <- 'full'
}

thresh_file <- paste0('./output/tables/thresholds_', suff, '.csv')

clean_df <- read_csv(thresh_file)

# ── 1. Threshold table — base, doy, and full Tb sweep (gdd_tb0..6) ──────────
thresh_all <- clean_df %>%
  mutate(
    metric_fam = case_when(
      grepl("omin", model_name) ~ "omin",
      grepl("max",  model_name) ~ "max",
      TRUE                      ~ "mean"
    ),
    model_type = case_when(
      grepl("_base_", model_name) ~ "base",
      grepl("_doy_",  model_name) ~ "doy",
      grepl("gddtb",  model_name) ~ paste0("gdd_tb", str_extract(model_name, "(?<=gddtb)\\d")),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(model_type)) %>%
  select(species, metric_fam, datarange, model_type,
         onset_p50_med, onset_p25_med, onset_rel50_med, onset_rel25_med,
         obs_max_prob, t_q05, t_q95)

# ── 1b. Monotonicity check across full Tb sweep (p50 and rel50) ────────────
thresh_sweep_all <- clean_df %>%
  filter(grepl("gddtb", model_name), grepl("_10$", model_name)) %>%
  mutate(
    metric_fam = case_when(
      grepl("omin", model_name) ~ "omin",
      grepl("max",  model_name) ~ "max",
      TRUE                      ~ "mean"
    ),
    tb = as.integer(str_extract(model_name, "(?<=gddtb)\\d"))
  ) %>%
  select(species, metric_fam, datarange, tb, onset_p50_med, onset_rel50_med) %>%
  arrange(species, metric_fam, datarange, tb) %>%
  group_by(species, metric_fam, datarange) %>%
  summarize(
    monotonic_p50   = all(diff(onset_p50_med)   >= -0.5, na.rm = TRUE),
    monotonic_rel50 = all(diff(onset_rel50_med) >= -0.5, na.rm = TRUE),
    .groups = "drop"
  )

# ── 2. Pivot wide — base, doy, gdd_tb0..6 side by side ──────────────────────
delta_gdd_all <- thresh_all %>%
  pivot_wider(
    names_from  = model_type,
    values_from = c(onset_p50_med, onset_p25_med,
                    onset_rel50_med, onset_rel25_med,
                    obs_max_prob, t_q05, t_q95)
  ) %>%
  # ── p50 deltas: GDD sweep vs base ──────────────────────────────────────────
  mutate(across(
    starts_with("onset_p50_med_gdd_tb"),
    ~ . - onset_p50_med_base,
    .names = "delta_p50_{.col}"
  )) %>%
  # ── rel50 deltas: GDD sweep vs base ────────────────────────────────────────
  mutate(across(
    starts_with("onset_rel50_med_gdd_tb"),
    ~ . - onset_rel50_med_base,
    .names = "delta_rel50_{.col}"
  )) %>%
  # ── doy deltas (both onset types) ──────────────────────────────────────────
  mutate(
    delta_p50_doy   = onset_p50_med_doy   - onset_p50_med_base,
    delta_rel50_doy = onset_rel50_med_doy - onset_rel50_med_base
  )

# ── 3. Clean up GDD-sweep column names ──────────────────────────────────────
delta_gdd_all <- delta_gdd_all %>%
  rename_with(~ str_replace(., "delta_p50_onset_p50_med_gdd_tb", "delta_p50_tb"),
              starts_with("delta_p50_onset_p50_med_gdd_tb")) %>%
  rename_with(~ str_replace(., "delta_rel50_onset_rel50_med_gdd_tb", "delta_rel50_tb"),
              starts_with("delta_rel50_onset_rel50_med_gdd_tb"))

# ── 3b. Join monotonicity flags ──────────────────────────────────────────────
delta_gdd_all <- delta_gdd_all %>%
  left_join(thresh_sweep_all, by = c("species", "metric_fam", "datarange"))

# ── 4. Worst-case shift flags across GDD sweep + doy + monotonicity ─────────
delta_gdd_all <- delta_gdd_all %>%
  mutate(
    max_abs_delta_p50   = pmax(!!!syms(paste0("delta_p50_tb", 0:6)), na.rm = TRUE) %>% abs(),
    max_abs_delta_rel50 = pmax(!!!syms(paste0("delta_rel50_tb", 0:6)), na.rm = TRUE) %>% abs(),
    flag_large_shift_gdd_p50   = if_else(!is.na(max_abs_delta_p50)   & max_abs_delta_p50   > 5, "large_shift_gdd_p50",   NA_character_),
    flag_large_shift_gdd_rel50 = if_else(!is.na(max_abs_delta_rel50) & max_abs_delta_rel50 > 5, "large_shift_gdd_rel50", NA_character_),
    flag_large_shift_doy_p50   = if_else(!is.na(delta_p50_doy)   & abs(delta_p50_doy)   > 5, "large_shift_doy_p50",   NA_character_),
    flag_large_shift_doy_rel50 = if_else(!is.na(delta_rel50_doy) & abs(delta_rel50_doy) > 5, "large_shift_doy_rel50", NA_character_),
    flag_nonmonotone_p50   = if_else(!monotonic_p50,   "nonmonotone_p50",   NA_character_),
    flag_nonmonotone_rel50 = if_else(!monotonic_rel50, "nonmonotone_rel50", NA_character_)
  ) %>%
  unite(
    col = "flag",
    flag_large_shift_gdd_p50, flag_large_shift_gdd_rel50,
    flag_large_shift_doy_p50, flag_large_shift_doy_rel50,
    flag_nonmonotone_p50, flag_nonmonotone_rel50,
    sep = ";", na.rm = TRUE
  ) %>%
  mutate(flag = na_if(flag, ""))

# ── 5. Write ─────────────────────────────────────────────────────────────────
delta_file <- paste0('./output/tables/delta_table_', suff, '.csv')
delta_fileDE <- paste0('./output/tables/german/delta_table_', suff, '.csv')
write_csv(delta_gdd_all, delta_file, na = "")
write_csv2(delta_gdd_all, delta_fileDE, na = "")
