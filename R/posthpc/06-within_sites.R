library(tidyverse)

# ── 1. Read compiled seg files ───────────────────────────────────────────────
seg_lad <- read_csv('./hpc/output/seg_Ladebow.csv') %>% mutate(site = 'Ladebow')
seg_dem <- read_csv('./hpc/output/seg_Demmin.csv')  %>% mutate(site = 'Demmin')

seg_all <- bind_rows(seg_lad, seg_dem)

# ── 2. Filter to primary models only, extract both metrics ──────────────────
primary_models <- c('gam_occ_mean_base_10', 'gam_occ_omin_base_10')

thresh <- seg_all %>%
  filter(model_name %in% primary_models) %>%
  mutate(
    metric = case_when(
      grepl('omin', model_name) ~ 'omin',
      grepl('max',  model_name) ~ 'max',
      TRUE                      ~ 'mean'
    ),
    # Flag boundary hits — threshold returned at temp_min or temp_max edge
    boundary_hit = (onset_rel50 == temp_min | onset_rel50 == temp_max | is.na(onset_rel50))
  ) %>%
  select(site, species, boot_id, metric, onset_rel50, obs_max_prob, boundary_hit)

# ── 3. Pivot wide by species, compute pairwise diffs per replicate ───────────
compute_pairwise_diffs <- function(thresh_site, site_name, species_list) {
  
  wide <- thresh_site %>%
    filter(site == site_name, species %in% species_list) %>%
    # Only use replicates where BOTH species have non-boundary fits
    select(boot_id, metric, species, onset_rel50, boundary_hit) %>%
    pivot_wider(
      names_from  = species,
      values_from = c(onset_rel50, boundary_hit)
    )
  
  # Generate all species pairs
  pairs <- combn(species_list, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    sp1 <- pair[1]; sp2 <- pair[2]
    col1 <- paste0('onset_rel50_', sp1)
    col2 <- paste0('onset_rel50_', sp2)
    bh1  <- paste0('boundary_hit_', sp1)
    bh2  <- paste0('boundary_hit_', sp2)
    
    wide %>%
      mutate(
        diff = .data[[col1]] - .data[[col2]],
        both_valid = !.data[[bh1]] & !.data[[bh2]]
      ) %>%
      filter(both_valid) %>%
      group_by(metric) %>%
      summarize(
        n_reps       = n(),
        diff_med     = median(diff, na.rm = TRUE),
        diff_q025    = quantile(diff, 0.025, na.rm = TRUE),
        diff_q975    = quantile(diff, 0.975, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(site = site_name, sp1 = sp1, sp2 = sp2)
  })
}

sps_lad <- c("Pinus sylvestris", "Quercus robur", "Betula pendula")
sps_dem <- c("Picea abies", "Larix decidua", "Pseudotsuga menziesii", "Fagus sylvatica", "Quercus robur")

diffs_lad <- compute_pairwise_diffs(thresh, 'Ladebow', sps_lad)
diffs_dem <- compute_pairwise_diffs(thresh, 'Demmin', sps_dem)

pairwise_diffs <- bind_rows(diffs_lad, diffs_dem)

# network_table2: species, metric, network_onset_rel50
# (you supply this — either hardcoded or read from your main thresholds table)
network_table2 <- read_csv('./output/tables/thresholds_full.csv') %>%
  filter(model_name %in% primary_models) %>%
  mutate(
    metric = case_when(
      grepl('omin', model_name) ~ 'omin',
      grepl('max',  model_name) ~ 'max',
      TRUE                      ~ 'mean'
    )) %>%
  filter(datarange == '60-172')

S3 <- pairwise_diffs %>%
  left_join(
    network_table2 %>% rename(sp1 = species, net_sp1 = onset_rel50_med),
    by = c('sp1', 'metric')
  ) %>%
  left_join(
    network_table2 %>% rename(sp2 = species, net_sp2 = onset_rel50_med),
    by = c('sp2', 'metric')
  ) %>%
  mutate(
    network_diff = net_sp1 - net_sp2,
    # Decision criterion: same sign AND network diff inside within-site CI
    same_sign   = sign(network_diff) == sign(diff_med),
    inside_ci   = network_diff >= diff_q025 & network_diff <= diff_q975,
    consistent  = same_sign & inside_ci
  ) %>%
  select(site, sp1, sp2, metric, n_reps,
         within_diff_med = diff_med, within_diff_q025 = diff_q025, within_diff_q975 = diff_q975,
         network_diff, consistent)

S3 %>% print(n=Inf)

write_csv(S3, './paper/tables/tableS4_thresholds_within_site_vs_network.csv')
