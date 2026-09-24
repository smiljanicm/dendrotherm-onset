library(tidyverse)
library(mgcv)
library(lubridate)

clean_data <- read_csv('./data/clean_data.csv', guess_max = Inf) 

clean_data <- 
  clean_data %>% 
  rename(Temperature = Temp_mean_GRO_period,
         TempMax = Temp_max_GRO_period,
         TempOverallMin = Temp_min_GRO_period)

clean_data <- 
  clean_data %>% 
  mutate(occ = occ_max) %>%
  mutate(doy = lubridate::yday(date))

# ── GDD columns (accumulated from DOY 31, Tb = 0–6°C) ───────────────────────
tb_seq <- 0:6

clean_data <- clean_data %>%
  filter(doy >= 31) %>%        # rows before DOY 31 get no GDD columns
  arrange(label, year, doy) %>%
  group_by(label, year)

for (tb in tb_seq) {
  clean_data <- clean_data %>%
    mutate(!!paste0("gdd_tb", tb) := cumsum(pmax(Temperature - tb, 0)))
}

clean_data <- clean_data %>% ungroup()

# ── Window datasets (GDD columns inherited automatically) ─────────────────────
window_data <- clean_data %>%
  filter(between(doy, 60, 172))

wide_data <- clean_data %>%
  filter(between(doy, 32, 172))

narrow_data <- clean_data %>%
  filter(between(doy, 91, 172))

write_csv(window_data, './data/gs_data.csv')

write_csv(wide_data, './data/wide_data.csv')

write_csv(narrow_data, './data/narrow_data.csv')
