library(tidyverse)

dat <- 
  read_csv('./data/initial_data.csv') %>% 
  filter(is.finite(GRO_max_rate)) %>%
  arrange(gro_window_date) %>%
  group_by(label.y, year=lubridate::year(gro_window_date)) %>%
  mutate(GRO_max_sum = cumsum(GRO_max_rate)) %>%
  mutate(above_Q90_max = ifelse(GRO_max_rate >= quantile(GRO_max_rate, 0.90, na.rm = TRUE), 1, 0),
         occ_max = as.numeric(GRO_max_rate > 0)) %>% 
  filter(year > 2013)

dat_meta <- dat %>% mutate(year = lubridate::year(gro_window_date)) %>% ungroup() %>%
  group_by(species.y, label.y, year) %>% 
  summarise(number_of_days_observed=n(),
            number_of_growth_days_above_zero = sum(GRO_max_rate>0, na.rm=T),
            ratio_of_growth_days_above_zero = mean(GRO_max_rate>0, na.rm=T),
            min_date = min(gro_window_date, na.rm=T),
            min_doy = min(lubridate::yday(min_date), na.rm=T),
            max_date = max(gro_window_date, na.rm=T),
            max_doy = max(lubridate::yday(max_date), na.rm=T))

dat_meta %>%
  ungroup() %>%
  # nrow()
  group_by(species.y) %>%
  summarize(tree_years_number = n()) %>%
  arrange(tree_years_number) %>%
  print()


tys <- dat_meta %>%
  filter(number_of_days_observed > 60) %>%
  filter(ratio_of_growth_days_above_zero > 0.1) %>% 
  filter(min_doy < 185)  %>%
  #  nrow()
  # ungroup() %>%
  # group_by(species) %>%
  # summarize(n = n()) %>%
  # arrange(n)
  ungroup() %>%
  mutate(ty = paste(label.y, year, sep=".")) %>%
  select(ty) %>% unlist()

clean_data <- 
  dat %>% filter(paste(label.y, year, sep=".") %in% tys) %>%
  mutate(tree_year = interaction(label.y, year, drop = TRUE))

cldat_meta <- 
  clean_data  %>% mutate(year = lubridate::year(gro_window_date)) %>%
  group_by(species.y, label.y, year) %>% 
  summarise(number_of_days_observed=n(),
            number_of_growth_days_above_zero = sum(GRO_max_rate>0, na.rm=T),
            ratio_of_growth_days_above_zero = mean(GRO_max_rate>0, na.rm=T),
            min_date = min(gro_window_date, na.rm=T),
            min_doy = min(lubridate::yday(min_date), na.rm=T),
            max_date = max(gro_window_date, na.rm=T),
            max_doy = max(lubridate::yday(max_date), na.rm=T))

cleaning_effect <- 
  dat_meta %>%
  ungroup() %>%
  group_by(species.y) %>%
  summarize(n_before = n(),
            min_before = min(year)) %>%
  left_join(cldat_meta %>%
              ungroup() %>%
              group_by(species.y) %>%
              summarize(n_after = n(),
                        min_after = min(year))) %>%
  mutate(diff = n_after - n_before) %>%
  arrange(n_after) %>%
  mutate(ratio = diff/n_before)

write_csv(cleaning_effect, file="./output/tables/cleaning_effect.csv")
write_csv2(cleaning_effect, file="./output/tables/german/cleaning_effect.csv")

write_csv(clean_data %>%
            rename(species = species.y,
                   label = label.y,
                   date = gro_window_date) %>%
            mutate(GRO = GRO_max_rate), 
          './data/clean_data.csv')
