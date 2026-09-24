library(tidyverse)

gs_data <- read_csv('./data/gs_data.csv')

gdmm <- gs_data %>% 
  mutate(GRO_nonzero = case_when(GRO_max_rate==0 ~ NA_real_,
                                 T ~ GRO_max_rate)) %>%
  group_by(species) %>%
  summarize(mean_gro_days = mean(GRO_max_rate>0, na.rm=T),
            mean_gro_magnitude = mean(GRO_nonzero, na.rm=T)) %>%
  left_join(
    gs_data %>% 
      group_by(species) %>%
      distinct(year,tree_year) %>%
      summarize(tree_years_num = n(),
                min_year = min(year),
                max_year = max(year)),
    by = join_by(species)
  ) %>%
  left_join(gs_data %>%
              group_by(species) %>%
              distinct(label) %>%
              summarize(number_of_trees = n())
  )

gdmm %>%
  write_csv("./output/tables/summary_stats.csv")

gdmm %>%
  write_csv2("./output/tables/german/summary_stats.csv")
