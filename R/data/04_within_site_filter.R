library(tidyverse)
library(dendrogreif)

con <- dendrogreif::con_init("../../../creds_local")

locs <- list_locations(con = con, site_ids = c(6,7))

gs_data <- read_csv('./data/gs_data.csv')

initial_data <- 
  gs_data %>% 
  filter(label %in% locs$label %>% unlist()) %>%
  left_join(locs %>% 
              filter(!is.na(label)) %>% 
              filter(height < 6) %>% 
              filter(grepl("Dendr", description)) %>% 
              select(site, label), by = join_by(label == label))

write_csv(initial_data, './data/gs_within_site.csv')
