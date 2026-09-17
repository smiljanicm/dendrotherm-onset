library(dendrogreif)
library(tidyverse)
library(treenetproc)
library(boot)

set.seed(42)
treenet <- function(loc, timestamp, datacol, what='twd') {
  df <- data.frame(series = loc, ts=timestamp,value=datacol) 
  try({
    df1 <- treenetproc:::calcmax(df)
    df2 <- treenetproc:::calctwdgro(df1)
  })
  if(what == 'twd') {
    return(df2$twd)
  } else if(what == 'max'){
    return(df2$max)
  } else if(what == 'gro_yr'){
    return(df2$gro_yr)
  }
}

con <- dendrogreif::con_init("../../../creds_local")

get_site_data <- function(site_id = 6, temp_ref = 255) {
  
  locations <- list_locations(con, site_ids=site_id)
  locs <-
    locations %>%
    filter(description == 'Dendrometer' | location_id == temp_ref) %>%
    filter(location_id != 278) %>%
    filter(location_id != 142) %>%
    filter(location_id != 125) %>%
    filter(location_id != 129) %>%
    filter(location_id != 132) %>%
    filter(location_id != 136) %>%
    dplyr::select(location_id) %>%
    unlist()
  vars = c(1,2)
  print(locs)
  if(!(temp_ref %in% locs)) {
    locs <- c(locs, temp_ref)
  }
  print(locs)
  all_data <- get_data(con, locs, vars,
                       start = '2013-01-01',
                       end = '2026-12-31',
                       toclean="clean")
  
  all_data <- all_data %>% filter(variable_id == 1) %>%
    group_by(location_id) %>%
    #  filter(all(!is.na(value))) %>%
    mutate(twd = treenet(location_id, timestamp, value)) %>%
    mutate(max = treenet(location_id, timestamp, value, 'max')) %>%
    mutate(gro = treenet(location_id, timestamp, value, 'gro_yr')) %>%
    dplyr::select(timestamp, gro, label, species, location_description, height, variable, variable_id, location_id) %>%
    rename("value" = gro) %>%
    mutate(variable = 'GRO',
           variable_id = 1000) %>%
    mutate(value = value - quantile(value, 0.025, na.rm=T)) %>%
    bind_rows(all_data)
  
  return(all_data)
  
}

thresh_window <- function(tree_ref, temp_ref = 255, all_data) {
  all_data$timestamp <- force_tz(all_data$timestamp, 'Etc/GMT-1')
  
  temp_data <-
    all_data %>%
    filter(location_id %in% temp_ref) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(date = as.Date(timestamp, tz=tz(timestamp))) %>%
    filter(Temperature < 37,
           Temperature > -25) %>%
    mutate(timestamp = force_tz(timestamp, tzone='Etc/GMT-1'))
  
  tree_data <- 
    all_data %>%
    filter(location_id %in% tree_ref) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(
      date = as.Date(timestamp, tz=tz(timestamp)),
      time = hms::as_hms(timestamp),
      gro_window_date = if_else(
        time > hms::as_hms("14:00:00"),
        date + 1,
        date
      )
    ) %>%
    mutate(ifelse(GRO >= 0, GRO, 0)) %>%
    group_by(species, label, gro_window_date) %>%
    summarise(
      GRO_max = max(GRO, na.rm = TRUE),
      GRO_max_ts = timestamp[match(GRO_max, GRO)],
      GRO_mean = mean(GRO, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = lubridate::year(gro_window_date)) %>%
    arrange(species, label, gro_window_date) %>%
    group_by(species, label, year) %>%
    mutate(GRO_mean_rate = GRO_mean - lag(GRO_mean)) %>%
    mutate(GRO_max_rate = GRO_max - lag(GRO_max)) %>%
    mutate(GRO_max_ts_start = GRO_max_ts - 60*60*24) %>%
    mutate(GRO_mean_rate = ifelse(GRO_mean_rate >= 0, GRO_mean_rate, 0)) %>%
    mutate(GRO_max_rate = ifelse(GRO_max_rate >= 0, GRO_max_rate, 0))
  
  tt_data <-
    temp_data %>%
    full_join(tree_data, by=join_by(timestamp < GRO_max_ts,
                                    timestamp >= GRO_max_ts_start)) %>%
    ungroup() %>%
    arrange(gro_window_date, timestamp)
  
  tt_data %>%
    filter(timestamp >= GRO_max_ts_start,
           timestamp <= GRO_max_ts) %>%
    group_by(species.y, label.y, gro_window_date, GRO_max_ts, GRO_max, GRO_mean, GRO_max_rate, GRO_mean_rate) %>%
    summarise(
      Temp_mean_GRO_period = mean(Temperature, na.rm = TRUE),
      Temp_max_GRO_period = max(Temperature, na.rm = TRUE),
      Temp_min_GRO_period = min(Temperature, na.rm=T),
      .groups = "drop"
    ) 
}

site_thresh <- function(all_data, temp_ref, thresh_func = thresh) {
  locs_dms <- all_data %>% filter(!is.na(label)) %>% distinct(location_id)
  
  tt <- 
    locs_dms %>% unlist() %>%
    map(thresh_func, temp_ref = temp_ref, all_data=all_data) %>% 
    bind_rows()
  
  tt  
}

sites <- list_sites(con)
site_ids = c(2,4,5,6,7,8,9,10,11,12,13,14,15,16,19,20,21,27,28,29,30,22,3,32,41,42)
temps = c(237,346,471,255,248,266,267,268,263,265,264,269,270,271,249,382,423,400,363,436,459,271,237,551,602,605)

data_map <- map2(site_ids, temps, get_site_data)
site_names <- sites[match(site_ids, sites$site_id), 'name'] %>% unlist()
names(data_map) <- site_names

thresh_map_window <- map2(data_map, temps, site_thresh, thresh_func = thresh_window)
thresh_df_window <- thresh_map_window %>%
  bind_rows()

write_csv(thresh_df_window, './data/initial_data.csv')
