# Get command-line arguments
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check that at least one argument was provided
  if (length(args) < 5) {
    stop("Usage: Rscript withinsite_analysis.R analysis_target start end datafile file_suffix
         analysis_target must be one of: DataCheck, Ladebow, Demmin, PostProcess")
  }
  
  # Use the first argument as the analysis target
  analysis_target <- args[1] 
  st <- as.numeric(args[2])
  en <- as.numeric(args[3])
  datafile <- args[4]
  suff <- args[5]
} else {
  # Defaults for interactive testing
  analysis_target <- 'Ladebow'
  st <- 1
  en <- 1
  datafile <- './data/gs_within_site.csv'
  suff <- 'Ladebow'
}

library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(lubridate)
library(mgcv)
library(tidyverse)

options(warn = 1)        # print warnings immediately as they occur
message("Logging started: ", Sys.time())

# Only load heavy modeling functions if we are actually bootstrapping
source("./hpc/R/fit_models.R")
source('./hpc/R/get_all_thresholds.R')
source("./hpc/R/boot_func.R")

print(datafile)

# ── 1. Read & Pre-process Data ──────────────────────────────────────────
initial_data <- read_csv(datafile)

# Define Focal Species per Site
sps_lad <- c("Pinus sylvestris", "Quercus robur", "Betula pendula")
sps_dem <- c("Picea abies", "Larix decidua", "Pseudotsuga menziesii", "Fagus sylvatica", "Quercus robur")

# Helper to find common springs for a specific site and list of species
get_common_springs <- function(dat, target_site, target_sps) {
  sub_dat <- dat %>% 
    filter(site == target_site) %>%
    filter(species %in% target_sps)
  
  combos <- sub_dat %>% distinct(site, species)
  if(nrow(combos) == 0) return(numeric(0))
  
  sp_years <- lapply(1:nrow(combos), function(i) {
    unique(sub_dat$year[sub_dat$site == combos$site[i] & sub_dat$species == combos$species[i]])
  })
  
  Reduce(intersect, sp_years)
}

lad_common <- get_common_springs(initial_data, "Ladebow", sps_lad)
dem_common <- get_common_springs(initial_data, "Demmin", sps_dem)

# ── 2. Synchronous Bootstrap Execution ─────────────────────────────────
run_site_boot <- function(start, end, dat, analysis_name) {
  set.seed(42) # Ensure identical year draws across parallel nodes
  n_boot <- 2000
  
  if (analysis_name == "Ladebow") {
    target_site <- "Ladebow"
    target_sps <- sps_lad
    common_years <- lad_common
  } else if (analysis_name == "Demmin") {
    target_site <- "Demmin"
    target_sps <- sps_dem
    common_years <- dem_common
  } else {
    stop("Invalid analysis target for bootstrap.")
  }
  
  # Filter to exact site, species, and common years for this run
  sp_dat <- dat %>%
    filter(site == target_site) %>%
    filter(species %in% target_sps) %>%
    filter(year %in% common_years)

  years_obs <- common_years
  n_years <- length(years_obs)
  
  # STAGE 1: Single global year draw per replicate (synchronous for all species at the site)
  boots_years <- map(seq_len(n_boot), ~ years_obs[sample(seq_along(years_obs), n_years, replace = TRUE)])
  
  combos <- sp_dat %>% distinct(site, species)
  
  for (b in start:end) {
    print(paste('Bootstrap:', b, '| Analysis:', analysis_name))
    yr_draw <- boots_years[[b]]
    rep_results <- list()
    
    # STAGE 2: Within each year draw, resample tree-years per species
    for (i in seq_len(nrow(combos))) {
      curr_site <- combos$site[i]
      curr_sp <- combos$species[i]
      spec <- curr_sp
      curr_sp_dat <- sp_dat %>% filter(site == curr_site, species == curr_sp)
      
      boot_dat_sp <- map_df(seq_along(yr_draw), function(yi) {
        yr <- yr_draw[yi]
        ty_obs <- unique(curr_sp_dat$tree_year[curr_sp_dat$year == yr])
        if (length(ty_obs) == 0) return(NULL) 
        
        ty_draw <- ty_obs[sample(seq_along(ty_obs), length(ty_obs), replace = TRUE)]
        map_df(ty_draw, function(ty) curr_sp_dat %>% filter(year == yr, tree_year == ty)) %>% 
          mutate(boot_yr_idx = yi)
      })

      out <- boot_func(boot_dat_sp, id = b, spec = curr_sp) #%>%
#          mutate(site = curr_site, species = curr_sp, replicate = b)
      filename =  paste0('./hpc/temp/bootseg_', curr_sp, '_', en, '_', suff, '.csv')
      print(filename)
      if(file.exists(filename)){
        out %>% 
          #bind_rows() %>% 
          write_csv(filename, append = T)
      } else {
        out %>% #bind_rows() %>% 
          write_csv(filename, append = T, col_names = T)
      }

    }
  }
}

if(analysis_target %in% c("Ladebow", "Demmin")) {
  run_site_boot(st, en, initial_data, analysis_target)
}
