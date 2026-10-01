# Get command-line arguments
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check that at least one argument was provided
  if (length(args) < 5) {
    stop("Usage: Rscript example.R species start end datafile file_suffix")
  }
  
  # Use the first argument
  spec <- args[1]
  st <- as.numeric(args[2])
  en<- as.numeric(args[3])
  datafile <- args[4]
  suff <- args[5]
}


library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(lubridate)
library(mgcv)
library(qgam)
library(tidyverse)

options(warn = 1)        # print warnings immediately as they occur, not at the end
message("Logging started: ", Sys.time())

source("./hpc/R/fit_models.R")
source('./hpc/R/get_all_thresholds.R')
source("./hpc/R/boot_func.R")

if(interactive()) {
  datafile <- './data/wide_data.csv'
  suff <- 'test'
}

print(datafile)
initial_data <- 
  read_csv(datafile)

initial_data %>% ungroup() %>% distinct(species) %>% arrange(species) %>% print()
initial_data %>% group_by(species) %>% distinct(year) %>% summarize(n())
species_an <- function(sp, start = 1, end = NULL) {
  set.seed(42)
  sp_dat <- initial_data %>% 
    filter(grepl(sp, species)) %>%
    mutate(
      year = lubridate::year(date),
      tree_year = interaction(label, year, drop = TRUE))
  n_boot = 10000
  # ── Stage 1: define year-level bootstrap draws ───────────────────────────
  years_obs <- unique(sp_dat$year)
  n_years   <- length(years_obs)

  # Each replicate: resample n_years years with replacement
  boots_years <- map(seq_len(n_boot), ~ years_obs[sample(seq_along(years_obs), n_years, replace = TRUE)])
  
  # ── Stage 2: within each year draw, resample tree-years ─────────────────
  # For each replicate, iterate over the resampled year vector (with duplicates),
  # pull the tree-years available in that year, resample with replacement,
  # then row-bind across years
  boots_sp <- map(boots_years, function(year_draw) {
    map(year_draw, function(yr) {
      ty_obs <- unique(sp_dat$tree_year[sp_dat$year == yr])
      ty_obs[sample(seq_along(ty_obs), length(ty_obs), replace = TRUE)]
    })
  })

  if(is.null(end)) {
    end <- length(boots_sp)
  }
  
  test_list_sp <- map(start:end, function(b) { #length(boots_sp$splits)
    print(datafile)
    print(paste('Bootstrap:',b))
    print(boots_sp[[b]])
    boot_dat <- map_df(seq_along(boots_sp[[b]]), function(i) {
      yr      <- boots_years[[b]][i]
      ty_draw <- boots_sp[[b]][[i]]
      map_df(ty_draw, function(ty) sp_dat %>% filter(year == yr, tree_year == ty)) %>% mutate(i = i)
    })
    out <- boot_func(boot_dat, id = b, spec) %>% mutate(datafile = datafile)
    filename = paste0('./hpc/temp/bootseg_',sp,'_',en,'_',suff,'.csv')
    print(filename)
    if(file.exists(filename)){
      out %>% 
        write_csv(filename, append = T)
    } else {
      out %>% 
        write_csv(filename, append = T, col_names = T)
    }
    
    out
  })
}

if(interactive()){
  spec = 'Ulm'
  en = 4
  st = 4
}

species_an(spec, st, en)
