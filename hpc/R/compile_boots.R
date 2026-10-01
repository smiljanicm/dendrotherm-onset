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
  suff <- 'Demmin'
}

dir_path <- './hpc/temp'
diag_path <- dir_path


diag_file <- paste0('./hpc/output/diag_', suff, '.csv')
seg_file <- paste0('./hpc/output/seg_', suff, '.csv')

print(diag_file)
print(seg_file)

diag_pattern <- paste0('^bootdiag_.*_', suff, '\\.csv$')
seg_pattern  <- paste0('^bootseg_.*_',  suff, '\\.csv$')

archive_path <- file.path(dir_path, "compiled_archive")
if (!dir.exists(archive_path)) dir.create(archive_path)

walk(list.files(path = dir_path, pattern=diag_pattern, full.names=T), 
     function(x) { 
       print(x)
       chunk <- read_csv(x, guess_max = Inf, show_col_types = F)
       if(file.exists(diag_file)) {
         chunk %>% write_csv(diag_file, append=T)
       } else {
         chunk %>% write_csv(diag_file)
       }
       
       gc(verbose=F)
       file.rename(x, file.path(archive_path, basename(x)))
       NULL
     }, .progress = T)

seg_path <- dir_path

walk(list.files(path = dir_path, pattern=seg_pattern, full.names=T), 
     function(x) { 
       print(x)
       chunk <- read_csv(x, guess_max = Inf, show_col_types = F)
       if(file.exists(seg_file)) {
         chunk %>% write_csv(seg_file, append=T)
       } else {
         chunk %>% write_csv(seg_file)
       }
       gc(verbose=F)
       file.rename(x, file.path(archive_path, basename(x)))
       NULL
     }, .progress = T) 
