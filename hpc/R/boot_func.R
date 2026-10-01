boot_func <- function(df_sp, id = 'test', spec = 'Fag') {
  species_list <- unique(df_sp$species)
  results <- map_df(species_list, function(sp) {
    
    message("Processing species: ", sp)
    # ── Loop over models ─────────────────────────────────────────────────────
    diagnostic <- list()
    thresh_out <- map_df(names(fit_models), function(mod_name) {
      message("Processing model: ", sp, ".", mod_name)
      
      tryCatch({
        model <- fit_models[[mod_name]](df_sp)
        
        diagnostic[[mod_name]] <<- extract_diagnostic(model, mod_name, sp, id, datafile)
        
        if(grepl('omin', mod_name)) {
          # print('Overallmin')
          temp_min <- floor(quantile(df_sp$TempOverallMin, 0.005, na.rm=TRUE))
          temp_max <- ceiling(quantile(df_sp$TempOverallMin, 0.995, na.rm=TRUE))
          low_05_perc <- quantile(df_sp$TempOverallMin, 0.05)
          high_05_perc <- quantile(df_sp$TempOverallMin, 0.95)
          low_10_perc <- quantile(df_sp$TempOverallMin, 0.1)
          high_10_perc <- quantile(df_sp$TempOverallMin, 0.9)
          
        } else if(grepl('max', mod_name)) {
          # print('Max')
          temp_min <- floor(quantile(df_sp$TempMax, 0.005, na.rm=TRUE))
          temp_max <- ceiling(quantile(df_sp$TempMax, 0.995, na.rm=TRUE))
          low_05_perc <- quantile(df_sp$TempMax, 0.05)
          high_05_perc <- quantile(df_sp$TempMax, 0.95)
          low_10_perc <- quantile(df_sp$TempMax, 0.1)
          high_10_perc <- quantile(df_sp$TempMax, 0.9)
          
        } else {
          # print('Mean')
          temp_min <- floor(quantile(df_sp$Temperature, 0.005, na.rm=TRUE))
          temp_max <- ceiling(quantile(df_sp$Temperature, 0.995, na.rm=TRUE))
          low_05_perc <- quantile(df_sp$Temperature, 0.05)
          high_05_perc <- quantile(df_sp$Temperature, 0.95)
          low_10_perc <- quantile(df_sp$Temperature, 0.1)
          high_10_perc <- quantile(df_sp$Temperature, 0.9)
        }        
        temp_seq <- seq(
          temp_min,
          temp_max,
          by=0.1
        )
        
        
        doy_min <- min(df_sp$doy, na.rm=T)
        doy_max <- max(df_sp$doy, na.rm=T)
        
        doy_seq <- seq(doy_min, doy_max, by = 1)
        
        gdd_tb <- if (grepl("gdd", mod_name)) {
          sub(".*gddtb(\\d+).*", "tb\\1", mod_name)   # extracts "tb0" .. "tb6"
        } else {
          NULL
        }
        gdd_col <- if (!is.null(gdd_tb)) paste0("gdd_", gdd_tb) else NULL  # "gdd_tb0" .. "gdd_tb6"

        # ── newdata construction (replace existing expand.grid block) ────────────────
        if(!grepl('doy',mod_name)) {
          doy_seq <- median(df_sp$doy, na.rm=T)
        }
        newdata <- expand.grid(Temperature = temp_seq,
                               doy = doy_seq)
        
        newdata <- newdata %>% 
          mutate(TempNightMin = Temperature,
                 TempDayMin = Temperature,
                 TempOverallMin = Temperature,
                 TempMax = Temperature,
                 temp_mid = Temperature)
        
        if (!is.null(gdd_col)) {
          newdata[[gdd_col]] <- median(df_sp[[gdd_col]], na.rm = TRUE)
        }
        
        pred    <- predict(model, newdata = newdata, type = "response")
        
        pred_df <- 
          newdata %>% 
          as.data.frame() %>%
          select(Temperature, doy, any_of(paste0("gdd_tb", 0:6))) %>%
          mutate(fit = as.numeric(pred),
                 model_name = mod_name,
                 boot_id = id,
                 t_min = as.vector(temp_min),
                 t_q05 = as.vector(low_05_perc),
                 t_q10 = as.vector(low_10_perc),
                 t_q90 = as.vector(high_10_perc),
                 t_q95 = as.vector(high_05_perc),
                 t_max = as.vector(temp_max))
        
        pred_df <- pred_df %>% filter(is.finite(fit))

        if (grepl('doyti', mod_name)) {
          return(NULL)
        }
        
        pred_df_1d <- if (grepl('doy', mod_name)) {
          pred_df %>%
            group_by(Temperature) %>%
            summarise(fit = mean(fit, na.rm = TRUE), .groups = "drop")
        } else {
          pred_df %>% select(Temperature, fit)
        }
        
        pred_df_1d <- 
          pred_df_1d %>%
          mutate(model_name = mod_name,
                 boot_id = id,
                 t_min = as.vector(temp_min),
                 t_q05 = as.vector(low_05_perc),
                 t_q10 = as.vector(low_10_perc),
                 t_q90 = as.vector(high_10_perc),
                 t_q95 = as.vector(high_05_perc),
                 t_max = as.vector(temp_max))
        # Extract all thresholds
        thresholds <- get_all_thresholds(mod_name, pred_df_1d)
        
        thresholds %>% mutate(species = sp, boot_id = id, temp_min = temp_min, temp_max = temp_max)
        
      }, error = function(e) {
        message("  Failed (error): ", mod_name, 
                " | boot=", id,
                " | species=", sp,
                " — ", e$message)
        NULL
      }, warning = function(w) {
        message("  Failed (warning): ", mod_name,
                " | boot=", id, 
                " | species=", sp,
                " — ", w$message,
                " | n=", nrow(df_sp),
                " | NA_temp=", sum(is.na(df_sp$Temperature)),
                " | unique_temp=", length(unique(na.omit(df_sp$Temperature))))
      })
    })
    
    filename_diag = paste0('./hpc/temp/bootdiag_',spec,'_',en,'_',suff,'.csv')
    # print(filename_pred)
    
    if (file.exists(filename_diag)) {
      bind_rows(diagnostic) %>%
        write_csv(filename_diag, append = TRUE, col_names = FALSE)
    } else {
      bind_rows(diagnostic) %>%
        write_csv(filename_diag, col_names = TRUE)
    }
    
    thresh_out
  })
  print(results)
  results
}
