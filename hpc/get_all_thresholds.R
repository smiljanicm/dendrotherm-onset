safe_threshold <- function(temp, fit, thr, fun = min) {
  idx <- fit >= thr
  if (!any(idx)) {
    return(NA_real_)
  }
  fun(temp[idx], na.rm = TRUE)
}

get_all_thresholds <- function(model_name, pred_df) {
  tryCatch({
    
    # Absolute thresholds (Rossi-comparable, max_fit = 1)
    abs_max <- 1.0
    
    # Relative thresholds (relative to observed maximum probability)
    obs_max <- max(pred_df$fit, na.rm = TRUE)
    
    tibble(
      model_name = model_name,
      # Data range metadata
      t_q05 = median(pred_df$t_q05, na.rm = TRUE),
      t_q10 = median(pred_df$t_q10, na.rm = TRUE),
      t_q90 = median(pred_df$t_q90, na.rm = TRUE),
      t_q95 = median(pred_df$t_q95, na.rm = TRUE),
      # Observed maximum fitted probability (useful diagnostic)
      obs_max_prob = obs_max,
      
      # ── Absolute thresholds (directly comparable to Rossi) ───────────
      onset_p50 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.50 * abs_max, min),
      onset_p25 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.25 * abs_max, min),
      onset_p10 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.10 * abs_max, min),
      onset_p05 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.05 * abs_max, min),
      
      # ── Relative thresholds (relative to observed max probability) ────
      onset_rel50 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.50 * obs_max, min),
      onset_rel25 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.25 * obs_max, min),
      onset_rel10 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.10 * obs_max, min),
      onset_rel05 = safe_threshold(pred_df$Temperature, pred_df$fit, 0.05 * obs_max, min)
    )
    
  }, error = function(e) {
    message("  Threshold extraction failed: ", e$message)
    tibble(model_name = model_name,
           obs_max_prob = NA,
           onset_p50 = NA, onset_p25 = NA, onset_p10 = NA, onset_p05 = NA,
           onset_rel50 = NA, onset_rel25 = NA, onset_rel10 = NA, onset_rel05 = NA)
  })
}

extract_diagnostic <- function(model, model_name, species, boot_id, datafile) {
  
  # ── AIC / BIC ─────────────────────────────────────────────────────────────
  aic_bic <- tibble(
    datafile = datafile,
    model_name = model_name,
    species    = species,
    boot_id    = boot_id,
    term       = "model",
    metric     = c("AIC", "BIC"),
    value      = c(AIC(model), BIC(model))
  )
  
  # ── k.check ───────────────────────────────────────────────────────────────
  kc <- tryCatch({
    k.check(model) %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      pivot_longer(-term, names_to = "metric") %>%
      mutate(
        datafile = datafile,
        model_name = model_name,
        species    = species,
        boot_id    = boot_id
      )
  }, error = function(e) {
    message("  k.check failed: ", e$message)
    tibble(term = NA, metric = "k.check", value = NA,
           model_name = model_name, species = species, boot_id = boot_id)
  })
  
  # ── Concurvity ────────────────────────────────────────────────────────────
  cc <- tryCatch({
    concurvity(model, full = FALSE) %>%
      map_df(function(mat) {
        mat %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          pivot_longer(-term, names_to = "by_term")
      }, .id = "metric") %>%
      mutate(
        datafile = datafile,
        metric     = paste(metric, by_term, sep = ":"),
        model_name = model_name,
        species    = species,
        boot_id    = boot_id
      ) %>%
      select(datafile, model_name, species, boot_id, term, metric, value)
  }, error = function(e) {
    message("  concurvity failed: ", e$message)
    tibble(term = NA, metric = "concurvity", value = NA,
           model_name = model_name, species = species, boot_id = boot_id)
  })
  
  # ── Combine ───────────────────────────────────────────────────────────────
  bind_rows(aic_bic, kc, cc)
}