library(tidyverse)
library(mgcv)

wide_data <- read_csv('./data/wide_data.csv')

# ── 1. Identify species with >6 distinct years ──────────────────────────────
species_years <- wide_data %>%
  group_by(species) %>%
  summarize(n_years = n_distinct(year), .groups = "drop")

eligible_species <- species_years %>%
  filter(n_years > 6) %>%
  pull(species)

message(Sys.time(), " | Eligible species (>6 years): ", length(eligible_species))
print(species_years %>% filter(species %in% eligible_species))

# ── 2. Model formula builder — both univariate and additive specs ──────────
build_formula <- function(model_type, tb = NULL) {
  full_term <- case_when(
    # univariate — single term alone
    model_type == "base"       ~ "s(Temperature, k = 10)",
    model_type == "doy"        ~ "s(doy, k = 5, bs = 'cc')",
    model_type == "gdd"        ~ paste0("s(gdd_tb", tb, ", k = 10)"),
    # additive — T plus confound
    model_type == "t_plus_doy" ~ "s(Temperature, k = 10) + s(doy, k = 5, bs = 'cc')",
    model_type == "t_plus_gdd" ~ paste0("s(Temperature, k = 10) + s(gdd_tb", tb, ", k = 5)")
  )
  as.formula(paste("occ ~", full_term))
}

# ── 3. LOYO-CV for one species × one model spec, with per-fold timing ───────
loyo_brier <- function(sp_data, model_type, tb = NULL, sp_name = "") {
  years <- sort(unique(sp_data$year))
  fmla  <- build_formula(model_type, tb)
  spec_label <- if (model_type %in% c("gdd", "t_plus_gdd")) paste0(model_type, "_tb", tb) else model_type
  
  briers <- map_dbl(years, function(y_test) {
    t0 <- Sys.time()
    
    train <- filter(sp_data, year != y_test)
    test  <- filter(sp_data, year == y_test)
    
    model <- tryCatch(
      gam(fmla, data = train, family = binomial(link = "logit"), method = "REML"),
      error = function(e) NULL
    )
    
    if (is.null(model)) {
      message(sprintf("  [%s | %s | year=%s] FAILED FIT (%.1fs)",
                      sp_name, spec_label, y_test, as.numeric(Sys.time() - t0, units = "secs")))
      return(NA_real_)
    }
    
    pred <- tryCatch(
      predict(model, newdata = test, type = "response"),
      error = function(e) rep(NA_real_, nrow(test))
    )
    brier <- mean((pred - test$occ)^2, na.rm = TRUE)
    
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    message(sprintf("  [%s | %s | year=%s] brier=%.4f (%.1fs)",
                    sp_name, spec_label, y_test, brier, elapsed))
    
    brier
  })
  
  tibble(year = years, brier = briers)
}

# ── 4. Run across all eligible species × model types × Tb sweep ────────────
model_specs <- bind_rows(
  tibble(model_type = "base",       tb = NA_character_),
  tibble(model_type = "doy",        tb = NA_character_),
  tibble(model_type = "gdd",        tb = as.character(0:6)),
  tibble(model_type = "t_plus_doy", tb = NA_character_),
  tibble(model_type = "t_plus_gdd", tb = as.character(0:6))
)

total_specs <- length(eligible_species) * nrow(model_specs)
spec_counter <- 0
overall_t0 <- Sys.time()

loyo_results <- map_dfr(eligible_species, function(sp) {
  sp_data <- wide_data %>% filter(species == sp)
  n_years_sp <- n_distinct(sp_data$year)
  
  message(sprintf("\n=== %s | species: %s (%d years) ===", Sys.time(), sp, n_years_sp))
  
  map_dfr(seq_len(nrow(model_specs)), function(i) {
    spec_counter <<- spec_counter + 1
    model_type <- model_specs$model_type[i]
    tb         <- model_specs$tb[i]
    
    spec_t0 <- Sys.time()
    result <- loyo_brier(sp_data, model_type, tb, sp_name = sp)
    spec_elapsed <- as.numeric(Sys.time() - spec_t0, units = "secs")
    
    total_elapsed <- as.numeric(Sys.time() - overall_t0, units = "secs")
    message(sprintf(">>> [%d/%d specs done] this spec: %.1fs | total elapsed: %.1f min",
                    spec_counter, total_specs, spec_elapsed, total_elapsed / 60))
    
    result %>% mutate(species = sp, model_type = model_type, tb = tb)
  })
})

message(sprintf("\nDONE. Total time: %.1f min", as.numeric(Sys.time() - overall_t0, units = "mins")))

loyo_file   <- './output/tables/loyo_brier_full.csv'
loyo_fileDE <- './output/tables/german/loyo_brier_full.csv'

write_csv(loyo_results, loyo_file)
write_csv2(loyo_results, loyo_fileDE)

# ── 5. Summarize: mean Brier per species × model ────────────────────────────
loyo_summary <- loyo_results %>%
  group_by(species, model_type, tb) %>%
  summarize(
    n_years_used = sum(!is.na(brier)),
    mean_brier   = mean(brier, na.rm = TRUE),
    .groups = "drop"
  )

summary_file   <- './output/tables/loyo_brier_summary.csv'
summary_fileDE <- './output/tables/german/loyo_brier_summary.csv'
write_csv(loyo_summary, summary_file)
write_csv2(loyo_summary, summary_fileDE)
