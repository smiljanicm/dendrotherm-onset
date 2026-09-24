library(tidyverse)

gs <- read_csv('./data/gs_data.csv')

# ── Parameters ───────────────────────────────────────────────────────────────
thaw_window <- 1  # days after sub-zero minimum to flag as freeze-thaw risk

# ── Classify days ────────────────────────────────────────────────────────────
frost_check <- gs %>% 
  arrange(species, tree_year, date) %>%
  group_by(species, tree_year) %>%
  mutate(
    sub_zero    = TempOverallMin < 0,
    freeze_thaw = map_lgl(row_number(), function(i) {
      any(sub_zero[pmax(1, i - thaw_window):(i - 1)], na.rm = TRUE)
    })
  ) %>%
  ungroup() %>%
  mutate(
    day_class = case_when(
      sub_zero                   ~ "sub_zero",
      freeze_thaw                ~ "freeze_thaw",
      TRUE                       ~ "clean"
    ),
    growth = as.integer(GRO > 0)
  )

# ── Summarise per species x day_class ────────────────────────────────────────
frost_summary <- frost_check %>%
  group_by(species, day_class) %>%
  summarise(
    n_days       = n(),
    n_growth     = sum(growth, na.rm = TRUE),
    growth_rate  = n_growth / n_days,
    .groups      = "drop"
  ) %>%
  pivot_wider(
    names_from  = day_class,
    values_from = c(n_days, n_growth, growth_rate)
  ) %>%
  mutate(
    # Rate ratio: freeze-thaw vs clean — key diagnostic
    rr_freeze_thaw = growth_rate_freeze_thaw / growth_rate_clean,
    rr_sub_zero    = growth_rate_sub_zero    / growth_rate_clean
  ) %>%
  arrange(desc(rr_freeze_thaw))

print(frost_summary, n = Inf)

# ── Fisher's exact test per species: freeze-thaw vs clean ────────────────────
frost_tests <- frost_check %>%
  filter(day_class %in% c("clean", "freeze_thaw")) %>%
  group_by(species) %>%
  summarise(
    p_fisher = fisher.test(
      matrix(c(
        sum(growth[day_class == "freeze_thaw"]),
        sum(day_class == "freeze_thaw") - sum(growth[day_class == "freeze_thaw"]),
        sum(growth[day_class == "clean"]),
        sum(day_class == "clean")       - sum(growth[day_class == "clean"])
      ), nrow = 2)
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_fisher, method = "BH"))

frost_summary <- frost_summary %>%
  left_join(frost_tests, by = "species")

frost_summary %>% 
  select(species, rr_freeze_thaw, rr_sub_zero, p_fisher, p_adj) %>%
  write_csv("./output/tables/freeze_thaw.csv")

frost_summary %>% 
  select(species, rr_freeze_thaw, rr_sub_zero, p_fisher, p_adj) %>%
  write_csv2("./output/tables/german/freeze_thaw_german.csv")
