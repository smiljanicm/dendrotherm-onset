fit_models <- list(
  # ── Base Models across k = 20 and k = 10 ──────────────────────────────────
  gam_occ_mean_base_20 = function(df) {
    gam(occ ~ s(Temperature, k = 20), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_base_20 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 20), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_base_20 = function(df) {
    gam(occ ~ s(TempMax, k = 20), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_mean_base_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_base_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_base_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── DOY-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_doy_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(doy, k = 5, bs = "cc"), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_doy_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(doy, k = 5, bs = "cc"), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_doy_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(doy, k = 5, bs = "cc"), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_0-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb0_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb0, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb0_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb0, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb0_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb0, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_1-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb1_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb1, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb1_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb1, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb1_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb1, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_2-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb2_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb2, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb2_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb2, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb2_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb2, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_3-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb3_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb3, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb3_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb3, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb3_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb3, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_4-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb4_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb4, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb4_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb4, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb4_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb4, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_5-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb5_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb5, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb5_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb5, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb5_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb5, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  # ── GDD_6-adjusted triplet (ΔT supplement) ──────────────────────────────────
  gam_occ_mean_gddtb6_10 = function(df) {
    gam(occ ~ s(Temperature, k = 10) + s(gdd_tb6, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_omin_gddtb6_10 = function(df) {
    gam(occ ~ s(TempOverallMin, k = 10) + s(gdd_tb6, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  },
  gam_occ_max_gddtb6_10 = function(df) {
    gam(occ ~ s(TempMax, k = 10) + s(gdd_tb6, k = 5), data = df, family = binomial(link = "logit"), method = "REML")
  }
)
