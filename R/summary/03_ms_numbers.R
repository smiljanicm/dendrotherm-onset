# ms_numbers.R
#
# Port of ms_numbers.py.
# Run from the project root:
# source("code/ms_numbers.R")
#
# Inputs:  data/output/
# Outputs: results/ms_numbers_report.txt
#          results/table3_cv.md
library(tidyverse)

run_ms_numbers <- function() {
  D <- file.path("output", "tables")
  RES <- "paper/tables"
  dir.create(RES, recursive = TRUE, showWarnings = FALSE)
  
  read_input <- function(filename) {
    read.csv(
      file.path(D, filename),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  ANAT <- c(
    "Abies alba" = "tracheid",
    "Picea abies" = "tracheid",
    "Pinus sylvestris" = "tracheid",
    "Larix decidua" = "tracheid",
    "Pseudotsuga menziesii" = "tracheid",
    "Quercus robur" = "ring",
    "Quercus rubra" = "ring",
    "Ulmus laevis" = "ring",
    "Fagus sylvatica" = "diffuse",
    "Carpinus betulus" = "diffuse",
    "Betula pendula" = "diffuse",
    "Acer pseudoplatanus" = "diffuse",
    "Alnus glutinosa" = "diffuse",
    "Tilia cordata" = "diffuse"
  )
  
  M <- c(
    min = "gam_occ_omin_base_10",
    mean = "gam_occ_mean_base_10",
    max = "gam_occ_max_base_10"
  )
  
  SIX <- c(
    "Quercus robur",
    "Fagus sylvatica",
    "Carpinus betulus",
    "Acer pseudoplatanus",
    "Pinus sylvestris",
    "Betula pendula"
  )
  
  out <- character()
  
  P <- function(x) {
    out <<- c(out, as.character(x))
  }
  
  show_table <- function(x, digits = NULL) {
    if (!is.null(digits)) {
      numeric_cols <- vapply(x, is.numeric, logical(1))
      x[numeric_cols] <- lapply(x[numeric_cols], round, digits)
    }
    P(capture.output(print(x, width = 200)))
  }
  
  t <- read_input("thresholds_full.csv")
  tr <- read_input("thresholds_restricted.csv")
  gd <- read_input("summary_stats.csv")
  
  stopifnot(!anyDuplicated(gd$species))
  rownames(gd) <- gd$species
  
  sel <- function(df, m, dr = "60-172") {
    x <- df[
      which(df$model_name == M[[m]] & df$datarange == dr),
      ,
      drop = FALSE
    ]
    stopifnot(!anyDuplicated(x$species))
    rownames(x) <- x$species
    x
  }
  
  # ------------------------------------------------------------
  # Table 1
  # ------------------------------------------------------------
  
  P("=== Table 1: period, springs, dendrometer-years (tree_years_num), number of trees, dbh range")
  
  # Python sorts anatomy labels alphabetically for Table 1.
  species_order <- names(ANAT)[order(unname(ANAT), names(ANAT))]
  dbh_df <- data.frame(species = c(
  "Larix decidua",
  "Abies alba",
  "Picea abies",
  "Pinus sylvestris",
  "Pseudotsuga menziesii",
  "Quercus robur",
  "Quercus rubra",
  "Ulmus laevis",
  "Fagus sylvatica",
  "Carpinus betulus",
  "Acer pseudoplatanus",
  "Tilia cordata",
  "Betula pendula",
  "Alnus glutinosa"
  ),
  dbh = c(
  "38.4 ± 9.5",   
  "73.0 ± 6.5",
  "39.8 ± 11.4",
  "43.3 ± 6.8",
  "51.1 ± 8.4",
  "59.7 ± 22.2",
  "39.8 ± 5.7",
  "49.4 ± 9.4",
  "58.7 ± 27.6",
  "41.8 ± 11.0",
  "51.9 ± 8.7",
  "48.9 ± 7.5",
  "36.6 ± 7.6",
  "32.0 ± 4.5"
  ))
  
  gd <- gd %>% left_join(dbh_df) %>% column_to_rownames("species")
  for (sp in species_order) {
    r <- gd[sp, , drop = FALSE]
    P(sprintf(
      "%-24s %d-%d  springs %2d  dy %3d trees %2d %-24s",
      sp,
      as.integer(r$min_year),
      as.integer(r$max_year),
      as.integer(r$max_year - r$min_year + 1),
      as.integer(r$tree_years_num),
      as.integer(r$number_of_trees),
      r$dbh
    ))
  }
  
  P(sprintf(
    "total dendrometer-years %d",
    as.integer(sum(gd$tree_years_num, na.rm = TRUE))
  ))
  
  # ------------------------------------------------------------
  # Table 2
  # ------------------------------------------------------------
  
  P(paste0(
    "\n=== Table 2: P = 0.50, DOY 60-172, full period: ",
    "median [95 % CI]; attainment; bootstrap replicates attained"
  ))
  
  T <- setNames(lapply(names(M), function(m) sel(t, m)), names(M))
  
  species_order <- names(ANAT)[order(
    match(unname(ANAT), c("tracheid", "ring", "diffuse")),
    T$mean[names(ANAT), "onset_rel50_med"]
  )]
  
  for (sp in species_order) {
    cells <- vapply(names(M), function(m) {
      r <- T[[m]][sp, , drop = FALSE]
      sprintf(
        "%.1f [%.1f, %.1f]",
        r$onset_rel50_med,
        r$onset_rel50_low,
        r$onset_rel50_high
      )
    }, character(1))
    
    n_boot <- vapply(names(M), function(m) {
      sprintf("%d", as.integer(T[[m]][sp, "n"]))
    }, character(1))
    
    P(paste0(
      sprintf("%-24s ", sp),
      paste(cells, collapse = " | "),
      sprintf(" | att(min) %.2f", T$min[sp, "obs_max_prob"]),
      sprintf(" | att(mean) %.2f", T$mean[sp, "obs_max_prob"]),
      sprintf(" | att(max) %.2f", T$max[sp, "obs_max_prob"]),
      " | n_boot ",
      paste(n_boot, collapse = ",")
    ))
  }
  
  P(sprintf(
    "§2.7 5th percentile of daily mean T: %.2f-%.2f °C",
    min(T$mean$t_q05, na.rm = TRUE),
    max(T$mean$t_q05, na.rm = TRUE)
  ))
  
  P(sprintf(
    "Tilia daily max: t_q95 %.1f °C",
    T$max["Tilia cordata", "t_q95"]
  ))
  
  # ------------------------------------------------------------
  # Table S2
  # ------------------------------------------------------------
  
  P(paste0(
    "\n=== Table S2: P = 0.50, DOY 60-172, full period: ",
    "median [95 % CI]; attainment; bootstrap replicates attained"
  ))
  
  T <- setNames(lapply(names(M), function(m) sel(t, m)), names(M))
  
  species_order <- names(ANAT)[order(
    match(unname(ANAT), c("tracheid", "ring", "diffuse")),
    T$mean[names(ANAT), "onset_p50_med"]
  )]
  
  for (sp in species_order) {
    cells <- vapply(names(M), function(m) {
      r <- T[[m]][sp, , drop = FALSE]
      sprintf(
        "%.1f [%.1f, %.1f]",
        r$onset_p50_med,
        r$onset_p50_low,
        r$onset_p50_high
      )
    }, character(1))
    
    n_boot <- vapply(names(M), function(m) {
      sprintf("%d", as.integer(T[[m]][sp, "n"]))
    }, character(1))
    
    P(paste0(
      sprintf("%-24s ", sp),
      paste(cells, collapse = " | "),
      sprintf(" | att(min) %.2f", T$min[sp, "obs_max_prob"]),
      sprintf(" | att(mean) %.2f", T$mean[sp, "obs_max_prob"]),
      sprintf(" | att(max) %.2f", T$max[sp, "obs_max_prob"]),
      " | n_boot ",
      paste(n_boot, collapse = ",")
    ))
  }
  
  P(sprintf(
    "§2.7 5th percentile of daily mean T: %.2f-%.2f °C",
    min(T$mean$t_q05, na.rm = TRUE),
    max(T$mean$t_q05, na.rm = TRUE)
  ))
  
  P(sprintf(
    "Tilia daily max: t_q95 %.1f °C",
    T$max["Tilia cordata", "t_q95"]
  ))
  
  
  # ------------------------------------------------------------
  # Frost
  # ------------------------------------------------------------
  
  P("\n=== §3.5 frost")
  
  f <- read_input("freeze_thaw.csv")
  
  P(sprintf(
    paste0(
      "RR freeze-thaw %.2f-%.2f (median %.2f); ",
      "sub-zero %.2f-%.2f (median %.2f); max adj p %.2g"
    ),
    min(f$rr_freeze_thaw, na.rm = TRUE),
    max(f$rr_freeze_thaw, na.rm = TRUE),
    median(f$rr_freeze_thaw, na.rm = TRUE),
    min(f$rr_sub_zero, na.rm = TRUE),
    max(f$rr_sub_zero, na.rm = TRUE),
    median(f$rr_sub_zero, na.rm = TRUE),
    max(f$p_adj, na.rm = TRUE)
  ))
  
  # ------------------------------------------------------------
  # Window sensitivity
  # ------------------------------------------------------------
  
  P("\n=== §3.5 window sensitivity (range across DOY 32/60/91-172)")
  
  for (m in names(M)) {
    x <- t[
      which(t$model_name == M[[m]]),
      c("species", "datarange", "onset_rel50_med"),
      drop = FALSE
    ]
    
    stopifnot(!anyDuplicated(x[c("species", "datarange")]))
    
    species <- sort(unique(x$species))
    windows <- sort(unique(x$datarange))
    
    w <- matrix(
      NA_real_,
      nrow = length(species),
      ncol = length(windows),
      dimnames = list(species, windows)
    )
    
    w[cbind(
      match(x$species, species),
      match(x$datarange, windows)
    )] <- x$onset_rel50_med
    
    rg <- apply(w, 1, function(z) {
      if (all(is.na(z))) return(NA_real_)
      max(z, na.rm = TRUE) - min(z, na.rm = TRUE)
    })
    
    P(sprintf(
      "%s: median range %.2f, max %.1f (%s)",
      m,
      median(rg, na.rm = TRUE),
      max(rg, na.rm = TRUE),
      names(rg)[which.max(rg)]
    ))
    
    if (m == "min") {
      show_table(as.data.frame(w[c(
        "Pseudotsuga menziesii",
        "Pinus sylvestris",
        "Abies alba",
        "Quercus rubra"
      ), , drop = FALSE]))
    }
  }
  
  # ------------------------------------------------------------
  # Period restriction
  # ------------------------------------------------------------
  
  P("\n=== §3.5 period restriction (restricted - full)")
  
  # Align explicitly by species, as pandas does when subtracting.
  species <- sort(unique(c(t$species, tr$species)))
  
  sh <- as.data.frame(
    matrix(
      NA_real_,
      nrow = length(species),
      ncol = length(M),
      dimnames = list(species, names(M))
    )
  )
  
  for (m in names(M)) {
    full <- sel(t, m)
    restricted <- sel(tr, m)
    
    sh[[m]] <- round(
      restricted[species, "onset_rel50_med"] -
        full[species, "onset_rel50_med"],
      2
    )
  }
  
  sh_matrix <- as.matrix(sh)
  
  P(sprintf(
    "within 1 °C: %d of %d; median %.1f",
    sum(abs(sh_matrix) <= 1, na.rm = TRUE),
    length(sh_matrix),
    median(as.vector(sh_matrix))
  ))
  
  show_table(sh[
    rowSums(abs(sh_matrix) > 0.45, na.rm = TRUE) > 0,
    ,
    drop = FALSE
  ])
  
  # ------------------------------------------------------------
  # Concurvity
  # ------------------------------------------------------------
  
  P("\n=== §3.5 concurvity (daily mean, 60-172)")
  
  cc <- read_input("concurvity_all_sweep_full.csv")
  cc <- cc[
    which(cc$datarange == "60-172" & cc$metric_fam == "mean"),
    ,
    drop = FALSE
  ]
  
  cd <- cc[which(cc$confound_type == "doy"), , drop = FALSE]
  cg <- cc[
    which(cc$confound_type == "gdd" & cc$tb == 0),
    ,
    drop = FALSE
  ]
  
  P(sprintf(
    paste0(
      "T vs DOY: observed median %.2f; worst median %.2f ",
      "(range %.2f-%.2f); T vs GDD(Tb=0) worst median %.2f, max %.3f"
    ),
    median(cd$observed_concurvity_med, na.rm = TRUE),
    median(cd$worst_concurvity_med, na.rm = TRUE),
    min(cd$worst_concurvity_med, na.rm = TRUE),
    max(cd$worst_concurvity_med, na.rm = TRUE),
    median(cg$worst_concurvity_med, na.rm = TRUE),
    max(cg$worst_concurvity_med, na.rm = TRUE)
  ))
  
  # ------------------------------------------------------------
  # Basis-dimension diagnostics
  # ------------------------------------------------------------
  
  P("\n=== §3.4 edf/k' (six cross-validated species, daily mean, 60-172)")
  
  k <- read_input("kcheck_summary_full.csv")
  k <- k[
    which(
      k$metric_fam == "mean" &
        k$datarange == "60-172" &
        k$species %in% SIX
    ),
    ,
    drop = FALSE
  ]
  
  terms <- c(
    gdd_tb0 = "s(gdd_tb0)",
    doy = "s(doy)",
    base = "s(Temperature)"
  )
  
  for (mt in names(terms)) {
    term <- terms[[mt]]
    x <- k$edf_ratio[which(k$model_type == mt & k$term == term)]
    
    P(sprintf(
      "%s: %.2f-%.2f",
      term,
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE)
    ))
  }
  
  # ------------------------------------------------------------
  # Leave-one-year-out paired comparisons
  # ------------------------------------------------------------
  
  P(paste0(
    "\n=== §3.4 leave-one-year-out paired differences ",
    "(mean [range]; folds; Wilcoxon p)"
  ))
  
  l <- read_input("loyo_brier_full.csv")
  
  get_brier <- function(mt, tb = NULL) {
    keep <- l$model_type == mt
    if (!is.null(tb)) {
      keep <- keep & !is.na(l$tb) & l$tb == tb
    }
    
    x <- l[which(keep), c("species", "year", "brier"), drop = FALSE]
    stopifnot(!anyDuplicated(x[c("species", "year")]))
    x
  }
  
  pair_brier <- function(a, b) {
    x <- merge(
      a, b,
      by = c("species", "year"),
      all = TRUE,
      suffixes = c("_a", "_b"),
      sort = TRUE
    )
    x$d <- x$brier_a - x$brier_b
    x
  }
  
  wilcoxon_p <- function(d) {
    # These conditions hold for the supplied archive.
    # Stop if future data require a different inference branch.
    stopifnot(
      length(d) > 0L,
      length(d) <= 50L,
      all(is.finite(d)),
      !any(d == 0),
      !anyDuplicated(abs(d))
    )
    
    stats::wilcox.test(
      d,
      alternative = "two.sided",
      exact = TRUE,
      correct = FALSE
    )$p.value
  }
  
  base <- get_brier("base")
  doy <- get_brier("doy")
  g0 <- get_brier("gdd", 0)
  tpd <- get_brier("t_plus_doy")
  tpg0 <- get_brier("t_plus_gdd", 0)
  
  COMP <- list(
    list(name = "M2 - M3", a = g0, b = doy, sign = ">"),
    list(name = "M1 - M2", a = base, b = g0, sign = "<"),
    list(name = "M4 - M3", a = tpd, b = doy, sign = "<"),
    list(name = "M1 - M3", a = base, b = doy, sign = ">"),
    list(name = "M5 - M2", a = tpg0, b = g0, sign = "<")
  )
  
  fmt <- function(v) {
    text <- if (abs(v) < 5e-4) "0.000" else sprintf("%+.3f", v)
    gsub("-", "\u2212", text, fixed = TRUE)
  }
  
  rows <- character()
  
  for (sp in SIX) {
    cells <- character()
    
    for (comparison in COMP) {
      paired <- pair_brier(comparison$a, comparison$b)
      
      d <- paired$d[
        paired$species == sp & !is.na(paired$d)
      ]
      
      n <- length(d)
      kk <- if (comparison$sign == ">") {
        sum(d > 0)
      } else {
        sum(d < 0)
      }
      
      p <- wilcoxon_p(d)
      
      P(sprintf(
        paste0(
          "%-22s %-13s mean %+.4f median %+.4f ",
          "[%+.3f, %+.3f] %d/%d p=%.4f"
        ),
        sp,
        comparison$name,
        mean(d),
        median(d),
        min(d),
        max(d),
        kk,
        n,
        p
      ))
      
      star <- if (p < 0.05) "**" else ""
      p_text <- if (p < 0.001) {
        sprintf("%.4f", p)
      } else {
        sprintf("%.3f", p)
      }
      
      cells <- c(cells, paste0(
        star, fmt(mean(d)), star,
        " [", fmt(min(d)), ", ", fmt(max(d)), "]; ",
        kk, "/", n, "; p = ", p_text
      ))
    }
    
    # Retains the Python script's Springs value from the last comparison.
    rows <- c(rows, paste0(
      "| *", sp, "* | ", n, " | ",
      paste(cells, collapse = " | "),
      " |"
    ))
  }
  
  # ------------------------------------------------------------
  # Largest forcing-model failures
  # ------------------------------------------------------------
  
  x <- pair_brier(base, g0)
  x <- x[
    complete.cases(x[c("brier_a", "brier_b")]),
    ,
    drop = FALSE
  ]
  
  names(x)[names(x) == "brier_a"] <- "base"
  names(x)[names(x) == "brier_b"] <- "gdd0"
  x$d <- x$gdd0 - x$base
  
  P("largest forcing-model failures (gdd0 - base):")
  
  top <- head(x[order(x$d, decreasing = TRUE), , drop = FALSE], 3)
  rownames(top) <- NULL
  show_table(top, digits = 3)
  
  P(sprintf(
    "max fold where temperature trails forcing: %.3f",
    max(x$base - x$gdd0, na.rm = TRUE)
  ))
  
  # ------------------------------------------------------------
  # Base-temperature sweep
  # ------------------------------------------------------------
  
  sm <- read_input("loyo_brier_summary.csv")
  
  P("Tb sweep, mean Brier(gdd) - Brier(doy), Tb = 0..6:")
  
  for (sp in sort(unique(sm$species))) {
    xx <- sm[sm$species == sp, , drop = FALSE]
    dv <- xx$mean_brier[which(xx$model_type == "doy")][1]
    
    gg <- xx[which(xx$model_type == "gdd"), , drop = FALSE]
    gg <- gg[order(gg$tb), , drop = FALSE]
    
    monotone <- isTRUE(all(diff(gg$mean_brier) < 0))
    
    P(paste0(
      sprintf("  %-22s ", sp),
      paste(sprintf("%+.4f", gg$mean_brier - dv), collapse = " "),
      "  monotone decrease: ",
      if (monotone) "True" else "False"
    ))
  }
  
  # ------------------------------------------------------------
  # Write outputs
  # ------------------------------------------------------------
  
  table3 <- c(
    paste0(
      "| Species | Springs | M2 \u2212 M4 | (M1 + M4) \u2212 M4 | ",
      "M1 \u2212 M2 | M1 \u2212 M4 |"
    ),
    "|---|---|---|---|---|---|",
    rows
  )
  
  # writeLines(
  #   enc2utf8(table3),
  #   file.path(RES, "old_table4_cv_rel.md"),
  #   useBytes = TRUE
  # )
  
  writeLines(
    enc2utf8(out),
    file.path(RES, "table1_table2_tableS2_ms_numbers_report.txt"),
    useBytes = TRUE
  )
  
  cat(paste(out, collapse = "\n"), "\n", sep = "")
  
  invisible(list(
    report = out,
    table3 = table3
  ))
}

ms_numbers_results <- run_ms_numbers()
