# anatomy_meta_analysis.R
# Source from the project root; sourcing defines functions without running analyses.
# Requires your existing R/summary/01_meta.R and output/tables/*.csv.
# For direct calls to fit_anatomy() or summarise_anatomy(), first source
# R/summary/01_meta.R, or pass an environment containing its functions as meta.
# run_anatomy_meta_analysis() loads that dependency automatically per run.
# Example wrapper (save these uncommented lines in R/run_all.R):
# source("code/anatomy_meta_analysis.R")
# anatomy_rel <- run_anatomy_meta_analysis(
#   onset = "onset_rel50",
#   files = c(groups = "main/table3_anatomy_meta_regression.csv",
#             contrasts = "supplement/tableS4_anatomy_contrasts.csv",
#             loso = "supplement/tableS5_anatomy_loso.csv",
#             sensitivity = "supplement/tableS6_anatomy_sensitivity.csv"))
# anatomy_p50 <- run_anatomy_meta_analysis(
#   onset = "onset_p50",
#   files = c(groups = "supplement/tableS7_anatomy_meta_regression.csv",
#             contrasts = "supplement/tableS8_anatomy_contrasts.csv",
#             loso = "supplement/tableS9_anatomy_loso.csv",
#             sensitivity = "supplement/tableS10_anatomy_sensitivity.csv"))
# Run the wrapper with: Rscript R/run_all.R
# Table names above are examples; choose your manuscript numbering.

# Constants and helper functions are available in the sourcing environment.
ANAT <- c(
  "Abies alba" = "tracheid", "Picea abies" = "tracheid",
  "Pinus sylvestris" = "tracheid", "Larix decidua" = "tracheid",
  "Pseudotsuga menziesii" = "tracheid",
  "Quercus robur" = "ring", "Quercus rubra" = "ring", "Ulmus laevis" = "ring",
  "Fagus sylvatica" = "diffuse", "Carpinus betulus" = "diffuse",
  "Betula pendula" = "diffuse", "Acer pseudoplatanus" = "diffuse",
  "Alnus glutinosa" = "diffuse", "Tilia cordata" = "diffuse"
)
# Conifers use the same label as the model's tracheid group.
GRP <- c("tracheid", "ring", "diffuse")
METRICS <- c(min = "gam_occ_omin_base_10", mean = "gam_occ_mean_base_10",
             max = "gam_occ_max_base_10")
Z <- 1.959964

load_thresholds <- function(metric, onset, data_dir = file.path("output", "tables"),
                            se_rule = "symmetric", min_dy = NULL) {
  stopifnot(metric %in% names(METRICS),
            onset %in% c("onset_rel50", "onset_p50"),
            se_rule %in% c("symmetric", "wider"))
  t <- read.csv(file.path(data_dir, "thresholds_full.csv"), stringsAsFactors = FALSE)
  onset_cols <- paste0(onset, c("_med", "_low", "_high"))
  required <- c("model_name", "datarange", "species", onset_cols)
  missing <- setdiff(required, names(t))
  if (length(missing)) stop("Missing threshold columns: ", paste(missing, collapse = ", "))
  keep <- !is.na(t$model_name) & !is.na(t$datarange) &
    t$model_name == METRICS[[metric]] & t$datarange == "60-172"
  s <- t[keep, c("species", onset_cols), drop = FALSE]
  names(s)[match(onset_cols, names(s))] <- c("onset_med", "onset_low", "onset_high")
  s$anat <- unname(ANAT[s$species])
  stopifnot(nrow(s) > 0L, !anyDuplicated(s$species), !anyNA(s$anat))
  onset_values <- s[c("onset_med", "onset_low", "onset_high")]
  if (!all(vapply(onset_values, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(onset_values)))) {
    stop("Onset values must be finite numeric values for metric ", metric)
  }
  if (any(s$onset_low > s$onset_med | s$onset_med > s$onset_high))
    stop("Invalid onset interval ordering for metric ", metric)
  s$se <- if (se_rule == "symmetric") {
    (s$onset_high - s$onset_low) / (2 * Z)
  } else {
    pmax(s$onset_high - s$onset_med, s$onset_med - s$onset_low) / Z
  }
  if (any(s$se <= 0)) stop("Standard errors must be positive for metric ", metric)
  if (!is.null(min_dy)) {
    dy <- read.csv(file.path(data_dir, "summary_stats.csv"), stringsAsFactors = FALSE)
    stopifnot(all(c("species", "tree_years_num") %in% names(dy)),
              !anyDuplicated(dy$species))
    n_dy <- dy$tree_years_num[match(s$species, dy$species)]
    s <- s[!is.na(n_dy) & n_dy >= min_dy, , drop = FALSE]
  }
  rownames(s) <- NULL
  s
}

fit_anatomy <- function(s, meta = parent.frame()) {
  rma_fun = meta$rma
  if (!is.function(rma_fun)) {
    stop(
      "`meta$rma` must be your meta-analysis fitting function. ",
      "Pass the environment or list containing it explicitly."
    )
  }
  if (nrow(s) <= length(GRP) || !all(GRP %in% s$anat))
    stop("Anatomy fitting requires all three groups and more than three species.")
  X <- vapply(GRP, function(g) as.numeric(s$anat == g), numeric(nrow(s)))
  list(f0 = meta$rma(s$onset_med, s$se^2, matrix(1, nrow(s), 1)),
       f1 = meta$rma(s$onset_med, s$se^2, X))
}

summarise_anatomy <- function(s, label, meta = parent.frame()) {
  fits <- fit_anatomy(s, meta)
  f0 <- fits$f0
  f1 <- fits$f1
  omnibus <- meta$ftest(f1, rbind(c(-1, 0, 1), c(-1, 1, 0)))
  R2 <- if (f0$tau2 > 0) max(0, (f0$tau2 - f1$tau2) / f0$tau2) else NA_real_
  make_row <- function(fit, c_i, group, n_species) {
    ci <- meta$contrast(fit, c_i)
    pi <- meta$prediction_interval(fit, c_i)
    data.frame(analysis = label, group = group, n_species = n_species,
               estimate = unname(ci["estimate"]), ci_low = unname(ci["ci_low"]),
               ci_high = unname(ci["ci_high"]), pi_low = unname(pi["pi_low"]),
               pi_high = unname(pi["pi_high"]), stringsAsFactors = FALSE)
  }
  rows <- lapply(seq_along(GRP), function(i)
    make_row(f1, diag(length(GRP))[i, ], GRP[i], sum(s$anat == GRP[i])))
  rows[[length(rows) + 1L]] <- make_row(f0, 1, "all", nrow(s))
  g <- do.call(rbind, rows)
  g$tau2_null <- f0$tau2
  g$tau2_anatomy <- f1$tau2
  g$tau_anatomy <- sqrt(f1$tau2)
  g$R2_anatomy <- R2
  g$F <- unname(omnibus["F"])
  g$df1 <- unname(omnibus["df1"])
  g$df2 <- unname(omnibus["df2"])
  g$p_anatomy <- unname(omnibus["p"])
  vectors <- list("diffuse-tracheid" = c(-1, 0, 1),
                  "ring-tracheid" = c(-1, 1, 0), "diffuse-ring" = c(0, -1, 1))
  cons <- lapply(names(vectors), function(nm) {
    ci <- meta$contrast(f1, vectors[[nm]])
    data.frame(analysis = label, contrast = nm, estimate = unname(ci["estimate"]),
               ci_low = unname(ci["ci_low"]), ci_high = unname(ci["ci_high"]),
               p = unname(ci["p"]), stringsAsFactors = FALSE)
  })
  list(groups = g, contrasts = do.call(rbind, cons))
}

loso <- function(s, metric) {
  # Preserve the original unweighted arithmetic-mean predictions.
  rows <- lapply(seq_len(nrow(s)), function(i) {
    rest <- s[-i, , drop = FALSE]
    grp <- rest$onset_med[rest$anat == s$anat[i]]
    data.frame(metric = metric, species = s$species[i], anatomy = s$anat[i],
               observed = s$onset_med[i], pred_overall = mean(rest$onset_med),
               pred_group = if (length(grp)) mean(grp) else NA_real_,
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$abs_err_overall <- abs(d$observed - d$pred_overall)
  d$abs_err_group <- abs(d$observed - d$pred_group)
  d$group_better <- d$abs_err_group < d$abs_err_overall
  d$group_better[is.na(d$group_better)] <- FALSE
  d
}

round_numeric <- function(d, digits) {
  cols <- vapply(d, is.numeric, logical(1))
  d[cols] <- lapply(d[cols], round, digits = digits)
  d
}

write_results <- function(d, filename, results_dir = "paper/tables") {
  d <- round_numeric(d, 4)
  cols <- vapply(d, is.logical, logical(1))
  d[cols] <- lapply(d[cols], function(x)
    ifelse(is.na(x), NA_character_, ifelse(x, "True", "False")))
  path <- file.path(results_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(d, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

run_anatomy_meta_analysis <- function(onset, files, data_dir = file.path("output", "tables"),
                                      results_dir = "paper/tables", meta_file = "R/summary/01_meta.R",
                                      print_results = TRUE) {
  stopifnot(is.character(onset), length(onset) == 1L, !is.na(onset),
            onset %in% c("onset_rel50", "onset_p50"))
  keys <- c("groups", "contrasts", "loso", "sensitivity")
  if (!is.character(files) || length(files) != length(keys) ||
      is.null(names(files)) || anyDuplicated(names(files)) ||
      !setequal(names(files), keys) || anyNA(files) ||
      any(!nzchar(files)) || anyDuplicated(files)) {
    stop("files must contain four unique filenames named groups, contrasts, loso, sensitivity.")
  }
  # Load the existing statistical implementation into a per-run environment.
  meta <- new.env(parent = environment())
  sys.source(meta_file, envir = meta)
  for (nm in c("rma", "ftest", "contrast", "prediction_interval")) {
    if (!exists(nm, envir = meta, inherits = FALSE) || !is.function(meta[[nm]]))
      stop("Missing statistical helper in meta_file: ", nm)
  }
  G <- C <- L <- S <- list()
  settings <- list(se_wider_half_interval = list(se_rule = "wider"),
                   min_20_dendrometer_years = list(min_dy = 20))
  for (m in names(METRICS)) {
    s <- load_thresholds(m, onset, data_dir)
    main <- summarise_anatomy(s, m, meta)
    G[[m]] <- data.frame(metric = m, main$groups, stringsAsFactors = FALSE)
    C[[m]] <- data.frame(metric = m, main$contrasts, stringsAsFactors = FALSE)
    L[[m]] <- loso(s, m)
    for (lab in names(settings)) {
      s2 <- do.call(load_thresholds,
                    c(list(metric = m, onset = onset, data_dir = data_dir), settings[[lab]]))
      sensitivity <- summarise_anatomy(s2, lab, meta)
      g2 <- data.frame(metric = m, sensitivity$groups, stringsAsFactors = FALSE)
      c2 <- sensitivity$contrasts
      dt <- c2[c2$contrast == "diffuse-tracheid", , drop = FALSE]
      g2$diffuse_minus_tracheid <- dt$estimate
      g2$diffuse_minus_tracheid_low <- dt$ci_low
      g2$diffuse_minus_tracheid_high <- dt$ci_high
      S[[paste(m, lab, sep = "_")]] <- g2
    }
  }
  G <- do.call(rbind, G)
  C <- do.call(rbind, C)
  L <- do.call(rbind, L)
  S <- do.call(rbind, S)
  rownames(G) <- rownames(C) <- rownames(L) <- rownames(S) <- NULL
  write_results(G, files[["groups"]], results_dir)
  if(files[["contrasts"]]!="con_none") {
    write_results(C, files[["contrasts"]], results_dir)
  }
  if(files[["loso"]]!="loso_none") {
    write_results(L, files[["loso"]], results_dir)
  }
  if(files[["sensitivity"]]!="sens_none") {
    write_results(S, files[["sensitivity"]], results_dir)
  }
  
  if (print_results) {
    print(round_numeric(G[, c("metric", "group", "n_species", "estimate", "ci_low",
                              "ci_high", "pi_low", "pi_high", "tau_anatomy",
                              "R2_anatomy", "F", "p_anatomy")], 3), row.names = FALSE)
    print(round_numeric(C, 3), row.names = FALSE)
    loso_summary <- do.call(rbind, lapply(sort(unique(L$metric)), function(m) {
      d <- L[L$metric == m, , drop = FALSE]
      data.frame(metric = m, mae_overall = mean(d$abs_err_overall, na.rm = TRUE),
                 mae_group = mean(d$abs_err_group, na.rm = TRUE),
                 group_better = sum(d$group_better))
    }))
    print(round_numeric(loso_summary, 2), row.names = FALSE)
    print(round_numeric(S[S$group != "all", c("metric", "analysis", "group", "estimate",
                                              "R2_anatomy", "p_anatomy", "diffuse_minus_tracheid",
                                              "diffuse_minus_tracheid_low", "diffuse_minus_tracheid_high")], 3),
          row.names = FALSE)
  }
  invisible(list(groups = G, contrasts = C, loso = L, sensitivity = S))
}

source("R/summary/01_meta.R")
anatomy_rel <- run_anatomy_meta_analysis(
  onset = "onset_rel50",
  files = c(groups = "table3_anatomy_meta_regression.csv",
            contrasts = "text3.3_anatomy_contrasts.csv",
            loso = "tableS1_anatomy_loso.csv",
            sensitivity = "sens_none"))

anatomy_p50 <- run_anatomy_meta_analysis(
  onset = "onset_p50",
  files = c(groups = "tableS3_anatomy_meta_regression.csv",
            contrasts = "con_none",
            loso = "loso_none",
            sensitivity = "sens_none"))

