# Run: Rscript brier_table.R loyo_brier_full.csv
# Base R only. Writes a formatted Markdown table and an unrounded CSV.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else "output/tables/loyo_brier_full.csv"
tb_value <- 6                  # This threshold reproduces the example's GDD columns.
n_boot <- 100000L
set.seed(42)

x <- read.csv(input, stringsAsFactors = FALSE)
required <- c("year", "brier", "species", "model_type", "tb")
stopifnot(all(required %in% names(x)))
models <- c("base", "gdd", "doy", "t_plus_gdd", "t_plus_doy")
labels <- c("Base", "GDD", "DOY", "T + GDD", "T + DOY")
species_order <- c("Fagus sylvatica", "Carpinus betulus", "Quercus robur",
                   "Acer pseudoplatanus", "Pinus sylvestris", "Betula pendula")
# Select an explicit threshold rather than optimizing it on held-out scores.
x <- x[x$model_type %in% models &
         (!(x$model_type %in% c("gdd", "t_plus_gdd")) |
            (!is.na(x$tb) & x$tb == tb_value)), ]
if (!setequal(unique(x$species), species_order)) {
  stop("Species differ from the example; update species_order.")
}
if (any(!is.finite(x$brier)) || any(x$brier < 0 | x$brier > 1) ||
    anyNA(x$year)) stop("Missing or invalid Brier scores/years.")
if (anyDuplicated(x[c("species", "year", "model_type")])) {
  stop("Duplicate species/year/model rows after selecting tb.")
}

summarise_species <- function(sp) {
  z <- x[x$species == sp, ]
  # Strict pairing: require every model to have exactly the same held-out years.
  years <- sort(unique(z$year))
  bs <- vapply(models, function(m) {
    q <- z[z$model_type == m, ]
    if (!setequal(q$year, years)) stop(paste("Unpaired years:", sp, m))
    q$brier[match(years, q$year)]
  }, numeric(length(years)))
  avg <- colMeans(bs)
  if (any(avg[c("base", "gdd", "doy")] <= 0)) stop("Zero reference mean.")
  # Ratio of means, not mean of annual ratios. Each year has equal weight.
  skill <- 100 * (1 - avg / avg["base"])
  added_gdd <- unname(100 * (1 - avg["t_plus_gdd"] / avg["gdd"]))
  added <- 100 * (1 - avg["t_plus_doy"] / avg["doy"])
  # Resample held-out years as paired blocks, preserving both model scores.
  draws <- replicate(n_boot, {
    i <- sample.int(nrow(bs), nrow(bs), replace = TRUE)
    denom <- mean(bs[i, "doy"])
    if (denom <= 0) stop("Zero DOY mean in bootstrap draw.")
    denom_gdd <- mean(bs[i, "gdd"])
    if (denom_gdd <= 0) stop("Zero GDD mean in bootstrap draw.")
    c(gdd = 100 * (1 - mean(bs[i, "t_plus_gdd"]) / denom_gdd),
      doy = 100 * (1 - mean(bs[i, "t_plus_doy"]) / denom))
  })
  ci <- unname(quantile(draws["doy", ], c(0.025, 0.975), type = 7))
  ci_gdd <- unname(quantile(draws["gdd", ], c(0.025, 0.975), type = 7))
  list(years = length(years), brier = avg, skill = skill, added = unname(added), ci = ci,
       added_gdd = added_gdd, ci_gdd = ci_gdd)
}
res <- lapply(species_order, summarise_species)
skill <- do.call(rbind, lapply(res, `[[`, "skill"))
brier <- do.call(rbind, lapply(res, `[[`, "brier"))
median_brier <- apply(brier, 2, median)
median_skill <- apply(skill, 2, median)
# Rank models within each species; 1 = lowest mean Brier / highest skill.
mean_rank <- colMeans(t(apply(skill, 1, function(v) rank(-v, ties.method = "average"))))
added <- vapply(res, `[[`, numeric(1), "added")
added_gdd <- vapply(res, `[[`, numeric(1), "added_gdd")

fmt <- function(v) gsub("-", "\u2212", sprintf("%.1f", v), fixed = TRUE)
format_scores <- function(bs, ss) {
  out <- paste0(sprintf("%.4f", bs), " (", fmt(ss), "%)")
  # Bold the lowest unrounded Brier score(s).
  winners <- which(abs(bs - min(bs)) < 1e-10)
  out[winners] <- paste0("**", out[winners], "**")
  out
}
rows <- lapply(seq_along(res), function(i) {
  r <- res[[i]]
  c(paste0("*", species_order[i], "*"), r$years, format_scores(r$brier, r$skill),
    paste0(fmt(r$added_gdd), " % (", fmt(r$ci_gdd[1]), "\u2013", fmt(r$ci_gdd[2]), ")"),
    paste0(fmt(r$added), " % (", fmt(r$ci[1]), "\u2013", fmt(r$ci[2]), ")"))
})
rows <- c(rows, list(
  c("**Median**", "", format_scores(median_brier, median_skill),
    paste0(fmt(median(added_gdd)), " %"),
    paste0(fmt(median(added)), " %")),
  c("**Mean rank**", "", {
    v <- fmt(mean_rank)
    v[which.min(mean_rank)] <- paste0("**", v[which.min(mean_rank)], "**")
    v
  }, "", "")
))
header <- c("Species", "Years", labels, "T added to GDD (95% CI)", "T added to DOY (95% CI)")
md_row <- function(v) paste0("| ", paste(v, collapse = " | "), " |")
md <- c(md_row(paste0("**", header, "**")),
        md_row(c(":---", rep("---:", length(header) - 1))),
        vapply(rows, md_row, character(1)),
        "",
        "Model cells show mean LOYO Brier score (Brier skill relative to Base, %). Each held-out year has equal weight. Lower Brier scores and higher skill are better; bold indicates the lowest Brier score in each score row.",
        "",
        "The Median row reports separate across-species medians of Brier scores and skill percentages; its percentages are not calculated from the median Brier scores. Mean rank uses within-species model ranks (1 = best). T added columns show relative Brier reductions against GDD or DOY, with paired-year percentile bootstrap 95% CIs.")
writeLines(md, "paper/tables/table4_brier_table.md", useBytes = TRUE)
cat(paste(md, collapse = "\n"), "\n")
raw <- data.frame(Species = species_order,
                  Years = vapply(res, `[[`, integer(1), "years"),
                  setNames(as.data.frame(brier), paste0(models, "_brier")),
                  setNames(as.data.frame(skill), paste0(models, "_skill_percent")),
                  T_added_to_GDD = added_gdd,
                  GDD_CI_lower = vapply(res, function(r) r$ci_gdd[1], numeric(1)),
                  GDD_CI_upper = vapply(res, function(r) r$ci_gdd[2], numeric(1)),
                  T_added_to_DOY = added,
                  CI_lower = vapply(res, function(r) r$ci[1], numeric(1)),
                  CI_upper = vapply(res, function(r) r$ci[2], numeric(1)),
                  check.names = FALSE)
# write.csv(raw, "R/summary/brier_table_values.csv", row.names = FALSE)

# CIs are percentile bootstrap intervals across held-out years. They do not
# refit models or account for dependence between overlapping LOYO training sets.
