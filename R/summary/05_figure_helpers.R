# Shared plotting helpers.
# Requires anatomy_meta_analysis.R to have been sourced.

if (!exists("load_thresholds", mode = "function") ||
    !exists("fit_anatomy", mode = "function")) {
  stop("First run source('R/summary/02_anatomy_meta_analysis.R').")
}

FIG_COL <- c(
  tracheid = "#2E6E8E",
  ring = "#C7743A",
  diffuse = "#5B8C3E"
)

FIG_PCH <- c(tracheid = 24, ring = 22, diffuse = 21)

FIG_LABEL <- c(
  tracheid = "Coniferous",
  ring = "Ring / semi-ring-porous",
  diffuse = "Diffuse-porous"
)

fig_abbrev <- function(species) {
  vapply(strsplit(species, " ", fixed = TRUE), function(x) {
    paste0(substr(x[1], 1, 1), ". ", x[2])
  }, character(1), USE.NAMES = FALSE)
}

fig_read <- function(filename) {
  read.csv(
    file.path("output", "tables", filename),
    stringsAsFactors = FALSE
  )
}

fig_png <- function(filename, width, height, draw) {
  dir.create("paper/figures", showWarnings = FALSE)
  
  args <- list(
    filename = file.path("paper", "figures", filename),
    width = width,
    height = height,
    units = "in",
    res = 300,
    pointsize = 8,
    bg = "white"
  )
  
  if (capabilities("cairo")) args$type <- "cairo"
  
  do.call(grDevices::png, args)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  par(
    family = "sans",
    las = 1,
    bty = "l",
    mgp = c(2, 0.5, 0),
    tcl = -0.2,
    lend = "butt"
  )
  
  draw()
  message("Written: ", file.path("paper", "figures", filename))
}

fig_alpha <- function(col, alpha) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

# Returns plotmath text for parse(text = fmt_p(p))[[1]]
fmt_p <- function(p) {
  if (p < 0.01) return("italic(P) < 0.01")
  if (p < 0.05 && round(p, 2) >= 0.05) {
    return("italic(P) < 0.05")
  }
  sprintf("italic(P) == %.2f", round(p, 2))
}

# Use a second decimal when one decimal would display a nonzero value as zero.
fmt_num <- function(x) {
  vapply(x, function(value) {
    digits <- 1L
    if (value != 0 && round(value, digits) == 0) {
      digits <- 2L
      while (round(value, digits) == 0) digits <- digits + 1L
    }
    gsub("-", "\u2212",
         sprintf(paste0("%.", digits, "f"), value),
         fixed = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

fmt_ci <- function(est, lo, hi) {
  sprintf("%s [%s, %s]", fmt_num(est), fmt_num(lo), fmt_num(hi))
}

# For contrasts displayed with an explicit positive sign.
fmt_signed <- function(x) {
  ifelse(x > 0, paste0("+", fmt_num(x)), fmt_num(x))
}
