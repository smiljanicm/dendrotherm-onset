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
