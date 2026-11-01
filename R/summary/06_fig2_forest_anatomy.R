source("R/summary/05_figure_helpers.R")

fig_diamond <- function(y, estimate, low, high, col, h = 0.32) {
  polygon(
    c(low, estimate, high, estimate),
    c(y, y + h, y, y - h),
    col = col,
    border = col,
    lwd = 0.6
  )
}

forest_panel <- function(metric, xlim, ref = NULL,
                         label_right = TRUE,
                         show_names = TRUE,
                         species_order = NULL,
                         meta = parent.frame()) {
  s <- load_thresholds(metric, onset = "onset_rel50")
  fits <- fit_anatomy(s, meta)
  f0 <- fits$f0
  f1 <- fits$f1
  
  s$w <- f1$w
  wmax <- max(s$w)
  
  y <- 0
  headers <- numeric(length(GRP))
  group_y <- numeric(length(GRP))
  parts <- vector("list", length(GRP))
  
  for (i in seq_along(GRP)) {
    g <- GRP[i]
    headers[i] <- y
    y <- y - 1
    
    sub <- s[s$anat == g, , drop = FALSE]
    
    if (is.null(species_order)) {
      sub <- sub[order(sub$onset_med), , drop = FALSE]
    } else {
      wanted <- species_order[species_order %in% sub$species]
      sub <- sub[match(wanted, sub$species), , drop = FALSE]
    }
    
    sub$y <- y - seq_len(nrow(sub)) + 1
    y <- y - nrow(sub)
    
    parts[[i]] <- sub
    group_y[i] <- y
    y <- y - 1.6
  }
  
  all_y <- y
  
  plot(
    NA,
    xlim = xlim,
    ylim = c(all_y - 1, if (is.null(ref)) 1 else 2.2),
    xaxs = "i", yaxs = "i",
    axes = FALSE,
    xlab = "", ylab = ""
  )
  
  axis(
    1,
    at = if (metric == "mean") seq(6, 18, 2) else 
      if (metric == "max") seq(10, 24, 2) else 
        if (metric == "min") seq(2, 14, 2) else NULL,
    lwd = 0.7,
    lwd.ticks = 0.7
  )
  
  overall <- contrast(f0, 1)
  
  abline(
    v = overall["estimate"],
    col = "grey60",
    lty = 3,
    lwd = 0.6
  )
  
  if (!is.null(ref)) {
    abline(v = ref, col = "#B03A2E", lty = 2, lwd = 0.7)
    text(
      ref, 1,
      sprintf("Rossi et al. (2008)\nconifers, %.1f °C", ref),
      adj = c(0.5, 0),
      col = "#B03A2E",
      cex = 0.8
    )
  }
  
  left_label <- function(y, label, italic = FALSE) {
    text(
      xlim[1] - 0.025 * diff(xlim),
      y,
      label,
      adj = 1,
      xpd = NA,
      font = if (italic) 3 else 1,
      cex = 0.9,
      col = if (italic) "black" else "grey35"
    )
  }
  
  right_label <- function(y, estimate, low, high,
                          col = "grey30", bold = FALSE) {
    text(
      xlim[2], y,
      fmt_ci(estimate, low, high),
      adj = 1,
      cex = 0.825,
      col = col,
      font = if (bold) 2 else 1
    )
  }
  
  for (i in seq_along(GRP)) {
    g <- GRP[i]
    sub <- parts[[i]]
    col <- FIG_COL[[g]]
    
    text(
      xlim[1], headers[i], FIG_LABEL[[g]],
      adj = 0, font = 2, col = col, cex = 0.975
    )
    
    segments(
      sub$onset_low, sub$y,
      sub$onset_high, sub$y,
      col = col, lwd = 1
    )
    
    points(
      sub$onset_med, sub$y,
      pch = FIG_PCH[[g]],
      bg = col,
      col = "white",
      lwd = 0.4,
      cex = sqrt(14 + 60 * sub$w / wmax) / 4
    )
    
    if (show_names) {
      left_label(sub$y, fig_abbrev(sub$species), italic = TRUE)
    }
    
    if (label_right) {
      right_label(
        sub$y,
        sub$onset_med,
        sub$onset_low,
        sub$onset_high
      )
    }
    
    ci <- contrast(f1, diag(length(GRP))[i, ])
    
    fig_diamond(
      group_y[i],
      ci["estimate"], ci["ci_low"], ci["ci_high"],
      col
    )
    
    if (show_names) {
      left_label(group_y[i], "group estimate")
    }
    
    if (label_right) {
      right_label(
        group_y[i],
        ci["estimate"], ci["ci_low"], ci["ci_high"],
        col = col, bold = TRUE
      )
    }
  }
  
  fig_diamond(
    all_y,
    overall["estimate"],
    overall["ci_low"],
    overall["ci_high"],
    "grey35"
  )
  
  if (show_names) left_label(all_y, "all species")
  
  if (label_right) {
    right_label(
      all_y,
      overall["estimate"],
      overall["ci_low"],
      overall["ci_high"],
      bold = TRUE
    )
  }
  
  omnibus <- ftest(f1, rbind(c(-1, 0, 1), c(-1, 1, 0)))
  dt <- contrast(f1, c(-1, 0, 1))
  dr <- contrast(f1, c(0, -1, 1))
  
  variance_text <- if (f0$tau2 > 0.1) {
    R2 <- max(0, (f0$tau2 - f1$tau2) / f0$tau2)
    sprintf(
      "explains %.0f %% of between-species variance",
      100 * R2
    )
  } else {
    "τ² = 0 with and without anatomy (variance explained not defined)"
  }
  
  p_label <- parse(text = fmt_p(omnibus["p"]))[[1]]
  
  list(
    test = bquote(
      "Anatomy: " *
        italic(F)(.(as.integer(omnibus["df1"])),
                  .(as.integer(omnibus["df2"]))) ==
        .(sprintf("%.2f", omnibus["F"])) *
        ", " * .(p_label)
    ),
    variance = variance_text,
    tau = if (f0$tau2 > 0.1) {
      sprintf("(τ² %.2f → %.2f)", f0$tau2, f1$tau2)
    } else {
      ""
    },
    dt = sprintf(
      "diffuse − coniferous %s °C [%s, %s]",
      fmt_signed(dt["estimate"]),
      fmt_num(dt["ci_low"]),
      fmt_num(dt["ci_high"])
    ),
    dr = sprintf(
      "diffuse − ring %s °C [%s, %s]",
      fmt_signed(dr["estimate"]),
      fmt_num(dr["ci_low"]),
      fmt_num(dr["ci_high"])
    )
  )
}

fig_png(
  "fig2_forest_anatomy_dailymean.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("mean", c(4, 24), ref = 8.4, meta = meta)
    title(xlab = "Onset threshold, daily mean (°C)")
    
    caption <- list(
      txt$test,
      trimws(paste(txt$variance, txt$tau)),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[[i]], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }  
  }
)

fig_png(
  "figS2_forest_anatomy_dailymin.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("min", c(-1, 18), ref = 4.1, meta = meta)
    title(xlab = "Onset threshold, daily minimum (°C)")
    
    caption <- list(
      txt$test,
      trimws(paste(txt$variance, txt$tau)),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[[i]], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }
  }
)

fig_png(
  "figS3_forest_anatomy_dailymax.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("max", c(6, 30), ref = NULL, meta = meta)
    title(xlab = "Onset threshold, daily maximum (°C)")
    
    caption <- list(
      txt$test,
      trimws(paste(txt$variance, txt$tau)),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[[i]], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }
  }
)
