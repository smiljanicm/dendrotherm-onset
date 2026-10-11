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
  
  # Establish the original row positions.
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
  
  axis(1, lwd = 0.7, lwd.ticks = 0.7)
  
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
      sprintf("%.1f [%.1f, %.1f]", estimate, low, high),
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
  R2 <- max(0, (f0$tau2 - f1$tau2) / f0$tau2)
  
  dt <- contrast(f1, c(-1, 0, 1))
  dr <- contrast(f1, c(0, -1, 1))
  if(f0$tau2 > 0.1) {
    variance_text = sprintf(
      "explains %.0f %% of between-species variance",
      100 * R2)
  } else{
    variance_text = sprintf(
      "\u03c4\u00b2 = 0 with and without anatomy (variance explained not defined)")
  }
  list(
    test = sprintf(
      "Anatomy: F(%d, %d) = %.2f, p = %.3f",
      as.integer(omnibus["df1"]),
      as.integer(omnibus["df2"]),
      omnibus["F"], omnibus["p"]
    ),
    
    variance = variance_text,
    tau = sprintf(
      "(\u03c4\u00b2 %.2f \u2192 %.2f)",
      f0$tau2, f1$tau2
    ),
    dt = sprintf(
      "diffuse \u2212 tracheid-bearing %+.1f °C [%.1f, %.1f]",
      dt["estimate"], dt["ci_low"], dt["ci_high"]
    ),
    dr = sprintf(
      "diffuse \u2212 ring %+.1f °C [%.1f, %.1f]",
      dr["estimate"], dr["ci_low"], dr["ci_high"]
    )
  )
}

fig_png(
  "fig2_forest_anatomy_dailymean_rel.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("mean", c(6, 23.5), ref = 8.4, meta = meta)
    title(xlab = "Onset threshold, daily mean (°C)")
    
    caption <- c(
      txt$test,
      paste(txt$variance, txt$tau),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[i], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }
  }
)

fig_png(
  "figS1a_forest_anatomy_dailymin_rel.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("min", c(1, 15.5), ref = 4.1, meta = meta)
    title(xlab = "Onset threshold, daily mean (°C)")
    
    caption <- c(
      txt$test,
      paste(txt$variance, txt$tau),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[i], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }
  }
)

fig_png(
  "figS1b_forest_anatomy_dailymax_rel.png",
  6.6, 6.2,
  function() {
    par(mar = c(3.4, 10, 1.5, 1), oma = c(3.5, 0, 0, 0))
    
    meta <- new.env(parent = environment())
    sys.source("R/summary/01_meta.R", envir = meta)
    
    txt <- forest_panel("max", c(5, 30.5), ref = NULL, meta = meta)
    title(xlab = "Onset threshold, daily max (°C)")
    
    caption <- c(
      txt$test,
      paste(txt$variance, txt$tau),
      paste(txt$dt, txt$dr, sep = ";  ")
    )
    
    for (i in seq_along(caption)) {
      mtext(
        caption[i], side = 1, outer = TRUE,
        line = i - 0.5, adj = 0.02,
        cex = 0.82, col = "grey25"
      )
    }
  }
)

fig_png(
  "figS1_forest_anatomy_3metrics_rel.png",
  12.5, 6.4,
  function() {
    layout(matrix(1:3, nrow = 1), widths = c(1.25, 1, 1))
    par(oma = c(1.5, 0, 0, 0))
    
    s <- load_thresholds("mean", onset = "onset_rel50")
    species_order <- s$species[order(s$onset_med)]
    
    metrics <- c("min", "mean", "max")
    titles <- c("Daily minimum", "Daily mean", "Daily maximum")
    limits <- list(c(0, 15), c(6, 21), c(9, 28))
    refs <- list(4.1, 8.4, NULL)
    
    for (i in seq_along(metrics)) {
      par(
        mar = c(8.3, if (i == 1) 10 else 1, 3.7, 1),
        cex = 1
      )
      
      meta <- new.env(parent = environment())
      sys.source("R/summary/01_meta.R", envir = meta)
      
      txt <- forest_panel(
        metrics[i],
        limits[[i]],
        ref = refs[[i]],
        label_right = FALSE,
        show_names = i == 1,
        species_order = species_order,
        meta = meta
      )
      
      title(main = titles[i], line = 2.1, cex.main = 1.125)
      title(xlab = paste0(
        "Onset threshold, ", tolower(titles[i]), " (°C)"
      ))
      
      caption <- unlist(txt, use.names = FALSE)
      
      for (j in seq_along(caption)) {
        mtext(
          caption[j], side = 1,
          line = 3.2 + (j - 1) * 0.9,
          cex = 0.78, col = "grey25"
        )
      }
    }
    
    mtext(
      "Species ordered within group by daily-mean threshold in all panels.",
      side = 1, outer = TRUE, line = 0.2,
      adj = 0.01, cex = 0.875, col = "grey35"
    )
  }
)
