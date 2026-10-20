source("R/summary/05_figure_helpers.R")

fig3_jitter <- c(
  -0.018255158965433727,
  0.048471388557274772,
  -0.10997483754018413,
  -0.043486834020995249,
  -0.077713704020235128,
  -0.089685509150864484,
  -0.069022753496912392,
  -0.0339766400505295,
  -0.022711155669252608,
  0.0085396814807385257,
  -0.01777720683127515,
  0.040748290087287101,
  -0.065020505059066175
)

fig_png(
  "figS1_LeaveOneYearOut_comparison.png",
  9.2, 4.3,
  function() {
    f <- fig_read("loyo_brier_full.csv")
    sm <- fig_read("loyo_brier_summary.csv")
    
    get_scores <- function(mt, tb = NULL) {
      keep <- f$model_type == mt
      
      if (!is.null(tb)) {
        keep <- keep & !is.na(f$tb) & f$tb == tb
      }
      
      x <- f[
        which(keep),
        c("species", "year", "brier"),
        drop = FALSE
      ]
      
      stopifnot(!anyDuplicated(x[c("species", "year")]))
      x
    }
    
    base <- get_scores("base")
    doy <- get_scores("doy")
    gdd0 <- get_scores("gdd", 0)
    
    names(base)[3] <- "base"
    names(doy)[3] <- "doy"
    names(gdd0)[3] <- "gdd0"
    
    piv <- merge(base, doy, by = c("species", "year"), all = TRUE)
    piv <- merge(
      piv, gdd0,
      by = c("species", "year"),
      all.x = TRUE
    )
    
    piv <- piv[order(piv$species, piv$year), , drop = FALSE]
    piv$d <- piv$gdd0 - piv$doy
    
    stopifnot(all(is.finite(piv$d)))
    
    means <- tapply(piv$d, piv$species, mean)
    species_order <- names(sort(means))
    
    pvalues <- vapply(species_order, function(sp) {
      d <- piv$d[piv$species == sp]
      
      stopifnot(
        length(d) <= 50,
        !any(d == 0),
        !anyDuplicated(abs(d))
      )
      
      wilcox.test(
        d,
        alternative = "two.sided",
        exact = TRUE,
        correct = FALSE
      )$p.value
    }, numeric(1))
    
    layout(matrix(1:2, nrow = 1), widths = c(1.15, 1))
    
    # ----------------------------------------------------------
    # Panel a: paired fold differences
    # ----------------------------------------------------------
    
    par(mar = c(5.2, 10, 2, 1.5), cex = 1)
    
    plot(
      NA,
      xlim = c(-0.07, max(piv$d) + 0.06),
      ylim = c(-0.5, length(species_order) - 0.5),
      xaxs = "i", yaxs = "i",
      axes = FALSE,
      xlab = "", ylab = ""
    )
    
    abline(v = 0, col = "grey40", lwd = 0.8)
    axis(1, lwd = 0.7)
    
    axis(
      2,
      at = seq_along(species_order) - 1,
      labels = fig_abbrev(species_order),
      tick = FALSE,
      font = 3,
      cex.axis = 0.925
    )
    
    for (i in seq_along(species_order)) {
      sp <- species_order[i]
      sub <- piv[piv$species == sp, , drop = FALSE]
      col <- FIG_COL[[ANAT[[sp]]]]
      
      stopifnot(nrow(sub) <= length(fig3_jitter))
      
      yy <- i - 1 + fig3_jitter[seq_len(nrow(sub))]
      is13 <- sub$year == 2013
      
      points(
        sub$d[!is13], yy[!is13],
        pch = 16,
        col = fig_alpha(col, 0.55),
        cex = 0.75
      )
      
      if (any(is13)) {
        points(
          sub$d[is13], yy[is13],
          pch = 1, col = col, lwd = 1.1, cex = 1.1
        )
        text(
          sub$d[is13], yy[is13] + 0.30,
          "2013", col = "grey30", cex = 0.75
        )
      }
      
      segments(
        mean(sub$d), i - 1 - 0.16,
        mean(sub$d), i - 1 + 0.16,
        col = col, lwd = 1.8
      )
      
      if (pvalues[sp] < 0.05) {
        text(
          mean(sub$d) + 0.012, i - 1 + 0.30,
          sprintf("p = %.2f", pvalues[sp]),
          adj = 0, cex = 0.8, col = "grey35"
        )
      }
    }
    
    title(xlab = "\u0394Brier  (forcing \u2212 seasonal position)")
    
    mtext("a", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.25)
    mtext(
      "\u2190 forcing better",
      side = 1, line = 3.5, adj = 0,
      cex = 0.85, col = "grey45"
    )
    mtext(
      "seasonal position better \u2192",
      side = 1, line = 3.5, adj = 1,
      cex = 0.85, col = "grey45"
    )
    
    # ----------------------------------------------------------
    # Panel b: base-temperature sweep
    # ----------------------------------------------------------
    
    curves <- lapply(species_order, function(sp) {
      gg <- sm[
        sm$species == sp & sm$model_type == "gdd",
        ,
        drop = FALSE
      ]
      gg <- gg[order(gg$tb), , drop = FALSE]
      
      dv <- sm$mean_brier[
        sm$species == sp & sm$model_type == "doy"
      ]
      
      stopifnot(length(dv) == 1L)
      
      data.frame(
        species = sp,
        tb = gg$tb,
        d = gg$mean_brier - dv
      )
    })
    
    ends <- do.call(rbind, lapply(curves, function(x) {
      data.frame(
        species = x$species[1],
        actual = tail(x$d, 1),
        label_y = tail(x$d, 1)
      )
    }))
    
    ends <- ends[order(ends$label_y), , drop = FALSE]
    
    if (nrow(ends) > 1L) {
      for (i in 2:nrow(ends)) {
        if (ends$label_y[i] - ends$label_y[i - 1] < 0.0035) {
          ends$label_y[i] <- ends$label_y[i - 1] + 0.0035
        }
      }
    }
    
    yrange <- range(
      0,
      unlist(lapply(curves, function(x) x$d)),
      ends$label_y
    )
    pad <- diff(yrange) * 0.06
    
    par(mar = c(5.2, 5.5, 2, 1), cex = 1)
    
    plot(
      NA,
      xlim = c(-0.3, 9.4),
      ylim = yrange + c(-pad, pad),
      xaxs = "i", yaxs = "i",
      axes = FALSE,
      xlab = "", ylab = ""
    )
    
    abline(h = 0, col = "grey40", lwd = 0.9)
    axis(1, at = 0:6, lwd = 0.7)
    axis(2, lwd = 0.7)
    
    for (x in curves) {
      col <- FIG_COL[[ANAT[[x$species[1]]]]]
      
      lines(x$tb, x$d, col = col, lwd = 1.2)
      points(
        x$tb, x$d,
        pch = 21, bg = col, col = "white",
        lwd = 0.6, cex = 0.7
      )
    }
    
    for (i in seq_len(nrow(ends))) {
      sp <- ends$species[i]
      col <- FIG_COL[[ANAT[[sp]]]]
      
      segments(
        6.05, ends$actual[i],
        6.45, ends$label_y[i],
        col = fig_alpha(col, 0.7), lwd = 0.6
      )
      
      text(
        6.55, ends$label_y[i],
        fig_abbrev(sp),
        adj = 0, font = 3, col = col, cex = 0.85
      )
    }
    
    title(xlab = expression("Base temperature " * T[b] * " (°C)"))
    
    mtext(
      "Mean \u0394Brier  (forcing \u2212 seasonal position)",
      side = 2, line = 3.5, las = 0
    )
    
    mtext("b", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.25)
    
    text(
      9.2, yrange[2] + pad * 0.45,
      "above the line: seasonal position predicts better",
      adj = c(1, 1), cex = 0.86, col = "grey35"
    )
  }
)
