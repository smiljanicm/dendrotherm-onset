source("R/summary/05_figure_helpers.R")

if (!requireNamespace("minpack.lm", quietly = TRUE)) {
  stop("Install 'minpack.lm' before running this script.")
}

fig_logistic <- function(T, k, x0, M) {
  M / (1 + exp(-k * (T - x0)))
}

mock_parameters <- function(r) {
  temperatures <- as.numeric(r[1, c(
    "onset_rel05_med", "onset_rel10_med",
    "onset_rel25_med", "onset_rel50_med"
  )])
  
  probabilities <- c(0.05, 0.10, 0.25, 0.50)
  M <- max(r$obs_max_prob, 0.51)
  start <- c(k = 0.4, x0 = r$onset_rel50_med)
  
  # Same least-squares objective, start, and fallback as Python.
  pars <- tryCatch({
    stopifnot(all(is.finite(temperatures)))
    
    fit <- minpack.lm::nls.lm(
      par = start,
      fn = function(par) {
        fig_logistic(temperatures, par["k"], par["x0"], M) -
          probabilities
      },
      control = minpack.lm::nls.lm.control(
        ftol = 1.49012e-8,
        ptol = 1.49012e-8,
        gtol = 0,
        maxfev = 20000,
        maxiter = 1000
      )
    )
    
    if (!fit$info %in% 1:4) stop("Nonlinear fit did not converge.")
    fit$par
  }, error = function(e) start)
  
  c(pars, M = M)
}

fig_png(
  "fig1_probability_curves_by_anatomy_MOCK.png",
  9.6, 6.8,
  function() {
    t <- fig_read("thresholds_full.csv")
    
    s <- t[
      t$model_name == "gam_occ_mean_base_10" &
        t$datarange == "60-172",
      ,
      drop = FALSE
    ]
    
    stopifnot(!anyDuplicated(s$species))
    rownames(s) <- s$species
    
    par(
      mfrow = c(3, 6),
      mar = c(1.8, 2.3, 1.1, 0.6),
      oma = c(3.2, 3.4, 2.2, 3),
      cex = 1
    )
    
    for (ri in seq_along(GRP)) {
      g <- GRP[ri]
      species <- sort(names(ANAT)[ANAT == g])
      col <- FIG_COL[[g]]
      
      for (ci in 1:6) {
        plot(
          NA,
          xlim = c(-2, 26), ylim = c(0, 1.06),
          xaxs = "i", yaxs = "i",
          axes = FALSE, xlab = "", ylab = ""
        )
        
        if (ci <= length(species)) {
          sp <- species[ci]
          r <- s[sp, , drop = FALSE]
          pars <- mock_parameters(r)
          
          rect(
            r$t_q05, 0, r$t_q95, 1.06,
            col = "grey93", border = NA
          )
          
          rect(
            r$onset_rel50_low, 0, r$onset_rel50_high, 1.06,
            col = fig_alpha(col, 0.18), border = NA
          )
          
          abline(h = 0.5, col = "grey55", lty = 2, lwd = 0.6)
          
          xs <- seq(r$t_min, r$t_max, length.out = 300)
          
          lines(
            xs,
            fig_logistic(xs, pars["k"], pars["x0"], pars["M"]),
            col = col, lwd = 1.6
          )
          
          segments(
            r$onset_rel50_med, 0,
            r$onset_rel50_med, 0.5,
            col = col, lty = 3, lwd = 0.8
          )
          
          points(
            r$onset_rel50_med, 0.5,
            pch = 21, bg = col, col = "white",
            lwd = 0.8, cex = 0.9
          )
          
          axis(
            1, at = c(0, 10, 20),
            labels = ri == 3,
            cex.axis = 0.85, lwd = 0.7
          )
          
          axis(
            2, at = c(0, 0.25, 0.5, 0.75, 1),
            labels = ci == 1,
            cex.axis = 0.85, lwd = 0.7
          )
          
          text(
            -2 + 0.04 * 28, 0.99 * 1.06,
            fig_abbrev(sp),
            adj = c(0, 1), font = 3, cex = 0.86
          )
          
          text(
            -2 + 0.96 * 28, 0.05 * 1.06,
            sprintf("%.1f °C", r$onset_rel50_med),
            adj = c(1, 0), col = col, cex = 0.86
          )
        }
        
        # Legend occupies unused cells of the middle row.
        if (ri == 2 && ci == 5) {
          legend(
            "center",
            legend = c(
              "central 90 % of observed temperatures",
              "95 % CI of the P = 0.50 crossing",
              "P = 0.50"
            ),
            pch = c(22, 22, NA),
            pt.bg = c("grey93", fig_alpha("grey55", 0.25), NA),
            col = c("grey93", fig_alpha("grey55", 0.25), "grey55"),
            lty = c(NA, NA, 2),
            pt.cex = 1.4,
            cex = 0.875,
            bty = "n",
            xpd = NA
          )
        }
        
        if (ci == 6) {
          mtext(
            tolower(FIG_LABEL[[g]]),
            side = 4, line = 0.6, las = 0,
            col = col, cex = 1.06
          )
        }
      }
    }
    
    mtext(
      "Daily mean air temperature (°C)",
      side = 1, outer = TRUE, line = 1.1, cex = 1.125
    )
    
    mtext(
      "Fitted probability of registered growth",
      side = 2, outer = TRUE, line = 1.8,
      las = 0, cex = 1.125
    )
    
    mtext(
      paste0(
        "LAYOUT MOCK \u2014 curves reconstructed from tabulated ",
        "crossings, not the fitted GAMs"
      ),
      side = 3, outer = TRUE, line = 0.6,
      col = "#B03A2E", cex = 0.94
    )
  }
)
