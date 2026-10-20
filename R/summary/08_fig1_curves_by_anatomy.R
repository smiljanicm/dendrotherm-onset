library(tidyverse)
library(mgcv)

gs_data <- read_csv('./data/gs_data.csv')

# ── 1. Fit point-estimate models per species ─────────────────────────────────
species_list <- unique(gs_data$species)

models <- map(species_list, function(sp) {
  sp_data <- gs_data %>% filter(species == sp)
  gam(occ ~ s(Temperature, k = 10),
      data = sp_data,
      family = binomial(link = 'logit'),
      method = 'REML')
}) %>% setNames(species_list)

# ── 2. Generate prediction curves ───────────────────────────────────────────
pred_curves <- map_dfr(species_list, function(sp) {
  sp_data <- gs_data %>% filter(species == sp)
  t_seq <- seq(
    floor(quantile(sp_data$Temperature, 0.005)),
    ceiling(quantile(sp_data$Temperature, 0.995)),
    by = 0.1
  )
  nd <- data.frame(Temperature = t_seq)
  nd$fit <- predict(models[[sp]], newdata = nd, type = 'response')
  nd$species <- sp
  nd
})

# ── 3. Pull CI and range from your existing thresholds table ─────────────────
# thresholds already has onset_rel50_med, onset_rel50_low, onset_rel50_high,
# t_q05, t_q95 per species from the bootstrap
thresh_plot <- read_csv('./output/tables/thresholds_full.csv') %>%
  filter(model_name == 'gam_occ_mean_base_10', datarange == '60-172') %>%
  select(species, onset_rel50_med, onset_rel50_low, onset_rel50_high,
         t_q05, t_q95, obs_max_prob)

# install.packages("ggh4x")  # Run once if needed

library(tidyverse)
library(mgcv)
library(patchwork)
library(ggh4x)

# ── Species order and metadata ───────────────────────────────────────────────
species_meta <- tribble(
  ~species,                ~label,             ~wood_type,               ~col,
  "Abies alba",            "A. alba",           "tracheid-bearing",       "#2E6D9E",
  "Larix decidua",         "L. decidua",        "tracheid-bearing",       "#2E6D9E",
  "Picea abies",           "P. abies",          "tracheid-bearing",       "#2E6D9E",
  "Pinus sylvestris",      "P. sylvestris",     "tracheid-bearing",       "#2E6D9E",
  "Pseudotsuga menziesii", "P. menziesii",      "tracheid-bearing",       "#2E6D9E",
  "Quercus robur",         "Q. robur",          "ring / semi-ring-porous", "#C4622D",
  "Quercus rubra",         "Q. rubra",          "ring / semi-ring-porous", "#C4622D",
  "Ulmus laevis",          "U. laevis",         "ring / semi-ring-porous", "#C4622D",
  "Acer pseudoplatanus",   "A. pseudoplatanus", "diffuse-porous",          "#3A7A3A",
  "Alnus glutinosa",       "A. glutinosa",      "diffuse-porous",          "#3A7A3A",
  "Betula pendula",        "B. pendula",        "diffuse-porous",          "#3A7A3A",
  "Carpinus betulus",      "C. betulus",        "diffuse-porous",          "#3A7A3A",
  "Fagus sylvatica",       "F. sylvatica",      "diffuse-porous",          "#3A7A3A",
  "Tilia cordata",         "T. cordata",        "diffuse-porous",          "#3A7A3A"
)

# Width of each individual species panel
panel_width_mm <- 35

# ── Helper: build one row plot ───────────────────────────────────────────────

build_row_plot <- function(species_in_row, show_x_axis = FALSE) {
  
  row_curves <- pred_curves %>%
    filter(species %in% species_in_row) %>%
    left_join(species_meta, by = "species") %>%
    mutate(species = factor(species, levels = species_in_row))
  
  row_thresh <- thresh_plot %>%
    filter(species %in% species_in_row) %>%
    left_join(species_meta, by = "species") %>%
    mutate(
      species = factor(species, levels = species_in_row),
      
      # Height on the fitted GAM at the bootstrap median threshold
      fit_at_threshold = map2_dbl(
        as.character(species),
        onset_rel50_med,
        function(sp, temperature) {
          if (!is.finite(temperature)) return(NA_real_)
          
          as.numeric(predict(
            models[[sp]],
            newdata = data.frame(Temperature = temperature),
            type = "response"
          ))
        }
      )
    )
  
  p <- ggplot(row_curves, aes(x = Temperature, y = fit)) +
    
    # Grey band: central 90% observed temperatures
    geom_rect(
      data = row_thresh,
      aes(xmin = t_q05, xmax = t_q95, ymin = -Inf, ymax = Inf),
      fill = "grey85", alpha = 0.4,
      inherit.aes = FALSE
    ) +
    
    # Coloured band: threshold confidence interval
    geom_rect(
      data = row_thresh,
      aes(
        xmin = onset_rel50_low, xmax = onset_rel50_high,
        ymin = -Inf, ymax = Inf, fill = col
      ),
      alpha = 0.20,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    
    # Horizontal reference at the GAM height of the threshold
    geom_hline(
      data = row_thresh,
      aes(yintercept = fit_at_threshold),
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.4,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    
    # GAM curve
    geom_line(aes(colour = col), linewidth = 1.0) +
    scale_colour_identity() +
    
    # Dotted vertical at the bootstrap median threshold
    geom_vline(
      data = row_thresh,
      aes(xintercept = onset_rel50_med, colour = col),
      linetype = "dotted", linewidth = 0.5,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    
    # Threshold point on the fitted GAM
    geom_point(
      data = row_thresh,
      aes(
        x = onset_rel50_med,
        y = fit_at_threshold,
        colour = col
      ),
      size = 2.5, shape = 21, fill = "white",
      stroke = 1.2,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    
    # Threshold label
    geom_text(
      data = row_thresh,
      aes(
        x = pmin(
          pmax(onset_rel50_med, t_q05 + 0.5),
          t_q95 - 0.5
        ) + 1,
        y = 0.04,
        label = paste0(round(onset_rel50_med, 1), " °C"),
        colour = col
      ),
      hjust = 0, size = 2.8,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    
    # Species label top-left
    geom_text(
      data = row_thresh,
      aes(label = label, colour = col),
      x = -Inf, y = 0.97,
      hjust = -0.1, vjust = 1,
      size = 3,
      fontface='italic',
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~ species, nrow = 1) +
    
    # Constant width for every species panel
    ggh4x::force_panelsizes(
      cols = grid::unit(panel_width_mm, "mm")
    ) +
    
    scale_x_continuous(
      breaks = c(0, 10, 20),
      limits = c(-5, 25)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00)
    ) +
    
    labs(
      x = if (show_x_axis) "Daily mean air temperature (°C)" else NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 10) +
    theme(
      strip.text       = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      panel.spacing    = grid::unit(0.4, "lines"),
      axis.text        = element_text(size = 8),
      axis.title.x     = if (show_x_axis) element_text() else element_blank(),
      axis.text.x      = if (show_x_axis) element_text() else element_blank(),
      axis.ticks.x     = if (show_x_axis) element_line() else element_blank(),
      plot.margin      = margin(t = 4, r = 4, b = 2, l = 4)
    )
  
  p
}

# ── Build three rows ─────────────────────────────────────────────────────────
row1_sps <- species_meta$species[
  species_meta$wood_type == "tracheid-bearing"
]

row2_sps <- species_meta$species[
  species_meta$wood_type == "ring / semi-ring-porous"
]

row3_sps <- species_meta$species[
  species_meta$wood_type == "diffuse-porous"
]

p_row1 <- build_row_plot(row1_sps, show_x_axis = FALSE)
p_row2 <- build_row_plot(row2_sps, show_x_axis = FALSE)
p_row3 <- build_row_plot(row3_sps, show_x_axis = TRUE)

# ── Explanatory key in the empty space beside the middle row ─────────────────
annotation_key <- grid::grobTree(
  # Grey temperature-range band
  grid::rectGrob(
    x = 0.06, y = 0.62,
    width = grid::unit(3, "mm"),
    height = grid::unit(3, "mm"),
    gp = grid::gpar(fill = "grey85", col = NA, alpha = 0.4)
  ),
  grid::textGrob(
    "Central 90% of observed temperatures",
    x = 0.12, y = 0.62, just = "left",
    gp = grid::gpar(fontsize = 7)
  ),
  
  # Threshold confidence band: neutral swatch represents all species colours
  grid::rectGrob(
    x = 0.06, y = 0.50,
    width = grid::unit(3, "mm"),
    height = grid::unit(3, "mm"),
    gp = grid::gpar(fill = "grey50", col = NA, alpha = 0.20)
  ),
  grid::textGrob(
    "95% bootstrap CI of threshold temperature",
    x = 0.12, y = 0.50, just = "left",
    gp = grid::gpar(fontsize = 7)
  ),
  
  # Horizontal reference line
  grid::segmentsGrob(
    x0 = 0.025, x1 = 0.095,
    y0 = 0.38, y1 = 0.38,
    gp = grid::gpar(col = "grey60", lty = "dashed", lwd = 1)
  ),
  grid::textGrob(
    "Fitted GAM probability at threshold",
    x = 0.12, y = 0.38, just = "left",
    gp = grid::gpar(fontsize = 7)
  )
)

# ── Left-align panels and add an anatomy label at the far right ──────────────
left_align_row <- function(p, anatomy_label, anatomy_colour,
                           key = NULL) {
  
  gt <- ggplotGrob(p)
  
  # Flexible empty space after the fixed-width species panels
  gt <- gtable::gtable_add_cols(
    gt,
    widths = grid::unit(1, "null"),
    pos = -1
  )
  space_col <- ncol(gt)
  
  # Find the vertical extent of the actual panels
  panel_cells <- gt$layout[
    grepl("^panel", gt$layout$name),
    ,
    drop = FALSE
  ]
  
  panel_top <- min(panel_cells$t)
  panel_bottom <- max(panel_cells$b)
  
  # Optional key, placed in the flexible empty space
  if (!is.null(key)) {
    gt <- gtable::gtable_add_grob(
      gt,
      grobs = key,
      t = panel_top,
      b = panel_bottom,
      l = space_col,
      r = space_col,
      clip = "off",
      name = "annotation_key"
    )
  }
  
  # Fixed-width column for the rotated anatomy label
  gt <- gtable::gtable_add_cols(
    gt,
    widths = grid::unit(7, "mm"),
    pos = -1
  )
  
  gt <- gtable::gtable_add_grob(
    gt,
    grobs = grid::textGrob(
      anatomy_label,
      rot = 90,
      gp = grid::gpar(
        col = anatomy_colour,
        fontsize = 9
      )
    ),
    t = panel_top,
    b = panel_bottom,
    l = ncol(gt),
    r = ncol(gt),
    clip = "off",
    name = "anatomy_label"
  )
  
  patchwork::wrap_elements(full = gt)
}

# ── Assemble annotated rows ─────────────────────────────────────────────────
p_combined <- (
  left_align_row(
    p_row1,
    anatomy_label = "coniferous",
    anatomy_colour = "#2E6D9E"
  ) /
    left_align_row(
      p_row2,
      anatomy_label = "ring / semi-ring-porous",
      anatomy_colour = "#C4622D",
      key = annotation_key
    ) /
    left_align_row(
      p_row3,
      anatomy_label = "diffuse-porous",
      anatomy_colour = "#3A7A3A"
    )
) +
  plot_layout(heights = c(1, 1, 1))

# ── Shared y-axis label ──────────────────────────────────────────────────────
p_final <- wrap_elements(p_combined) +
  labs(tag = "Fitted probability of registered growth") +
  theme(
    plot.tag          = element_text(angle = 90, size = 10),
    plot.tag.position = "left",
    plot.margin       = margin(2, 4, 2, 4)
  )

ggsave(
  "paper/figures/fig1_probability_curves_by_anatomy.png",
  p_final,
  width = 250,
  height = 180,
  units = "mm"
)
