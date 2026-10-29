#-------------------------------------------------------------------------------
# Summary: shared style (sizes, fonts, colours, theme, save) for the main and supplementary figures
#-------------------------------------------------------------------------------
# Author: Kuldeep Kumar; 
# Date: 21 October, 2026 
# Source after project_dir is set. Font block as in the Fig. 2 and Fig. 5 scripts.
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(ggprism) })

# ------------------------------- CONSTANTS ------------------------------------
FIG_WIDTH_IN <- 7.20   # ~183 mm, Nature-family full width
PNG_DPI      <- 600
FONT_FAMILY  <- "sans"

# ------------------------------- FONTS (points at final size) -----------------
buffer_PT <- 3
FONT_MIN_PT            <- 5 + buffer_PT
FONT_BASE_PT           <- 7 + buffer_PT
FONT_AXIS_TITLE_PT     <- 7 + buffer_PT
FONT_AXIS_TEXT_PT      <- 6 + buffer_PT
FONT_LEGEND_TITLE_PT   <- 7 + buffer_PT
FONT_LEGEND_TEXT_PT    <- 6 + buffer_PT
FONT_ANNOTATION_PT     <- 6 + buffer_PT
FONT_INTERNAL_LABEL_PT <- 7 + buffer_PT
FONT_PANEL_TAG_PT      <- 12 + buffer_PT
FONT_BRAIN_PT          <- 7 + buffer_PT
FONT_ROI_LABEL_PT      <- 5.5
# geom_text()/annotate() sizes are in mm
PT_PER_MM <- 72.27 / 25.4
pt_to_geom <- function(x) x / PT_PER_MM

# ------------------------------- COLOURS --------------------------------------
COL_DEL <- "#B7222E"; COL_DUP <- "#2878B8"                       # CNV identity
COL_SIG_DEL <- "#D95F02"; COL_SIG_DUP <- "#7570B3"               # significance classes only (never identity)
COL_SIG_BOTH <- "#1B9E77"; COL_SIG_NS <- "#D3D3D3"
COL_SENSORIMOTOR <- "#51127c"; COL_ASSOCIATION <- "#fc8961"; COL_NEUTRAL <- "#F2F2F2"
COL_DELTA_BAR <- "#8E4C97"                                       # cDEC bars
COL_POINT <- "grey60"
COL_LINE <- "grey70"
# band tints: the effective v_oct12 tints (#59A14F / #F28E2B at alpha 0.13 on white), as opaque fills
COL_BAND_ENERGY <- "#E9F3E8"; COL_BAND_EXPANSION <- "#FDF0E3"

# ------------------------------- THEME ----------------------------------------
theme_ng <- function(grid_x = FALSE) {
  t <- theme_prism(base_size = FONT_BASE_PT, base_family = FONT_FAMILY, base_fontface = "plain") +
    theme(axis.title = element_text(size = FONT_AXIS_TITLE_PT, face = "plain"),
          axis.text = element_text(size = FONT_AXIS_TEXT_PT, face = "plain"),
          legend.title = element_text(size = FONT_LEGEND_TITLE_PT, face = "plain"),
          legend.text = element_text(size = FONT_LEGEND_TEXT_PT, face = "plain"),
          plot.title = element_text(size = FONT_AXIS_TITLE_PT, face = "plain", hjust = 0.5),
          strip.text = element_text(size = FONT_AXIS_TEXT_PT, face = "plain"),
          plot.tag = element_text(size = FONT_PANEL_TAG_PT, face = "bold"),
          panel.grid.minor = element_blank())
  if (grid_x) t <- t + theme(panel.grid.major.x = element_line(colour = "grey80", linetype = "dotted", linewidth = 0.3))
  t
}

# ------------------------------- TEXT HELPERS ---------------------------------
format_pspin <- function(p) {   # as Figure2_r3_12oct2026.R
  if (is.na(p)) return("NA")
  if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
}
stat_line <- function(rho, p, prefix = "") paste0(prefix, "ρ = ", sprintf("%.2f", rho), ", p-spin = ", format_pspin(p))

# ------------------------------- SAVE -----------------------------------------
# flag_hires = 1 (set in the calling script): PNG at 600 dpi + SVG; otherwise PNG at 300 dpi only. No PDF.
fn_hires <- function() exists("flag_hires") && isTRUE(flag_hires == 1)
save_figure <- function(plot_obj, stem, width, height, out_dir = OUT_DIR) {
  if (fn_hires()) { svglite::svglite(file.path(out_dir, paste0(stem, ".svg")), width = width, height = height, bg = "white"); print(plot_obj); invisible(grDevices::dev.off()) }
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, units = "in", dpi = if (fn_hires()) PNG_DPI else 300, bg = "white")
  for (e in c(if (fn_hires()) "svg", "png")) if (!file.exists(file.path(out_dir, paste0(stem, ".", e)))) stop("save_figure: ", stem, ".", e, " not written")
}
