#-------------------------------------------------------------------------------
# Summary: R-script to plot Fig. 2 (cortical localization of deletion and duplication effects on cognitive ability)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Plots Fig. 2 a-d and writes the plotted values to the Source Data sheets
#             a) Del / Dup effect-size maps (180 Glasser regions) + pies of Bonferroni-significant regions
#             b) Del vs Dup effect sizes, all 180 regions (Spearman rho, two-sided p-spin) + significance-class inset map
#             c) cortical dosage-effect contrast map, cDEC = z(Del) - z(Dup)
#             d) cDEC vs cortical hierarchy rank (H = -T1w/T2w)
# Key steps:
#       1) PART A: Parameters and input files
#       2) PART B: Load effect sizes (ST4), H, profile statistics (ST7)
#       3) PART C: Panels a-d
#       4) PART D: Assemble; save PNG + SVG; write Source Data
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 29 October, 2026
# Version: 1
#
#-----
# UPDATE log:
#       v1 (29 Oct 2026): release version of Figure2_r3_12oct2026.R; inputs from the final Supplementary Tables
#                         (ST4 effect sizes, ST7 profile statistics), regions joined by name; Source Data sheets added
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (sheets ST4, ST7)
#       2) data/inputs/inputs_main_neuromaps20_Glasser180.csv (T1w/T2w)
#
#-----
# Output:
#       1) Figs/Fig2/Fig2.png, Fig2.svg
#       2) data/source_data_sheets/Fig2a-d.csv, Fig2_stats.csv
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

## ---- fig2_setup ----
#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Kuldeep's local system; Else: run from the repository root
flag_hires = 0  # 0: PNG at 300 dpi only (GitHub); 1: PNG at 600 dpi + SVG (publication files, Zenodo)
if( flag_local == 1){
  project_dir <- "/Users/kuldeepkumar/Documents/Projects/Temp_CNV_GWAS_oct2022/Release_NG_R1_29oct2026/GitHub"
} else {
  project_dir <- getwd()
}

data_dir <- paste0(project_dir,"/data")
inputs_dir <- paste0(data_dir,"/inputs")
source_data_dir <- paste0(data_dir,"/source_data_sheets")
OUT_DIR <- paste0(project_dir,"/Figs/Fig2")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(ggrepel)
  library(tibble)
  library(dplyr)
  library(readxl)
  library(ggseg)
  library(ggsegGlasser)
  library(cowplot)
  library(patchwork)
})
source(file.path(project_dir, "functions/fn_write_SourceData_sheet.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

FIG_WIDTH_IN  <- 7.20       # ~183 mm, Nature-family full-width figure
FIG_HEIGHT_IN <- 5.6
PNG_DPI       <- 600

filename_ST_C <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
filename_maps_C <- file.path(inputs_dir, "inputs_main_neuromaps20_Glasser180.csv")

# Global font-size parameters (points at final figure size).
buffer_PT = 3
FONT_MIN_PT <- 5 + buffer_PT
FONT_BASE_PT <- 7 + buffer_PT
FONT_AXIS_TITLE_PT <- 7 + buffer_PT
FONT_AXIS_TEXT_PT <- 6 + buffer_PT
FONT_LEGEND_TITLE_PT <- 7 + buffer_PT
FONT_LEGEND_TEXT_PT <- 6 + buffer_PT
FONT_ANNOTATION_PT <- 6 + buffer_PT
FONT_ROI_LABEL_PT <- 5.5
FONT_INTERNAL_LABEL_PT <- 7 + buffer_PT
FONT_PANEL_TAG_PT <- 12 + buffer_PT ## Panel Tag label
FONT_BRAIN_PT <- 7 + buffer_PT

# geom_text()/annotate() use mm-based text-size units, not points.
PT_PER_MM <- 72.27 / 25.4
pt_to_geom <- function(x) x / PT_PER_MM
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Colours: significance classes (Del / Dup / both / ns) and the sensorimotor-association scale
COL_DEL <- "#D95F02"
COL_DUP <- "#7570B3"
COL_BOTH <- "#1B9E77"
COL_NS <- "#D3D3D3"
COL_SENSORIMOTOR <- "#51127c"
COL_ASSOCIATION  <- "#fc8961"
COL_NEUTRAL <- "#F2F2F2"

## ---- fig2_functions ----
assert_file <- function(x) {
  if (!file.exists(x)) stop("Missing required input file: ", x, call. = FALSE)
}

zscore <- function(x) as.numeric(scale(as.numeric(x)))

format_pspin <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
}

# Bonferroni over the 360 regional tests (180 Del + 180 Dup)
bonferroni_del_dup <- function(p_del, p_dup) {
  pmat <- cbind(p_del, p_dup)
  padj <- matrix(p.adjust(as.vector(pmat), method = "bonferroni"), ncol = 2)
  list(del = padj[, 1], dup = padj[, 2])
}

make_pie <- function(n_sig, n_total = 180) {
  d <- data.frame(
    class = factor(c("ns", "sig"), levels = c("ns", "sig")),
    n = c(n_total - n_sig, n_sig),
    label = c("", n_sig)
  )
  ggplot(d, aes(x = "", y = n, fill = class)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = pt_to_geom(FONT_INTERNAL_LABEL_PT)) +
    scale_fill_manual(values = c("white", "#999999")) +
    theme_void() +
    theme(legend.position = "none")
}

make_effect_map <- function(effect, p_adj, cbar_lim = 0.13, show_legend = FALSE) {
  roi <- names(glasser$palette)[1:180]
  dat <- tibble(
    region = roi,
    effect = pmax(pmin(effect, cbar_lim), -cbar_lim),
    edge_col = ifelse(p_adj < 0.05, "black", "grey70"),
    edge_size = ifelse(p_adj < 0.05, 0.15, 0.03)
  )
  p <- ggseg(
    dat, atlas = "glasser", position = "dispersed", hemisphere = "left",
    mapping = aes(fill = effect, colour = I(edge_col), size = I(edge_size)),
    adapt_scales = TRUE
  ) +
    theme_brain(text.family = "sans", text.size = FONT_BRAIN_PT) +
    scale_fill_gradientn(
      limits = c(-cbar_lim, cbar_lim), breaks = c(-cbar_lim, 0, cbar_lim),
      colours = colorRampPalette(c("blue4", "blue", "white", "red", "red4"))(300),
      oob = scales::squish
    ) +
    guides(fill = guide_colourbar(title = "Effect size", title.position = "top",
                                  barwidth = 7.5, barheight = 0.55, title.hjust = 0.5)) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = FONT_LEGEND_TITLE_PT), legend.text = element_text(size = FONT_LEGEND_TEXT_PT),
      plot.margin = margin(-7, 0, -7, 0)
    )
  p <- p + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))
  p
}

make_significance_map <- function(sig_class) {
  roi <- names(glasser$palette)[1:180]
  code <- c("ns" = 0, "Del" = 1, "Dup" = 2, "Del-Dup" = 3)[sig_class]
  dat <- tibble(region = roi, effect = as.numeric(code))
  p = ggseg(dat, atlas = "glasser", position = "dispersed", hemisphere = "left",
        mapping = aes(fill = effect), colour = "grey65", size = 0.01,
        adapt_scales = TRUE) +
    theme_brain(text.family = "sans", text.size = FONT_BRAIN_PT) +
    scale_fill_gradientn(
      limits = c(0, 3), breaks = c(0, 1, 2, 3),
      colours = c(COL_NS, COL_DEL, COL_DUP, COL_BOTH)
    ) +
    theme(legend.position = "none", plot.margin = margin(-8, -8, -8, -8))
  p <- p + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))
  p
}

make_delta_map <- function(delta_z, cbar_lim) {
  roi <- names(glasser$palette)[1:180]
  dat <- tibble(region = roi, delta_z = pmax(pmin(delta_z, cbar_lim), -cbar_lim))
  p = ggseg(dat, atlas = "glasser", position = "dispersed", hemisphere = "left",
        mapping = aes(fill = delta_z), colour = "grey70", size = 0.03,
        adapt_scales = TRUE) +
    theme_brain(text.family = "sans", text.size = FONT_BRAIN_PT) +
    scale_fill_gradient2(
      low = COL_SENSORIMOTOR, mid = COL_NEUTRAL, high = COL_ASSOCIATION,   # cDEC = z(Del) - z(Dup)
      midpoint = 0, limits = c(-cbar_lim, cbar_lim), oob = scales::squish,
      breaks = c(-cbar_lim, 0, cbar_lim),
      labels = c("relatively Del", "0", "relatively Dup")
    ) +
    guides(fill = guide_colourbar(
      title = expression("Cortical dosage-effect contrast (cDEC)"),
      title.position = "top",
      barwidth = 7.0, barheight = 0.55, title.hjust = 0.5
    )) +
    theme(legend.position = "bottom", legend.title = element_text(size = FONT_LEGEND_TITLE_PT),
          legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0),
          legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
          legend.spacing.y = grid::unit(0, "cm"),
          legend.text = element_text(size = FONT_LEGEND_TEXT_PT), plot.margin = margin(-5, 0, -5, 0))
  p <- p + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))
  p
}

save_figure <- function(plot_obj, stem, width, height) {
  svg_file <- file.path(OUT_DIR, paste0(stem, ".svg"))
  png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
  if (flag_hires == 1) {   # SVG only for the publication files
    svglite::svglite(svg_file, width = width, height = height, bg = "white")
    print(plot_obj); invisible(grDevices::dev.off())
  }
  ggsave(png_file, plot_obj, width = width, height = height, units = "in",
         dpi = if (flag_hires == 1) PNG_DPI else 300, bg = "white")
}

#-----------------------------------------------------------------------------
# PART B: Load effect sizes, H and profile statistics
#-----------------------------------------------------------------------------

## ---- fig2_data ----
assert_file(filename_ST_C); assert_file(filename_maps_C)
roi_names <- names(glasser$palette)[1:180]                   # ggseg Glasser left hemisphere, 180 regions
roi_ST_names <- paste0("L_", gsub("-", ".", roi_names, fixed = TRUE))   # gene_set_name in the Supplementary Tables

df_ST4 <- read_excel(filename_ST_C, sheet = "ST4") %>% filter(gene_set_exp == "All_data_analysis_Fig2")
get_profile <- function(cnv, col) {
  x <- df_ST4[df_ST4$cnv_type == cnv, ]
  x <- x[match(roi_ST_names, x$gene_set_name), ]              # join by region name
  if (anyNA(x$gene_set_name)) stop("Region join failed for ", cnv)
  as.numeric(x[[col]])
}

array_es_DEL <- get_profile("DEL", "Estimate")
array_es_DUP <- get_profile("DUP", "Estimate")
array_pval_DEL <- get_profile("DEL", "pvalue")
array_pval_DUP <- get_profile("DUP", "pvalue")
padj <- bonferroni_del_dup(array_pval_DEL, array_pval_DUP)
array_pval_DEL_adj <- padj$del
array_pval_DUP_adj <- padj$dup

sig_class <- rep("ns", 180)
sig_class[array_pval_DEL_adj < 0.05 & array_pval_DUP_adj >= 0.05] <- "Del"
sig_class[array_pval_DEL_adj >= 0.05 & array_pval_DUP_adj < 0.05] <- "Dup"
sig_class[array_pval_DEL_adj < 0.05 & array_pval_DUP_adj < 0.05] <- "Del-Dup"
sig_class <- factor(sig_class, levels = c("ns", "Del", "Dup", "Del-Dup"))

# Cortical hierarchy H = -z(T1w/T2w), oriented sensorimotor -> association
neuromaps <- read.csv(filename_maps_C, check.names = FALSE)
neuromaps <- neuromaps[match(roi_names, neuromaps$glasser_region), ]
if (anyNA(neuromaps$glasser_region)) stop("Region join failed for the normative maps")
hierarchy_z <- -zscore(neuromaps$T1T2ratio)
hierarchy_rank <- rank(hierarchy_z, ties.method = "average")

del_z <- zscore(array_es_DEL)
dup_z <- zscore(array_es_DUP)
delta_z <- del_z - dup_z   # cDEC = z(Del) - z(Dup); effects are negative, so negative cDEC = relatively larger Del effect

# Profile statistics from Supplementary Table 7 (cleaned rotations, n_rot = 9,999, two-sided); rho recomputed only as a guard
st7 <- read_excel(filename_ST_C, sheet = "ST7")
iq <- st7[st7$gene_set_exp == "All_data_analysis_Fig2", , drop = FALSE]
if (nrow(iq) != 1) stop("ST7 must contain exactly one All_data_analysis_Fig2 row.")
rho_del_dup <- iq$rho_DEL_DUP_all; pspin_del_dup <- iq$pspin_DEL_DUP_all
rho_cdec_h <- iq$rho_cDEC_H; pspin_cdec_h <- iq$pspin_cDEC_H
if (abs(cor(array_es_DEL, array_es_DUP, method = "spearman") - rho_del_dup) > 1e-9 ||
    abs(cor(delta_z, hierarchy_z, method = "spearman") - rho_cdec_h) > 1e-9) stop("Recomputed rho does not match ST7.")
if (abs(rho_cdec_h - 0.368) > 0.001 || abs(pspin_cdec_h - 0.0064) > 1e-9) stop("Guard failed: rho(cDEC,H) = +0.368, p_spin = 0.0064.")
if (sum(array_pval_DEL_adj < 0.05) != iq$n_sig_DEL || sum(array_pval_DUP_adj < 0.05) != iq$n_sig_DUP) stop("Pie counts do not match ST7 n_sig_DEL / n_sig_DUP.")

#-----------------------------------------------------------------------------
# PART C: Panels
#-----------------------------------------------------------------------------

## ---- fig2_source_data ----
# Source Data: one row per region with every value plotted in a-d; printed statistics in Fig2_stats
df_source_Fig2 <- data.frame(
  gene_set_name = roi_ST_names,
  glasser_region = roi_names,
  Estimate_DEL = array_es_DEL,
  pvalue_DEL = array_pval_DEL,
  pvalue_bonferroni_360_DEL = array_pval_DEL_adj,
  Estimate_DUP = array_es_DUP,
  pvalue_DUP = array_pval_DUP,
  pvalue_bonferroni_360_DUP = array_pval_DUP_adj,
  significance_class = as.character(sig_class),
  z_DEL = del_z,
  z_DUP = dup_z,
  cDEC = delta_z,
  T1T2ratio = neuromaps$T1T2ratio,
  H = hierarchy_z,
  H_rank = hierarchy_rank
)
fn_write_SourceData_sheet(df_source_Fig2, "Fig2a-d", source_data_dir)
df_stats_Fig2 <- data.frame(
  panel = c("a", "a", "b", "b", "d", "d"),
  statistic = c("n_significant_regions_DEL (Bonferroni 0.05/360)", "n_significant_regions_DUP (Bonferroni 0.05/360)",
                "rho_DEL_DUP_all (Spearman)", "pspin_DEL_DUP_all (two-sided)", "rho_cDEC_H (Spearman)", "pspin_cDEC_H (two-sided)"),
  value = c(sum(array_pval_DEL_adj < 0.05), sum(array_pval_DUP_adj < 0.05), rho_del_dup, pspin_del_dup, rho_cdec_h, pspin_cdec_h),
  printed = c(sum(array_pval_DEL_adj < 0.05), sum(array_pval_DUP_adj < 0.05), sprintf("%.2f", rho_del_dup), format_pspin(pspin_del_dup),
              sprintf("%.2f", rho_cdec_h), format_pspin(pspin_cdec_h)),
  source = "Supplementary Table 7 (row All_data_analysis_Fig2)"
)
fn_write_SourceData_sheet(df_stats_Fig2, "Fig2_stats", source_data_dir)

## ---- fig2_panel_a ----
#----- Panel a: Del / Dup effect-size maps + pies of significant regions (Fig. 2a)
p_del <- make_effect_map(array_es_DEL, array_pval_DEL_adj, cbar_lim = 0.13, show_legend = TRUE)
legend_eff <- cowplot::get_legend(p_del)
p_del <- p_del + theme(legend.position = "none")
p_dup <- make_effect_map(array_es_DUP, array_pval_DUP_adj, cbar_lim = 0.13, show_legend = FALSE)

p_del_pie <- ggdraw(p_del) + draw_plot(make_pie(sum(array_pval_DEL_adj < 0.05)), x = 0.33, y = 0.41, width = 0.33, height = 0.42)
p_dup_pie <- ggdraw(p_dup) + draw_plot(make_pie(sum(array_pval_DUP_adj < 0.05)), x = 0.33, y = 0.41, width = 0.33, height = 0.42)

panel_A_body <- plot_grid(p_del_pie, p_dup_pie, ncol = 1, labels = c("Del", "Dup"),
                          label_fontface = "bold", label_size = FONT_INTERNAL_LABEL_PT, label_x = 0.44,
                          rel_heights = c(1, 1))
legend_eff_shifted <- ggdraw() + draw_plot(legend_eff, x = 0, y = 0.23, width = 1, height = 1)
panel_A <- plot_grid(panel_A_body, legend_eff_shifted, ncol = 1, rel_heights = c(1, 0.13))

## ---- fig2_panel_b ----
#----- Panel b: Del vs Dup effect sizes, all 180 regions (Fig. 2b)
df_b <- data.frame(
  region = roi_names,
  del = array_es_DEL,
  dup = array_es_DUP,
  significance = sig_class
)
cols_sig <- c("ns" = COL_NS, "Del" = COL_DEL, "Dup" = COL_DUP, "Del-Dup" = COL_BOTH)
stat_label_b <- paste0("ρ = ", sprintf("%.2f", rho_del_dup),
                       "\np-spin = ", format_pspin(pspin_del_dup))

p_b <- ggplot(df_b, aes(x = dup, y = del)) +
  geom_point(aes(color = significance), size = 1.8, alpha = 0.95) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotdash", linewidth = 0.45,
              color = "black", fullrange = TRUE) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.35) +
  scale_color_manual(values = cols_sig, breaks = c("Del", "Dup", "Del-Dup", "ns"),
                     labels = c("Del", "Dup", "Del-Dup", "ns"), name = "Signif.") +
  annotate("text", x = -Inf, y = -Inf, label = stat_label_b,
           hjust = -0.05, vjust = -0.35, size = pt_to_geom(FONT_ANNOTATION_PT)) +
  labs(x = "Effect size Dup", y = "Effect size Del") +
  theme_minimal(base_size = FONT_BASE_PT) +
  theme(legend.position = "bottom", legend.title = element_text(size = FONT_LEGEND_TITLE_PT),
        legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0),
        legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
        legend.spacing.y = grid::unit(0, "cm"),
        legend.text = element_text(size = FONT_LEGEND_TEXT_PT), legend.key.width = grid::unit(0.7, "lines"),
        axis.title = element_text(size = FONT_AXIS_TITLE_PT),
        axis.text = element_text(size = FONT_AXIS_TEXT_PT),
        panel.grid.minor = element_blank())

p_sig_map <- make_significance_map(as.character(sig_class))
panel_B <- ggdraw(p_b) + draw_plot(p_sig_map, x = 0.68, y = 0.79, width = 0.31, height = 0.28)

## ---- fig2_panel_c ----
#----- Panel c: cDEC map (Fig. 2c)
delta_lim <- ceiling(max(abs(delta_z), na.rm = TRUE) * 2) / 2
if (!is.finite(delta_lim) || delta_lim < 1) delta_lim <- 1
panel_C <- make_delta_map(delta_z, delta_lim)

## ---- fig2_panel_d ----
#----- Panel d: cDEC vs cortical hierarchy rank (Fig. 2d)
df_d <- data.frame(
  region = roi_names,
  hierarchy = hierarchy_z,
  hierarchy_rank = hierarchy_rank,
  delta_z = delta_z
)
stat_label_d <- paste0("ρ = ", sprintf("%.2f", rho_cdec_h),   # ST7 rho_cDEC_H, cDEC = z(Del) - z(Dup)
                       "\np-spin = ", format_pspin(pspin_cdec_h))

p_d <- ggplot(df_d, aes(x = hierarchy_rank, y = delta_z, color = delta_z)) +
  # Shade bottom/top 30% of cortical hierarchy
  annotate("rect",xmin = 1, xmax = 54,ymin = -Inf, ymax = Inf,fill = "grey70", alpha = 0.12) +
  annotate("rect",xmin = 127, xmax = 180,ymin = -Inf, ymax = Inf,fill = "grey70", alpha = 0.12) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, inherit.aes = TRUE,
              linetype = "dashed", linewidth = 0.45, color = "grey30") +
  geom_point(size = 1.7, alpha = 0.95) +
  scale_color_gradient2(low = COL_SENSORIMOTOR, mid = COL_NEUTRAL, high = COL_ASSOCIATION,
                        midpoint = 0, limits = c(-delta_lim, delta_lim), oob = scales::squish,
                        guide = "none") +
  annotate("text", x = -Inf, y = Inf, label = stat_label_d,
           hjust = -0.05, vjust = 1.15, size = pt_to_geom(FONT_ANNOTATION_PT)) +
  scale_x_continuous(breaks = c(1,180),
                     labels = NULL,   # hide default labels
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous( breaks = c(-2, 0, 2)) +
  labs(x = "Cortical hierarchy rank", y = expression("Cortical dosage-effect contrast (cDEC)")) +
  theme_minimal(base_size = FONT_BASE_PT) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = FONT_AXIS_TITLE_PT),
        axis.text = element_text(size = FONT_AXIS_TEXT_PT),
        axis.text.x = element_text(size = FONT_AXIS_TEXT_PT)) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(8, 10, 22, 10)   # extra bottom room for labels/arrows
  ) +
  annotate("text",x = 14, y = -4.7, label = "Sensorimotor",hjust = 0, size = pt_to_geom(FONT_AXIS_TEXT_PT)) +
  annotate("text", x = 167, y = -4.7, label = "Association", hjust = 1, size = pt_to_geom(FONT_AXIS_TEXT_PT)) +
  annotate("segment", x = 16, xend = 3, y = -4.4, yend = -4.4,
           arrow = arrow(length = unit(0.10, "cm"), type = "closed"),linewidth = 0.35) +
  annotate("segment", x = 165, xend = 178, y = -4.4, yend = -4.4,
           arrow = arrow(length = unit(0.10, "cm"), type = "closed"),linewidth = 0.35)
panel_D = ggdraw(p_d)

#-----------------------------------------------------------------------------
# PART D: Assemble and save
#-----------------------------------------------------------------------------

## ---- fig2_assemble ----
design <- "
AB
CD
#D
"
final_fig <- panel_A + panel_B + panel_C + panel_D +
  plot_layout(design = design, widths = c(1, 1), heights = c(1, 0.5, 0.5)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(face = "bold", size = FONT_PANEL_TAG_PT))) &
  theme(plot.tag = element_text(face = "bold", size = FONT_PANEL_TAG_PT))

save_figure(final_fig, "Fig2", FIG_WIDTH_IN, FIG_HEIGHT_IN)
message("Saved Fig. 2 to: ", OUT_DIR)
