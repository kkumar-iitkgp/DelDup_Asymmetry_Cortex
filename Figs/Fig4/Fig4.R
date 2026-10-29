#-------------------------------------------------------------------------------
# Summary: R-script to plot Fig. 4 (gene dosage effects on cognitive ability across levels of regional specificity)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Plots Fig. 4 a-e and writes the plotted values to the Source Data sheets
#             a,c) Del (a) and Dup (c) effect-size maps for gene sets at Z > 0.2, 0.5, 0.8, 1.0 + pies of significant regions
#             b,d) mean effect size (95% CI) of sensorimotor and association poles; * = exact KS test, FDR q < 0.05 (94 tests, ST11)
#             e)   Spearman rho(cDEC, H) per threshold (ST7); * = two-sided p-spin, FDR q < 0.05 across the 4 thresholds
# Key steps:
#       1) PART A: Parameters and input files
#       2) PART B: Load effect sizes (ST4), poles, KS statistics (ST11), profile statistics (ST7)
#       3) PART C: Panels a-e
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
#       v1 (29 Oct 2026): release version of Figure4_r7_12oct2026.R; inputs from the final Supplementary Tables
#                         (ST4, ST7, ST11), regions joined by name; PNG + SVG only; Source Data sheets added.
#                         The published figure adds column headers (Del / Dup / cDEC) and moves the legends by hand.
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (sheets ST4, ST7, ST11)
#       2) data/inputs/inputs_main_neuromaps20_Glasser180.csv (T1w/T2w); data/inputs/inputs_pole_membership_Glasser180.csv
#
#-----
# Output:
#       1) Figs/Fig4/Fig4.png, Fig4.svg
#       2) data/source_data_sheets/Fig4a-d.csv, Fig4b-e_stats.csv
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

## ---- fig4_setup ----
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
OUT_DIR <- paste0(project_dir,"/Figs/Fig4")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(dplyr)
  library(tibble)
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

FIG_WIDTH_IN <- 7.20 # ~183 mm, Nature-family full-width figure
FIG_HEIGHT_IN <- 3.30
PNG_DPI <- 600

filename_ST_C <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
filename_maps_C <- file.path(inputs_dir, "inputs_main_neuromaps20_Glasser180.csv")
filename_poles_C <- file.path(inputs_dir, "inputs_pole_membership_Glasser180.csv")

# Global font-size parameters (points at final figure size).
FONT_MIN_PT <- 6
FONT_BASE_PT <- 8.5
FONT_AXIS_TITLE_PT <- 9
FONT_AXIS_TEXT_PT <- 8
FONT_LEGEND_TITLE_PT <- 8.5
FONT_LEGEND_TEXT_PT <- 7.5
FONT_ANNOTATION_PT <- 7.5
FONT_ROI_LABEL_PT <- 6
FONT_INTERNAL_LABEL_PT <- 8.5
FONT_PANEL_TAG_PT <- 13  ## Panel Tag label
FONT_BRAIN_PT <- 7

PT_PER_MM <- 72.27 / 25.4
pt_to_geom <- function(x) x / PT_PER_MM
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COL_SENSORIMOTOR <- "#51127c"
COL_ASSOCIATION <- "#fc8961"
COL_DELTA_BAR <- "#8E4C97"

## ---- fig4_functions ----
assert_file <- function(x) if (!file.exists(x)) stop("Missing required input file: ", x, call. = FALSE)
zscore <- function(x) as.numeric(scale(as.numeric(x)))

save_figure <- function(plot_obj, stem, width, height) {
  svg_file <- file.path(OUT_DIR, paste0(stem, ".svg")); png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
  if (flag_hires == 1) {   # SVG only for the publication files
    svglite::svglite(svg_file, width = width, height = height, bg = "white")
    print(plot_obj); invisible(grDevices::dev.off())
  }
  ggsave(png_file, plot_obj, width = width, height = height, units = "in", dpi = if (flag_hires == 1) PNG_DPI else 300, bg = "white")
}

add_panel_tag <- function(plot_obj, tag) {
  ggdraw(plot_obj) +
    draw_plot_label(tag, x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = FONT_PANEL_TAG_PT)
}

# Bonferroni within each regional-specificity threshold: 180 Del + 180 Dup = 360 tests per threshold.
bonferroni_del_dup <- function(p_del, p_dup) {
  z <- matrix(p.adjust(c(p_del, p_dup), method = "bonferroni"), ncol = 2)
  list(del = z[, 1], dup = z[, 2])
}

make_pie <- function(n_sig, n_total = 180) {
  d <- data.frame(x = factor(c("ns", "sig"), levels = c("ns", "sig")), n = c(n_total - n_sig, n_sig), label = c("", n_sig))
  ggplot(d, aes(x = "", y = n, fill = x)) + geom_col(width = 1) + coord_polar("y") +
    geom_text(aes(label = label), position = position_stack(vjust = .5), size = pt_to_geom(FONT_INTERNAL_LABEL_PT)) +
    scale_fill_manual(values = c("white", "#999999")) + theme_void() + theme(legend.position = "none")
}

make_effect_map <- function(effect, p_adj, cbar_lim = 0.30, show_legend = FALSE) {
  roi <- names(glasser$palette)[1:180]
  dat <- tibble(region = roi, effect = pmax(pmin(effect, cbar_lim), -cbar_lim),
                edge_col = ifelse(p_adj < .05, "black", "grey70"),
                edge_size = ifelse(p_adj < .05, 0.15, 0.03))
  p <- ggseg(dat, atlas = "glasser", position = "dispersed", hemisphere = "left",
        mapping = aes(fill = effect, colour = I(edge_col), size = I(edge_size)), adapt_scales = TRUE) +
    theme_brain(text.family = "sans", text.size = FONT_BRAIN_PT) +
    scale_fill_gradientn(limits = c(-cbar_lim, cbar_lim), breaks = c(-cbar_lim, 0, cbar_lim),
                         colours = colorRampPalette(c("blue4", "blue", "white", "red", "red4"))(300),
                         oob = scales::squish) +
    guides(fill = guide_colourbar(title = "Effect size", title.position = "top", barwidth = 5.25, barheight = .5, title.hjust = .5)) +
    theme(legend.position = if(show_legend) "bottom" else "none", legend.title = element_text(size = FONT_LEGEND_TITLE_PT),
          legend.text = element_text(size = FONT_LEGEND_TEXT_PT), plot.margin = margin(-8, 0, -8, 0))
  p <- p + theme(axis.title = element_blank(), axis.text.x = element_blank()) +
    theme(plot.margin = grid::unit(c(-28, 0, -28, 0), "mm"))
  p
}

# Sensorimotor / association poles: bottom and top 30% of H = -z(T1w/T2w) (54 regions each)
assign_poles <- function(h) {
  q <- quantile(h, probs = c(0.30, 0.70), na.rm = TRUE)
  case_when(h <= q[1] ~ "sensorimotor", h >= q[2] ~ "association", TRUE ~ "middle")
}

mean_ci95 <- function(x) {
  x <- x[is.finite(x)]; m <- mean(x); se <- sd(x) / sqrt(length(x))
  c(mean = m, lo = m - qt(.975, df = length(x)-1) * se, hi = m + qt(.975, df = length(x)-1) * se, n = length(x))
}

#-----------------------------------------------------------------------------
# PART B: Load effect sizes, poles and statistics
#-----------------------------------------------------------------------------

## ---- fig4_data ----
assert_file(filename_ST_C); assert_file(filename_maps_C); assert_file(filename_poles_C)
df_ST4 <- read_excel(filename_ST_C, sheet = "ST4")
df_ST11 <- read_excel(filename_ST_C, sheet = "ST11")   # sensorimotor/association summaries and exact KS tests (FDR over 94 tests)
df_ST7 <- read_excel(filename_ST_C, sheet = "ST7")     # profile statistics (rho_cDEC_H, pspin_cDEC_H)

thresholds <- c("Z0p2", "Z0p5", "Z0p8", "Z1p0")
threshold_labels <- c("0.2", "0.5", "0.8", "1.0")
d4 <- df_ST4 %>% filter(gene_set_exp_category == "regional_specificity_threshold", gene_set_exp %in% thresholds)

roi_names <- names(glasser$palette)[1:180]
roi_ST_names <- paste0("L_", gsub("-", ".", roi_names, fixed = TRUE))
neuromaps <- read.csv(filename_maps_C, check.names = FALSE)
neuromaps <- neuromaps[match(roi_names, neuromaps$glasser_region), ]
hierarchy_z <- -zscore(neuromaps$T1T2ratio)
poles <- read.csv(filename_poles_C)
tercile <- poles$pole[match(roi_names, poles$glasser_name)]
tercile[tercile == "transition"] <- "middle"
stopifnot(!anyNA(tercile), identical(tercile, assign_poles(hierarchy_z)), sum(tercile == "sensorimotor") == 54, sum(tercile == "association") == 54)

#-----------------------------------------------------------------------------
# PART C: Panels
#-----------------------------------------------------------------------------

## ---- fig4_panels_a_c ----
#----- Panels a, c: Del and Dup maps per threshold (Fig. 4a,c); pole summaries for b, d
maps_del <- maps_dup <- vector("list", length(thresholds))
source_maps <- list(); source_sa <- list()

for (i in seq_along(thresholds)) {
  th <- thresholds[i]
  xdel <- d4 %>% filter(gene_set_exp == th, cnv_type == "DEL"); xdel <- xdel[match(roi_ST_names, xdel$gene_set_name), ]
  xdup <- d4 %>% filter(gene_set_exp == th, cnv_type == "DUP"); xdup <- xdup[match(roi_ST_names, xdup$gene_set_name), ]
  if (anyNA(xdel$gene_set_name) || anyNA(xdup$gene_set_name)) stop("Expected 180 DEL and DUP regions at ", th)
  es_del <- as.numeric(xdel$Estimate); es_dup <- as.numeric(xdup$Estimate)
  p_del <- as.numeric(xdel$pvalue); p_dup <- as.numeric(xdup$pvalue)
  pa <- bonferroni_del_dup(p_del, p_dup)

  mdel0 <- make_effect_map(es_del, pa$del, 0.30, show_legend = (i == length(thresholds)))
  if (i == length(thresholds)) effect_legend <- cowplot::get_legend(mdel0)
  mdel <- mdel0 + theme(legend.position = "none")
  mdup <- make_effect_map(es_dup, pa$dup, 0.30, show_legend = FALSE)
  maps_del[[i]] <- ggdraw(mdel) + draw_plot(make_pie(sum(pa$del < .05)), x=.36, y=.59, width=.28, height=.37)
  maps_dup[[i]] <- ggdraw(mdup) + draw_plot(make_pie(sum(pa$dup < .05)), x=.36, y=.59, width=.28, height=.37)

  source_maps[[i]] <- bind_rows(
    data.frame(gene_set_exp = th, threshold_z = threshold_labels[i], cnv_type = "DEL", gene_set_name = roi_ST_names, glasser_region = roi_names,
               Estimate = es_del, pvalue = p_del, pvalue_bonferroni_360 = pa$del, significant = pa$del < .05,
               H = hierarchy_z, pole = tercile),
    data.frame(gene_set_exp = th, threshold_z = threshold_labels[i], cnv_type = "DUP", gene_set_name = roi_ST_names, glasser_region = roi_names,
               Estimate = es_dup, pvalue = p_dup, pvalue_bonferroni_360 = pa$dup, significant = pa$dup < .05,
               H = hierarchy_z, pole = tercile)
  )

  for (cnv in c("DEL", "DUP")) {
    es <- if (cnv == "DEL") es_del else es_dup
    for (grp in c("sensorimotor", "association")) {
      z <- mean_ci95(es[tercile == grp])
      source_sa[[length(source_sa)+1]] <- data.frame(threshold=th, threshold_z=threshold_labels[i], CNV=cnv,
                                                      group=grp, mean=z["mean"], ci95_lo=z["lo"], ci95_hi=z["hi"], n=z["n"])
    }
  }
}
source_maps_df <- bind_rows(source_maps)
source_sa_df <- bind_rows(source_sa)
fn_write_SourceData_sheet(source_maps_df, "Fig4a-d", source_data_dir)

# KS statistics from Supplementary Table 11: exact two-sided KS, BH FDR across all 94 sensorimotor-association tests
ks_stats <- df_ST11 %>%
  filter(gene_set_exp_category == "regional_specificity_threshold", gene_set_exp %in% thresholds) %>%
  transmute(threshold = gene_set_exp, threshold_z = threshold_labels[match(gene_set_exp, thresholds)], CNV = cnv_type,
            KS_D = KolmogorovSmirnov_D_stat_sensori_vs_assoc, KS_p = KolmogorovSmirnov_pvalue_sensori_vs_assoc,
            KS_q_FDR_94 = KolmogorovSmirnov_pvalue_sensori_vs_assoc_FDR, KS_star = ifelse(KS_q_FDR_94 < 0.05, "*", ""),
            mean_ST11 = mean_effectsize_sensorimotor, mean_assoc_ST11 = mean_effectsize_association)
stopifnot(nrow(ks_stats) == 8, nrow(df_ST11) == 94,
          max(abs(p.adjust(df_ST11$KolmogorovSmirnov_pvalue_sensori_vs_assoc, "BH") - df_ST11$KolmogorovSmirnov_pvalue_sensori_vs_assoc_FDR)) < 1e-12)
# guard: pole means reproduce ST11
chk <- source_sa_df %>% left_join(ks_stats %>% select(threshold, CNV, mean_ST11, mean_assoc_ST11), by = c("threshold", "CNV")) %>%
  mutate(m = ifelse(group == "sensorimotor", mean_ST11, mean_assoc_ST11))
stopifnot(max(abs(chk$mean - chk$m)) < 1e-10)

source_sa_df <- source_sa_df %>%
  left_join(ks_stats %>% select(threshold, threshold_z, CNV, KS_D, KS_p, KS_q_FDR_94, KS_star), by = c("threshold", "threshold_z", "CNV"))

## ---- fig4_layout_a_c ----
stack_A <- plot_grid(plotlist = maps_del, ncol = 1, rel_heights = c(1, 1, 1, 1))
panel_A_body <- ggdraw() + draw_plot(stack_A, x = 0, y = 0.05, width = 1, height = 0.90)
panel_A <- plot_grid(panel_A_body, effect_legend, ncol = 1, rel_heights = c(1, .08))
stack_C <- plot_grid(plotlist = maps_dup, ncol = 1, rel_heights = c(1, 1, 1, 1))
panel_C_body <- ggdraw() + draw_plot(stack_C, x = 0, y = 0.05, width = 1, height = 0.90)
panel_C <- plot_grid(panel_C_body, effect_legend, ncol = 1, rel_heights = c(1, .08))

## ---- fig4_panels_b_d ----
#----- Panels b, d: pole means + 95% CI, KS stars (Fig. 4b,d)
source_sa_df$threshold_z <- factor(source_sa_df$threshold_z, levels = rev(threshold_labels))
source_sa_df$group <- factor(source_sa_df$group, levels = c("sensorimotor", "association"))
custom_colors <- c("sensorimotor" = COL_SENSORIMOTOR, "association" = COL_ASSOCIATION)

make_sa_plot <- function(cnv, show_y = FALSE, show_legend = FALSE) {
  dd <- source_sa_df %>% filter(CNV == cnv)
  rng <- diff(range(c(dd$ci95_lo, dd$ci95_hi), na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- 0.01
  star_df <- dd %>%
    group_by(threshold_z) %>%
    summarise(y_star = max(ci95_hi, na.rm = TRUE) + 0.1 * rng, label = dplyr::first(KS_star), .groups = "drop") %>%
    filter(label == "*")
  ggplot(dd, aes(x = threshold_z, y = mean, ymin = ci95_lo, ymax = ci95_hi, color = group, group = group)) +
    geom_pointrange(position = position_dodge(width = .32), linewidth = .7, size = .35) +
    geom_text(data = star_df, aes(x = threshold_z, y = y_star, label = label), inherit.aes = FALSE,
              color = "black", fontface = "bold", size = pt_to_geom((FONT_INTERNAL_LABEL_PT + 1) * 2)) +
    coord_flip() +
    scale_color_manual(values = custom_colors, name = NULL, labels = c("sensorimotor", "association")) +
    labs(x = NULL, y = "Effect size") +
    theme_prism(base_size = FONT_BASE_PT) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = if(show_y) element_line() else element_blank(),
      axis.text.y = if(show_y) element_text(size = FONT_AXIS_TEXT_PT) else element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", color = "grey70"),
      panel.grid.minor = element_blank(),
      legend.position = if(show_legend) "bottom" else "none",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.text = element_text(size = FONT_LEGEND_TEXT_PT),
      axis.title = element_text(size = FONT_AXIS_TITLE_PT),
      axis.text.x = element_text(size = FONT_AXIS_TEXT_PT),
      plot.margin = margin(4, 2, 2, 2)
    )
}

panel_B_base <- make_sa_plot("DEL", FALSE, TRUE)
legend_B <- cowplot::get_legend(
  panel_B_base + guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.box = "vertical", legend.box.spacing = grid::unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))
)
panel_B <- plot_grid(panel_B_base + theme(legend.position = "none"), legend_B, ncol = 1, rel_heights = c(1, 0.10))

panel_D_base <- make_sa_plot("DUP", FALSE, TRUE)
legend_D <- cowplot::get_legend(
  panel_D_base + guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.box = "vertical", legend.box.spacing = grid::unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))
)
panel_D <- plot_grid(panel_D_base + theme(legend.position = "none"), legend_D, ncol = 1, rel_heights = c(1, 0.10))

## ---- fig4_panel_e ----
#----- Panel e: rho(cDEC, H) per threshold, FDR across the 4 thresholds (Fig. 4e)
stats <- df_ST7 %>%
  filter(gene_set_exp_category == "regional_specificity_threshold", gene_set_exp %in% thresholds) %>%
  mutate(threshold_z = threshold_labels[match(gene_set_exp, thresholds)]) %>%
  arrange(match(gene_set_exp, thresholds))
if (nrow(stats) != 4) stop("Expected four regional-specificity rows in ST7.")
stats$rho_DeltaZ_H <- stats$rho_cDEC_H                      # cDEC = z(Del) - z(Dup)
stats$q_spin_DeltaZ_H <- p.adjust(stats$pspin_cDEC_H, "BH")  # FDR across the 4 thresholds (as in the legend)
stats$star <- ifelse(stats$q_spin_DeltaZ_H < .05, "*", "")
stats$threshold_z <- factor(stats$threshold_z, levels = rev(threshold_labels))

df_stats_Fig4 <- bind_rows(
  source_sa_df %>% transmute(panel = ifelse(CNV == "DEL", "b", "d"), gene_set_exp = threshold, threshold_z = as.character(threshold_z), cnv_type = CNV,
                             pole = as.character(group), mean_effect_size = mean, ci95_lo, ci95_hi, n_regions = n,
                             KS_D, KS_pvalue_exact = KS_p, KS_q_FDR_94tests = KS_q_FDR_94, star = KS_star),
  data.frame(panel = "e", gene_set_exp = stats$gene_set_exp, threshold_z = as.character(stats$threshold_z), cnv_type = "cDEC",
             rho_cDEC_H = stats$rho_cDEC_H, pspin_cDEC_H = stats$pspin_cDEC_H, q_FDR_4_thresholds = stats$q_spin_DeltaZ_H, star = stats$star)
)
fn_write_SourceData_sheet(df_stats_Fig4, "Fig4b-e_stats", source_data_dir)

stats$value_label <- sprintf("%.2f", stats$rho_DeltaZ_H)
panel_E_base <- ggplot(stats, aes(y = threshold_z, x = rho_DeltaZ_H)) +
  geom_vline(xintercept = 0, linewidth = .4, color = "grey35") +
  geom_col(width = .26, fill = COL_DELTA_BAR) +
  geom_text(aes(x = rho_DeltaZ_H + 0.015, label = value_label), hjust = 0, color = "black", size = pt_to_geom(FONT_ANNOTATION_PT)) +
  geom_text(aes(x = rho_DeltaZ_H + 0.15, label = star), hjust = 0, color = "black", fontface = "bold", size = pt_to_geom(FONT_ANNOTATION_PT * 2)) +
  scale_x_continuous(limits = c(-0.05, max(0.65, max(stats$rho_DeltaZ_H) + 0.14)), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Correlation", y = NULL) +
  theme_prism(base_size = FONT_BASE_PT) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey75"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = FONT_AXIS_TITLE_PT),
    axis.text.x = element_text(size = FONT_AXIS_TEXT_PT),
    plot.margin = margin(4, 2, 2, 2)
  )
panel_E_footer <- ggdraw() +
  draw_label("*", x = 0.43, y = 0.5, hjust = 1, vjust = 0.5, size = 2 * (FONT_MIN_PT + 1), fontface = "bold") +
  draw_label("FDR q < 0.05", x = 0.45, y = 0.5, hjust = 0, vjust = 0.5, size = (FONT_MIN_PT + 1))
panel_E <- plot_grid(panel_E_base, panel_E_footer, ncol = 1, rel_heights = c(1, 0.10))

#-----------------------------------------------------------------------------
# PART D: Assemble and save
#-----------------------------------------------------------------------------

## ---- fig4_assemble ----
# Left-side guide shown once for all five panels: vertical arrow + regional-specificity thresholds
y_pos <- c(0.83, 0.63, 0.43, 0.23)
panel_guide <- ggdraw() +
  draw_line(x = c(0.88, 0.88), y = c(0.90, 0.12),
            arrow = grid::arrow(type = "closed", length = grid::unit(0.08, "in")), color = "grey35", size = 0.4) +
  draw_label("Increasing regional specificity (Z)", x = 0.18, y = 0.51, angle = 90,
             hjust = 0.5, vjust = 0.5, size = FONT_AXIS_TITLE_PT, fontface = "plain")
for (i in seq_along(threshold_labels)) {
  panel_guide <- panel_guide + draw_label(threshold_labels[i], x = 0.59, y = y_pos[i], hjust = 0.5, vjust = 0.5, size = FONT_AXIS_TEXT_PT)
}

panel_A_t <- add_panel_tag(panel_A, "a")
panel_B_t <- add_panel_tag(panel_B, "b")
panel_C_t <- add_panel_tag(panel_C, "c")
panel_D_t <- add_panel_tag(panel_D, "d")
panel_E_t <- add_panel_tag(panel_E, "e")

final_core <- plot_grid(panel_guide, panel_A_t, panel_B_t, panel_C_t, panel_D_t, panel_E_t,
                        nrow = 1, align = "h", rel_widths = c(0.22, 1.08, 0.84, 1.08, 0.84, 0.76))
final_fig <- ggdraw() + draw_plot(final_core, x = 0, y = 0.03, width = 1, height = 0.95)

save_figure(final_fig, "Fig4", FIG_WIDTH_IN, FIG_HEIGHT_IN)
message("Saved Fig. 4 to: ", OUT_DIR)
