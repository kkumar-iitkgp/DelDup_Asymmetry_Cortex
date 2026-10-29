#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 18 (genetic constraint vs similarity to the cortical gradient)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Plots Supplementary Fig. 18 a-e and writes the plotted values to the Source Data sheets
#          -> y axis (all panels): similarity of the gene's cortical expression to the cortical gradient H = -T1w/T2w,
#             i.e. -1 x TMC (TMC = similarity to T1w/T2w); low = sensorimotor, high = association.
#          -> Spearman rho (pairwise complete, all protein-coding genes with both scores); no p-values.
# Key steps:
#       1) PART A: Parameters and input files
#       2) PART B: Load gene-level scores; LOEUF strata
#       3) PART C: Panels a-b (LOEUF)
#       4) PART D: Panels c-e (pHaplo, pTriplo, missense Z)
#       5) PART E: Assemble; save PNG + SVG; write Source Data
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: 1
#
#
#-----
# Input files:
#       1) data/inputs/inputs_SupFig18_gene_scores.csv (from Final_scores_Collins_GeneExpPC_protein_for_scatter_dec2022.RData;
#          source_data_build/06_build_SupFig18_inputs.R in the Zenodo archive)
#
#-----
# Output:
#       1) SupFigs/output/SupFig18.png, SupFig18.svg
#       2) data/source_data_sheets/SupFig18a-e.csv, SupFig18_stats.csv
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

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
output_dir <- paste0(project_dir,"/SupFigs/output")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggExtra)
  library(cowplot)
  library(svglite)
})
source(file.path(project_dir, "functions", "00_figure_style.R"))
source(file.path(project_dir, "functions", "fn_write_SourceData_sheet.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

flag_marginal_type <- "histogram"   # "histogram" (as the previous version) or "density"
flag_fitline_row2 <- TRUE           # linear fit line in panels c-e (as panel a)
fig_width_N <- 7.2; fig_height_N <- 6.5
col_contour_a <- "#3366FF"          # ggplot2 default contour blue (previous version)
col_contour_row2 <- c(pHaplo = "#1B9E77", pTriplo = "#D95F02", missense_Z = "#7570B3")   # Dark2
col_density_b <- "#4C72B0"
xlim_misZ <- c(-5, 7)               # view window for missense Z (coord_cartesian; no data dropped)
loeuf_levels <- c("Loeuf1: 0.0-0.35","Loeuf2: 0.35-1.0","Loeuf3: 1.0-2.0")
ylab_C <- "Similarity to cortical gradient\n(sensorimotor − / association +)"

#-----------------------------------------------------------------------------
# PART B: Load gene-level scores; LOEUF strata
#-----------------------------------------------------------------------------

df_scores <- read.csv(file.path(inputs_dir, "inputs_SupFig18_gene_scores.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df_scores) == 19689, !anyDuplicated(df_scores$ensembl_gene_id))
df_scores$similarity_H <- -1 * df_scores$TMC      # H = -T1w/T2w (Kuldeep's decision, 29 Oct 2026)

# LOEUF strata: same rule as the LOEUF gene-set strata (0-0.35 inclusive; 0.35-1.0 and 1.0-2.0 right-closed)
df_scores <- df_scores %>%
  mutate(LOEUF_bin = case_when(
    LOEUF >= 0 & LOEUF <= 0.35 ~ "Loeuf1: 0.0-0.35",
    LOEUF > 0.35 & LOEUF <= 1.0 ~ "Loeuf2: 0.35-1.0",
    LOEUF > 1.0 & LOEUF <= 2.0 ~ "Loeuf3: 1.0-2.0",
    TRUE ~ NA_character_))
df_scores$LOEUF_bin <- factor(df_scores$LOEUF_bin, levels = loeuf_levels)
df_plot <- df_scores[!is.na(df_scores$similarity_H), ]
ylim_N <- range(df_plot$similarity_H)

fn_rho <- function(x, y) cor(x, y, method = "spearman", use = "pairwise.complete.obs")
fn_n <- function(x, y) sum(!is.na(x) & !is.na(y))
fn_rho_label <- function(r) sprintf("ρ = %.2f", r)

df_stats <- data.frame(panel = c("a", "c", "d", "e"),
                       pair = c("similarity_H ~ LOEUF", "similarity_H ~ pHaplo", "similarity_H ~ pTriplo", "similarity_H ~ missense_Z"),
                       n = c(fn_n(df_plot$LOEUF, df_plot$similarity_H), fn_n(df_plot$pHaplo, df_plot$similarity_H),
                             fn_n(df_plot$pTriplo, df_plot$similarity_H), fn_n(df_plot$missense_Z, df_plot$similarity_H)),
                       spearman_rho = c(fn_rho(df_plot$LOEUF, df_plot$similarity_H), fn_rho(df_plot$pHaplo, df_plot$similarity_H),
                                        fn_rho(df_plot$pTriplo, df_plot$similarity_H), fn_rho(df_plot$missense_Z, df_plot$similarity_H)),
                       gene_universe = "all protein-coding genes with both scores (pairwise complete); 19,689 genes in the scores file")
df_stats_b <- df_plot %>% filter(!is.na(LOEUF_bin)) %>% group_by(LOEUF_bin) %>%
  summarise(n = n(), median_similarity_H = median(similarity_H), .groups = "drop") %>%
  transmute(panel = "b", pair = paste0("similarity_H | ", LOEUF_bin), n, spearman_rho = NA_real_,
            gene_universe = paste0("median similarity_H = ", sprintf("%.3f", median_similarity_H)))
df_stats <- rbind(df_stats, as.data.frame(df_stats_b))

#-----------------------------------------------------------------------------
# PART C: Panels a-b (LOEUF)
#-----------------------------------------------------------------------------

#----- Panel a: LOEUF vs similarity to cortical gradient (Supplementary Fig. 18a)
df_a <- df_plot[!is.na(df_plot$LOEUF), ]
p1 <- ggplot(df_a, aes(x = LOEUF, y = similarity_H)) +
  geom_point(size = 0.2, alpha = 0.3, colour = "grey60") +
  geom_density_2d(colour = col_contour_a, linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "solid", fullrange = TRUE, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.35, linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 1.0, linewidth = 0.4, linetype = "dotted") +
  scale_x_continuous(breaks = c(0, 0.35, 1, 2), labels = c(0, 0.35, 1, 2)) +
  coord_cartesian(xlim = c(0, 2), ylim = ylim_N) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, label = fn_rho_label(df_stats$spearman_rho[1]), size = pt_to_geom(FONT_ANNOTATION_PT)) +
  labs(x = "Intolerance to loss-of-function (LOEUF)", y = ylab_C) +
  theme_ng() + theme(legend.position = "none")
p_a <- ggMarginal(p1, type = flag_marginal_type, margins = "both", fill = "skyblue", colour = "black", size = 5)

#----- Panel b: density of similarity to cortical gradient by LOEUF stratum (Supplementary Fig. 18b)
df_b <- df_plot[!is.na(df_plot$LOEUF_bin), ]
p_b <- ggplot(df_b, aes(x = similarity_H)) +
  geom_density(fill = col_density_b, colour = "black", alpha = 0.6, linewidth = 0.4) +
  facet_wrap(~ LOEUF_bin, ncol = 1, scales = "fixed") +
  labs(x = "Similarity to cortical gradient", y = "Density") +
  theme_bw(base_size = FONT_BASE_PT, base_family = FONT_FAMILY) +
  theme(axis.title = element_text(size = FONT_AXIS_TITLE_PT), axis.text = element_text(size = FONT_AXIS_TEXT_PT),
        strip.background = element_rect(fill = "white"), strip.text = element_text(size = FONT_AXIS_TEXT_PT),
        panel.spacing = unit(0.4, "lines"), panel.grid.minor = element_blank())

#-----------------------------------------------------------------------------
# PART D: Panels c-e (pHaplo, pTriplo, missense Z)
#-----------------------------------------------------------------------------

fn_row2 <- function(df_in, x_col, x_lab, col_contour, rho_N, xlim = NULL) {
  df_in <- df_in[!is.na(df_in[[x_col]]), ]
  p <- ggplot(df_in, aes(x = .data[[x_col]], y = similarity_H)) +
    geom_point(size = 0.2, alpha = 0.3, colour = "grey60") +
    geom_density_2d(colour = col_contour, linewidth = 0.3)
  if (flag_fitline_row2) p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, fullrange = TRUE, colour = "black", linewidth = 0.5)
  p + coord_cartesian(xlim = xlim, ylim = ylim_N) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, label = fn_rho_label(rho_N), size = pt_to_geom(FONT_ANNOTATION_PT)) +
    labs(x = x_lab, y = ylab_C) + theme_ng() + theme(legend.position = "none")
}
#----- Panel c: pHaplo vs similarity to cortical gradient (Supplementary Fig. 18c)
p_c <- fn_row2(df_plot, "pHaplo", "Haploinsufficiency (pHaplo)", col_contour_row2[["pHaplo"]], df_stats$spearman_rho[2])
#----- Panel d: pTriplo vs similarity to cortical gradient (Supplementary Fig. 18d)
p_d <- fn_row2(df_plot, "pTriplo", "Triplosensitivity (pTriplo)", col_contour_row2[["pTriplo"]], df_stats$spearman_rho[3]) + theme(axis.title.y = element_blank())
#----- Panel e: missense Z vs similarity to cortical gradient (Supplementary Fig. 18e)
n_outside_misZ <- sum(df_plot$missense_Z < xlim_misZ[1] | df_plot$missense_Z > xlim_misZ[2], na.rm = TRUE)
message(sprintf("missense Z view %g to %g (coord_cartesian): %d of %d genes outside the view; all genes used for contours, fit and rho",
                xlim_misZ[1], xlim_misZ[2], n_outside_misZ, sum(!is.na(df_plot$missense_Z))))
p_e <- fn_row2(df_plot, "missense_Z", "Missense Z-score", col_contour_row2[["missense_Z"]], df_stats$spearman_rho[4], xlim = xlim_misZ) + theme(axis.title.y = element_blank())

#-----------------------------------------------------------------------------
# PART E: Assemble; save PNG + SVG; write Source Data
#-----------------------------------------------------------------------------

lab_args <- list(label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
row1 <- do.call(plot_grid, c(list(p_a, p_b, nrow = 1, rel_widths = c(1.45, 1), labels = c("a", "b")), lab_args))
row2 <- do.call(plot_grid, c(list(p_c, p_d, p_e, nrow = 1, rel_widths = c(1.18, 1, 1), labels = c("c", "d", "e"), align = "h", axis = "tb"), lab_args))
fig <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1.25, 1)) + theme(plot.background = element_rect(fill = "white", colour = NA))

if (flag_hires == 1) { svglite::svglite(file.path(output_dir, "SupFig18.svg"), width = fig_width_N, height = fig_height_N, bg = "white"); print(fig); invisible(dev.off()) }
ggsave(file.path(output_dir, "SupFig18.png"), fig, width = fig_width_N, height = fig_height_N, units = "in", dpi = if (flag_hires == 1) PNG_DPI else 300, bg = "white")

fn_write_SourceData_sheet(df_plot[, c("ensembl_gene_id", "gene_symbol", "LOEUF", "LOEUF_bin", "pHaplo", "pTriplo", "missense_Z", "TMC", "similarity_H")],
                          "SupFig18a-e", source_data_dir)
fn_write_SourceData_sheet(df_stats, "SupFig18_stats", source_data_dir)
print(df_stats)
