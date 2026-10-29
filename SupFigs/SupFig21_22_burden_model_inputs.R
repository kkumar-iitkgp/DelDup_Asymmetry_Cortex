#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Figs. 21 and 22 (Gene-Set Burden X1 per individual; burden coefficients and CNV characteristics)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Supplementary Fig. 21: number of genes of the regional gene set carried per individual (X1), 6 regions
#             (3 sensorimotor, 3 association), Del and Dup (individuals per X1 bin, log scale).
#          -> Supplementary Fig. 22: a) beta1 (inside set) and beta2 (outside set) across the 180 regions, Genome Burden G1 (ST16);
#             b) genes encompassed per unique genic CNV; c) unique carriers of exactly that CNV (1, 2, ..., 20, > 20).
#          -> The counts were computed from individual-level data (source_data_build/05_build_SupFig21_22_aggregates.R in the Zenodo archive);
#             this script plots the aggregate counts. Plotting code as in the 14 Oct 2026 scripts.
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
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (ST16, G1 beta1 / beta2)
#       2) data/inputs/inputs_SupFig21_X1_histograms.csv, inputs_SupFig22_genes_per_CNV_bins.csv, inputs_SupFig22_carriers_per_CNV_bins.csv
#
#-----
# Output:
#       1) SupFigs/output/SupFig21.png/.svg, SupFig22.png/.svg
#       2) data/source_data_sheets/SupFig21.csv, SupFig22a.csv, SupFig22b.csv, SupFig22c.csv
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
inputs_dir <- file.path(project_dir, "data", "inputs")
source_data_dir <- file.path(project_dir, "data", "source_data_sheets")
OUT_DIR <- file.path(project_dir, "SupFigs", "output")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(data.table); library(readxl)
  library(cowplot); library(patchwork); library(scales)
})
source(file.path(project_dir, "functions", "00_figure_style.R"))
source(file.path(project_dir, "functions", "fn_write_SourceData_sheet.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

XL <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
TYPES <- c("DEL", "DUP"); TLAB <- c(DEL = "Del", DUP = "Dup"); TCOL <- c(Del = COL_DEL, Dup = COL_DUP)
N <- 258292
BINS <- c("0", "1", "2", "3", "4", "5–9", "≥10"); GB <- c(as.character(1:9), "≥10"); RB2 <- c(as.character(1:20), ">20")
SEL <- c(L_V1 = "V1 (primary visual)", L_4 = "4 (primary motor)", L_3b = "3b (primary somatosensory)", L_46 = "46 (dorsolateral prefrontal)", L_PGi = "PGi (inferior parietal)", L_TE1a = "TE1a (lateral temporal)")

#-----------------------------------------------------------------------------
# PART B: Data
#-----------------------------------------------------------------------------

f2 <- fread(file.path(inputs_dir, "inputs_SupFig21_X1_histograms.csv"), encoding = "UTF-8")[, .(TYPE = cnv_type, region, bin = factor(X1_bin, levels = BINS), n_individuals)]
stopifnot(all(f2[, sum(n_individuals), by = .(TYPE, region)]$V1 == N), setequal(f2$region, names(SEL)))
st16 <- as.data.table(read_excel(XL, "ST16"))[gene_set_exp == "GB_G1_primary"]
pb <- st16[, .(TYPE = cnv_type, region = gene_set_name, beta1 = Estimate, beta2 = Estimate_outside)]
stopifnot(nrow(pb) == 360, !anyNA(pb))
pc <- fread(file.path(inputs_dir, "inputs_SupFig22_genes_per_CNV_bins.csv"), encoding = "UTF-8")[, .(TYPE = cnv_type, bin = factor(genes_per_CNV_bin, levels = GB), N = n_unique_cnvs, pct)]
pr <- fread(file.path(inputs_dir, "inputs_SupFig22_carriers_per_CNV_bins.csv"), encoding = "UTF-8")[, .(TYPE = cnv_type, bin = factor(unique_carriers_bin, levels = RB2), N = n_unique_cnvs, pct)]

#-----------------------------------------------------------------------------
# PART C: Supplementary Fig. 21 (X1 per individual, 6 regions)
#-----------------------------------------------------------------------------

## ---- supfig21 ----
f2[, `:=`(cnv = TLAB[TYPE], panel = factor(SEL[region], levels = SEL))]
fn_write_SourceData_sheet(f2[, .(cnv_type = cnv, region, region_label = as.character(panel), X1_bin = as.character(bin), n_individuals)], "SupFig21", source_data_dir)
fig2 <- ggplot(f2, aes(bin, n_individuals, fill = cnv)) + geom_col(width = .75) + facet_grid(cnv ~ panel, labeller = label_wrap_gen(14)) +
  scale_fill_manual(values = TCOL, guide = "none") + scale_y_log10(labels = comma, breaks = 10^(0:5)) +
  labs(x = "Genes in the regional set carried per individual (X1)", y = "Individuals") + theme_ng() +
  theme(strip.background = element_blank(), strip.text = element_text(size = FONT_AXIS_TEXT_PT), axis.text.x = element_text(size = FONT_MIN_PT, angle = 45, hjust = 1))
save_figure(fig2 + theme(plot.background = element_rect(fill = "white", colour = NA)), "SupFig21", FIG_WIDTH_IN * 1.15, 4.2, out_dir = OUT_DIR)

#-----------------------------------------------------------------------------
# PART D: Supplementary Fig. 22 (beta1 / beta2; genes per CNV; carriers per CNV)
#-----------------------------------------------------------------------------

## ---- supfig22 ----
pbl <- melt(pb, id.vars = c("TYPE", "region"), variable.name = "term", value.name = "beta")
pbl[, term := factor(c(beta1 = "β1\ninside set", beta2 = "β2\noutside set")[as.character(term)], levels = c("β1\ninside set", "β2\noutside set"))]
pbl[, cnv := TLAB[TYPE]]
fn_write_SourceData_sheet(pb[, .(cnv_type = TYPE, region, beta1_inside = beta1, beta2_outside = beta2)], "SupFig22a", source_data_dir)
fb <- ggplot(pbl, aes(term, beta, fill = cnv)) + geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_violin(alpha = .35, colour = NA, scale = "width", width = .8) + geom_jitter(aes(colour = cnv), width = .12, height = 0, size = .5, alpha = .55) +
  geom_boxplot(width = .18, outlier.shape = NA, fill = "white", alpha = .8, linewidth = .35) +
  facet_wrap(~ cnv, nrow = 1) + scale_fill_manual(values = TCOL, guide = "none") + scale_colour_manual(values = TCOL, guide = "none") +
  labs(x = NULL, y = "Effect size per gene") + theme_ng() + theme(strip.background = element_blank(), strip.text = element_text(size = FONT_AXIS_TITLE_PT))
fn_back_to_back_bars <- function(d, xlab, leg = c(.88, .88), lab_pt = FONT_MIN_PT) { d <- copy(d)[, cnv := TLAB[TYPE]][, y := ifelse(TYPE == "DUP", N, -N)]
  lim <- max(d$N) * 1.25
  ggplot(d, aes(bin, y, fill = cnv)) + geom_col(width = .75) + geom_hline(yintercept = 0, linewidth = .35) +
    geom_text(aes(label = sprintf("%.1f%%", pct), vjust = ifelse(y > 0, -.3, 1.3)), size = pt_to_geom(lab_pt)) +
    scale_fill_manual(values = TCOL, name = NULL) + scale_y_continuous(labels = function(x) comma(abs(x)), limits = c(-lim, lim)) +
    labs(x = xlab, y = "Unique CNVs") + theme_ng() + theme(legend.position = "inside", legend.position.inside = leg) }
fn_write_SourceData_sheet(pc[, .(cnv_type = TYPE, genes_per_CNV_bin = as.character(bin), n_unique_cnvs = N, pct)], "SupFig22b", source_data_dir)
fc <- fn_back_to_back_bars(pc, "Genes encompassed per CNV")
fn_write_SourceData_sheet(pr[, .(cnv_type = TYPE, unique_carriers_bin = as.character(bin), n_unique_cnvs = N, pct)], "SupFig22c", source_data_dir)
fr <- local({ d <- copy(pr)[, cnv := TLAB[TYPE]][, y := ifelse(TYPE == "DUP", N, -N)]   # pseudo-log axis so that bins 3-20 are visible; labels = counts
  lim <- max(d$N) * 3
  ggplot(d, aes(bin, y, fill = cnv)) + geom_col(width = .75) + geom_hline(yintercept = 0, linewidth = .35) +
    geom_text(aes(label = comma(N), vjust = ifelse(y > 0, -.35, 1.35)), size = pt_to_geom(FONT_MIN_PT - 1.5)) +
    scale_fill_manual(values = TCOL, name = NULL) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10), breaks = c(-10000, -1000, -100, -10, 0, 10, 100, 1000, 10000),
                       labels = function(x) comma(abs(x)), limits = c(-lim, lim)) +
    labs(x = "Unique carriers per CNV", y = "Unique CNVs") + theme_ng() + theme(legend.position = "inside", legend.position.inside = c(.88, .85)) })
lab_args <- list(label_size = FONT_PANEL_TAG_PT, label_fontface = "bold", hjust = -0.2, vjust = 1.1)
top <- do.call(plot_grid, c(list(fb, fc, ncol = 2, labels = c("a", "b"), align = "h", axis = "tb"), lab_args))
fig1 <- plot_grid(top, do.call(plot_grid, c(list(fr, labels = "c"), lab_args)), ncol = 1, rel_heights = c(1, 0.95)) + theme(plot.background = element_rect(fill = "white", colour = NA))
save_figure(fig1, "SupFig22", FIG_WIDTH_IN * 1.15, 7.6, out_dir = OUT_DIR)
