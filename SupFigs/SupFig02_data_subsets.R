#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 2 (replication in data subsets)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a-d: Del and Dup effect-size maps for UK Biobank, five other cohorts, discovery and replication halves (outline: Bonferroni 0.05/360)
#          -> e-h: region-wise Del vs Del and Dup vs Dup effect sizes; rho and p-spin from ST6
#          -> Statistics are read from the Supplementary Tables (final numbering); the data frame passed to each
#             ggplot call is written to the Source Data sheet(s) listed in SHEET_NAMES.
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: 1
#
#-----
# UPDATE log:
#       v1 (29 Oct 2026): release version of code_plots_NG_Aug2026/SuppFig_DataSubsets_v2_12oct2026.R
#                         (paths, final Supplementary Table numbering, output names; plotting code unchanged)
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig02.png, SupFig02.svg (+ checks/SupFig02_checks.txt)
#       2) data/source_data_sheets/<sheet>.csv
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

#-----------------------------------------------------------------------------
# Packages, style and helpers
#-----------------------------------------------------------------------------

source(file.path(project_dir, "functions", "00_figure_style.R"))
source(file.path(project_dir, "functions", "00_supfig_helpers.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

FIG_ID <- "SupFig02"
SHEET_NAMES <- c(maps = "SupFig2a-d", scatters = "SupFig2e-h")

STEM <- "SuppFig_DataSubsets"; LIM <- 0.15
SUB <- c(UKBB = "UK Biobank", Other_5cohort = "Five cohorts", Discovery = "Discovery", Replication = "Replication")
P <- lapply(names(SUB), profile); names(P) <- names(SUB)

# ---- maps
src_maps <- list(); txt <- c()
cols <- lapply(names(SUB), function(e) { p <- P[[e]]; sD <- pmin(p$pD * 360, 1) < .05; sU <- pmin(p$pU * 360, 1) < .05
  src_maps[[e]] <<- data.table(subset = SUB[[e]], region = p$region, effect_size_DEL = p$D, sig_bonferroni_DEL = sD, effect_size_DUP = p$U, sig_bonferroni_DUP = sU)
  plot_grid(ggdraw() + draw_label(SUB[[e]], size = FONT_INTERNAL_LABEL_PT),
            plot_grid(plot_grid(vlabel("Del"), effect_map(p$D, sD, LIM), rel_widths = c(.09, 1)),
                      plot_grid(vlabel("Dup"), effect_map(p$U, sU, LIM), rel_widths = c(.09, 1)), ncol = 1),
            ncol = 1, rel_heights = c(.12, 1)) })
txt <- c(txt, SUB, "Del", "Dup", "Effect size")

# ---- scatters (ST6 values; recomputed rho checked)
PAIRS <- tribble(~cnv, ~e1, ~e2, ~x, ~y,
  "del-del", "UKBB", "Other_5cohort", "UKBB", "Other_5cohort", "dup-dup", "UKBB", "Other_5cohort", "UKBB", "Other_5cohort",
  "del-del", "Replication", "Discovery", "Discovery", "Replication", "dup-dup", "Replication", "Discovery", "Discovery", "Replication")
src_sc <- list()
sc <- lapply(seq_len(nrow(PAIRS)), function(k) { q <- PAIRS[k, ]; r <- st6_row(q$cnv, q$e1, q$e2); t <- if (q$cnv == "del-del") "D" else "U"
  xv <- P[[q$x]][[t]]; yv <- P[[q$y]][[t]]; lab <- if (t == "D") "Del" else "Dup"
  chk(abs(cor(xv, yv, method = "spearman") - r$cor_input1_vs_input2) < 1e-10, sprintf("%s %s vs %s: rho from plotted data %.4f == ST6; p-spin %s", lab, q$x, q$y, r$cor_input1_vs_input2, format_pspin(r$pspin_input1_vs_input2)))
  lines <- c(fmt_rho(r$cor_input1_vs_input2), fmt_p(r$pspin_input1_vs_input2)); txt <<- c(txt, lines)
  src_sc[[k]] <<- data.table(cnv_type = lab, x_subset = SUB[[q$x]], y_subset = SUB[[q$y]], region = ROI, x_effect_size = xv, y_effect_size = yv,
                             rho = r$cor_input1_vs_input2, p_spin = r$pspin_input1_vs_input2)
  g <- pvp_scatter(xv, yv, if (t == "D") COL_DEL else COL_DUP, paste0("Effect size ", lab, "\n", SUB[[q$x]]), paste0("Effect size ", lab, "\n", SUB[[q$y]]), lines)
  g <- g + scale_x_continuous(n.breaks = 3); no_smooth(g, paste("scatter", lab, q$x, q$y)); g })
txt <- c(txt, unlist(lapply(seq_len(nrow(PAIRS)), function(k) { t <- if (PAIRS$cnv[k] == "del-del") "Del" else "Dup"; paste0("Effect size ", t, "\n", SUB[c(PAIRS$x[k], PAIRS$y[k])]) })))

# ---- assemble
row1 <- plot_grid(plotlist = cols, nrow = 1, labels = letters[1:4], label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
row1 <- plot_grid(row1, plot_grid(NULL, colourbar_only(LIM, "Effect size"), NULL, nrow = 1, rel_widths = c(1, 1, 1)), ncol = 1, rel_heights = c(1, .22))
row2 <- plot_grid(plotlist = sc, nrow = 1, align = "h", labels = letters[5:8], label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
fig <- plot_grid(row1, row2, ncol = 1, rel_heights = c(.95, 1))
audit_text(txt)
save_v2(fig, STEM, FIG_WIDTH_IN, 5.4)
write_sd(rbindlist(src_maps), STEM, "maps"); write_sd(rbindlist(src_sc), STEM, "scatters"); write_checks(STEM)
