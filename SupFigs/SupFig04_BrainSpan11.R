#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 4 (replication with BrainSpan adult RNA-seq, 11 regions)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: Del vs Dup effect sizes, 11 BrainSpan regions (Bonferroni 0.05/22)
#          -> b: box plots per pole (sensorimotor V1, A1, S1, M1; association IPC, ITC, DLPFC, VLPFC, OFC)
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
#       v1 (29 Oct 2026): release version of code_plots_NG_Aug2026/SuppFig_BrainSpan11_v3_12oct2026.R
#                         (paths, final Supplementary Table numbering, output names; plotting code unchanged)
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig04.png, SupFig04.svg (+ checks/SupFig04_checks.txt)
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

FIG_ID <- "SupFig04"
SHEET_NAMES <- c(main = "SupFig4a-b")

STEM <- "SuppFig_BrainSpan11"
FULL <- c(V1 = "V1 primary visual cortex", A1 = "A1 primary auditory cortex", S1 = "S1 primary somatosensory cortex", M1 = "M1 primary motor cortex",
          DLPFC = "DLPFC dorsolateral prefrontal cortex", VLPC = "VLPFC ventrolateral prefrontal cortex", OFC = "OFC orbitofrontal cortex",
          ACC = "ACC anterior cingulate cortex", IPC = "IPC posterior inferior parietal cortex", ITC = "ITC inferior temporal cortex",
          STC = "STC posterior superior temporal cortex")
SM <- c("V1", "A1", "S1", "M1"); AS <- c("IPC", "ITC", "DLPFC", "VLPC", "OFC")
p <- profile("BrainSpan_adult_RNAseq_11regions", ord = names(FULL))
chk(nrow(p) == 11, "BrainSpan: 11 regions joined by code")
p[, sig := sig_class(pD, pU, 22)][, lab := fifelse(sig != "ns", FULL[region], NA_character_)]
p[, pole := fifelse(region %in% SM, "sensorimotor", fifelse(region %in% AS, "association", "intermediate"))]
chk(sum(p$pole == "sensorimotor") == 4 && sum(p$pole == "association") == 5, "poles (Option B): 4 sensorimotor, 5 association; ACC, STC intermediate")
pa <- dvd_scatter(p$D, p$U, p$sig, p$lab, size = 2.2); no_smooth(pa, "BrainSpan Del-vs-Dup scatter")
chk(abs(cor(p$D, p$U, method = "spearman") - st7_row("BrainSpan_adult_RNAseq_11regions")$rho_DEL_DUP_all) < 1e-10, "Del-Dup rho (legend only) == ST7 rho_DEL_DUP_all")
kb <- KS_ALT[region_set == "BrainSpan11_poles_OptionB"]; chk(nrow(kb) == 2, "KS_ALT: BrainSpan KS rows for DEL and DUP")
pp <- p[pole != "intermediate"]
bx <- function(t) { v <- if (t == "DEL") pp$D else pp$U; k <- kb[cnv_type == t]
  chk(abs(unname(ks.test(v[pp$pole == "sensorimotor"], v[pp$pole == "association"])$statistic) - k$KS_D) < 1e-12,
      sprintf("%s: KS D from plotted values == KS_ALT (D = %.3f, p = %.3f, min attainable %.3f)", t, k$KS_D, k$KS_p_exact, k$min_attainable_p))
  g <- ggplot(tibble(v = v, pole = factor(pp$pole, levels = names(COL_POLE))), aes(pole, v, fill = pole)) +
    geom_hline(yintercept = 0, colour = COL_LINE, linewidth = .3) +
    geom_boxplot(alpha = .45, outlier.shape = NA, width = .55, linewidth = .35) +
    geom_jitter(aes(colour = pole), width = .08, height = 0, size = 1.6, seed = 1) +
    scale_fill_manual(values = COL_POLE) + scale_colour_manual(values = COL_POLE) +
    scale_x_discrete(labels = c("Sensori-\nmotor", "Asso-\nciation")) +
    labs(x = NULL, y = paste("Effect size", if (t == "DEL") "Del" else "Dup")) + theme_ng() + theme(legend.position = "none")
  if (k$KS_p_exact < .05) g <- g + annotate("text", x = 1.5, y = Inf, label = "*", vjust = 1.2, size = pt_to_geom(FONT_PANEL_TAG_PT * 1.4))
  g }
pb <- bx("DEL") | bx("DUP")
fig <- plot_grid(pa, pb, ncol = 1, rel_heights = c(0.7, 0.3), labels = c("a", "b"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(c(na.omit(p$lab), "Effect size Del", "Effect size Dup", "Sensori-\nmotor", "Asso-\nciation"))
save_v2(fig, STEM, FIG_WIDTH_IN, 7.0)
write_sd(p[, .(region_code = region, region_name = FULL[region], pole, effect_size_DEL = D, effect_size_DUP = U, p_DEL = pD, p_DUP = pU, significance_class = sig, label = lab)], STEM)
write_checks(STEM)
