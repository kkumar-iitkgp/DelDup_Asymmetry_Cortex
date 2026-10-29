#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 3 (replication with Human Protein Atlas RNA-seq, 73 regions)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: Del vs Dup effect sizes, 73 HPA regions (significance classes, Bonferroni 0.05/146)
#          -> b: sensorimotor vs association core densities; KS from inputs_KS_alternative_parcellations.csv
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
#       v1 (29 Oct 2026): release version of code_plots_NG_Aug2026/SuppFig_HPA73_v3_12oct2026.R
#                         (paths, final Supplementary Table numbering, output names; plotting code unchanged)
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig03.png, SupFig03.svg (+ checks/SupFig03_checks.txt)
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

FIG_ID <- "SupFig03"
SHEET_NAMES <- c(main = "SupFig3a-b")

STEM <- "SuppFig_HPA73"
hn <- fread(file.path(inputs_dir, "inputs_HPA73_regions.csv"))$region_name
ids <- paste0("HPA_nTPM_region", seq_along(hn))
p <- profile("hpa_RNAseq_73regions", ord = ids)[, name := hn]
poles <- fread(file.path(inputs_dir, "inputs_HPA73_core_poles.csv")); p[, pole := poles$pole[match(name, poles$region)]]
chk(nrow(p) == 73 && !anyNA(p$pole) && sum(p$pole == "sensorimotor core") == 10 && sum(p$pole == "association core") == 9,
    "HPA-73: 73 regions joined by id -> name -> core pole (10 sensorimotor, 9 association)")
p[, sig := sig_class(pD, pU, 146)]
p[, lab := fifelse(sig != "ns" & pole != "neither", name, NA_character_)]
chk(all(is.na(p$lab[p$sig == "ns"])) && all(p$pole[!is.na(p$lab)] != "neither"), sprintf("labels: %d significant core regions; n.s. and non-core points unlabelled", sum(!is.na(p$lab))))
pa <- dvd_scatter(p$D, p$U, p$sig, p$lab, size = 1.8); no_smooth(pa, "HPA Del-vs-Dup scatter")
chk(abs(cor(p$D, p$U, method = "spearman") - st7_row("hpa_RNAseq_73regions")$rho_DEL_DUP_all) < 1e-10, "Del-Dup rho (legend only) == ST7 rho_DEL_DUP_all")

core <- p[pole != "neither"][, pole := factor(sub(" core", "", pole), levels = names(COL_POLE))]
kb <- KS_ALT[region_set == "HPA73_anatomical_cores"]; chk(nrow(kb) == 2, "KS_ALT: HPA core KS rows for DEL and DUP")
dens <- function(t) { v <- if (t == "DEL") core$D else core$U; k <- kb[cnv_type == t]
  chk(abs(unname(ks.test(v[core$pole == "sensorimotor"], v[core$pole == "association"])$statistic) - k$KS_D) < 1e-12, sprintf("%s: KS D from plotted values == KS_ALT (D = %.3f, p = %.3f)", t, k$KS_D, k$KS_p_exact))
  g <- ggplot(tibble(v = v, pole = core$pole), aes(v, fill = pole, colour = pole)) +
    geom_density(alpha = .45, linewidth = .35, adjust = 1.2) + geom_rug(length = unit(.06, "npc"), linewidth = .4) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .3) +
    scale_fill_manual(values = COL_POLE, labels = c("Sensorimotor core", "Association core"), name = NULL) +
    scale_colour_manual(values = COL_POLE, labels = c("Sensorimotor core", "Association core"), name = NULL) +
    scale_x_continuous(limits = range(v) + c(-1, 1) * 0.6 * diff(range(v)), n.breaks = 4) + labs(x = paste("Effect size", if (t == "DEL") "Del" else "Dup"), y = "Density") + theme_ng() + theme(legend.position = "bottom")
  if (k$KS_p_exact < .05) g <- g + annotate("text", x = Inf, y = Inf, label = "*", hjust = 1.5, vjust = 1.2, size = pt_to_geom(FONT_PANEL_TAG_PT * 1.4))
  g }
pb <- (dens("DEL") | dens("DUP")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
fig <- plot_grid(pa, pb, ncol = 1, rel_heights = c(0.7, 0.3), labels = c("a", "b"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(c(na.omit(p$lab), "Effect size Del", "Effect size Dup", "Density", "Sensorimotor core", "Association core"))
save_v2(fig, STEM, FIG_WIDTH_IN, 7.4)
write_sd(p[, .(region_id = region, region_name = name, core_pole = pole, effect_size_DEL = D, effect_size_DUP = U, p_DEL = pD, p_DUP = pU, significance_class = sig, label = lab)], STEM)
write_checks(STEM)
