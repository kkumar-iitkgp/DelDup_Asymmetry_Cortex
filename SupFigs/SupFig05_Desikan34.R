#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 5 (alternative AHBA gene sets: Desikan-34 parcellation)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: Del and Dup maps (Desikan-34, outline Bonferroni 0.05/68)
#          -> b: Del vs Dup scatter
#          -> c: sensorimotor vs association densities (H poles, 30th/70th percentiles)
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
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig05.png, SupFig05.svg (+ checks/SupFig05_checks.txt)
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

FIG_ID <- "SupFig05"
SHEET_NAMES <- c(main = "SupFig5a-c")

STEM <- "SuppFig_Desikan34"
d3 <- sheet("ST4")[gene_set_exp == "AHBA_PausLab_Desikan"]; regs <- sort(unique(d3$gene_set_name))
p <- profile("AHBA_PausLab_Desikan", ord = regs)[, `:=`(label = sub("^ctx_", "", region), short = sub("^ctx_lh_", "", region))]
chk(nrow(p) == 34 && all(p$label %in% dk$data$label), "Desikan: 34 regions; all joined to ggseg dk by label (ctx_ prefix removed)")
p[, sig := sig_class(pD, pU, 68)][, `:=`(sD = pmin(pD * 68, 1) < .05, sU = pmin(pU * 68, 1) < .05)]
dk_map <- function(v, s, lim, ttl) {
  dat <- tibble(label = p$label, effect = pmax(pmin(v, lim), -lim), edge = ifelse(s, "black", "grey75"), lw = ifelse(s, .35, .08))
  ggseg(dat, atlas = "dk", hemisphere = "left", mapping = aes(fill = effect, colour = I(edge), size = I(lw))) +
    scale_fill_gradientn(limits = c(-lim, lim), breaks = c(-lim, 0, lim), labels = function(x) sprintf("%g", x),
                         colours = colorRampPalette(c("blue4", "blue", "white", "red", "red4"))(300), oob = squish, name = ttl, na.value = "grey85") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = .5, barwidth = 5, barheight = .35)) +
    theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
    theme(legend.position = "bottom", legend.title = element_text(size = FONT_LEGEND_TITLE_PT), legend.text = element_text(size = FONT_LEGEND_TEXT_PT),
          axis.text = element_blank(), axis.title = element_blank()) }
pa <- plot_grid(dk_map(p$D, p$sD, .3, "Effect size Del"), dk_map(p$U, p$sU, .15, "Effect size Dup"), nrow = 1)
p[, lab := fifelse(sig != "ns", short, NA_character_)]
pb <- dvd_scatter(p$D, p$U, p$sig, p$lab, size = 2); no_smooth(pb, "Desikan Del-vs-Dup scatter")
r6 <- st7_row("AHBA_PausLab_Desikan")
chk(abs(cor(p$D, p$U, method = "spearman") - r6$rho_DEL_DUP_all) < 1e-10, sprintf("Del-Dup rho (legend) %.3f == ST7; p-spin %.4f", r6$rho_DEL_DUP_all, r6$pspin_DEL_DUP_all))
chk(!any(grepl("ctx_", na.omit(p$lab))), sprintf("labels without ctx_lh_ prefix; %d significant regions labelled", sum(!is.na(p$lab))))
# c: pole densities (Desikan H poles) with KS from ST7b_v3
dp <- fread(file.path(inputs_dir, "inputs_Desikan34_poles.csv")); p[, pole := dp$pole[match(short, dp$region)]]
chk(!anyNA(p$pole) && sum(p$pole == "sensorimotor") == 10 && sum(p$pole == "association") == 10, "c: Desikan poles joined by name (10 sensorimotor, 10 association)")
K <- KS_ALT[region_set == "Desikan34_H_poles_30_70"]; chk(nrow(K) == 2, "c: Desikan KS rows for DEL and DUP (inputs_KS_alternative_parcellations.csv)")
pp <- p[pole != "transition"][, pole := factor(pole, levels = names(COL_POLE))]
dens <- function(t) { v <- if (t == "DEL") pp$D else pp$U; k <- K[cnv_type == t]
  chk(abs(unname(ks.test(v[pp$pole == "sensorimotor"], v[pp$pole == "association"])$statistic) - k$KS_D) < 1e-12, sprintf("c %s: KS D from plotted values == ST7b_v3 (D = %.2f, p = %.3f)", t, k$KS_D, k$KS_p_exact))
  g <- ggplot(tibble(v = v, pole = pp$pole), aes(v, fill = pole, colour = pole)) +
    geom_density(alpha = .45, linewidth = .35) + geom_rug(length = unit(.06, "npc"), linewidth = .4) + geom_vline(xintercept = 0, linetype = "dotted", linewidth = .3) +
    scale_fill_manual(values = COL_POLE, labels = c("Sensorimotor", "Association"), name = NULL) + scale_colour_manual(values = COL_POLE, labels = c("Sensorimotor", "Association"), name = NULL) +
    scale_x_continuous(limits = range(v) + c(-1, 1) * 0.6 * diff(range(v)), n.breaks = 4) +
    labs(x = paste("Effect size", if (t == "DEL") "Del" else "Dup"), y = "Density") + theme_ng()
  if (k$KS_p_exact < .05) g <- g + annotate("text", x = -Inf, y = Inf, label = "*", hjust = -0.6, vjust = 1.1, size = pt_to_geom(FONT_PANEL_TAG_PT * 1.4))
  g }
pc <- (dens("DEL") | dens("DUP")) + plot_layout(guides = "collect") & theme(legend.position = "right")
fig <- plot_grid(pa, pb, pc, ncol = 1, rel_heights = c(0.2, 0.6, 0.2), labels = c("a", "b", "c"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(c(na.omit(p$lab), "Effect size Del", "Effect size Dup", "Density", "Sensorimotor", "Association"))
save_v2(fig, STEM, FIG_WIDTH_IN, 9.0)
write_sd(p[, .(region = short, pole, effect_size_DEL = D, effect_size_DUP = U, p_DEL = pD, p_DUP = pU, sig_bonferroni68_DEL = sD, sig_bonferroni68_DUP = sU, significance_class = sig, label = lab)], STEM)
write_checks(STEM)
