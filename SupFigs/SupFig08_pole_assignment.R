#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 8 (sensorimotor and association poles)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: Glasser-180 map of the poles (<= 30th / >= 70th percentile of H = -T1w/T2w)
#          -> b: H by rank
#          -> Statistics are read from the Supplementary Tables (final numbering); the data frame passed to each
#             ggplot call is written to the Source Data sheet(s) listed in SHEET_NAMES.
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
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig08.png, SupFig08.svg (+ checks/SupFig08_checks.txt)
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

FIG_ID <- "SupFig08"
SHEET_NAMES <- c(main = "SupFig8a-b")

STEM <- "SuppFig_PoleAssignment"
pm <- .pm[match(ROI, region)]; q <- quantile(pm$H, c(.30, .70), type = 7)
chk(sum(pm$H <= q[1]) == 54 && sum(pm$H >= q[2]) == 54 && all((pm$pole == "sensorimotor") == (pm$H <= q[1])) && all((pm$pole == "association") == (pm$H >= q[2])),
    sprintf("poles == bottom / top 30%% of H (cut-offs %.3f / %.3f): 54 sensorimotor, 54 association, 72 transition", q[1], q[2]))
pc <- c(sensorimotor = COL_SENSORIMOTOR, transition = COL_SIG_NS, association = COL_ASSOCIATION)
pl <- c(sensorimotor = "Sensorimotor", transition = "Transition", association = "Association")
pa <- ggseg(tibble(region = GLASSER_REGION, pole = factor(pm$pole, levels = names(pc))), atlas = "glasser", position = "dispersed", hemisphere = "left",
            mapping = aes(fill = pole), colour = "grey40", size = .05, adapt_scales = TRUE) +
  scale_fill_manual(values = pc, labels = pl, name = NULL, drop = FALSE, na.translate = FALSE) +
  theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
  theme(legend.position = "bottom", legend.text = element_text(size = FONT_LEGEND_TEXT_PT), axis.title = element_blank(), axis.text = element_blank())
dd <- pm[, .(region = glasser_name, H, T1T2ratio, percentile, pole)][, rank := rank(H, ties.method = "first")]
pb <- ggplot(dd, aes(rank, H, colour = factor(pole, levels = names(pc)))) +
  geom_hline(yintercept = q, linetype = "dashed", linewidth = .35, colour = "grey35") + geom_point(size = 1.1) +
  scale_colour_manual(values = pc, guide = "none") + labs(x = "Hierarchy rank", y = "H (−T1w/T2w)") + theme_ng()
fig <- plot_grid(pa, pb, nrow = 1, labels = c("a", "b"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(c(pl, "Hierarchy rank", "H (−T1w/T2w)"))
save_v2(fig, STEM, FIG_WIDTH_IN, 3.0)
write_sd(dd[, .(region, H, T1T2ratio, percentile, rank, pole)], STEM); write_checks(STEM)
