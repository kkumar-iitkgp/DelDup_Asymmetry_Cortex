#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 6 (alternative AHBA gene sets: spatial co-expression modules)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: Del vs Dup effect sizes for the 16 co-expression modules (Bonferroni 0.05/32)
#          -> b: module PC1 maps (Glasser-180, z-scored)
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
#       1) SupFigs/output/SupFig06.png, SupFig06.svg (+ checks/SupFig06_checks.txt)
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

FIG_ID <- "SupFig06"
SHEET_NAMES <- c(scatter = "SupFig6a", modulePC1_z = "SupFig6b")

STEM <- "SuppFig_WagstylModules"
M <- paste0("M", 1:16)
p <- profile("gene_expression_AHBA_Wagstyl_16modules", ord = M)
chk(nrow(p) == 16, "Wagstyl: 16 modules joined by id")
p[, sig := sig_class(pD, pU, 32)][, lab := region]
S_NS <- 2.2
pb <- ggplot(p[order(sig != "ns")], aes(U, D)) +
  geom_hline(yintercept = 0, colour = COL_LINE, linewidth = .3) + geom_vline(xintercept = 0, colour = COL_LINE, linewidth = .3) +
  geom_point(aes(fill = sig, size = sig != "ns"), shape = 21, colour = "white", stroke = .15) +
  geom_text_repel(aes(label = lab, fontface = ifelse(sig != "ns", "bold", "plain")), size = pt_to_geom(FONT_ANNOTATION_PT), min.segment.length = 0,
                  segment.colour = "grey50", segment.size = .2, max.overlaps = Inf, seed = 1, box.padding = .35) +
  scale_fill_manual(values = COL_SIG, drop = FALSE, guide = "none") + scale_size_manual(values = c(`FALSE` = S_NS, `TRUE` = 1.2 * S_NS), guide = "none") +
  labs(x = "Effect size Dup", y = "Effect size Del") + theme_ng()
no_smooth(pb, "Wagstyl Del-vs-Dup scatter")
chk(all(!is.na(p$lab)) && nrow(p) == 16, sprintf("all 16 modules labelled; %d significant (bold, 1.2x point size)", sum(p$sig != "ns")))
X <- as.data.frame(fread(file.path(inputs_dir, "inputs_Wagstyl_module_PC1_Glasser180.csv"))); rownames(X) <- X$glasser_region; X <- as.matrix(X[, -1])
chk(identical(sub(" g=.*", "", colnames(X)), M), "module PC1 object: 16 columns M1..M16")
ix <- match(GLASSER_REGION, gsub("[._]", "-", rownames(X)))
chk(!anyNA(ix), "module PC1 maps joined to the 180 Glasser regions by name")
X <- X[ix, ]
mod_map <- function(k) { v <- as.numeric(scale(X[, k]))
  ggseg(tibble(region = GLASSER_REGION, v = pmax(pmin(v, 2.5), -2.5)), atlas = "glasser", hemisphere = "left", view = "lateral",
        mapping = aes(fill = v), colour = NA, adapt_scales = TRUE) +
    scale_fill_gradientn(limits = c(-2.5, 2.5), colours = colorRampPalette(c("blue4", "blue", "white", "red", "red4"))(300), oob = squish, name = "Module PC1 (z)") +
    labs(title = M[k]) + theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
    theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0),
          plot.title = element_text(size = FONT_AXIS_TEXT_PT, hjust = .5, face = "plain")) }
maps <- plot_grid(plotlist = lapply(1:16, mod_map), ncol = 8)
leg <- get_legend(mod_map(1) + theme(legend.position = "bottom") + guides(fill = guide_colourbar(title.position = "top", title.hjust = .5, barwidth = 5, barheight = .35)))
pa <- plot_grid(maps, leg, ncol = 1, rel_heights = c(1, .22))
fig <- plot_grid(pb, pa, ncol = 1, rel_heights = c(0.5, 0.5), labels = c("a", "b"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(c(M, na.omit(p$lab), "Module PC1 (z)", "Effect size Del", "Effect size Dup"))
save_v2(fig, STEM, FIG_WIDTH_IN, 6.4)
write_sd(p[, .(module = region, effect_size_DEL = D, effect_size_DUP = U, p_DEL = pD, p_DUP = pU, significance_class = sig, label = lab)], STEM, "scatter")
write_sd(data.table(region = ROI, setNames(as.data.table(apply(X, 2, function(v) as.numeric(scale(v)))), M)), STEM, "modulePC1_z")
write_checks(STEM)
