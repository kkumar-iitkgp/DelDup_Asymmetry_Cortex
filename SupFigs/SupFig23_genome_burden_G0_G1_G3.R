#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 23 (Genome Burden specifications G0, G1 and G3)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: regional effect sizes under G0, G1, G3
#          -> b: G0 vs G1 and G3 vs G1 (rho, p-spin from ST6)
#          -> c: cDEC vs hierarchy rank per specification (rho, p-spin from ST7)
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
#       1) SupFigs/output/SupFig23.png, SupFig23.svg (+ checks/SupFig23_checks.txt)
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

FIG_ID <- "SupFig23"
SHEET_NAMES <- c(a = "SupFig23a", b = "SupFig23b", c = "SupFig23c")

STEM <- "SuppFig_GenomeBurden_G0G1G3"
SPEC <- c(GB_G0_noGenomeBurden = "G0", GB_G1_primary = "G1", GB_G3_totalCNVgenes = "G3")
P <- lapply(names(SPEC), profile, src = "ST16"); names(P) <- names(SPEC)
H <- .pm$H[match(ROI, .pm$region)]; chk(!anyNA(H), "H joined to the 180 regions by name")
txt <- c()
# a
la <- rbindlist(lapply(names(SPEC), function(e) data.table(spec = SPEC[[e]], region = ROI, Del = P[[e]]$D, Dup = P[[e]]$U))) %>%
  pivot_longer(c(Del, Dup), names_to = "cnv", values_to = "es") %>% mutate(spec = factor(spec, levels = SPEC))
pa <- ggplot(la, aes(spec, es)) + geom_line(aes(group = region), colour = "grey80", linewidth = .2, alpha = .6) +
  geom_violin(aes(fill = cnv), alpha = .25, colour = NA, width = .7) + geom_boxplot(fill = "white", width = .12, outlier.shape = NA, linewidth = .3) +
  facet_wrap(~cnv, scales = "free_y", labeller = as_labeller(c(Del = "Del", Dup = "Dup"))) +
  scale_fill_manual(values = c(Del = COL_DEL, Dup = COL_DUP), guide = "none") +
  labs(x = NULL, y = "Effect size") + theme_ng() + theme(strip.background = element_blank(), strip.text = element_text(size = FONT_AXIS_TITLE_PT))
txt <- c(txt, SPEC, "Effect size", "Del", "Dup")
# b
sc <- function(e, t) { r <- st6_row(if (t == "D") "del-del" else "dup-dup", "GB_G1_primary", e)
  x <- P$GB_G1_primary[[t]]; y <- P[[e]][[t]]; cn <- if (t == "D") "Del" else "Dup"
  chk(abs(cor(x, y, method = "spearman") - r$cor_input1_vs_input2) < 1e-10, sprintf("%s %s vs G1: rho %.3f from plotted data == ST6; p-spin %s", cn, SPEC[[e]], r$cor_input1_vs_input2, format_pspin(r$pspin_input1_vs_input2)))
  lines <- c(fmt_rho(r$cor_input1_vs_input2), fmt_p(r$pspin_input1_vs_input2)); txt <<- c(txt, lines, paste0("Effect size ", cn, "\nG1"), paste0("Effect size ", cn, "\n", SPEC[[e]]))
  g <- pvp_scatter(x, y, if (t == "D") COL_DEL else COL_DUP, paste0("Effect size ", cn, "\nG1"), paste0("Effect size ", cn, "\n", SPEC[[e]]), lines, identity = TRUE) + scale_x_continuous(n.breaks = 3)
  no_smooth(g, paste(cn, SPEC[[e]], "vs G1")); list(g = g, sd = data.table(cnv_type = cn, spec = SPEC[[e]], region = ROI, G1 = x, other = y, rho = r$cor_input1_vs_input2, p_spin = r$pspin_input1_vs_input2)) }
B <- list(sc("GB_G0_noGenomeBurden", "D"), sc("GB_G3_totalCNVgenes", "D"), sc("GB_G0_noGenomeBurden", "U"), sc("GB_G3_totalCNVgenes", "U"))
pb <- wrap_plots(lapply(B, `[[`, "g"), nrow = 1)
# c (Fig. 2 style): x = hierarchy rank (1 = most sensorimotor, 180 = most association), points coloured by cDEC, lm line
z <- function(v) as.numeric(scale(v))
hr <- rank(H, ties.method = "average")
CD <- lapply(names(SPEC), function(e) z(P[[e]]$D) - z(P[[e]]$U)); names(CD) <- names(SPEC)
dlim <- ceiling(max(abs(unlist(CD))) * 2) / 2   # shared colour limit (Fig. 2 rule)
C <- lapply(names(SPEC), function(e) { cd <- CD[[e]]; r <- st7_row(e)
  chk(abs(cor(cd, hr, method = "spearman") - r$rho_cDEC_H) < 1e-10, sprintf("%s: rho(cDEC, hierarchy rank) %.3f from plotted data == ST7; p-spin %s", SPEC[[e]], r$rho_cDEC_H, format_pspin(r$pspin_cDEC_H)))
  lines <- c(fmt_rho(r$rho_cDEC_H), fmt_p(r$pspin_cDEC_H)); txt <<- c(txt, lines)
  g <- ggplot(tibble(hr = hr, cd = cd), aes(hr, cd, colour = cd)) +
    annotate("rect", xmin = 1, xmax = 54, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.12) +
    annotate("rect", xmin = 127, xmax = 180, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.12) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = .45, colour = "grey30") +
    geom_point(size = 1.1, alpha = .95) +
    scale_colour_gradient2(low = COL_SENSORIMOTOR, mid = COL_NEUTRAL, high = COL_ASSOCIATION, midpoint = 0, limits = c(-dlim, dlim), oob = squish, guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -.05, vjust = 1.2, label = paste(lines, collapse = "\n"), size = pt_to_geom(FONT_ANNOTATION_PT), lineheight = .95) +
    scale_x_continuous(breaks = c(1, 180), labels = c("← Sensorimotor", "Association →"), expand = expansion(mult = .01)) +
    coord_cartesian(ylim = c(min(cd), max(cd) + diff(range(cd)) * .4), clip = "off") +
    labs(title = SPEC[[e]], x = "Cortical hierarchy rank", y = "cDEC") +
    theme_ng() + theme(plot.title = element_text(size = FONT_AXIS_TITLE_PT, hjust = .5), plot.margin = margin(4, 8, 4, 4), axis.text.x = element_text(hjust = c(0, 1), size = FONT_AXIS_TEXT_PT - 1))
  list(g = g, sd = data.table(spec = SPEC[[e]], region = ROI, hierarchy_rank = hr, cDEC = cd, rho = r$rho_cDEC_H, p_spin = r$pspin_cDEC_H)) })
pc <- wrap_plots(lapply(C, `[[`, "g"), nrow = 1)
txt <- c(txt, "Cortical hierarchy rank", "cDEC", "← Sensorimotor", "Association →")
fig <- plot_grid(pa, pb, pc, ncol = 1, rel_heights = c(1, 1, 1), labels = c("a", "b", "c"), label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
audit_text(txt)
save_v2(fig, STEM, FIG_WIDTH_IN, 8.4)
write_sd(as.data.table(la), STEM, "a"); write_sd(rbindlist(lapply(B, `[[`, "sd")), STEM, "b"); write_sd(rbindlist(lapply(C, `[[`, "sd")), STEM, "c")
write_checks(STEM)
