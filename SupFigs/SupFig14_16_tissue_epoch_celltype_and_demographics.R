#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Figs. 14 and 16 (tissue, epoch, cell type; age, sex, ancestry, assessment)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Rows: name
#          -> Del map
#          -> Del density
#          -> Dup map
#          -> Dup density; maps outline Bonferroni 0.05/360; stars = ST11 FDR q < 0.05
#          -> Column 1 (Supplementary Fig. 14): median regional gene-set size and unique Del / Dup genes (ST4)
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
#       1) SupFigs/output/SupFig14.png, SupFig14.svg (+ checks/SupFig14_checks.txt)
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

FIG_ID <- "SupFig14"
OUT_MAP <- c(SuppFig_TissueTemporalCellType = "SupFig14", SuppFig_Demographics_AssessmentType = "SupFig16")
SHEET_NAMES <- c("SuppFig_TissueTemporalCellType:main" = "SupFig14a-c", "SuppFig_TissueTemporalCellType:counts" = "SupFig14_counts", "SuppFig_Demographics_AssessmentType:main" = "SupFig16a-b")

LAB_PT <- FONT_INTERNAL_LABEL_PT - 1; TICKS_DEL <- c(-.3, -.2, -.1, 0); TICKS_DUP <- c(-.1, -.05, 0, .05)
vlabel <- function(txt) ggdraw() + draw_label(txt, angle = 90, size = LAB_PT, fontfamily = FONT_FAMILY)
LIM_DEL <- .30; LIM_DUP <- .15; XL_DEL <- c(-.30, .08); XL_DUP <- c(-.13, .08)
gc <- function(x) paste0("genetic_constraint_", x)

row_plot <- function(exp, label, last = FALSE) { p <- profile(exp); sD <- pmin(p$pD * 360, 1) < .05; sU <- pmin(p$pU * 360, 1) < .05
  stD <- ks_star(exp, "DEL"); stU <- ks_star(exp, "DUP")
  list(g = plot_grid(vlabel(label), effect_map(p$D, sD, LIM_DEL), density_old(p$D, XL_DEL, stD, if (last) TICKS_DEL, FONT_AXIS_TEXT_PT - 1), effect_map(p$U, sU, LIM_DUP), density_old(p$U, XL_DUP, stU, if (last) TICKS_DUP, FONT_AXIS_TEXT_PT - 1),
                     nrow = 1, rel_widths = c(.13, 1, .75, 1, .75)),
       sd = data.table(row = label, gene_set_exp = exp, region = ROI, pole = POLE, effect_size_DEL = p$D, sig_bonferroni360_DEL = sD,
                       effect_size_DUP = p$U, sig_bonferroni360_DUP = sU, KS_star_DEL = stD, KS_star_DUP = stU)) }
count_bar <- function(exps, labels) {
  z <- sheet("ST4")[gene_set_exp %in% exps, .(n_unique = median(n_unique_genes), set = median(n_gene_set_size)), by = .(gene_set_exp, cnv_type)]
  dd <- rbind(z[cnv_type == "DEL", .(gene_set_exp, Type = "Gene set", value = set)], z[, .(gene_set_exp, Type = c(DEL = "Del", DUP = "Dup")[cnv_type], value = n_unique)]) %>%
    mutate(label = factor(labels[match(gene_set_exp, exps)], levels = rev(labels)), Type = factor(Type, levels = c("Dup", "Del", "Gene set")))
  list(g = ggplot(dd, aes(value, label, fill = Type)) + geom_col(position = position_dodge(.8), width = .75) +
         scale_fill_manual(values = c(`Gene set` = "grey65", Del = COL_DEL, Dup = COL_DUP), breaks = c("Gene set", "Del", "Dup"), name = NULL) +
         scale_x_continuous(n.breaks = 3, expand = expansion(mult = c(0, .05))) + labs(x = "Genes", y = NULL) +
         theme_ng() + theme(legend.position = "none", plot.margin = margin(14, 4, 2, 12), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                            axis.text.x = element_text(size = FONT_MIN_PT), axis.title.x = element_text(size = FONT_MIN_PT), legend.text = element_text(size = FONT_MIN_PT),
                            legend.key.size = unit(3, "mm")),
       sd = as.data.table(dd)[, label := as.character(label)]) }
header <- function(with_col1) { h <- plot_grid(NULL, ggdraw() + draw_label("Del", size = LAB_PT), ggdraw() + draw_label("Dup", size = LAB_PT),
                                             nrow = 1, rel_widths = c(.13, 1.75, 1.75)); if (with_col1) plot_grid(NULL, h, nrow = 1, rel_widths = c(.5, 3.63)) else plot_grid(NULL, h, nrow = 1, rel_widths = c(.04, 1)) }
count_key <- function() get_legend(ggplot(tibble(x = 1:3, Type = factor(c("Gene set", "Del", "Dup"), levels = c("Gene set", "Del", "Dup"))), aes(x, x, fill = Type)) + geom_col() +
  scale_fill_manual(values = c(`Gene set` = "grey65", Del = COL_DEL, Dup = COL_DUP), name = NULL) + theme_ng() + theme(legend.position = "bottom", legend.text = element_text(size = FONT_LEGEND_TEXT_PT)))
legend_row <- function(with_col1) plot_grid(plotlist = c(if (with_col1) list(count_key()), list(colourbar_only(LIM_DEL, "Effect size Del"), colourbar_only(LIM_DUP, "Effect size Dup"), pole_key())), nrow = 1,
                                           rel_widths = c(if (with_col1) 1.2, .85, .85, 1.6))

build <- function(stem, blocks, with_col1, h) {
  .chk$log <- character(); SD <- list(); CB <- list(); txt <- c("Del", "Dup", "Effect size Del", "Effect size Dup", "Sensorimotor", "Association")
  B <- lapply(seq_along(blocks), function(b) { bl <- blocks[[b]]
    R <- lapply(seq_along(bl$exp), function(i) row_plot(bl$exp[i], bl$lab[i], last = i == length(bl$exp))); SD[[b]] <<- rbindlist(lapply(R, `[[`, "sd")); txt <<- c(txt, bl$lab)
    rows <- plot_grid(plotlist = lapply(R, `[[`, "g"), ncol = 1)
    if (with_col1) { cb <- count_bar(bl$exp, bl$lab); CB[[b]] <<- cb$sd[, block := bl$name]; txt <<- c(txt, "Genes", "Gene set")
      rows <- plot_grid(cb$g, rows, nrow = 1, rel_widths = c(.5, 3.63)) } else rows <- plot_grid(NULL, rows, nrow = 1, rel_widths = c(.04, 1))
    list(g = rows, n = length(bl$exp)) })
  GAP <- 0.45; nb <- length(B); gl <- list(); hh <- c(); lb <- c()
  for (b in seq_len(nb)) { gl <- c(gl, list(B[[b]]$g)); hh <- c(hh, B[[b]]$n + 0.25 + if (with_col1) 1.2 else 0); lb <- c(lb, letters[b])
    if (b < nb) { gl <- c(gl, list(NULL)); hh <- c(hh, GAP); lb <- c(lb, "") } }
  body <- plot_grid(plotlist = gl, ncol = 1, rel_heights = hh, labels = lb, label_size = FONT_PANEL_TAG_PT, label_fontface = "bold")
  fig <- plot_grid(header(with_col1), body, legend_row(with_col1), ncol = 1, rel_heights = c(.25, sum(hh), 1.1))
  S <- rbindlist(SD); st <- unique(S[, .(row, gene_set_exp, KS_star_DEL, KS_star_DUP)])
  chk(TRUE, sprintf("stars (ST11 FDR q < 0.05): %s", paste(st[, sprintf("%s %s/%s", row, ifelse(KS_star_DEL, "*", "-"), ifelse(KS_star_DUP, "*", "-"))], collapse = "; ")))
  audit_text(txt)
  save_v2(fig, stem, FIG_WIDTH_IN, h); write_sd(S, stem); if (length(CB)) write_sd(rbindlist(CB), stem, "counts")
  write_checks(stem); invisible(st) }

# ---- Supplementary Fig. 16 (old S6)
st6 <- build("SuppFig_Demographics_AssessmentType", list(
  list(name = "demographics", exp = c("Ancestry_EUR", "Sex_Male", "Sex_Female", "Age_lessthan_60", "Age_greaterthan_60"),
       lab = c("EUR", "Male", "Female", "Age ≤ 60", "Age > 60")),
  list(name = "assessment", exp = c("UKBB_FI_center_visit", "UKBB_FI_online"), lab = c("FI centre", "FI online"))), FALSE, 7.0)
# ---- Supplementary Fig. 14 (old S8)
st8 <- build("SuppFig_TissueTemporalCellType", list(
  list(name = "tissue", exp = c("hpa_brain_elevated", "hpa_other_tissue_elevated", "hpa_low_tissue_specificity"), lab = c("Brain", "Other tissue", "Low spec.")),
  list(name = "epoch", exp = paste0("BrainVar_PeakEpoch", 1:3), lab = paste("Epoch", 1:3)),
  list(name = "cell", exp = c("Excitatory", "Inhibitory", "Astrocytes", "Microglial", "Oligodendrocyte", "OligodendrocytePrecursorCells_OPC"),
       lab = c("Excitatory", "Inhibitory", "Astrocytes", "Microglia", "Oligodend.", "OPC"))), TRUE, 10.4)
# Supplementary Fig. 19 (old S9) is drawn by SupFig19_trajectory_pHaplo_pTriplo_missense.R

