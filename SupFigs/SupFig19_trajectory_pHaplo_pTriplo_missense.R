#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 19 (temporal trajectory, haploinsufficiency, triplosensitivity and missense constraint)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Rows: BrainVar trajectory, pHaplo, pTriplo, missense Z (low / mid / high)
#          -> Layout as Supplementary Fig. 14; missense rows: Dup density window -0.30 to 0.08
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
#       1) SupFigs/output/SupFig19.png, SupFig19.svg (+ checks/SupFig19_checks.txt)
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

FIG_ID <- "SupFig19"
SHEET_NAMES <- c(main = "SupFig19a-d", counts = "SupFig19_counts")

# v3b (12 Oct 2026): Fig. S9 only = v3 with ONE change: Missense-Z rows, Dup density view -0.30..0.08 (coord_cartesian, all 108 pole
#   regions), ticks -0.3, -0.2, -0.1, 0 on the block's last row (as v3). Gaps, ticks, fonts and everything else as v3.
XL_DUP_MIS <- c(-.30, .08); TICKS_DUP_MIS <- c(-.3, -.2, -.1, 0); DENS <- list()
LAB_PT <- FONT_INTERNAL_LABEL_PT - 1; TICKS_DEL <- c(-.3, -.2, -.1, 0); TICKS_DUP <- c(-.1, -.05, 0, .05)
vlabel <- function(txt) ggdraw() + draw_label(txt, angle = 90, size = LAB_PT, fontfamily = FONT_FAMILY)
LIM_DEL <- .30; LIM_DUP <- .15; XL_DEL <- c(-.30, .08); XL_DUP <- c(-.13, .08)
gc <- function(x) paste0("genetic_constraint_", x)

row_plot <- function(exp, label, last = FALSE) { mis <- grepl("missense_Z", exp); p <- profile(exp); sD <- pmin(p$pD * 360, 1) < .05; sU <- pmin(p$pU * 360, 1) < .05
  stD <- ks_star(exp, "DEL"); stU <- ks_star(exp, "DUP")
  list(g = plot_grid(vlabel(label), effect_map(p$D, sD, LIM_DEL), density_old(p$D, XL_DEL, stD, if (last) TICKS_DEL, FONT_AXIS_TEXT_PT - 1), effect_map(p$U, sU, LIM_DUP), { du <- density_old(p$U, if (mis) XL_DUP_MIS else XL_DUP, stU, if (last) (if (mis) TICKS_DUP_MIS else TICKS_DUP), FONT_AXIS_TEXT_PT - 1); DENS[[exp]] <<- list(n = nrow(du$data), xlim = du$coordinates$limits$x); du },
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

# ---- old S9
st9 <- build("SuppFig_TrajectoryHaploTriploMissense", list(
  list(name = "trajectory", exp = paste0("BrainVar_trajectory_", c("Falling", "NonTransitional", "Rising")), lab = c("Falling", "Non-trans.", "Rising")),
  list(name = "pHaplo", exp = gc(paste0("pHaplo_cat_", c("low", "mid", "high"))), lab = paste("pHaplo", c("low", "mid", "high"))),
  list(name = "pTriplo", exp = gc(paste0("pTriplo_cat_", c("low", "mid", "high"))), lab = paste("pTriplo", c("low", "mid", "high"))),
  list(name = "missense", exp = gc(paste0("missense_Z_cat_", c("low", "mid", "high"))), lab = paste("Mis-Z", c("low", "mid", "high")))), TRUE, 10.8)

# ---- v3b checks (S9 only)
.chk$log <- character()
m <- DENS[grepl("missense_Z", names(DENS))]; o <- DENS[!grepl("missense_Z", names(DENS))]
chk(length(m) == 3 && all(sapply(m, `[[`, "n") == 108) && all(sapply(m, function(z) identical(z$xlim, XL_DUP_MIS))),
    "Missense-Z low / mid / high Dup densities: 108 pole regions each (54 + 54), view -0.30..0.08 (coord_cartesian, no data dropped)")
chk(length(o) == 9 && all(sapply(o, function(z) identical(z$xlim, XL_DUP))), "all other Dup densities unchanged from v3: view -0.13..0.08 (9 rows)")
mh <- profile(gc("missense_Z_cat_high"))$U[POLE != "transition"]
chk(min(mh) > XL_DUP_MIS[1], sprintf("Missense-Z high Dup: full left tail inside the view (min %.3f > -0.30)", min(mh)))
writeLines(.chk$log, file.path(CHECKS_DIR, paste0(FIG_ID, "_extra_checks.txt")))
