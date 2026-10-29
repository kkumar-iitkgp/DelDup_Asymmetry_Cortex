#-------------------------------------------------------------------------------
# Summary: shared data access, checks and panel builders for the Supplementary Figures
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Source after functions/00_figure_style.R. Figures READ statistics from the Supplementary Tables
#             (data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx); nothing a table holds is recomputed
#             here (recomputations are checks only).
#          -> Tables (final numbering): ST4 regional effect sizes; ST6 profile-pair correlations; ST7 profile statistics
#             (SD/SE, I2, rho with H, p-spin); ST8 other phenotypes and All of Us; ST11 sensorimotor vs association KS tests
#             (exact, FDR across the 94 tests); ST12 CNV-level robustness; ST16 Genome Burden G0/G1/G3.
#          -> Small inputs not held in the Supplementary Tables (data/inputs/, see data/README.md):
#             pole membership (Glasser-180), All of Us pair statistics (sign concordance, GATK-SV threshold pairs),
#             KS tests for the alternative parcellations, and two patches (ST6 data-subset pairs; ST7 Schaefer-400 p-spin).
#          -> Rules: Spearman; two-sided p-spin (10,000 rotations); KS stars = ST11 FDR column (q < 0.05);
#             significance classes / outlines = Bonferroni 0.05 / (2 x n regions) within the analysis.
#          -> Output: PNG + SVG only (SupFigs/output/); Source Data sheets in data/source_data_sheets/.
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: Revision 1 (R1)
#
#-----
# 
#-----

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(data.table); library(readxl)
  library(ggseg); library(ggsegGlasser); library(cowplot); library(patchwork); library(ggrepel); library(scales)
})
source(file.path(project_dir, "functions", "fn_write_SourceData_sheet.R"))

data_dir <- file.path(project_dir, "data")
inputs_dir <- file.path(data_dir, "inputs")
source_data_dir <- file.path(data_dir, "source_data_sheets")
OUT_DIR <- file.path(project_dir, "SupFigs", "output")
CHECKS_DIR <- file.path(OUT_DIR, "checks")
dir.create(CHECKS_DIR, showWarnings = FALSE, recursive = TRUE)
XL <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
POLES_CSV <- file.path(inputs_dir, "inputs_pole_membership_Glasser180.csv")
for (f in c(XL, POLES_CSV)) if (!file.exists(f)) stop("Missing required input file: ", f, call. = FALSE)

COL_SIG <- c(Del = COL_SIG_DEL, Dup = COL_SIG_DUP, `Del-Dup` = COL_SIG_BOTH, ns = COL_SIG_NS)
COL_POLE <- c(sensorimotor = COL_SENSORIMOTOR, association = COL_ASSOCIATION)

# ------------------------------- CHECKS ---------------------------------------
.chk <- new.env(); .chk$log <- character()
chk <- function(ok, msg) { s <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", msg); message(s); .chk$log <- c(.chk$log, s); if (!isTRUE(ok)) stop(msg, call. = FALSE) }
# on-plot text audit: no x-sign, beta, arrows, FDR wording, line/method notes
audit_text <- function(txt) {
  bad <- grepl("×|β|↓|↑|FDR|q-spin|pending|dashed|dotted|line =|fit|visual aid|Bonferroni|n_rot|rotation", txt, ignore.case = TRUE)
  chk(!any(bad), sprintf("on-plot text audit: %d strings, none with x-sign / beta / arrows / FDR / line or method notes%s", length(txt),
                         if (any(bad)) paste0(" (offending: ", paste(txt[bad], collapse = " | "), ")") else ""))
}
no_smooth <- function(p, what) chk(!any(vapply(p$layers, function(l) inherits(l$stat, "StatSmooth") || inherits(l$geom, "GeomSmooth"), logical(1))),
                                   paste0(what, ": no fit line (no geom_smooth layer)"))

# ------------------------------- OUTPUT ---------------------------------------
# Each script sets FIG_ID (e.g. "SupFig07") and SHEET_NAMES (part -> Source Data sheet name; "main" when no part).
# OUT_MAP (optional): stem -> figure file name, for scripts that draw two figures.
save_v2 <- function(plot_obj, stem, width, height) {
  s <- if (exists("OUT_MAP") && stem %in% names(OUT_MAP)) OUT_MAP[[stem]] else FIG_ID
  save_figure(plot_obj, s, width, height, out_dir = OUT_DIR)
  chk(!file.exists(file.path(OUT_DIR, paste0(s, ".pdf"))), paste0(s, ": PNG + SVG written, no PDF"))
}
write_sd <- function(df, stem, part = "") {
  key <- if (nzchar(part)) part else "main"
  key2 <- paste(stem, key, sep = ":")
  nm <- if (key2 %in% names(SHEET_NAMES)) SHEET_NAMES[[key2]] else SHEET_NAMES[[key]]
  fn_write_SourceData_sheet(as.data.frame(df), nm, source_data_dir)
}
write_checks <- function(stem) {
  s <- if (exists("OUT_MAP") && stem %in% names(OUT_MAP)) OUT_MAP[[stem]] else FIG_ID
  writeLines(.chk$log, file.path(CHECKS_DIR, paste0(s, "_checks.txt")))
}

# ------------------------------- TABLES ---------------------------------------
.tab <- new.env()
sheet <- function(s) { if (is.null(.tab[[s]])) .tab[[s]] <- as.data.table(read_excel(XL, s, guess_max = 50000)); .tab[[s]] }
ST7 <- copy(sheet("ST7"))     # profile statistics (SD/SE, I2, rho with H, p-spin, gate columns)
ST6 <- copy(sheet("ST6"))     # profile-pair correlations (rho, 95% CI, p-spin)
ST11 <- copy(sheet("ST11"))   # sensorimotor vs association: means, 95% CI, KS D, exact p, FDR (94 tests)
for (cc in c("I2_pass_Del", "I2_pass_Dup", "Pair_interpretable")) if (is.character(ST7[[cc]])) ST7[[cc]] <- toupper(ST7[[cc]]) == "TRUE"

# patch 1: ST7 Schaefer-400 p-spin (pending in the 29 Oct workbook; values from the 13 Oct 2026 statistics run)
.p7 <- fread(file.path(inputs_dir, "inputs_patch_ST7_Schaefer400_pspin.csv"))
for (i in seq_len(nrow(.p7))) { ix <- which(ST7$gene_set_exp == .p7$gene_set_exp[i])
  for (cc in setdiff(names(.p7), c("gene_set_exp_category", "gene_set_exp"))) if (is.na(ST7[[cc]][ix])) set(ST7, ix, cc, .p7[[cc]][i]) }
# patch 2: ST6 data-subset pairs (UK Biobank vs five cohorts; Discovery vs Replication), plotted in Supplementary Fig. 2e-h
.p6 <- fread(file.path(inputs_dir, "inputs_patch_ST6_data_subset_pairs.csv"))
ST6 <- rbind(ST6, .p6[!paste(cor_pair, gene_set_exp_input1, gene_set_exp_input2) %in% ST6[, paste(cor_pair, gene_set_exp_input1, gene_set_exp_input2)]], fill = TRUE)
KS_ALT <- fread(file.path(inputs_dir, "inputs_KS_alternative_parcellations.csv"))   # HPA-73 cores, BrainSpan-11 poles, Desikan-34 poles (nominal)
AOU_PAIRS <- fread(file.path(inputs_dir, "inputs_AoU_pair_statistics.csv"))          # All of Us pairs: rho, p-spin, sign concordance

chk(nrow(ST11) == 94 && max(abs(p.adjust(ST11$KolmogorovSmirnov_pvalue_sensori_vs_assoc, "BH") - ST11$KolmogorovSmirnov_pvalue_sensori_vs_assoc_FDR)) < 1e-12,
    "ST11: 94 KS tests; FDR column == BH over the 94")
st7_row <- function(exp) { r <- ST7[gene_set_exp == exp]; if (nrow(r) != 1) stop("ST7 row missing/duplicated: ", exp, call. = FALSE); r }
st6_row <- function(pair, e1, e2) { r <- ST6[cor_pair == pair & gene_set_exp_input1 == e1 & gene_set_exp_input2 == e2]
  if (nrow(r) != 1) stop(sprintf("ST6 row missing: %s %s vs %s", pair, e1, e2), call. = FALSE); r }
ks_star <- function(exp, cnv) { r <- ST11[gene_set_exp == exp & toupper(cnv_type) == toupper(cnv)]
  if (nrow(r) != 1) stop("ST11 row missing for ", exp, " / ", cnv, call. = FALSE); r$KolmogorovSmirnov_pvalue_sensori_vs_assoc_FDR < 0.05 }

# Glasser region order = ST4 main profile (name-based everywhere)
ROI <- sheet("ST4")[gene_set_exp == "All_data_analysis_Fig2" & cnv_type == "DEL", gene_set_name]
GLASSER_REGION <- sub("^L_", "", gsub(".", "-", ROI, fixed = TRUE))
chk(length(ROI) == 180 && all(GLASSER_REGION %in% names(glasser$palette)), "Glasser: 180 ST4 regions, all found in the ggseg atlas by name")

# wide Del/Dup profile for one analysis, joined by region name to `ord` (default: the 180 Glasser regions)
profile <- function(exp, src = "ST4", ord = ROI, phen = NULL) {
  d <- if (src == "ST8") sheet("ST8")[Phenotype == phen & Experiment == exp] else sheet(src)[gene_set_exp == exp]
  a <- d[cnv_type == "DEL"]; a <- a[match(ord, a$gene_set_name)]; b <- d[cnv_type == "DUP"]; b <- b[match(ord, b$gene_set_name)]
  if (anyNA(a$Estimate) || anyNA(b$Estimate)) stop("profile join failed: ", exp, call. = FALSE)
  data.table(region = ord, D = a$Estimate, U = b$Estimate, pD = a$pvalue, pU = b$pvalue)
}
sig_class <- function(pD, pU, denom) {
  d <- pmin(pD * denom, 1) < .05; u <- pmin(pU * denom, 1) < .05
  factor(ifelse(d & u, "Del-Dup", ifelse(d, "Del", ifelse(u, "Dup", "ns"))), levels = names(COL_SIG))
}

# poles (Glasser, by name)
.pm <- fread(POLES_CSV)
POLE <- .pm$pole[match(ROI, .pm$region)]
chk(!anyNA(POLE) && sum(POLE == "sensorimotor") == 54 && sum(POLE == "association") == 54, "poles: pole membership joined by name, 54 sensorimotor / 54 association")

# ------------------------------- MAPS -----------------------------------------
# Glasser LH dispersed (lateral + medial); outlines = Bonferroni-significant regions
effect_map <- function(effect, sig, lim, legend_title = "Effect size", show_legend = FALSE) {
  dat <- tibble(region = GLASSER_REGION, effect = pmax(pmin(effect, lim), -lim),
                edge = ifelse(sig, "black", "grey75"), lw = ifelse(sig, .18, .03))
  ggseg(dat, atlas = "glasser", hemisphere = "left", position = "dispersed",
        mapping = aes(fill = effect, colour = I(edge), size = I(lw)), adapt_scales = TRUE) +
    scale_fill_gradientn(limits = c(-lim, lim), breaks = c(-lim, 0, lim), labels = function(x) sprintf("%g", x),
                         colours = colorRampPalette(c("blue4", "blue", "white", "red", "red4"))(300), oob = squish, name = legend_title) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = .5, barwidth = 5, barheight = .35)) +
    theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
    theme(legend.position = if (show_legend) "bottom" else "none", legend.title = element_text(size = FONT_LEGEND_TITLE_PT),
          legend.text = element_text(size = FONT_LEGEND_TEXT_PT), axis.title = element_blank(), axis.text = element_blank(),
          plot.margin = margin(-6, 0, -6, 0))
}
colourbar_only <- function(lim, title) get_legend(effect_map(rep(0, 180), rep(FALSE, 180), lim, title, show_legend = TRUE))

# ------------------------------- SCATTERS -------------------------------------
# Del (y) vs Dup (x), all regions, significance-class colours, labels only where `lab` is non-NA, no line
dvd_scatter <- function(D, U, sig, lab = NA_character_, size = 1.3) {
  dd <- tibble(D = D, U = U, sig = sig, lab = lab) %>% arrange(sig != "ns")
  ggplot(dd, aes(U, D)) +
    geom_hline(yintercept = 0, colour = COL_LINE, linewidth = .3) + geom_vline(xintercept = 0, colour = COL_LINE, linewidth = .3) +
    geom_point(aes(fill = sig), shape = 21, colour = "white", stroke = .15, size = size) +
    geom_text_repel(data = filter(dd, !is.na(lab)), aes(label = lab), size = pt_to_geom(FONT_MIN_PT), min.segment.length = 0,
                    segment.colour = "grey50", segment.size = .2, max.overlaps = Inf, seed = 1, box.padding = .3) +
    scale_fill_manual(values = COL_SIG, drop = FALSE, name = NULL, labels = c(Del = "Del sig.", Dup = "Dup sig.", `Del-Dup` = "Del & Dup sig.", ns = "n.s.")) +
    labs(x = "Effect size Dup", y = "Effect size Del") + theme_ng() + theme(legend.position = "none")
}
# profile-vs-profile scatter (same CNV type), identity colour, two short stats lines, no fit line
pvp_scatter <- function(x, y, col, xlab, ylab, lines, identity = FALSE) {
  dd <- tibble(x = x, y = y); r <- range(c(x, y)); pad <- diff(r) * .30
  p <- ggplot(dd, aes(x, y))
  if (identity) p <- p + annotate("segment", x = r[1], y = r[1], xend = r[2], yend = r[2], linetype = "dashed", colour = "grey45", linewidth = .35)
  p + geom_hline(yintercept = 0, colour = COL_LINE, linewidth = .3) + geom_vline(xintercept = 0, colour = COL_LINE, linewidth = .3) +
    geom_point(colour = col, alpha = .7, size = .9) +
    annotate("text", x = -Inf, y = Inf, hjust = -.05, vjust = 1.2, label = paste(lines, collapse = "\n"), size = pt_to_geom(FONT_ANNOTATION_PT), lineheight = .95) +
    { hr <- .35 + .20 * (length(lines) - 2); if (identity) coord_cartesian(xlim = r, ylim = c(r[1], r[2] + diff(r) * hr), clip = "off")
      else coord_cartesian(ylim = c(min(y), max(y) + diff(range(y)) * hr), clip = "off") } +
    scale_y_continuous(n.breaks = 4) + labs(x = xlab, y = ylab) + theme_ng() + theme(plot.margin = margin(4, 10, 2, 2))
}
fmt_rho <- function(r) sprintf("ρ = %.2f", r)
fmt_p <- function(p, min_le = FALSE) if (min_le && p <= 1e-4 + 1e-12) "p-spin ≤ 10⁻⁴" else paste0("p-spin = ", format_pspin(p))

# ------------------------------- DENSITIES (2023 layout) ------------------------
# sensorimotor vs association densities; no axes; dotted zero line; star on the left (ST11 FDR q < 0.05)
density_old <- function(effect, xlim, star = FALSE, ticks = NULL, tick_pt = FONT_AXIS_TEXT_PT) {   # ticks: x breaks to show (NULL = no axis)
  dd <- tibble(effect = effect[POLE != "transition"], pole = factor(POLE[POLE != "transition"], levels = names(COL_POLE)))
  p <- ggplot(dd, aes(effect, fill = pole, colour = pole)) +
    geom_density(alpha = .5, linewidth = .35) + geom_vline(xintercept = 0, linetype = "dotted", linewidth = .35) +
    scale_fill_manual(values = COL_POLE) + scale_colour_manual(values = COL_POLE) +
    coord_cartesian(xlim = xlim, expand = FALSE, clip = "on") + theme_classic(base_family = FONT_FAMILY) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = "none", plot.margin = margin(4, 2, 2, 2))
  if (isTRUE(star)) p <- p + annotate("text", x = xlim[1], y = Inf, label = "*", hjust = -.2, vjust = 1.1, fontface = "bold", size = pt_to_geom(FONT_PANEL_TAG_PT * 1.6))
  if (!is.null(ticks)) p <- p + scale_x_continuous(breaks = ticks, labels = function(x) sub("^(-?)0\\.", "\\1.", sprintf("%g", x))) +
    theme(axis.text.x = element_text(size = tick_pt, colour = "black"), axis.ticks.x = element_line(linewidth = .3), axis.ticks.length.x = unit(2, "pt"),
          axis.line.x = element_line(linewidth = .3))
  p
}
pole_key <- function() get_legend(ggplot(tibble(x = 1:2, pole = factor(names(COL_POLE), levels = names(COL_POLE))), aes(x, x, fill = pole)) +
  geom_tile(alpha = .5) + scale_fill_manual(values = COL_POLE, labels = c("Sensorimotor", "Association"), name = NULL) + theme_ng() +
  theme(legend.position = "bottom", legend.text = element_text(size = FONT_LEGEND_TEXT_PT)))
vlabel <- function(txt) ggdraw() + draw_label(txt, angle = 90, size = FONT_INTERNAL_LABEL_PT, fontfamily = FONT_FAMILY)

# ------------------------------- PROFILE SUMMARY ----------------------------------
# rows: data.frame(label, block, + ST7 columns); four columns: a SD/SE (+ I2 size), b rho(Del,H)/rho(Dup,H), c rho(cDEC,H), d Del-Dup rho
# star = two-sided p-spin < 0.05 (nominal); open = negative-control row (no stars)
# gate = TRUE: heterogeneity gate from the ST7 columns I2_pass_Del / I2_pass_Dup / Pair_interpretable.
#   R1: profile with I2 < 25% -> its rho with H drawn faded (alpha .35), no asterisk (b).
#   R2: Del or Dup I2 < 25% -> light grey band (#EFEFEF) behind the row in every panel, dagger after the label, no asterisk on c/d.
#   R3: `open` rows (negative controls) -> open symbols, light cDEC bar, no asterisks.
COL_PAIR_BAND <- "#EFEFEF"
profile_summary <- function(rows, extra = NULL, xlab_a = "SD / mean(SE)", gate = FALSE, title_a = NULL) {
  rows <- rows %>% mutate(open = if ("open" %in% names(.)) open else FALSE)
  if (gate) {
    stopifnot(all(c("I2_pass_Del", "I2_pass_Dup", "Pair_interpretable") %in% names(rows)))
    rows <- rows %>% mutate(fadeD = I2_pass_Del %in% FALSE, fadeU = I2_pass_Dup %in% FALSE, band = Pair_interpretable %in% FALSE,
                            label = ifelse(band, paste0(label, " †"), label))
  } else rows <- rows %>% mutate(fadeD = FALSE, fadeU = FALSE, band = FALSE)
  rows <- rows %>% mutate(label = factor(label, levels = rev(label)))
  brk <- which(diff(as.integer(factor(rows$block, levels = unique(rows$block)))) != 0); yb <- nrow(rows) - brk + .5
  yband <- nrow(rows) - which(rows$band) + 1
  hl <- function(p) { if (length(yband)) p <- p + scale_y_discrete() +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = yband - .5, ymax = yband + .5, fill = COL_PAIR_BAND, colour = NA)
    p + geom_hline(yintercept = yb, colour = "grey72", linewidth = .4) }
  st <- function(p, open) ifelse(!open & !is.na(p) & p < .05, "*", "")
  th <- theme_ng(grid_x = TRUE) + theme(axis.title.y = element_blank(), plot.title = element_text(size = FONT_AXIS_TITLE_PT, hjust = 0))
  nolab <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  la <- rows %>% select(label, open, SDbySE_DEL, SDbySE_DUP, I2_DEL, I2_DUP) %>%
    pivot_longer(-c(label, open), names_to = c(".value", "cnv"), names_pattern = "(SDbySE|I2)_(DEL|DUP)") %>% mutate(cnv = c(DEL = "Del", DUP = "Dup")[cnv])
  pa <- hl(ggplot(la, aes(SDbySE, label, colour = cnv, size = I2))) + labs(title = title_a) + geom_point(aes(shape = open), alpha = .85, stroke = .7, na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    scale_colour_manual(values = c(Del = COL_DEL, Dup = COL_DUP), name = NULL) + scale_size(range = c(.8, 3.4), name = expression(I^2 ~ "(%)"), breaks = c(0, 50, 100), limits = c(0, 100)) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40", linewidth = .35) +
    labs(x = xlab_a) + th + theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), axis.ticks.y = element_blank())
  lb <- rows %>% transmute(label, open, Del = rho_DEL_H, Dup = rho_DUP_H, sD = st(pspin_DEL_H, open | fadeD), sU = st(pspin_DUP_H, open | fadeU),
                           aD = ifelse(fadeD, .35, 1), aU = ifelse(fadeU, .35, 1))
  pb <- hl(ggplot(lb, aes(y = label))) + geom_vline(xintercept = 0, linewidth = .35) +
    geom_segment(aes(x = Del, xend = Dup, yend = label), colour = "grey75", linewidth = .6, na.rm = TRUE) +
    geom_point(aes(x = Del, shape = open, alpha = I(aD)), colour = COL_DEL, size = 1.9, stroke = .7, na.rm = TRUE) + geom_point(aes(x = Dup, shape = open, alpha = I(aU)), colour = COL_DUP, size = 1.9, stroke = .7, na.rm = TRUE) +
    geom_text(aes(x = Del, label = sD), colour = COL_DEL, vjust = .1, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3), na.rm = TRUE) +
    geom_text(aes(x = Dup, label = sU), colour = COL_DUP, vjust = .1, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3), na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") + scale_x_continuous(n.breaks = 4, expand = expansion(mult = .08)) + labs(x = "ρ with hierarchy (H)") + th + nolab
  lc <- rows %>% transmute(label, open, r = rho_cDEC_H, s = st(pspin_cDEC_H, open | band))
  pc <- hl(ggplot(lc, aes(r, label))) + geom_vline(xintercept = 0, linewidth = .35) +
    geom_col(aes(alpha = open), fill = COL_DELTA_BAR, width = .62, orientation = "y", na.rm = TRUE) + scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = .35), guide = "none") +
    geom_text(aes(x = r + sign(r) * .02, label = s, hjust = ifelse(r >= 0, 0, 1)), vjust = .75, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3), na.rm = TRUE) +
    scale_x_continuous(n.breaks = 3, expand = expansion(mult = .15)) + labs(x = "ρ(cDEC, H)") + th + nolab
  ld <- rows %>% transmute(label, open, all = rho_DEL_DUP_all, sig = rho_DEL_DUP_sig)
  pd <- hl(ggplot(ld, aes(y = label))) + geom_vline(xintercept = 0, linewidth = .35) +
    geom_segment(aes(x = all, xend = sig, yend = label), colour = "grey75", linewidth = .6, na.rm = TRUE) +
    geom_point(aes(x = all, shape = "All regions"), colour = "grey50", size = 1.7, na.rm = TRUE) +
    geom_point(aes(x = sig, shape = "Sig. regions"), colour = "black", size = 1.7, na.rm = TRUE) +
    scale_shape_manual(values = c(`All regions` = 16, `Sig. regions` = 18), name = NULL) +
    labs(x = "ρ(Del, Dup)") + th + nolab + theme(legend.position = "bottom", legend.direction = "vertical")
  pl <- list(pa, pb, pc, pd); if (!is.null(extra)) pl <- c(pl, list(extra(rows, th + nolab, hl)))
  list(plot = wrap_plots(pl, nrow = 1, widths = c(1.15, 1, .8, .9, if (!is.null(extra)) .8)) + plot_annotation(tag_levels = "a") &
         theme(plot.tag = element_text(size = FONT_PANEL_TAG_PT, face = "bold")),
       source = rows %>% mutate(label = as.character(label)))
}

# ----------------------------------------------------------------------
# For profile_summary(gate = TRUE) output: asterisk counts in b and c equal the table rule; gate columns equal I2 >= 25%.
check_gate <- function(ps) {
  s <- ps$source
  chk(all(s$I2_pass_Del == (s$I2_DEL >= 25)) && all(s$I2_pass_Dup == (s$I2_DUP >= 25)) && all(s$Pair_interpretable == (s$I2_pass_Del & s$I2_pass_Dup)),
      "gate columns (ST7) == I2 >= 25% recomputed from the ST7 I2 values (R1/R2)")
  stars <- function(k) { d <- ggplot_build(ps$plot[[k]])$data; x <- unlist(lapply(d, function(l) if ("label" %in% names(l) && is.character(l$label)) l$label)); sum(x == "*") }
  e <- function(p, ok) sum(!s$open & ok & !is.na(p) & p < .05)
  chk(stars(2) == e(s$pspin_DEL_H, s$I2_pass_Del) + e(s$pspin_DUP_H, s$I2_pass_Dup), sprintf("b: %d asterisks == rule (p-spin < 0.05 and I2 >= 25%% for that profile)", stars(2)))
  chk(stars(3) == e(s$pspin_cDEC_H, s$Pair_interpretable), sprintf("c: %d asterisks == rule (p-spin < 0.05 and Pair_interpretable)", stars(3)))
  g <- s[!s$Pair_interpretable, ]
  chk(TRUE, sprintf("R2 band + dagger (%d rows): %s", nrow(g), if (nrow(g)) paste(sprintf("%s [I2 %.1f/%.1f]", g$label, g$I2_DEL, g$I2_DUP), collapse = "; ") else "none"))
  r1 <- s[(!s$I2_pass_Del | !s$I2_pass_Dup), ]
  chk(TRUE, sprintf("R1 faded rho with H: %s", if (nrow(r1)) paste(sprintf("%s (%s)", r1$label, ifelse(!r1$I2_pass_Del & !r1$I2_pass_Dup, "Del, Dup", ifelse(!r1$I2_pass_Del, "Del", "Dup"))), collapse = "; ") else "none"))
}
