#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 17 (robustness to genes, loci, CNV size and frequency)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Left: profile concordance (Spearman rho) with the primary Genome Burden model G1 (ST6)
#          -> Right: rho(cDEC,H) per specification (ST7); dashed line at the G1 value
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
#       1) SupFigs/output/SupFig17.png, SupFig17.svg (+ checks/SupFig17_checks.txt)
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

FIG_ID <- "SupFig17"
SHEET_NAMES <- c(main = "SupFig17")

STEM <- "SuppFig_CNVRobustness_Forest"; PRIMARY <- "GB_G1_primary"
rows <- tribble(~fam, ~exp, ~label,
  "Genome burden", "GB_G0_noGenomeBurden", "G0 no genome burden", "Genome burden", "GB_G1_primary", "G1 primary", "Genome burden", "GB_G3_totalCNVgenes", "G3 total CNV genes",
  "Gene-count caps", "genecount_exclude_ge5", "exclude ≥5 genes", "Gene-count caps", "genecount_exclude_ge6", "exclude ≥6 genes",
  "Gene-count caps", "genecount_exclude_ge8", "exclude ≥8 genes", "Gene-count caps", "genecount_exclude_ge10", "exclude ≥10 genes",
  "Gene-count caps", "genecount_keep_ge5", "≥5 genes only", "Gene-count caps", "genecount_keep_ge10", "≥10 genes only",
  "GD loci", "GD_exclude_16major", "remove 16 major GD loci", "GD loci", "GD_exclude_41all", "remove all 41 GD regions",
  "Frequency caps", "freq_le_1in1000", "frequency ≤ 1/1,000", "Frequency caps", "freq_le_1in5000", "frequency ≤ 1/5,000",
  "Frequency caps", "freq_le_1in10000", "frequency ≤ 1/10,000", "Frequency caps", "freq_le_1in25000", "frequency ≤ 1/25,000",
  "Frequency caps", "freq_gt_1in10000", "frequency > 1/10,000",
  "Top-k genes", "topk_drop_top1", "drop top 1 gene", "Top-k genes", "topk_drop_top5", "drop top 5 genes",
  "Top-k genes", "topk_drop_top10", "drop top 10 genes", "Top-k genes", "topk_drop_top25", "drop top 25 genes")
chk(!any(grepl("G2|G4|G5|G6|G7|LOCO|leaveout", rows$exp)), "genome-burden rows G0, G1, G3 only; no LOCO / leave-one-locus rows")
v <- rbindlist(lapply(rows$exp, st7_row)); REF <- st7_row(PRIMARY)$rho_cDEC_H
conc <- function(e, pair) if (e == PRIMARY) c(1, NA) else { r <- st6_row(pair, PRIMARY, e); c(r$cor_input1_vs_input2, r$pspin_input1_vs_input2) }
d <- rows %>% mutate(rho_DEL = sapply(exp, function(e) conc(e, "del-del")[1]), rho_DUP = sapply(exp, function(e) conc(e, "dup-dup")[1]),
                     disp_DEL = v$SDbySE_DEL, disp_DUP = v$SDbySE_DUP, I2_DEL = v$I2_DEL, I2_DUP = v$I2_DUP, pass_DEL = v$I2_pass_Del, pass_DUP = v$I2_pass_Dup,
                     pair = v$Pair_interpretable, rho_cDEC_H = v$rho_cDEC_H, pspin_cDEC_H = v$pspin_cDEC_H) %>%
  mutate(label = ifelse(pair, label, paste0(label, " †")))
chk(all(d$pass_DEL == (d$I2_DEL >= 25)) && all(d$pass_DUP == (d$I2_DUP >= 25)) && all(d$pair == (d$pass_DEL & d$pass_DUP)), "gate columns (ST7) == I2 >= 25% (R1/R2)")
chk(nrow(d) == 20 && !anyNA(d$rho_DEL) && !anyNA(d$rho_cDEC_H), sprintf("20 rows (3 genome burden, 6 gene-count, 2 GD, 5 frequency, 4 top-k): concordance from ST6 (vs %s; G1 itself = 1 by definition), rho(cDEC,H) / p-spin / SD-SE from ST7", PRIMARY))
lv <- rev(d$label); fams <- unique(d$fam)
lc <- d %>% select(fam, label, rho_DEL, rho_DUP, pass_DEL, pass_DUP) %>%
  pivot_longer(-c(fam, label), names_to = c(".value", "cnv"), names_pattern = "(rho|pass)_(DEL|DUP)") %>%
  mutate(cnv = c(DEL = "Del", DUP = "Dup")[cnv], a = ifelse(pass, 1, .35), y = as.numeric(factor(label, levels = lv)) + ifelse(cnv == "Del", .15, -.15),
         panel = "Concordance with G1 (ρ)", fam = factor(fam, levels = fams))
rc <- d %>% transmute(fam = factor(fam, levels = fams), label, rho = rho_cDEC_H, pair,
                      star = ifelse(pair & pspin_cDEC_H < .05, "*", ""), y = as.numeric(factor(label, levels = lv)), panel = "ρ(cDEC, H)")
PAN <- c("Concordance with G1 (ρ)", "ρ(cDEC, H)"); lc$panel <- factor(lc$panel, PAN); rc$panel <- factor(rc$panel, PAN)
bd <- d %>% filter(!pair) %>% transmute(fam = factor(fam, levels = fams), y = as.numeric(factor(label, levels = lv))) %>%
  tidyr::crossing(panel = factor(PAN, PAN))
g <- ggplot() +
  geom_rect(data = bd, aes(xmin = -Inf, xmax = Inf, ymin = y - .5, ymax = y + .5), fill = COL_PAIR_BAND, colour = NA) +
  geom_blank(data = rc %>% tidyr::crossing(dy = c(-.5, .5)) %>% mutate(y = y + dy), aes(x = rho, y = y)) +   # same y range (rows +/- 0.5) in every facet
  geom_vline(data = tibble(panel = factor(PAN, PAN), x = c(1, 0)), aes(xintercept = x), colour = "grey60", linewidth = .3) +
  geom_vline(data = tibble(panel = factor(PAN[2], PAN), x = REF), aes(xintercept = x), linetype = "dashed", colour = "grey30", linewidth = .35) +
  geom_segment(data = rc, aes(x = 0, xend = rho, y = y, yend = y), colour = "grey75", linewidth = .35) +
  geom_point(data = lc, aes(rho, y, colour = cnv, alpha = I(a)), size = 1.8, stroke = .7) +
  geom_point(data = rc, aes(rho, y), size = 1.9, stroke = .7, colour = "black") +
  geom_text(data = rc, aes(rho, y, label = star), nudge_x = .035, vjust = .75, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3)) +
  scale_y_continuous(breaks = seq_along(lv), labels = lv, expand = expansion(add = .1)) +
  scale_colour_manual(values = c(Del = COL_DEL, Dup = COL_DUP), name = NULL) +
  facet_grid(fam ~ panel, scales = "free", space = "free_y", switch = "x") +
  labs(x = NULL, y = NULL) + theme_ng(grid_x = TRUE) +
  theme(strip.placement = "outside", strip.text.y = element_blank(), strip.text.x = element_text(size = FONT_AXIS_TITLE_PT),
        panel.spacing.y = unit(2.5, "mm"), panel.spacing.x = unit(6, "mm"), legend.position = "bottom")
no_smooth(g, "forest")
audit_text(c(d$label, PAN, "Del", "Dup"))
chk(sum(rc$star == "*") == sum(d$pair & d$pspin_cDEC_H < .05), sprintf("rho(cDEC,H): %d asterisks == rule (p-spin < 0.05 and Pair_interpretable)", sum(rc$star == "*")))
g2 <- d %>% filter(!pair); chk(TRUE, sprintf("R2 band + dagger (%d rows): %s", nrow(g2), paste(sprintf("%s [I2 %.1f/%.1f]", g2$label, g2$I2_DEL, g2$I2_DUP), collapse = "; ")))
chk(TRUE, "y axis: each facet shows only its own rows (range rows +/- 0.5, no stray breaks)"); chk(!any(grepl("1\\.2|reliab", c(d$label, PAN))), "no SD/SE 1.2 threshold or reliability wording")
save_v2(g, STEM, FIG_WIDTH_IN, 6.2)
write_sd(d %>% transmute(family = fam, specification = label, gene_set_exp = exp, concordance_rho_DEL_vs_G1 = rho_DEL, concordance_rho_DUP_vs_G1 = rho_DUP,
                         SDbySE_DEL = disp_DEL, SDbySE_DUP = disp_DUP, I2_DEL, I2_DUP, I2_pass_Del = pass_DEL, I2_pass_Dup = pass_DUP, Pair_interpretable = pair,
                         rho_cDEC_H, pspin_cDEC_H, primary_rho_cDEC_H = REF), STEM)
write_checks(STEM)
