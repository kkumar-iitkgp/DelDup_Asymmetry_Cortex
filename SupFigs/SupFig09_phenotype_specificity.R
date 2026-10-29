#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 9 (phenotypic specificity and generalization: profile summary)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Rows: cognitive ability, educational attainment, household income, birthplace north/east, calcium, LDL and total cholesterol (negative controls), All of Us income and educational attainment
#          -> Columns a-d as in Supplementary Fig. 7
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
#       1) SupFigs/output/SupFig09.png, SupFig09.svg (+ checks/SupFig09_checks.txt)
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

FIG_ID <- "SupFig09"
SHEET_NAMES <- c(main = "SupFig9a-d")

STEM <- "SuppFig11_PhenotypeSpecificity_Summary_v2"
rows <- tribble(~block, ~gene_set_exp, ~label, ~ref_pair,
  "main", "All_data_analysis_Fig2", "Cognitive ability (Fig. 2)", NA,
  "ukbb", "EA_1|UKBB", "Educational attainment", "IQ", "ukbb", "Income|UKBB", "Household income", "IQ",
  "neg", "BirthPlace_North|UKBB", "Birthplace north", "IQ", "neg", "BirthPlace_East|UKBB", "Birthplace east", "IQ", "neg", "Calcium|UKBB", "Calcium", "IQ",
  "neg", "LDL|UKBB", "LDL cholesterol", "IQ", "neg", "Cholesterol|UKBB", "Total cholesterol", "IQ",
  "aou", "Income|AoU_ArrayWGS_50kb", "All of Us household income", "Income|UKBB", "aou", "EA|AoU_ArrayWGS_50kb", "All of Us educational attainment", "EA_1|UKBB")
chk(nrow(rows) == 10 && !any(grepl("GATKSV|EA_advanced", rows$gene_set_exp)), "rows: the 10 phenotypes of the figure; no GATK-SV or June EA_advanced run")
v <- rbindlist(lapply(rows$gene_set_exp, st7_row))
chk(nrow(v) == 10 && identical(v$gene_set_exp, rows$gene_set_exp), "10 rows joined to ST7 by gene_set_exp")
d <- bind_cols(rows, v[, .(n_regions, SDbySE_DEL, SDbySE_DUP, I2_DEL, I2_DUP, rho_DEL_H, pspin_DEL_H, rho_DUP_H, pspin_DUP_H, rho_cDEC_H, pspin_cDEC_H,
                          rho_DEL_DUP_all, pspin_DEL_DUP_all, rho_DEL_DUP_sig, pspin_DEL_DUP_sig, I2_pass_Del, I2_pass_Dup, Pair_interpretable)]) %>% mutate(open = block == "neg")
conc <- function(e, ref, pair) { if (is.na(ref)) return(c(NA_real_, NA_real_))
  r <- if (ref == "IQ") st6_row(pair, "All_data_analysis_Fig2", e) else st6_row(pair, ref, e); c(r$cor_input1_vs_input2, r$pspin_input1_vs_input2) }
cc <- t(mapply(function(e, ref) c(conc(e, ref, "del-del"), conc(e, ref, "dup-dup")), d$gene_set_exp, d$ref_pair))
d <- d %>% mutate(conc_rho_DEL = cc[, 1], conc_pspin_DEL = cc[, 2], conc_rho_DUP = cc[, 3], conc_pspin_DUP = cc[, 4])
EXP <- c(`EA_1|UKBB` = .88, `Income|UKBB` = .83, `BirthPlace_North|UKBB` = -.74, `Income|AoU_ArrayWGS_50kb` = .54, `EA|AoU_ArrayWGS_50kb` = .64)
EXPU <- c(`EA_1|UKBB` = .25, `Income|UKBB` = .38, `BirthPlace_North|UKBB` = -.15, `Income|AoU_ArrayWGS_50kb` = .44, `EA|AoU_ArrayWGS_50kb` = .26)
ix <- match(names(EXP), d$gene_set_exp)
chk(all(abs(round(d$conc_rho_DEL[ix], 2) - EXP) < .006) && all(abs(round(d$conc_rho_DUP[ix], 2) - EXPU) < .006),
    sprintf("e: concordance as expected (Del/Dup): %s", paste(sprintf("%s %.2f/%.2f", d$label[ix], d$conc_rho_DEL[ix], d$conc_rho_DUP[ix]), collapse = "; ")))
chk(sum(!is.na(d$conc_rho_DEL)) == 9, "e: 9 concordance rows from ST6 (7 vs cognitive ability, 2 All of Us vs UK Biobank); cognitive-ability row empty")
extra <- function(rows, th, hl) {   # panel e; asterisk only if p-spin < 0.05, not a negative control and the profile passes the I2 gate
  st <- function(p, ok) ifelse(ok & !is.na(p) & p < .05, "*", "")
  le <- rows %>% transmute(label, open, Del = conc_rho_DEL, Dup = conc_rho_DUP, sD = st(conc_pspin_DEL, !open & fadeD == FALSE), sU = st(conc_pspin_DUP, !open & fadeU == FALSE),
                           aD = ifelse(fadeD, .35, 1), aU = ifelse(fadeU, .35, 1))
  hl(ggplot(le, aes(y = label))) + geom_vline(xintercept = 0, linewidth = .35) +
    geom_segment(aes(x = Del, xend = Dup, yend = label), colour = "grey75", linewidth = .6, na.rm = TRUE) +
    geom_point(aes(x = Del, shape = open, alpha = I(aD)), colour = COL_DEL, size = 1.9, stroke = .7, na.rm = TRUE) +
    geom_point(aes(x = Dup, shape = open, alpha = I(aU)), colour = COL_DUP, size = 1.9, stroke = .7, na.rm = TRUE) +
    geom_text(aes(x = Del, label = sD), colour = COL_DEL, vjust = .1, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3), na.rm = TRUE) +
    geom_text(aes(x = Dup, label = sU), colour = COL_DUP, vjust = .1, size = pt_to_geom(FONT_ANNOTATION_PT * 1.3), na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") + scale_x_continuous(breaks = c(-.5, 0, .5), labels = c("-.5", "0", ".5"), expand = expansion(mult = .1)) +
    labs(x = "ρ with reference") + th }
ps <- profile_summary(d, xlab_a = "Dispersion ratio\n(SD/SE)", gate = TRUE)   # v2: no panel e, no panel-a title
ps$plot[[1]] <- ps$plot[[1]] + theme(legend.box = "horizontal", legend.direction = "horizontal", legend.spacing.x = unit(10, "pt"))
ps$plot[[4]] <- ps$plot[[4]] + theme(legend.direction = "vertical", legend.margin = margin(-4, 0, 0, 0), legend.spacing.y = unit(0, "pt"), legend.key.height = unit(9, "pt"))
check_gate(ps)
s <- ps$source
txt <- c(as.character(s$label), "Dispersion ratio\n(SD/SE)", "ρ with hierarchy (H)", "ρ(cDEC, H)", "ρ(Del, Dup)", "All regions", "Sig. regions", "Del", "Dup")
audit_text(txt); chk(!any(grepl("adverse|reliab", txt, ignore.case = TRUE)), "no 'higher = more adverse', no 'reliability' wording")
save_v2(ps$plot, STEM, FIG_WIDTH_IN, 4.0)
sd <- s %>% transmute(label, gene_set_exp, block, negative_control = open, SDbySE_DEL, SDbySE_DUP, I2_DEL, I2_DUP, I2_pass_Del, I2_pass_Dup, Pair_interpretable,
                      rho_DEL_H, pspin_DEL_H, rho_DUP_H, pspin_DUP_H, rho_cDEC_H, pspin_cDEC_H, rho_DEL_DUP_all, rho_DEL_DUP_sig)   # v2: plotted columns only (no panel e)
write_sd(sd, STEM); write_checks(STEM); message("Saved ", FIG_ID)
