#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 7 (sensitivity analyses: cortical profile summary)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: dispersion ratio SD/SE (point size = I2)
#          -> b: rho(Del,H), rho(Dup,H)
#          -> c: rho(cDEC,H)
#          -> d: rho(Del,Dup), all and significant regions; heterogeneity gate (I2 < 25%)
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
#       1) SupFigs/output/SupFig07.png, SupFig07.svg (+ checks/SupFig07_checks.txt)
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

FIG_ID <- "SupFig07"
SHEET_NAMES <- c(main = "SupFig7a-d")

STEM <- "SuppFig_Sensitivity_ProfileSummary"

rows <- tribble(~block, ~gene_set_exp, ~label,
  "main", "All_data_analysis_Fig2", "Main analysis (Fig. 2)",
  "Z", "Z0p2", "Regional specificity Z > 0.2", "Z", "Z0p5", "Regional specificity Z > 0.5",
  "Z", "Z0p8", "Regional specificity Z > 0.8", "Z", "Z1p0", "Regional specificity Z > 1.0",
  "data", "UKBB", "UK Biobank", "data", "Other_5cohort", "Five cohorts", "data", "Discovery", "Discovery", "data", "Replication", "Replication",
  "demo", "UKBB_FI_center_visit", "FI centre visit", "demo", "UKBB_FI_online", "FI online", "demo", "Ancestry_EUR", "Ancestry: EUR",
  "demo", "Sex_Male", "Male", "demo", "Sex_Female", "Female", "demo", "Age_lessthan_60", "Age ≤ 60", "demo", "Age_greaterthan_60", "Age > 60",
  "parc", "AHBA_PausLab_Desikan", "Desikan-34", "parc", "AHBA_schaefer400x7MNI", "Schaefer-400 (LH)",
  "DS", "gene_exp_AHBA_quality_DScentile10", "DS 10th centile", "DS", "gene_exp_AHBA_quality_DScentile20", "DS 20th centile",
  "DS", "gene_exp_AHBA_quality_DScentile30", "DS 30th centile", "DS", "gene_exp_AHBA_quality_DScentile40", "DS 40th centile")
v <- rbindlist(lapply(rows$gene_set_exp, st7_row))
chk(nrow(v) == nrow(rows) && identical(v$gene_set_exp, rows$gene_set_exp), sprintf("%d rows joined to ST7 by gene_set_exp", nrow(rows)))
chk(!any(grepl("Z1p2|DScentile50", rows$gene_set_exp)), "no Z1p2 and no DS centile above 40th")
d <- bind_cols(rows, v[, .(n_regions, SDbySE_DEL, SDbySE_DUP, I2_DEL, I2_DUP, rho_DEL_H, pspin_DEL_H, rho_DUP_H, pspin_DUP_H,
                          rho_cDEC_H, pspin_cDEC_H, rho_DEL_DUP_all, pspin_DEL_DUP_all, rho_DEL_DUP_sig, pspin_DEL_DUP_sig, I2_pass_Del, I2_pass_Dup, Pair_interpretable)])
chk(all(!is.na(d$pspin_cDEC_H)), "cDEC column present: every row has rho(cDEC,H) and its p-spin (Schaefer400 filled in ST7)")
ps <- profile_summary(d, xlab_a = "Dispersion ratio\n(SD/SE)", gate = TRUE, title_a = "Regional\nheterogeneity")
check_gate(ps)
audit_text(c(as.character(ps$source$label), "Regional\nheterogeneity", "Dispersion ratio\n(SD/SE)", "ρ with hierarchy (H)", "ρ(cDEC, H)", "ρ(Del, Dup)", "All regions", "Sig. regions", "Del", "Dup"))
save_v2(ps$plot, STEM, FIG_WIDTH_IN, 7.4)
write_sd(ps$source %>% select(-block), STEM); write_checks(STEM)
