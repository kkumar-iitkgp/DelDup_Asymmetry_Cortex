#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 15 (biological strata: profile summary)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Rows: main analysis + 27 strata (LOEUF, pHaplo, pTriplo, missense Z, cell types, peak epoch, trajectory, HPA tissue specificity)
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
#       1) SupFigs/output/SupFig15.png, SupFig15.svg (+ checks/SupFig15_checks.txt)
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

FIG_ID <- "SupFig15"
SHEET_NAMES <- c(main = "SupFig15a-d")

STEM <- "SuppFig_BiologicalStrata_ProfileSummary"
gc <- function(x) paste0("genetic_constraint_", x)
rows <- tribble(~block, ~gene_set_exp, ~label,
  "main", "All_data_analysis_Fig2", "Main analysis (Fig. 2)",
  "LOEUF", gc("loeuf_0p0_0p35"), "LOEUF < 0.35", "LOEUF", gc("loeuf_0p35_1p0"), "LOEUF 0.35-1", "LOEUF", gc("loeuf_1p0_2p0"), "LOEUF > 1",
  "pHaplo", gc("pHaplo_cat_low"), "pHaplo low", "pHaplo", gc("pHaplo_cat_mid"), "pHaplo mid", "pHaplo", gc("pHaplo_cat_high"), "pHaplo high",
  "pTriplo", gc("pTriplo_cat_low"), "pTriplo low", "pTriplo", gc("pTriplo_cat_mid"), "pTriplo mid", "pTriplo", gc("pTriplo_cat_high"), "pTriplo high",
  "mis", gc("missense_Z_cat_low"), "Missense-Z low", "mis", gc("missense_Z_cat_mid"), "Missense-Z mid", "mis", gc("missense_Z_cat_high"), "Missense-Z high",
  "cell", "Excitatory", "Excitatory neurons", "cell", "Inhibitory", "Inhibitory neurons", "cell", "Astrocytes", "Astrocytes",
  "cell", "Microglial", "Microglia", "cell", "Oligodendrocyte", "Oligodendrocytes", "cell", "OligodendrocytePrecursorCells_OPC", "OPC",
  "epoch", "BrainVar_PeakEpoch1", "Peak epoch 1", "epoch", "BrainVar_PeakEpoch2", "Peak epoch 2", "epoch", "BrainVar_PeakEpoch3", "Peak epoch 3",
  "traj", "BrainVar_trajectory_Falling", "Falling", "traj", "BrainVar_trajectory_NonTransitional", "Non-transitional", "traj", "BrainVar_trajectory_Rising", "Rising",
  "hpa", "hpa_brain_elevated", "Brain elevated", "hpa", "hpa_other_tissue_elevated", "Other tissue elevated", "hpa", "hpa_low_tissue_specificity", "Low tissue specificity")
v <- rbindlist(lapply(rows$gene_set_exp, st7_row))
chk(nrow(v) == 28 && identical(v$gene_set_exp, rows$gene_set_exp), "28 rows (main + 27 strata: LOEUF 3, pHaplo 3, pTriplo 3, missense-Z 3, cell 6, epoch 3, trajectory 3, HPA 3) joined to ST7")
d <- bind_cols(rows, v[, .(n_regions, SDbySE_DEL, SDbySE_DUP, I2_DEL, I2_DUP, rho_DEL_H, pspin_DEL_H, rho_DUP_H, pspin_DUP_H,
                          rho_cDEC_H, pspin_cDEC_H, rho_DEL_DUP_all, pspin_DEL_DUP_all, rho_DEL_DUP_sig, pspin_DEL_DUP_sig, I2_pass_Del, I2_pass_Dup, Pair_interpretable)])
chk(all(!is.na(d$pspin_cDEC_H)), "cDEC column present for every row")
ps <- profile_summary(d, xlab_a = "Dispersion ratio\n(SD/SE)", gate = TRUE, title_a = "Regional\nheterogeneity")
check_gate(ps)
audit_text(c(as.character(ps$source$label), "Regional\nheterogeneity", "Dispersion ratio\n(SD/SE)", "ρ with hierarchy (H)", "ρ(cDEC, H)", "ρ(Del, Dup)", "All regions", "Sig. regions", "Del", "Dup"))
save_v2(ps$plot, STEM, FIG_WIDTH_IN, 9.2)
write_sd(ps$source %>% select(-block), STEM); write_checks(STEM)
