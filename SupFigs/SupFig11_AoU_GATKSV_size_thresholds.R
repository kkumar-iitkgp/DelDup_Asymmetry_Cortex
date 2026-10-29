#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 11 (CNV calling and minimum CNV size in All of Us)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Household income: UK Biobank vs GATK-SV 50 kb, then GATK-SV 50 kb vs 40, 30, 20, 10 kb; Del and Dup
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
#       1) SupFigs/output/SupFig11.png, SupFig11.svg (+ checks/SupFig11_checks.txt)
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

FIG_ID <- "SupFig11"
SHEET_NAMES <- c(main = "SupFig11a-j")

lab_of <- c(UKBB = "UKBB", AoU_ArrayWGS_50kb = "AoU", AoU_GATKSV_50kb = "AoU GATK-SV 50 kb", AoU_GATKSV_40kb = "AoU GATK-SV 40 kb",
            AoU_GATKSV_30kb = "AoU GATK-SV 30 kb", AoU_GATKSV_20kb = "AoU GATK-SV 20 kb", AoU_GATKSV_10kb = "AoU GATK-SV 10 kb")
one <- function(i1, i2, t, title, show_title = TRUE) {
  a <- AOU_PAIRS[input1 == i1 & input2 == i2 & cnv_type == t]; stopifnot(nrow(a) == 1)
  s1 <- strsplit(i1, "|", fixed = TRUE)[[1]]; s2 <- strsplit(i2, "|", fixed = TRUE)[[1]]
  x <- profile(s1[2], "ST8", phen = s1[1])[[if (t == "DEL") "D" else "U"]]; y <- profile(s2[2], "ST8", phen = s2[1])[[if (t == "DEL") "D" else "U"]]
  chk(abs(cor(x, y, method = "spearman") - a$rho) < 1e-10 && abs(100 * mean(sign(x) == sign(y)) - a$sign_concordance_pct) < 1e-10,
      sprintf("%s vs %s %s: rho %.3f and sign concordance %.0f%% from plotted data == AOU_PAIRS (p-spin %s, legend)", i1, i2, t, a$rho, a$sign_concordance_pct, format_pspin(a$pspin)))
  cn <- if (t == "DEL") "Del" else "Dup"; lines <- c(fmt_rho(a$rho), fmt_p(a$pspin), sprintf("sign conc. = %.0f%%", a$sign_concordance_pct))
  g <- pvp_scatter(x, y, if (t == "DEL") COL_DEL else COL_DUP, paste0("Effect size ", cn, "\n", lab_of[[s1[2]]]), paste0("Effect size ", cn, "\n", lab_of[[s2[2]]]), lines, identity = TRUE) +
    scale_x_continuous(n.breaks = 4) + labs(title = if (show_title) title else NULL) + theme(plot.title = element_text(size = FONT_AXIS_TITLE_PT, hjust = .5))
  no_smooth(g, paste(i1, i2, t))
  list(g = g, txt = c(lines, if (show_title) title, g$labels$x, g$labels$y),
       sd = data.table(panel = title, cnv_type = cn, x_profile = i1, y_profile = i2, region = ROI, x_effect_size = x, y_effect_size = y, rho = a$rho, p_spin = a$pspin, sign_concordance_pct = a$sign_concordance_pct)) }
build <- function(stem, specs, ncol, h, show_title = TRUE) {
  .chk$log <- character()
  L <- lapply(seq_len(nrow(specs)), function(i) one(specs$i1[i], specs$i2[i], specs$t[i], specs$title[i], show_title))
  fig <- wrap_plots(lapply(L, `[[`, "g"), ncol = ncol) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = FONT_PANEL_TAG_PT, face = "bold"))
  audit_text(unlist(lapply(L, `[[`, "txt")))
  save_v2(fig, stem, FIG_WIDTH_IN, h); write_sd(rbindlist(lapply(L, `[[`, "sd")), stem); write_checks(stem) }

# The cross-cohort figure of the original script is Supplementary Fig. 10 (SupFig10_AoU_cross_cohort.R).
thr <- c(40, 30, 20, 10)
sp <- rbind(tibble(i1 = "Income|UKBB", i2 = "Income|AoU_GATKSV_50kb", lab = "UKBB vs 50 kb"),
            tibble(i1 = "Income|AoU_GATKSV_50kb", i2 = paste0("Income|AoU_GATKSV_", thr, "kb"), lab = paste0("50 kb vs ", thr, " kb")))
sp <- rbind(mutate(sp, t = "DEL", title = paste0(lab, ", Del")), mutate(sp, t = "DUP", title = paste0(lab, ", Dup")))
sp <- sp[c(rbind(1:5, 6:10)), ]   # Del | Dup columns
build("SuppFig_AoU_GATKSV_SizeThresholds", sp, 2, 10.2, show_title = FALSE)
