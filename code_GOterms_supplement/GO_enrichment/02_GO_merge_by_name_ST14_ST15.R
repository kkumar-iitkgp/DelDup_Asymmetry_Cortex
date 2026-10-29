#-------------------------------------------------------------------------------
# Summary: R-script to merge the per-region fgsea results by GO-term and region NAME, and build Supplementary Tables 14 and 15
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> Reads the 180 x 3 per-region fgsea outputs (01_GO_fgsea_per_region.R), joins them by (pathway, region) names,
#             applies the term-size rule and the prevalence filter, and correlates the NES map of each retained term with the
#             deletion and duplication effect-size maps (Spearman; spin test with the GO NES map rotated and the effect-size
#             profile held fixed; cleaned rotation set, 9,999 rotations; two-sided; BH FDR per CNV type across the retained terms).
#             Correct merging: the 2023 summary step combined the per-region outputs by row position (cbind), but each region's
#             rows are sorted by its own NES and hold a different number of terms; only region 1 (L_V1) was name-correct.
#             Here every value is joined by term name and region name.
# Key steps:
#       1) PART A: Parameters and input files
#       2) PART B: Long table (all tested GO terms x 180 regions), global BH FDR
#       3) PART C: Prevalence filter and NES maps
#       4) PART D: Spearman correlation with Del/Dup effect-size maps (spin test), ST15
#       5) PART E: Write the long table, ST14 and ST15
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: 2 (Revision 1)
#
#-----
# Input files:
#       1) results_fgsea_pos_pathway_<go_mf|go_cc|go_bp>_flagSubset_0_roi_<1..180>.RData (2023 fgsea run)
#       2) results_GO_GeneSets_fgsea_pos_flagSubset_0.RData (GO gene sets of that run)
#       3) df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData (region order of the run = its columns)
#       4) SupplementaryTables_NatGenRevision_29oct2026.xlsx (ST4: Fig. 2 effect sizes)
#       5) perm_clean_GlasserLH180.rds (10,000 spherical rotations, Glasser LH 180; identity/duplicate rotations removed: 9,999)
#       6) fn_spin_functions.R
#
#-----
# Output (out_dir):
#       1) GO_regional_enrichment_fgsea_all_terms.csv.gz (all tested terms x 180 regions; pathway_cat, pathway, size, region, ES, NES, p, p_fdr)
#       2) GO_term_catalogue.csv (one row per tested term: size, regions tested, significant regions, retained)
#       3) ST14_GO_enrichment_retained_terms.csv, ST15_GO_NESmap_vs_DelDup_correlations.csv
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Kuldeep's local system; Else: run from the folder that holds data/ and results/
if( flag_local == 1){
  project_dir <- "/Users/kuldeepkumar/Documents/Projects/Temp_CNV_GWAS_oct2022"
  code_dir <- file.path(project_dir, "Release_NG_R1_29oct2026/GitHub/code_GOterms_supplement/GO_enrichment")
  out_dir <- file.path(project_dir, "Release_NG_R1_29oct2026/Zenodo/results/GO_fgsea")
  filename_ST_C <- file.path(project_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
  filename_perm_C <- file.path(project_dir, "outputs_FinalStats_12oct2026/perm_clean_GlasserLH180.rds")
} else {
  project_dir <- getwd()
  code_dir <- file.path(project_dir, "code_GOterms_supplement/GO_enrichment")
  out_dir <- file.path(project_dir, "results/GO_fgsea")
  filename_ST_C <- file.path(project_dir, "data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx")
  filename_perm_C <- file.path(project_dir, "data/perm_clean_GlasserLH180.rds")
}
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results/res_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(readxl) })
source(file.path(code_dir, "fn_spin_functions.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

array_sources <- c("go_mf", "go_cc", "go_bp")
size_min_N <- 50; size_max_N <- 700   # Methods: gene-set size (unique symbols) > 50 and < 700
n_regions_min_N <- 5                  # Methods: significant in more than 5 regions
q_thr_N <- 0.05

fn_gate <- function(ok, msg) { cat(sprintf("[GATE %s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", msg)); if (!isTRUE(ok)) stop("gate failed: ", msg) }
fn_normalize_glasser_name <- function(x) tolower(gsub("\\.", "-", sub("^L_", "", x)))

# region order of the fgsea run = columns of the expression matrix (roi_index 1..180)
e_exp <- new.env(); load(file.path(data_dir, "df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData"), envir = e_exp)
roi_names <- colnames(get(ls(e_exp)[1], envir = e_exp)); rm(e_exp)

# Fig. 2 effect sizes (ST4), joined by region name (raw effect sizes, as in ST15)
df_ST4 <- as.data.table(read_excel(filename_ST_C, "ST4"))[gene_set_exp == "All_data_analysis_Fig2"]
mat_beta <- sapply(c(del = "DEL", dup = "DUP"), function(t) df_ST4[cnv_type == t][match(fn_normalize_glasser_name(roi_names), fn_normalize_glasser_name(gene_set_name)), Estimate])
fn_gate(length(roi_names) == 180 && !anyNA(mat_beta), "ST4 Fig. 2 effect sizes joined to the 180 fgsea regions by name")

#-----------------------------------------------------------------------------
# PART B: Long table (all tested GO terms x 180 regions)
#-----------------------------------------------------------------------------

list_GO <- local({ e <- new.env(); load(file.path(results_dir, "results_GO_GeneSets_fgsea_pos_flagSubset_0.RData"), envir = e); e$myGO_min_max_list })
df_terms <- rbindlist(lapply(array_sources, function(s) { n <- lengths(lapply(list_GO[[s]], unique)); data.table(pathway_cat = s, pathway = names(n), pathway_size = n) }))
df_terms <- df_terms[pathway_size > size_min_N & pathway_size < size_max_N]
print(df_terms[, .N, by = pathway_cat])

df_long_all <- rbindlist(lapply(array_sources, function(s) rbindlist(lapply(1:180, function(i) {
  load(file.path(results_dir, sprintf("results_fgsea_pos_pathway_%s_flagSubset_0_roi_%d.RData", s, i)))
  df_fgRes[, .(pathway_cat = s, pathway, region = roi_names[i], size = as.numeric(size), ES, NES, p = pval)] }))))
fn_gate(!anyDuplicated(df_long_all, by = c("pathway", "region")), "(pathway, region) unique in the fgsea output")
fn_gate(all(df_terms$pathway %in% df_long_all$pathway), sprintf("all %d terms of the size rule present in the fgsea output", nrow(df_terms)))

df_long <- df_long_all[df_terms[, .(pathway, pathway_size)], on = "pathway", nomatch = NULL]
df_long[, p_fdr := p.adjust(p, "BH")]   # global BH over all tested term x region pairs
cat(sprintf("global BH family: %d term x region tests\n", nrow(df_long)))

#-----------------------------------------------------------------------------
# PART C: Prevalence filter and NES maps
#-----------------------------------------------------------------------------

df_term_info <- df_long[, .(n_regions_tested = .N, n_sig_regions = sum(p_fdr < q_thr_N)), by = .(pathway_cat, pathway, pathway_size)]
df_term_info[, retained := n_regions_tested == 180 & n_sig_regions > n_regions_min_N]
df_filt <- df_term_info[retained == TRUE]
df_filt[, pathway_cat := factor(pathway_cat, array_sources)]; setorder(df_filt, pathway_cat, pathway); df_filt[, pathway_cat := as.character(pathway_cat)]
print(df_filt[, .N, by = pathway_cat])

mat_NES <- as.matrix(dcast(df_long[pathway %in% df_filt$pathway], pathway ~ region, value.var = "NES"), rownames = "pathway")[df_filt$pathway, roi_names, drop = FALSE]
fn_gate(!anyNA(mat_NES) && nrow(mat_NES) == nrow(df_filt), "NES maps complete (180 regions per retained term)")

#-----------------------------------------------------------------------------
# PART D: Spearman correlation with Del/Dup maps (spin test), ST15
#-----------------------------------------------------------------------------

perm_id <- as.matrix(readRDS(filename_perm_C)); if (min(perm_id) == 0) perm_id <- perm_id + 1L; storage.mode(perm_id) <- "integer"
fn_gate(nrow(perm_id) == 180, sprintf("spin rotations: 180 regions x %d rotations", ncol(perm_id)))
list_sp <- lapply(c(del = "del", dup = "dup"), function(t) rbindlist(lapply(df_filt$pathway, function(k) {
  r <- spin_cor_one_map(mat_beta[, t], mat_NES[k, ], perm_id, "spearman", TRUE); data.table(r = r$rho, p = r$p_spin) })))   # NES map rotated
df_ST15 <- data.table(pathway_cat = df_filt$pathway_cat, pathway_name = df_filt$pathway, pathway_size = df_filt$pathway_size, n_sig_regions = df_filt$n_sig_regions,
                      rho_nes_del = list_sp$del$r, pval_cor_es_del = list_sp$del$p, rho_nes_dup = list_sp$dup$r, pval_cor_es_dup = list_sp$dup$p)
df_ST15[, `:=`(pval_cor_es_del_FDR = p.adjust(pval_cor_es_del, "BH"), pval_cor_es_dup_FDR = p.adjust(pval_cor_es_dup, "BH"))]
setcolorder(df_ST15, c("pathway_cat", "pathway_name", "pathway_size", "n_sig_regions", "rho_nes_del", "pval_cor_es_del", "pval_cor_es_del_FDR",
                       "rho_nes_dup", "pval_cor_es_dup", "pval_cor_es_dup_FDR"))

#-----------------------------------------------------------------------------
# PART E: Write
#-----------------------------------------------------------------------------

df_ST14 <- df_long[pathway %in% df_filt$pathway, .(pathway_cat, pathway, size, region, ES, NES, p, p_fdr)]
df_ST14[, `:=`(pathway = factor(pathway, df_filt$pathway), region = factor(region, roi_names))]; setorder(df_ST14, pathway, region)
df_ST14[, `:=`(pathway = as.character(pathway), region = as.character(region))]

df_long_out <- df_long[, .(pathway_cat, pathway, size, region, ES, NES, p, p_fdr)]
df_long_out[, region := factor(region, roi_names)]; setorder(df_long_out, pathway_cat, pathway, region); df_long_out[, region := as.character(region)]
f_long <- file.path(out_dir, "GO_regional_enrichment_fgsea_all_terms.csv")
fwrite(df_long_out, f_long); system2("gzip", c("-f", "-n", shQuote(f_long)))
fwrite(df_term_info, file.path(out_dir, "GO_term_catalogue.csv"))
fwrite(df_ST14, file.path(out_dir, "ST14_GO_enrichment_retained_terms.csv"))
fwrite(df_ST15, file.path(out_dir, "ST15_GO_NESmap_vs_DelDup_correlations.csv"))
print(df_ST15)
cat("done:", out_dir, "\n")
