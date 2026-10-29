#-------------------------------------------------------------------------------
# Summary: R-script to run fgsea (GO terms) on the ranked gene list of each of the 180 cortical regions
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> For each Glasser-180 region, ranks the genes with AHBA regional-specificity Z-score >= 0.5 (gene symbols,
#             15,633 AHBA genes) by their Z-score and runs fgsea against MSigDB C5 GO gene sets (BP, CC, MF; gene symbols).
#             This is the analysis step of the original 2023 run (RscriptMac_Fig5_fgsea_180geneset_final_plot_22June2023.R);
#             the summary of the per-region results is in 02_GO_merge_by_name_ST14_ST15.R.
# Key steps:
#       1) PART A: Parameters and input files
#       2) PART B: Ranked gene lists (Z > 0.5) per region
#       3) PART C: GO gene sets (msigdbr) and fgsea per region; one RData file per GO category x region
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026 (analysis code of 22 June 2023)
# Version: 2 (Revision 1)
#
#-----
# Input files:
#       1) df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData (AHBA expression, gene symbols, 15,633 genes x 180 regions).
#          Check: re-running fgsea on these ranked lists with the saved GO gene sets reproduces the 2023 ES and size
#          exactly (go_cc, region 5: 882/882 terms); the ENSEMBL / HPA-brain-detected matrix (13,678 genes) does not.
#
#-----
# Output:
#       1) <results_dir>/results_fgsea_pos_pathway_<go_mf|go_cc|go_bp>_flagSubset_0_roi_<1..180>.RData (data.table df_fgRes)
#       2) <results_dir>/results_GO_GeneSets_fgsea_pos_flagSubset_0.RData (GO gene sets used)
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Kuldeep's local system; Else: run from the folder that holds data/ and results/
if( flag_local == 1){
  project_dir <- "/Users/kuldeepkumar/Documents/Projects/Temp_CNV_GWAS_oct2022"
} else {
  project_dir <- getwd()
}
data_dir <- paste0(project_dir,"/data")
results_dir <- paste0(project_dir,"/results/res_enrichment")

flag_run_fgsea = 0   # 0: do not re-run (results of the 2023 run are used); 1: run fgsea for every region

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(msigdbr)
library(fgsea)

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

input_zscore_threshold = 0.5
flag_subset_myGO = 0                      # 0: all GO gene sets (size rule applied later, in 02_)
array_sources = c("go_mf","go_cc","go_bp")

#-----------------------------------------------------------------------------
# PART B: Ranked gene lists (Z > 0.5) per region
#-----------------------------------------------------------------------------

e_exp = new.env(); load(file=paste0(data_dir,"/df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData"), envir = e_exp)
df_GeneExp = get(ls(e_exp)[1], envir = e_exp)
# 15633 genes (gene symbols) x 180 regions

df_GeneExp_Zscore = as.data.frame(t(scale(t(df_GeneExp))))
df_GeneExp_Zscore_thres = df_GeneExp_Zscore
df_GeneExp_Zscore_thres[df_GeneExp_Zscore_thres < input_zscore_threshold] = 0

fn_ranked_list <- function(roi_index){
  input_gene_list = df_GeneExp_Zscore_thres[,roi_index]
  input_gene_names = rownames(df_GeneExp)
  remove_indx = which(input_gene_list == 0)
  input_gene_list = input_gene_list[-c(remove_indx)]
  input_gene_names = input_gene_names[-c(remove_indx)]
  # order the gene-list; remove duplicates if any
  names(input_gene_list) = input_gene_names
  input_ranked_gene_list = sort(input_gene_list,decreasing = TRUE)
  input_ranked_gene_list[!duplicated(names(input_ranked_gene_list))]
}

#-----------------------------------------------------------------------------
# PART C: GO gene sets and fgsea per region
#-----------------------------------------------------------------------------

if (flag_run_fgsea == 1) {
  list_msigdbr = list(go_bp = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
                      go_cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"),
                      go_mf = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"))
  myGO_min_max_list = vector(mode = "list", length = 0)

  for (flag_source_gene_set_suffix in array_sources) {
    print(flag_source_gene_set_suffix)
    in_msigdbr_df = list_msigdbr[[flag_source_gene_set_suffix]]
    # the ranked lists use gene symbols (rownames of df_GeneExp)
    myGO_subset = split(x = in_msigdbr_df$gene_symbol, f = in_msigdbr_df$gs_name)
    myGO_min_max_list[[flag_source_gene_set_suffix]] = myGO_subset

    for (roi_index in c(1:180)) {
      fgRes <- fgsea::fgsea(pathways = myGO_subset,
                            stats = fn_ranked_list(roi_index),
                            eps = 0.0,
                            scoreType = "pos") %>%
        as.data.frame() %>%
        arrange(desc(NES))
      df_fgRes = data.table::as.data.table(fgRes)
      save(df_fgRes,file=paste0(results_dir,"/results_fgsea_pos_pathway_",flag_source_gene_set_suffix,"_flagSubset_",flag_subset_myGO,"_roi_",roi_index,".RData"))
    }
  }
  save(myGO_min_max_list,file=paste0(results_dir,"/results_GO_GeneSets_fgsea_pos_flagSubset_",flag_subset_myGO,".RData"))
} else {
  cat("flag_run_fgsea = 0: per-region fgsea results of the 2023 run are used (", results_dir, ")\n")
}
