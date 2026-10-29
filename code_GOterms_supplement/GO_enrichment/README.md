# GO-term enrichment of the 180 cortical gene sets

Code for the regional GO-term enrichment analysis (Supplementary Tables 14 and 15; Supplementary Fig. 20).

## Scripts

| Script | What it does |
|---|---|
| `01_GO_fgsea_per_region.R` | For each Glasser-180 region, ranks the genes with AHBA regional-specificity Z-score >= 0.5 by Z-score and runs `fgsea(pathways, stats, eps = 0, scoreType = "pos")` against MSigDB C5 GO gene sets (BP, CC, MF). Writes one result file per GO category x region. |
| `02_GO_merge_by_name_ST14_ST15.R` | Joins the per-region results by GO-term and region **name**. Applies the term-size rule (>50 and <700 genes) and a global Benjamini-Hochberg FDR over all term x region tests. Keeps the terms tested in all 180 regions that are significant (q < 0.05) in more than 5 regions. Correlates each retained term's NES map with the deletion and duplication effect-size maps (Spearman; spin test with the NES map rotated; BH FDR per CNV type). |
| `fn_spin_functions.R` | Spin-permutation helper functions (`spin_cor_one_map()`). |

## Merging by name

The 2023 summary step combined the per-region outputs by row position (`cbind`). That is not correct, for three reasons:
- each region's rows are sorted by that region's own NES;
- each region keeps a different number of terms;
- the term labels were taken from region 1.

As a result, only region 1 (L_V1) carried the right values. `02_` joins every value by (term, region).

## Results (21 Oct 2026)

`02_` reproduces the final Supplementary Tables exactly:
- ST14: 9 terms x 180 regions; size, ES, NES, p and p_fdr identical.
- ST15: rho, p-spin and FDR identical.
- The long table of all 2,519 tested terms (452,989 term x region rows)

## Running

The scripts need inputs that are not in the GitHub repository: the AHBA regional expression object (`df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData`), the per-region fgsea result files of the 2023 run (540 files, about 100 MB; `02_` reads them from `results/res_enrichment/`) and the Glasser spin rotations. The GO gene sets used and all results are in the Zenodo archive (`results/GO_fgsea/`). The code is here to document the analysis behind Supplementary Tables 14 and 15.

- `01_` is set to `flag_run_fgsea = 0`: the paper uses the per-region results of the 2023 run (about 100 MB;). With these inputs, re-running fgsea reproduces the 2023 enrichment scores and sizes exactly (checked for GO CC, region 5).
- `01_` uses the GO gene sets saved with that run (`results/GO_fgsea/inputs/results_GO_GeneSets_fgsea_pos_flagSubset_0.RData` on Zenodo; msigdbr 7.5.1-era release). A newer msigdbr release gives slightly different gene sets.
- fgsea p-values come from random sampling, so a re-run gives slightly different p-values.
- Inputs: `df_AHBA_geneExp_GlasserLH_180regions_gene_x_roi.RData` (AHBA expression; gene symbols; genes x regions); Supplementary Table 4 (Fig. 2 effect sizes); `perm_clean_GlasserLH180.rds` (spin rotations; see `results/spin_null_rotations/README.md` on Zenodo).
- Outputs are written to `results/GO_fgsea/` on Zenodo:
  - `GO_regional_enrichment_fgsea_all_terms.csv.gz`
  - `GO_term_catalogue.csv`
  - `ST14_GO_enrichment_retained_terms.csv`
  - `ST15_GO_NESmap_vs_DelDup_correlations.csv`
