# Cortical gene dosage effects on cognitive ability: code and source data for the figures

This repository contains the code and the Source Data that reproduce the figures of the manuscript (Nature Genetics, revision 1). We aggregated rare deletions and duplications (copy number variants, CNVs) carried by 258,292 individuals into gene sets defined by relative high expression in each of 180 cortical regions (Glasser parcellation, Allen Human Brain Atlas). We then estimated the effect of each regional gene set on cognitive ability with cortical-gene-set burden analysis (CC-GSBA). The figures map the deletion and duplication effects across the cortex, relate them to the sensorimotor-association hierarchy (H = −T1w/T2w) and to 20 normative cortical maps, and test their robustness across cohorts, phenotypes, gene strata and model specifications.

- Zenodo archive (code, Source Data, gene sets, gene scores, GO enrichment results, CC-GSBA toy example): DOI to be added
- Citation: to be added

## Repository layout

| Folder / file | Content |
|---|---|
| `Figs/Fig2` … `Figs/Fig5` | Main figures: `FigN.R` (code), `FigN.Rmd` / `FigN.md` (code with legend, knitted), `FigN.png` (300 dpi preview) |
| `SupFigs/` | One R script per Supplementary Figure; `SupFigs/output/` holds the PNG files (300 dpi previews) and `checks/` (logs of the built-in checks) |
| `functions/` | Shared style (`00_figure_style.R`), helpers for the Supplementary Figures (`00_supfig_helpers.R`) and the Source Data writer (`fn_write_SourceData_sheet.R`) |
| `data/inputs/` | Inputs read by the figure scripts: the Supplementary Tables workbook and small region- or gene-level files (see `data/README.md`) |
| `data/source_data_sheets/` | One CSV per Source Data sheet, written by the figure scripts (the data frame passed to each plot) |
| `data/SourceData_R1.xlsx` | Source Data workbook: README sheet first, then one sheet per panel group (Figs. 2-5, then Supplementary Figs. 1-23) |
| `code_GOterms_supplement/GO_enrichment/` | Gene Ontology enrichment (fgsea) per regional gene set and the tables behind Supplementary Tables 14 and 15 (inputs and results on Zenodo) |
| `session_info.txt` | R session used to produce the figures |
| `LICENSE` | MIT |

## Figures, scripts and Source Data

| Figure | Script | Source Data sheet(s) |
|---|---|---|
| Fig. 1 | illustration (workflow) | none |
| Fig. 2 | `Figs/Fig2/Fig2.R` | `Fig2a-d`, `Fig2_stats` |
| Fig. 3 | `Figs/Fig3/Fig3.R` | `Fig3a-b`, `Fig3c-d`, `Fig3c-d_stats`, `Fig3e` |
| Fig. 4 | `Figs/Fig4/Fig4.R` | `Fig4a-d`, `Fig4b-e_stats` |
| Fig. 5 | `Figs/Fig5/Fig5.R` | `Fig5a-b`, `Fig5c-e`, `Fig5e_stats` |
| Supplementary Fig. 1 | `SupFigs/SupFig01_CNV_summary_genesets.R` (a-c, g, h; d-f illustration) | `SupFig1a-c`, `SupFig1g-h` |
| Supplementary Fig. 2 | `SupFigs/SupFig02_data_subsets.R` | `SupFig2a-d`, `SupFig2e-h` |
| Supplementary Fig. 3 | `SupFigs/SupFig03_HPA73.R` | `SupFig3a-b` |
| Supplementary Fig. 4 | `SupFigs/SupFig04_BrainSpan11.R` | `SupFig4a-b` |
| Supplementary Fig. 5 | `SupFigs/SupFig05_Desikan34.R` | `SupFig5a-c` |
| Supplementary Fig. 6 | `SupFigs/SupFig06_coexpression_modules.R` | `SupFig6a`, `SupFig6b` |
| Supplementary Fig. 7 | `SupFigs/SupFig07_sensitivity_profile_summary.R` | `SupFig7a-d` |
| Supplementary Fig. 8 | `SupFigs/SupFig08_pole_assignment.R` | `SupFig8a-b` |
| Supplementary Fig. 9 | `SupFigs/SupFig09_phenotype_specificity.R` | `SupFig9a-d` |
| Supplementary Fig. 10 | `SupFigs/SupFig10_AoU_cross_cohort.R` | `SupFig10a-d` |
| Supplementary Fig. 11 | `SupFigs/SupFig11_AoU_GATKSV_size_thresholds.R` | `SupFig11a-j` |
| Supplementary Fig. 12 | `SupFigs/SupFig12_normative_maps.R` | `SupFig12` |
| Supplementary Fig. 13 | `SupFigs/SupFig13_energy_expansion_unadjusted_replication.R` | `SupFig13a`, `SupFig13a_stats`, `SupFig13b`, `SupFig13c` |
| Supplementary Figs. 14, 16 | `SupFigs/SupFig14_16_tissue_epoch_celltype_and_demographics.R` | `SupFig14a-c`, `SupFig14_counts`, `SupFig16a-b` |
| Supplementary Fig. 15 | `SupFigs/SupFig15_biological_strata_profile_summary.R` | `SupFig15a-d` |
| Supplementary Fig. 17 | `SupFigs/SupFig17_genes_loci_size_frequency.R` | `SupFig17` |
| Supplementary Fig. 18 | `SupFigs/SupFig18_constraint_vs_cortical_gradient.R` | `SupFig18a-e`, `SupFig18_stats` |
| Supplementary Fig. 19 | `SupFigs/SupFig19_trajectory_pHaplo_pTriplo_missense.R` | `SupFig19a-d`, `SupFig19_counts` |
| Supplementary Fig. 20 | illustration (GO enrichment workflow); code in `code_GOterms_supplement/GO_enrichment/`, results on Zenodo | none |
| Supplementary Figs. 21, 22 | `SupFigs/SupFig21_22_burden_model_inputs.R` | `SupFig21`, `SupFig22a`, `SupFig22b`, `SupFig22c` |
| Supplementary Fig. 23 | `SupFigs/SupFig23_genome_burden_G0_G1_G3.R` | `SupFig23a`, `SupFig23b`, `SupFig23c` |

Some published figures had text or legend placement adjusted by hand after export (Figs. 3-5 and several Supplementary Figures). The code reproduces the plotted data, statistics and layout of the exported figures; hand-placed labels may differ.

## How to run

1. Install R (tested with R 4.2.2) and the packages: ggplot2, ggprism, ggrepel, ggpubr, ggseg, ggsegGlasser, ggsegYeo2011, ggExtra, cowplot, patchwork, dplyr, tidyr, tibble, data.table, readxl, openxlsx, scales, svglite, knitr (versions in `session_info.txt`).
2. Clone the repository and run each script from the repository root, for example `Rscript Figs/Fig2/Fig2.R` or `Rscript SupFigs/SupFig07_sensitivity_profile_summary.R`. Each script sets `project_dir` in its header (`flag_local = 2`: the working directory).
3. Figures are written as PNG at 300 dpi (previews). Set `flag_hires = 1` in a script's header to write the publication files instead (PNG at 600 dpi + SVG; these are in the Zenodo archive). The plotted data are written to `data/source_data_sheets/`. The Source Data workbook `data/SourceData_R1.xlsx` collects these sheets.

The scripts that built `data/inputs/` from internal result files, and the script that assembles the Source Data workbook, are in the Zenodo archive (`code_and_source_data/source_data_build/`); they are not needed to reproduce the figures.

## Data availability

No individual-level data are included. All statistics are region-level or gene-level aggregates, and All of Us counts follow the All of Us dissemination policy. Individual-level data are available from the cohorts (UK Biobank, the five other cohorts and All of Us) through their access procedures, as described in the manuscript.

## Conventions

Spearman correlations; two-sided spin-permutation p-values (p-spin, 10,000 rotations); Bonferroni-corrected significance within each analysis (0.05 / (2 × number of regions)); H = −T1w/T2w; cDEC = z(Del) − z(Dup); deletions in red, duplications in blue.

## License

MIT (see `LICENSE`).
