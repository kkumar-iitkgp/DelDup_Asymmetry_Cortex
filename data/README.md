# data/: inputs and Source Data

## Source Data workbook

`SourceData_R1.xlsx` (Figs. 2-5 and Supplementary Figs. 1-23) was assembled by `source_data_build/03_assemble_SourceData.R` (Zenodo archive, `code_and_source_data/source_data_build/`).
It starts with a README sheet (sheet, figure, panels, description, rows, columns). Each other sheet is one CSV in `source_data_sheets/`.
Each CSV is the data frame passed to the plotting call, written by the figure script, at full precision.

| Sheet | Figure | Panels | Description | Rows | Columns |
|---|---|---|---|---|---|
| `Fig2a-d` | Fig. 2 | a-d | 180 regions: Del and Dup effect sizes, SE, p, Bonferroni flags, significance class, z(Del), z(Dup), cDEC, H, hierarchy rank, pole | 180 | `gene_set_name`, `glasser_region`, `Estimate_DEL`, `pvalue_DEL`, `pvalue_bonferroni_360_DEL`, `Estimate_DUP`, `pvalue_DUP`, `pvalue_bonferroni_360_DUP`, `significance_class`, `z_DEL`, `z_DUP`, `cDEC`, `T1T2ratio`, `H`, `H_rank` |
| `Fig2_stats` | Fig. 2 | a-d | Printed statistics: significant regions (pies), rho(Del,Dup) and rho(cDEC,H) with p-spin | 6 | `panel`, `statistic`, `value`, `printed`, `source` |
| `Fig3a-b` | Fig. 3 | a-b | 20 normative maps: domain, rho(cDEC, map), rho(Del, map), rho(Dup, map), p-spin | 20 | `map`, `source_column`, `domain`, `rho_cDEC`, `pspin_cDEC`, `star_a`, `rho_del`, `pspin_del`, `rho_dup`, `pspin_dup`, `outline_del`, `outline_dup`, `delta_rho_dup_minus_del`, `pspin_delta`, `n_rot` |
| `Fig3c-d` | Fig. 3 | c-d | 180 regions: composite ranks, H rank, H-adjusted ranks (Del, Dup, composites) | 180 | `gene_set_name`, `H`, `rank_H`, `energy_composite`, `energy_composite_rank_percentile`, `energy_rank_resid_H`, `expansion_composite`, `expansion_composite_rank_percentile`, `expansion_rank_resid_H`, `del_rank_resid_H`, `dup_rank_resid_H` |
| `Fig3c-d_stats` | Fig. 3 | c-d | Partial rho | H and p-spin | 10 | `panel`, `statistic`, `value`, `pspin`, `n_rot` |
| `Fig3e` | Fig. 3 | e | Yeo-17 networks: mean and s.e.m. of Del and Dup effect sizes; rho | 17 | `network_index`, `network`, `colour_hex`, `n_vertices`, `n_glasser_regions`, `mean_del`, `sem_del_vertices`, `sem_del_regions`, `mean_dup`, `sem_dup_vertices`, `sem_dup_regions`, `rho_networks`, `p_parametric`, `pspin_profile_pair`, `n_rot` |
| `Fig4a-d` | Fig. 4 | a-d | Regional-specificity thresholds x region x CNV type: effect size, SE, p, significance, pole | 1440 | `gene_set_exp`, `threshold_z`, `cnv_type`, `gene_set_name`, `glasser_region`, `Estimate`, `pvalue`, `pvalue_bonferroni_360`, `significant`, `H`, `pole` |
| `Fig4b-e_stats` | Fig. 4 | b-e | Pole means, 95% CI, KS D, exact p, FDR (94 tests, ST11); rho(cDEC,H), p-spin, FDR (4 thresholds) | 20 | `panel`, `gene_set_exp`, `threshold_z`, `cnv_type`, `pole`, `mean_effect_size`, `ci95_lo`, `ci95_hi`, `n_regions`, `KS_D`, `KS_pvalue_exact`, `KS_q_FDR_94tests`, `star`, `rho_cDEC_H`, `pspin_cDEC_H`, `q_FDR_4_thresholds` |
| `Fig5a-b` | Fig. 5 | a-b | LOEUF bins (aggregated gene counts) and genes per LOEUF stratum | 209 | `panel`, `cnv_type`, `LOEUF_bin_lower`, `LOEUF_bin_upper`, `n_unique_genes`, `LOEUF_stratum` |
| `Fig5c-e` | Fig. 5 | c-e | LOEUF stratum x region x CNV type: effect size, SE, p, significance, cDEC, hierarchy rank | 540 | `gene_set_exp`, `LOEUF_category`, `gene_set_name`, `glasser_region`, `Estimate_DEL`, `pvalue_DEL`, `pvalue_bonferroni_1080_DEL`, `Estimate_DUP`, `pvalue_DUP`, `pvalue_bonferroni_1080_DUP`, `significance_class`, `z_DEL`, `z_DUP`, `cDEC`, `T1T2ratio`, `H`, `H_rank` |
| `Fig5e_stats` | Fig. 5 | c-e | Per stratum: rho(cDEC,H), p-spin, significant regions | 3 | `gene_set_exp`, `LOEUF_category`, `n_sig_regions_DEL_bonferroni_1080`, `n_sig_regions_DUP_bonferroni_1080`, `rho_DEL_DUP_all`, `pspin_DEL_DUP_all`, `rho_DEL_DUP_sig_descriptive`, `n_sig_union`, `rho_cDEC_H`, `pspin_cDEC_H`, `printed_e` |
| `SupFig1a-c` | Supplementary Fig. 1 | a-c | Unique genes (% of 13,678), total genes, CNV carriers (ST3) | 6 | `panel`, `quantity`, `cnv_type`, `count`, `percent_of_13678_genes` |
| `SupFig1g-h` | Supplementary Fig. 1 | g-h | Unique deleted / duplicated genes per region (ST4) | 180 | `region`, `n_unique_genes_DEL`, `n_unique_genes_DUP` |
| `SupFig2a-d` | Supplementary Fig. 2 | a-d | Effect-size maps per data subset (180 regions), Bonferroni flags | 720 | `subset`, `region`, `effect_size_DEL`, `sig_bonferroni_DEL`, `effect_size_DUP`, `sig_bonferroni_DUP` |
| `SupFig2e-h` | Supplementary Fig. 2 | e-h | Region-wise effect sizes between data subsets; rho and p-spin | 720 | `cnv_type`, `x_subset`, `y_subset`, `region`, `x_effect_size`, `y_effect_size`, `rho`, `p_spin` |
| `SupFig3a-b` | Supplementary Fig. 3 | a-b | HPA 73 regions: Del and Dup effect sizes, p, significance class, core pole | 73 | `region_id`, `region_name`, `core_pole`, `effect_size_DEL`, `effect_size_DUP`, `p_DEL`, `p_DUP`, `significance_class`, `label` |
| `SupFig4a-b` | Supplementary Fig. 4 | a-b | BrainSpan 11 regions: Del and Dup effect sizes, p, significance class, pole | 11 | `region_code`, `region_name`, `pole`, `effect_size_DEL`, `effect_size_DUP`, `p_DEL`, `p_DUP`, `significance_class`, `label` |
| `SupFig5a-c` | Supplementary Fig. 5 | a-c | Desikan-34: Del and Dup effect sizes, p, Bonferroni flags, pole | 34 | `region`, `pole`, `effect_size_DEL`, `effect_size_DUP`, `p_DEL`, `p_DUP`, `sig_bonferroni68_DEL`, `sig_bonferroni68_DUP`, `significance_class`, `label` |
| `SupFig6a` | Supplementary Fig. 6 | a | 16 co-expression modules: Del and Dup effect sizes, p, significance class | 16 | `module`, `effect_size_DEL`, `effect_size_DUP`, `p_DEL`, `p_DUP`, `significance_class`, `label` |
| `SupFig6b` | Supplementary Fig. 6 | b | Module PC1 maps (z-scored), 180 regions x 16 modules | 180 | `region`, `M1`, `M2`, `M3`, `M4`, `M5`, `M6`, `M7`, `M8`, `M9`, `M10`, `M11`, `M12`, `M13`, `M14`, `M15`, `M16` |
| `SupFig7a-d` | Supplementary Fig. 7 | a-d | Profile statistics per sensitivity analysis (ST7) | 22 | `gene_set_exp`, `label`, `n_regions`, `SDbySE_DEL`, `SDbySE_DUP`, `I2_DEL`, `I2_DUP`, `rho_DEL_H`, `pspin_DEL_H`, `rho_DUP_H`, `pspin_DUP_H`, `rho_cDEC_H`, `pspin_cDEC_H`, `rho_DEL_DUP_all`, `pspin_DEL_DUP_all`, `rho_DEL_DUP_sig`, `pspin_DEL_DUP_sig`, `I2_pass_Del`, `I2_pass_Dup`, `Pair_interpretable`, `open`, `fadeD`, `fadeU`, `band` |
| `SupFig8a-b` | Supplementary Fig. 8 | a-b | 180 regions: H, T1w/T2w, percentile, rank, pole | 180 | `region`, `H`, `T1T2ratio`, `percentile`, `rank`, `pole` |
| `SupFig9a-d` | Supplementary Fig. 9 | a-d | Profile statistics per phenotype (ST7) | 10 | `label`, `gene_set_exp`, `block`, `negative_control`, `SDbySE_DEL`, `SDbySE_DUP`, `I2_DEL`, `I2_DUP`, `I2_pass_Del`, `I2_pass_Dup`, `Pair_interpretable`, `rho_DEL_H`, `pspin_DEL_H`, `rho_DUP_H`, `pspin_DUP_H`, `rho_cDEC_H`, `pspin_cDEC_H`, `rho_DEL_DUP_all`, `rho_DEL_DUP_sig` |
| `SupFig10a-d` | Supplementary Fig. 10 | a-d | UK Biobank vs All of Us, 180 regions; rho, p-spin, sign concordance | 720 | `panel`, `cnv_type`, `x_profile`, `y_profile`, `region`, `x_effect_size`, `y_effect_size`, `rho`, `p_spin`, `sign_concordance_pct` |
| `SupFig11a-j` | Supplementary Fig. 11 | a-j | All of Us GATK-SV size thresholds, 180 regions; rho, p-spin, sign concordance | 1800 | `panel`, `cnv_type`, `x_profile`, `y_profile`, `region`, `x_effect_size`, `y_effect_size`, `rho`, `p_spin`, `sign_concordance_pct` |
| `SupFig12` | Supplementary Fig. 12 | all | 20 normative maps x 180 regions: value and rank percentile | 3600 | `order`, `map_group`, `map_name`, `source_column`, `region`, `value`, `sign_reversed`, `rank_percentile` |
| `SupFig13a` | Supplementary Fig. 13 | a | 180 regions: composite rank percentiles, Del and Dup effect sizes | 180 | `region`, `energy_composite_percentile`, `expansion_composite_percentile`, `effect_size_DEL`, `effect_size_DUP` |
| `SupFig13a_stats` | Supplementary Fig. 13 | a | Unadjusted rho and p-spin | 4 | `dimension`, `cnv_type`, `rho`, `p_spin` |
| `SupFig13b` | Supplementary Fig. 13 | b | 19 maps: partial rho | H and p-spin, Del and Dup | 38 | `map_group`, `map_name`, `cnv_type`, `rho_Hadj`, `p_spin_Hadj` |
| `SupFig13c` | Supplementary Fig. 13 | c | Replication: partial rho | H and I_H per sample; p-spin | 30 | `sample`, `quantity`, `value`, `p_spin` |
| `SupFig14a-c` | Supplementary Fig. 14 | a-c | Tissue, epoch, cell type x 180 regions: effect sizes, Bonferroni flags, pole, KS star (ST11) | 2160 | `row`, `gene_set_exp`, `region`, `pole`, `effect_size_DEL`, `sig_bonferroni360_DEL`, `effect_size_DUP`, `sig_bonferroni360_DUP`, `KS_star_DEL`, `KS_star_DUP` |
| `SupFig14_counts` | Supplementary Fig. 14 | column 1 | Median gene-set size and unique Del / Dup genes per stratum (ST4) | 36 | `gene_set_exp`, `Type`, `value`, `label`, `block` |
| `SupFig15a-d` | Supplementary Fig. 15 | a-d | Profile statistics per biological stratum (ST7) | 28 | `gene_set_exp`, `label`, `n_regions`, `SDbySE_DEL`, `SDbySE_DUP`, `I2_DEL`, `I2_DUP`, `rho_DEL_H`, `pspin_DEL_H`, `rho_DUP_H`, `pspin_DUP_H`, `rho_cDEC_H`, `pspin_cDEC_H`, `rho_DEL_DUP_all`, `pspin_DEL_DUP_all`, `rho_DEL_DUP_sig`, `pspin_DEL_DUP_sig`, `I2_pass_Del`, `I2_pass_Dup`, `Pair_interpretable`, `open`, `fadeD`, `fadeU`, `band` |
| `SupFig16a-b` | Supplementary Fig. 16 | a-b | Age, sex, ancestry, assessment x 180 regions: effect sizes, Bonferroni flags, pole, KS star | 1260 | `row`, `gene_set_exp`, `region`, `pole`, `effect_size_DEL`, `sig_bonferroni360_DEL`, `effect_size_DUP`, `sig_bonferroni360_DUP`, `KS_star_DEL`, `KS_star_DUP` |
| `SupFig17` | Supplementary Fig. 17 | left/right | Robustness specifications: concordance with G1 (ST6), rho(cDEC,H) and p-spin (ST7) | 20 | `family`, `specification`, `gene_set_exp`, `concordance_rho_DEL_vs_G1`, `concordance_rho_DUP_vs_G1`, `SDbySE_DEL`, `SDbySE_DUP`, `I2_DEL`, `I2_DUP`, `I2_pass_Del`, `I2_pass_Dup`, `Pair_interpretable`, `rho_cDEC_H`, `pspin_cDEC_H`, `primary_rho_cDEC_H` |
| `SupFig18a-e` | Supplementary Fig. 18 | a-e | Genes: LOEUF, LOEUF stratum, pHaplo, pTriplo, missense Z, TMC, similarity to H (= -TMC) | 14500 | `ensembl_gene_id`, `gene_symbol`, `LOEUF`, `LOEUF_bin`, `pHaplo`, `pTriplo`, `missense_Z`, `TMC`, `similarity_H` |
| `SupFig18_stats` | Supplementary Fig. 18 | a-e | Spearman rho and n per pair; genes per LOEUF stratum | 7 | `panel`, `pair`, `n`, `spearman_rho`, `gene_universe` |
| `SupFig19a-d` | Supplementary Fig. 19 | a-d | Trajectory, pHaplo, pTriplo, missense Z x 180 regions: effect sizes, Bonferroni flags, pole, KS star | 2160 | `row`, `gene_set_exp`, `region`, `pole`, `effect_size_DEL`, `sig_bonferroni360_DEL`, `effect_size_DUP`, `sig_bonferroni360_DUP`, `KS_star_DEL`, `KS_star_DUP` |
| `SupFig19_counts` | Supplementary Fig. 19 | column 1 | Median gene-set size and unique Del / Dup genes per stratum (ST4) | 36 | `gene_set_exp`, `Type`, `value`, `label`, `block` |
| `SupFig21` | Supplementary Fig. 21 | all | Individuals per X1 bin, 6 regions, Del and Dup (aggregate counts) | 84 | `cnv_type`, `region`, `region_label`, `X1_bin`, `n_individuals` |
| `SupFig22a` | Supplementary Fig. 22 | a | beta1 (inside set) and beta2 (outside set), G1, 180 regions (ST16) | 360 | `cnv_type`, `region`, `beta1_inside`, `beta2_outside` |
| `SupFig22b` | Supplementary Fig. 22 | b | Unique genic CNVs per bin of genes encompassed (aggregate counts) | 20 | `cnv_type`, `genes_per_CNV_bin`, `n_unique_cnvs`, `pct` |
| `SupFig22c` | Supplementary Fig. 22 | c | Unique genic CNVs per bin of unique carriers of the exact CNV (aggregate counts) | 42 | `cnv_type`, `unique_carriers_bin`, `n_unique_cnvs`, `pct` |
| `SupFig23a` | Supplementary Fig. 23 | a | Regional effect sizes under G0, G1, G3 (ST16) | 1080 | `spec`, `region`, `cnv`, `es` |
| `SupFig23b` | Supplementary Fig. 23 | b | G0 vs G1 and G3 vs G1, 180 regions; rho and p-spin (ST6) | 720 | `cnv_type`, `spec`, `region`, `G1`, `other`, `rho`, `p_spin` |
| `SupFig23c` | Supplementary Fig. 23 | c | cDEC vs hierarchy rank per specification; rho and p-spin (ST7) | 540 | `spec`, `region`, `hierarchy_rank`, `cDEC`, `rho`, `p_spin` |

Common column names:

- `region`: Glasser region (gene_set_name, e.g. L_V1) unless stated
- `gene_set_name`: region / gene-set name as in Supplementary Table 4
- `gene_set_exp`: analysis name (Supplementary Tables, Experiments_key)
- `cnv_type`: DEL (deletions) or DUP (duplications)
- `Estimate`: CC-GSBA effect size (z-scored cognitive ability per gene)
- `SE`: standard error
- `pvalue`: two-sided p-value
- `rho`: Spearman correlation
- `p_spin`: two-sided spin-permutation p-value (10,000 rotations)
- `pole`: sensorimotor / transition / association (H percentiles 30 / 70)
- `H`: cortical hierarchy, -T1w/T2w
- `cDEC`: cortical dosage-effect contrast z(Del) - z(Dup)
- `similarity_H`: similarity of the gene's cortical expression to H (= -TMC)
- `TMC`: similarity of the gene's cortical expression to the T1w/T2w map
- `ensembl_gene_id`: ENSEMBL gene ID (version suffix removed)
- `gene_symbol`: HGNC symbol

## inputs/

Region- or gene-level inputs read by the figure scripts. No individual-level data. The build scripts named below are in the Zenodo archive (`code_and_source_data/source_data_build/`).

| File | Content | Built by |
|---|---|---|
| `SupplementaryTables_NatGenRevision_29oct2026.xlsx` | Supplementary Tables 1-16 (final numbering); main input of all figures | copy |
| `inputs_pole_membership_Glasser180.csv` | 180 regions: H, T1w/T2w, percentile, pole (54 sensorimotor / 72 transition / 54 association) | `01_build_inputs.R` |
| `inputs_main_neuromaps20_Glasser180.csv` | 20 normative cortical maps parcellated to the 180 regions | `01b_build_inputs_main_figs.R` |
| `inputs_main_spin_rotations_Glasser180.rds` | Spin rotations, Glasser-180 (cleaned set, 9,999 rotations) | `01b_build_inputs_main_figs.R` |
| `inputs_Fig3_*` | fsaverage5 region / Yeo-17 network labels (Fig. 3 brain plots) and reference statistics (checks) | `01b_build_inputs_main_figs.R` |
| `inputs_Fig5_*` | Unique deleted / duplicated genes with LOEUF (ENSEMBL ID, symbol, LOEUF) and gene counts per LOEUF stratum | `01b_build_inputs_main_figs.R` |
| `inputs_SupFig12_13_composite_maps_Glasser180.csv` | Energy-metabolism and cortical-expansion composite maps (Fig. 3, Supplementary Figs. 12-13) | copy |
| `inputs_SupFig13_Fig3_dimension_statistics.csv`, `inputs_SupFig13_partial_rho_adjH_neuromaps20.csv` | Partial rho given H, interaction I_H and replication statistics (not held in the Supplementary Tables) | copy |
| `inputs_HPA73_regions.csv`, `inputs_HPA73_core_poles.csv` | HPA 73 regions: names and sensorimotor / association cores | `01_build_inputs.R` |
| `inputs_Desikan34_poles.csv` | Desikan-34 poles on H (30th / 70th percentiles) | `01_build_inputs.R` |
| `inputs_KS_alternative_parcellations.csv` | KS tests for HPA-73 cores, BrainSpan-11 poles, Desikan-34 poles (nominal; not in the 94-test family of ST11) | `01_build_inputs.R` |
| `inputs_AoU_pair_statistics.csv` | All of Us profile pairs: rho, p-spin, sign concordance, GATK-SV threshold pairs | `01_build_inputs.R` |
| `inputs_Wagstyl_module_PC1_Glasser180.csv` | PC1 of the 16 co-expression modules, 180 regions | `01_build_inputs.R` |
| `inputs_patch_ST6_data_subset_pairs.csv` | ST6 rows for UK Biobank vs five cohorts and Discovery vs Replication (Supplementary Fig. 2e-h); missing from the 29 Oct 2026 tables | `01_build_inputs.R` |
| `inputs_patch_ST7_Schaefer400_pspin.csv` | ST7 p-spin values for the Schaefer-400 row (pending in the 29 Oct 2026 tables) | `01_build_inputs.R` |
| `inputs_SupFig18_gene_scores.csv` | 19,689 protein-coding genes: LOEUF, missense Z, pHaplo, pTriplo, TMC, main gene-set flag | `06_build_SupFig18_inputs.R` |
| `inputs_SupFig21_X1_histograms.csv`, `inputs_SupFig22_*` | Aggregate counts (individuals per X1 bin; unique CNVs per bin) | `05_build_SupFig21_22_aggregates.R` |
