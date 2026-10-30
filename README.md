## Project: DelDup_Asymmetry_Cortex

## Project Title 
Genomic deletion and duplication effects on cognitive ability are asymmetrically distributed along the cortical sensorimotor–association hierarchy

## Citation
If you use this project or code in your research, please cite the following: 

- **Papers:**
1) Core Method: **Kuldeep Kumar**, **Sayeh Kazem**, et al. "Mirror effect of genomic deletions and duplications on cognitive ability across the human cerebral cortex." bioRxiv 2025.01.06.631492; 2025.
2) PGC-CNV application across 6 Psychiatric disorders: **Worrawat Engchuan**, et al. "Psychiatric disorders converge on common pathways but diverge in cellular context, spatial distribution, and directionality of genetic effects." medRxiv 2025.07.11.25331381; 2025.


## Repository Contents
Code for the analysis and figures for the Mirror effect of genomic deletions and duplications across the cortex. The repository is organized into three main components: 
* CC-GSBA: cerebral-cortex gene-set burden analysis (CC-GSBA) pipeline, with simulated toy phenotype data and run;
* Figs: R Markdown (.Rmd) files containing the steps and code to generate the main figures and their associated statistics;
* SupFigs: R scripts containing the steps and code for generating the supplementary figures. 

In addition, the gene-sets used in the paper (due to size), along with a version of this code, are hosted on **Zenodo** at: [https://doi.org/10.5281/zenodo.23776818](https://doi.org/10.5281/zenodo.23776818)

## CC-GSBA Workflow

Our methodology for associating copy number variants (CNVs) with cognitive ability leverages gene expression data across cortical regions. The analysis integrates genotypic and phenotypic data from large population cohorts with gene expression data from the Allen Human Brain Atlas (AHBA).

- Gene Assignment: Each gene is mapped to cortical regions where it shows preferential expression (z-score > 0.5). This process is repeated for all genes, resulting in 180 distinct regional gene sets.
- CNV Aggregation: For each individual, we calculate two burden metrics for each cortical region:
- - Gene-Set Burden: The count of CNVs overlapping genes within a specific regional gene set.
- - Genome Burden: The count of CNVs outside that specific regional gene set, serving as a genome-wide control.
- Regression Analysis: A regression model is applied to test for a region-specific association between CNV burden and cognitive ability, while controlling for the genome-wide burden and other covariates.  Cognitive ability (Y) is modeled as a function of the Gene-Set Burden and the Genome Burden. The model also adjusts for standard covariates such as age, sex, and ancestry.

<p align="center">
 <img src="CC-GSBA/CC-GSBA_Workflow.png" alt="CC-GSBA_Workflow" width="1000"/>
</p>

#### The core script:
The CC-GSBA script is the core pipeline for testing the association between CNVs aggregated in cortex-specific gene sets and cognitive ability.

#### Data Used in This Repository : 

This repository's pipeline runs on synthetic and non-sensitive toy data to allow for development, testing, and demonstration without the need for access to real human genetic data. The datasets included are:

- Toy Dataset of CNV by Individual: A simulated matrix representing CNV status across individuals.

- Toy Dataset of Gene by CNV: A mapping of CNVs to the genes they affect.

- Simulated Cognitive Ability Measure: A synthetic phenotype/score for the individuals.

- Quality Control Individuals: A subset of individuals who pass a simulated quality control process.

- Real Gene Sets: The actual gene sets (derived from CC_GSBA, assigning section) that are used to aggregate CNV burden.
  
#### File Types
The main CC-GSBA pipeline is implemented in Python Notebooks (*.ipynb).
The toy datasets and real gene sets are all in (*.tsv) formats. See **Zenodo** at: [https://doi.org/10.5281/zenodo.17593039](https://doi.org/10.5281/zenodo.17593039)
for a full set of gene sets. 

#### 💻 Usage

1. Prerequisites: Ensure you have Python 3.11 (for the CC-GSBA notebook) and R version 4.2.2 (for the figures) installed, along with the necessary packages (see requirements.txt).
   
2. Clone the repository:

```bash


git clone https://github.com/SayehKazem/CC_GSBA.git


cd CC_GSBA


```

3. Run the Core Pipeline: Open and run the CC-GSBA/core_pipeline.ipynb notebook to see the full workflow demonstrated on the synthetic data.


##  Code and source data for the figures

This repository contains the code and the Source Data that reproduce the figures of the manuscript (Nature Genetics, revision 1). We aggregated rare deletions and duplications (copy number variants, CNVs) carried by 258,292 individuals into gene sets defined by relative high expression in each of 180 cortical regions (Glasser parcellation, Allen Human Brain Atlas). We then estimated the effect of each regional gene set on cognitive ability with cortical-gene-set burden analysis (CC-GSBA). The figures map the deletion and duplication effects across the cortex, relate them to the sensorimotor-association hierarchy (H = −T1w/T2w) and to 20 normative cortical maps, and test their robustness across cohorts, phenotypes, gene strata and model specifications.

- **Zenodo** archive (code, Source Data, gene sets, gene scores, GO enrichment results, CC-GSBA toy example):[https://doi.org/10.5281/zenodo.23776818](https://doi.org/10.5281/zenodo.23776818)
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
See **README_Figs_SupFigs** for details on scripts and source data for each figure and supplement figure. 

## 📜 Copyright and Licensing

This repository contains both software (code) and non-software content (data, figures, documentation).

* **Code License:** All source code (including Python and R scripts) is licensed under the **MIT License**.
* **Data and Documentation License:** All generated data (including gene sets), figures, and documentation (including the README) are licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0) License**.

**Copyright (c) 2025** Kuldeep Kumar, Sayeh Kazem, Guillaume Dumas, and Sebastien Jacquemont.

