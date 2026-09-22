# GDM-DNAm-Biomarker
# DNA Methylation-Based Biomarker Discovery and Prediction of Gestational Diabetes Mellitus

## Overview

This repository contains the R code used for DNA methylation-based biomarker discovery, predictive modelling, and external validation for gestational diabetes mellitus (GDM).

The analytical workflow includes DNA methylation data preprocessing and quality control, handling of missing methylation values, epigenome-wide association analysis (EWAS), feature selection, biomarker panel development using machine-learning classifiers, and evaluation of the selected CpG biomarker panel in three external cohorts.

The repository has been organized to provide a transparent and reproducible record of the computational analyses performed in the study. Analysis scripts are numbered according to the order of the analytical workflow.

## Analysis Workflow

The overall analysis was conducted in the following stages:

1. **DNA methylation preprocessing and quality control**
   - Import of DNA methylation data
   - Sample-level quality control
   - Probe-level filtering
   - Removal of probes on sex chromosomes
   - Removal of SNP-associated and cross-reactive probes
   - Generation of processed methylation beta values

2. **Missingness filtering and imputation**
   - CpGs were filtered according to the proportion of observed methylation values.
   - CpGs with less than or equal to 80% observed values were excluded.
   - Missing beta values among the retained CpGs were imputed using `methyLImp`.
   - The resulting dataset was used as the input for subsequent EWAS.

3. **Epigenome-wide association analysis (EWAS)**
   - Quantile normalization was applied to the methylation data.
   - Differential methylation analysis was performed using `ChAMP`.
   - Multiple-testing correction was performed using the Benjamini-Hochberg procedure.
   - CpGs meeting the study-defined significance threshold (`adjusted P < 1 × 10^-6`) and having gene annotation were retained.
   - This resulted in 527 CpGs for downstream feature-selection analyses.

4. **Feature selection**
   - Three feature-selection approaches were applied to the 527 EWAS-selected CpGs:
     - Random Forest-based feature ranking
     - LASSO regression
     - Random Forest Recursive Feature Elimination (RF-RFE/backward selection)
   - The selected CpG subsets were subsequently used for biomarker panel development.

5. **Random Forest biomarker panel development**
   - Random Forest classification was applied to CpG subsets obtained from the three feature-selection approaches.
   - The analysis included repeated training/test partitioning, forward addition of candidate CpGs, cross-validation, model evaluation, and resampling-based performance assessment.

6. **Support Vector Machine biomarker panel development**
   - Support Vector Machine (SVM) classification with a radial kernel was applied to the CpG subsets obtained from the three feature-selection approaches.
   - Candidate CpGs were evaluated through an iterative panel-building procedure.
   - Model performance was assessed using multiple discrimination and classification metrics.

7. **External validation**
   - The final 11-CpG biomarker panel was evaluated in three external cohorts:
     - GARBH-INi
     - PRiDE
     - REVAMP
   - Cohort-specific logistic regression and ROC analyses were performed to evaluate the selected CpG panel.

## Repository Structure

```text
GDM-DNAm-Biomarker/
│
├── R/
│   ├── 01_Preprocessing.R
│   ├── 02_missingness_imputation.R
│   ├── 03_EWAS.R
│   ├── 04_feature_selection.R
│   ├── 05_RF_panel_building.R
│   ├── 06_SVM_panel_building.R
│   ├── 07_external_validation_GARBH_INi.R
│   ├── 08_external_validation_PRiDE.R
│   ├── 09_external_validation_REVAMP.R
│   └── README.md
│
├── data/
│   └── README.md
│
├── results/
│   └── README.md

## R Scripts

The analysis scripts are located in the `R/` directory and are numbered according to the sequence of the analytical workflow.

| Script | Analysis |
|---|---|
| `01_Preprocessing.R` | DNA methylation preprocessing and quality control |
| `02_missingness_imputation.R` | Missingness filtering and methylation-value imputation |
| `03_EWAS.R` | Epigenome-wide association analysis and identification of significant CpGs |
| `04_feature_selection.R` | Feature selection using Random Forest, LASSO, and RF-RFE/backward selection |
| `05_RF_panel_building.R` | Random Forest-based biomarker panel development and performance evaluation |
| `06_SVM_panel_building.R` | SVM-based biomarker panel development and performance evaluation |
| `07_external_validation_GARBH_INi.R` | External validation in the GARBH-INi cohort |
| `08_external_validation_PRiDE.R` | External validation in the PRiDE cohort |
| `09_external_validation_REVAMP.R` | External validation in the REVAMP cohort |

Additional information about each script and its position in the analytical workflow is provided in [`R/README.md`](R/README.md).

## Final CpG Biomarker Panel

The final biomarker panel consisted of 11 CpGs:

| CpG | Gene |
|---|---|
| cg10139015 | PSMC3 |
| cg11001216 | A1BG |
| cg04539775 | CLSTN3 |
| cg04985016 | TRIP12 |
| cg24745753 | CTNNA2 |
| cg15980170 | ZNF778 |
| cg09518293 | TTN |
| cg16311883 | ZNF664 |
| cg17967426 | PANK3 |
| cg09217522 | RECQL |
| cg00078968 | WNT3A |

These CpGs were subsequently evaluated in the external validation cohorts.

## Data Availability

The datasets required to reproduce the analyses will be made available separately. Because some of the methylation datasets
and intermediate analysis objects are too large for standard GitHub storage, large data files will be deposited in an appropriate
research data repository.

The corresponding repository links, file descriptions, and instructions for reproducing the analysis will be added to the
`data/README.md` file.

**Data repository/DOI:** To be added.

## Results

Large intermediate and final analysis outputs are not stored directly in this GitHub repository.
Relevant result files will be deposited separately,and links and file descriptions will be provided in `results/README.md`.

**Results repository/DOI:** To be added.

## Software

The analyses were conducted primarily in R. Major packages used across the workflow include packages for DNA methylation
preprocessing, EWAS, feature selection, machine learning, ROC analysis, and visualization.

Package-specific dependencies are loaded within the corresponding analysis scripts. Detailed software and package-version
information will be added to the repository to facilitate reproducibility.

## Reproducibility

The scripts should be executed sequentially where applicable, as outputs generated at earlier stages are used as inputs
to subsequent analyses.

The general workflow is:

```text
Raw DNA methylation data
        ↓
01_Preprocessing.R
        ↓
Processed methylation data
        ↓
02_missingness_imputation.R
        ↓
Imputed methylation dataset
        ↓
03_EWAS.R
        ↓
527 EWAS-selected CpGs
        ↓
04_feature_selection.R
        ↓
RF / LASSO / RF-RFE selected CpG subsets
        ↓
05_RF_panel_building.R
06_SVM_panel_building.R
        ↓
Final CpG biomarker panel
        ↓
External validation
        ├── 07_external_validation_GARBH_INi.R
        ├── 08_external_validation_PRiDE.R
        └── 09_external_validation_REVAMP.R
```

Users should refer to `R/README.md` for script-specific information and to the forthcoming `data/README.md` for the required input datasets.

## Citation

Citation information for the associated manuscript will be added following publication.

## Contact

For questions regarding the analysis or code, please contact the corresponding study authors.
│
├── docs/
│
└── README.md
```
