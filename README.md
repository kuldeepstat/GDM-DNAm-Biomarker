# DNA Methylation-Based Biomarker Discovery and Validation for Gestational Diabetes Mellitus Prediction

## Overview

This repository contains the R code used for DNA methylation-based biomarker discovery, predictive modelling, and external validation for gestational diabetes mellitus (GDM).

The analytical workflow includes DNA methylation data preprocessing and quality control, handling of missing methylation values, epigenome-wide association analysis (EWAS), feature selection, biomarker panel development using machine-learning classifiers, and evaluation of the selected CpG biomarker panel in three external cohorts.

The repository is organized to provide a transparent, reproducible record of the computational analyses performed in the study. Analysis scripts are numbered according to the order of the analytical workflow.

## Analysis Workflow

### 1. DNA methylation preprocessing and quality control

- Import of DNA methylation data
- Sample-level quality control
- Probe-level filtering
- Removal of probes on sex chromosomes
- Removal of SNP-associated and cross-reactive probes
- Generation of processed methylation beta values

### 2. Missingness filtering and imputation

- CpGs were filtered according to the proportion of observed methylation values.
- CpGs with less than or equal to 80% observed values were excluded.
- Missing beta values among the retained CpGs were imputed using `methyLImp`.
- The resulting dataset was used as the input for subsequent EWAS.

### 3. Epigenome-wide association analysis (EWAS)

- Quantile normalization was applied to the methylation data.
- Differential methylation analysis was performed using `ChAMP`.
- Multiple-testing correction was performed using the Benjamini-Hochberg procedure.
- CpGs meeting the study-defined significance threshold (adjusted P < 1 × 10^-6) and having gene annotation were retained.
- This resulted in 527 CpGs for downstream feature-selection analyses.

### 4. Feature selection

- Three feature-selection approaches were applied to the 527 EWAS-selected CpGs: Random Forest Forward Selection (RFFS)-based feature ranking, LASSO regression, and Random Forest Recursive Feature Elimination (RF-RFE/backward selection).
- The selected CpG subsets (40 CpGs by RFFS, 40 by RFE, and 21 by LASSO from the 527 CpGs) were subsequently used for biomarker panel development.

### 5. Random Forest biomarker panel development

- Random Forest classification was applied to CpG subsets obtained from the three feature-selection approaches.
- The analysis included repeated training/test partitioning, forward addition of candidate CpGs, cross-validation, model evaluation, and resampling-based performance assessment.

### 6. Support Vector Machine biomarker panel development

- Support Vector Machine (SVM) classification with a radial kernel was applied to the CpG subsets obtained from the three feature-selection approaches.
- Candidate CpGs were evaluated through an iterative panel-building procedure.
- Model performance was assessed using multiple discrimination and classification metrics.

### 7. External validation

- The final 11-CpG biomarker panel was evaluated in three external cohorts: GARBH-INi, PRiDE, and REVAMP.
- Cohort-specific logistic regression and ROC analyses were performed to evaluate the selected CpG panel.

## Repository Structure

```text
GDM-DNAm-Biomarker/
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
├── data/
│   └── README.md
├── results/
│   └── README.md
├── docs/
└── README.md
```

## R Scripts

| Script | Analysis |
|---|---|
| `01_Preprocessing.R` | DNA methylation preprocessing and quality control |
| `02_missingness_imputation.R` | Missingness filtering and methylation-value imputation |
| `03_EWAS.R` | EWAS and identification of significant CpGs |
| `04_feature_selection.R` | Feature selection using RFFS, LASSO, and RF-RFE/backward selection |
| `05_RF_panel_building.R` | Random Forest-based biomarker panel development and performance evaluation |
| `06_SVM_panel_building.R` | SVM-based biomarker panel development and performance evaluation |
| `07_external_validation_GARBH_INi.R` | External validation in the GARBH-INi cohort |
| `08_external_validation_PRiDE.R` | External validation in the PRiDE cohort |
| `09_external_validation_REVAMP.R` | External validation in the REVAMP cohort |

`R/README.md` provides additional information about each script and its position in the analytical workflow.

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

The datasets required to reproduce the analyses are available separately.

Descriptions of the datasets used at different stages of the analytical workflow are provided in `data/README.md`.

## Results

Large intermediate and final analysis outputs are not stored directly in this GitHub repository.

Descriptions of the relevant result files are provided in `results/README.md`.

## Software and Dependencies

All analyses were conducted in R. The workflow uses R packages for DNA methylation preprocessing and quality control, missing-value imputation, epigenome-wide association analysis (EWAS), feature selection, machine-learning model development, performance evaluation, ROC analysis, and visualization.

The required packages are specified within the corresponding analysis scripts. To facilitate computational reproducibility, the R version, package versions, and session information used for the analyses are provided in this repository.

## Reproducibility

The scripts should be executed sequentially where applicable, as outputs generated at earlier stages are used as inputs to subsequent analyses.

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
RFFS / LASSO / RF-RFE selected CpG subsets
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

Users should refer to `R/README.md` for script-specific information and to `data/README.md` for descriptions of the required input datasets.

## Contact

For questions regarding the analysis or code, please contact the corresponding author.
