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
│
├── docs/
│
└── README.md
```
