# R Analysis Scripts

## Overview

This directory contains the R scripts used for DNA methylation preprocessing, missing-data handling, epigenome-wide association analysis (EWAS), feature selection, biomarker panel development, and external validation.

The scripts are numbered according to the approximate order of the analytical workflow. The output generated at one stage is used as the input for subsequent stages where applicable.

## Analysis Sequence

```text
01_Preprocessing.R
        ↓
02_missingness_imputation.R
        ↓
03_EWAS.R
        ↓
04_feature_selection.R
        ↓
05_RF_panel_building.R    06_SVM_panel_building.R
        ↓                       ↓
        └──── Final CpG biomarker panel ────┘
                         ↓
        07 GARBH-INi | 08 PRiDE | 09 REVAMP
```

## 01_Preprocessing.R

### Purpose

Performs preprocessing and quality control of the DNA methylation data.

### Main analytical steps

- Import of methylation array data
- Calculation and assessment of detection P-values
- Sample-level quality control
- Normalization of methylation data
- Probe-level quality-control filtering
- Removal of probes located on sex chromosomes
- Removal of SNP-associated probes
- Removal of cross-reactive probes
- Extraction of methylation beta and M-values after preprocessing

The processed beta-value dataset generated from this stage provides the basis for the subsequent missingness and imputation analysis.

**Next step:** `02_missingness_imputation.R`

## 02_missingness_imputation.R

### Purpose

Performs CpG-level missingness filtering followed by imputation of missing methylation beta values.

### Main analytical steps

- Processed DNA methylation data generated during the preprocessing stage.
- CpGs are retained when more than 80% of their observations are non-missing: `colMeans(!is.na(DataProcessed)) > 0.8`.
- Thus, CpGs with 20% or greater missingness are excluded by this implementation.
- Missing beta values among the retained CpGs are imputed using the `methyLImp` package.
- Because of the high dimensionality of the methylation dataset, imputation is performed in blocks of CpGs before the blocks are combined into the final dataset.
- The binary outcome variable `Gr` is subsequently added to the final dataset.

**Main output:** `FDATAV1`, containing the imputed CpG methylation values together with the outcome variable.

**Next step:** `03_EWAS.R`

## 03_EWAS.R

### Purpose

Performs the epigenome-wide association analysis to identify CpGs associated with GDM.

### Main analytical steps

- Input: `FDATAV1`
- Extraction of the methylation beta-value matrix
- Quantile normalization
- Use of the binary GDM group variable as phenotype
- Differential methylation analysis using `ChAMP`
- Benjamini-Hochberg multiple-testing correction
- Retention of CpGs with adjusted P-value < 1 × 10^-6 and available gene annotation.

A total of 527 CpGs meeting the filtering criteria are retained for downstream feature selection.

**Next step:** `04_feature_selection.R`

## 04_feature_selection.R

### Purpose

Applies multiple feature-selection approaches to the 527 CpGs identified through EWAS.

### Main analytical steps

- Random Forest Forward Selection (RFFS)-based feature ranking
- LASSO regression
- Random Forest Recursive Feature Elimination (RF-RFE/backward selection)
- The workflow produces 40 CpGs from RFFS, 21 CpGs from LASSO, and 40 CpGs from RF-RFE/backward selection.
- The selected datasets are stored in the `DataRF`, `DataLS`, and `DataBC` worksheets.

These feature-selected datasets are subsequently evaluated using Random Forest and Support Vector Machine classifiers.

**Next step:** `05_RF_panel_building.R` and `06_SVM_panel_building.R`

## 05_RF_panel_building.R

### Purpose

Develops and evaluates CpG biomarker panels using Random Forest classification.

### Main analytical steps

- Inputs: `DataRF`, `DataLS`, and `DataBC`
- Repeated analysis iterations
- Stratified 67:33 training/test partitioning
- Forward addition of candidate CpG predictors
- Repeated 10-fold cross-validation with 5 repeats
- Random Forest model fitting
- Evaluation of candidate CpG panels
- Selection of the best-performing panel within each iteration
- Resampling-based assessment of model performance
- Performance measures include AUC, sensitivity, specificity, F-score, accuracy, Brier score, logarithmic score, AUPRC, divergence, and Kolmogorov-Smirnov statistic.

Results are generated separately for CpG subsets originating from Random Forest, LASSO, and RF-RFE/backward feature selection.

## 06_SVM_panel_building.R

### Purpose

Develops and evaluates CpG biomarker panels using Support Vector Machine classification.

### Main analytical steps

- Inputs: `DataRF`, `DataLS`, and `DataBC`
- 10 analysis iterations
- Stratified 67:33 training/test partitioning
- Forward addition of candidate CpGs
- Repeated 10-fold cross-validation with 5 repeats
- Radial-kernel Support Vector Machine modelling
- Evaluation of candidate panels
- Identification of the best-performing panel within each iteration
- Repeated resampling for assessment of model performance
- Performance measures include AUC, sensitivity, specificity, F-score, accuracy, Brier score, logarithmic score, AUPRC, divergence, and Kolmogorov-Smirnov statistic.

The results from the Random Forest and SVM analyses were used in the biomarker panel-development workflow.

## External Validation

The final 11-CpG biomarker panel was subsequently evaluated in three independent cohorts. The external-validation analyses are provided as separate scripts.

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

## 07_external_validation_GARBH_INi.R

Evaluates the final 11-CpG biomarker panel in the GARBH-INi external validation cohort.

**Outcome variable:** `GroupCode`

**Model 1a:** 11 CpGs.

**Model 1b:** 11 CpGs + Age + BMI + family history of diabetes.

Predicted probabilities are obtained from the fitted logistic regression models. ROC analysis is performed using `pROC`, including AUC, 95% confidence interval for AUC, and sensitivity and specificity at a probability threshold of 0.5.

## 08_external_validation_PRiDE.R

Evaluates the final 11-CpG biomarker panel in the PRiDE external validation cohort.

**Outcome variable:** `GDMgroup_Code`

**Model 1a:** 11 CpGs.

**Model 1b:** 11 CpGs + age + BMI + family history of diabetes.

Predicted probabilities are obtained from the fitted logistic regression models. ROC analysis is performed using `pROC`, including AUC, 95% confidence interval for AUC, and sensitivity and specificity at a probability threshold of 0.5.

## 09_external_validation_REVAMP.R

Evaluates the final 11-CpG biomarker panel in the REVAMP external validation cohort.

**Outcome variable:** `GroupCode`

**Model 1a:** 11 CpGs.

**Model 1b:** 11 CpGs + Age + BMI + family history of diabetes.

Predicted probabilities are obtained from the fitted logistic regression models. ROC analysis is performed using `pROC`, including AUC, 95% confidence interval for AUC, and sensitivity and specificity at a probability threshold of 0.5.

## Data and Large Analysis Objects

Some input datasets and intermediate analysis objects are too large to store directly in this GitHub repository. The corresponding file names and instructions are documented separately in `data/README.md`. Large analysis outputs are similarly documented in `results/README.md`.

## Reproducibility Notes

The scripts retain the analytical workflow used for the study. Repository-relative file paths are used where possible to facilitate execution on systems other than the original analysis computer.

Users reproducing the complete workflow should execute the scripts in numerical order and ensure that the required input files are available before running each stage.

All analyses were conducted in R. The required packages are specified within the corresponding analysis scripts. To facilitate computational reproducibility, the R version, package versions, and session information used for the analyses are provided in the repository.
