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
 ┌──────────────────────────────┐
 │                              │
05_RF_panel_building.R    06_SVM_panel_building.R
 │                              │
 └──────────────┬───────────────┘
                ↓
       Final CpG biomarker panel
                ↓
 ┌──────────────┼──────────────┐
 ↓              ↓              ↓
07             08             09
GARBH-INi      PRiDE          REVAMP
```

---

## 01_Preprocessing.R

### Purpose

Performs preprocessing and quality control of the DNA methylation data.

### Main analytical steps

The script includes:

- import of methylation array data;
- calculation and assessment of detection P-values;
- sample-level quality control;
- normalization of methylation data;
- probe-level quality-control filtering;
- removal of probes located on sex chromosomes;
- removal of SNP-associated probes;
- removal of cross-reactive probes; and
- extraction of methylation beta and M-values after preprocessing.

The processed beta-value dataset generated from this stage provides the basis for the subsequent missingness and imputation analysis.

### Next step

The processed methylation data are used by:

`02_missingness_imputation.R`

---

## 02_missingness_imputation.R

### Purpose

Performs CpG-level missingness filtering followed by imputation of missing methylation beta values.

### Input

Processed DNA methylation data generated during the preprocessing stage.

### Main analytical steps

CpGs are retained when more than 80% of their observations are non-missing:

```r
colMeans(!is.na(DataProcessed)) > 0.8
```

Thus, CpGs with 20% or greater missingness are excluded by this implementation.

Missing beta values among the retained CpGs are imputed using the `methyLImp` package.

Because of the high dimensionality of the methylation dataset, imputation is performed in blocks of CpGs before the 
blocks are combined into the final dataset.

The binary outcome variable `Gr` is subsequently added to the final dataset.

### Main output

`FDATAV1`

The final dataset contains the imputed CpG methylation values together with the outcome variable.

### Next step

`FDATAV1` is used as the input for:

`03_EWAS.R`

---

## 03_EWAS.R

### Purpose

Performs the epigenome-wide association analysis to identify CpGs associated with GDM.

### Input

`FDATAV1`

### Main analytical steps

The script:

- extracts the methylation beta-value matrix;
- performs quantile normalization;
- uses the binary GDM group variable as the phenotype;
- performs differential methylation analysis using `ChAMP`;
- applies Benjamini-Hochberg multiple-testing correction; and
- filters the EWAS results according to the study-defined significance criterion.

CpGs are retained when:

```text
Adjusted P-value < 1 × 10^-6
```

and a gene annotation is available.

### Main output

A total of **527 CpGs** meeting the filtering criteria are retained for downstream feature selection.

A dataset containing these CpGs together with the outcome variable is generated for the next stage.

### Next step

The 527-CpG dataset is used by:

`04_feature_selection.R`

---

## 04_feature_selection.R

### Purpose

Applies multiple feature-selection approaches to the 527 CpGs identified through EWAS.

### Input

Dataset containing the 527 EWAS-selected CpGs and the binary GDM outcome.

### Feature-selection approaches

Three approaches are implemented:

1. Random Forest-based feature ranking
2. LASSO regression
3. Random Forest Recursive Feature Elimination (RF-RFE/backward selection)

The resulting feature subsets are retained separately for subsequent classifier-based biomarker panel development.

### Selected feature sets

The workflow produces:

- 40 CpGs from Random Forest-based selection;
- 21 CpGs from LASSO; and
- 40 CpGs from RF-RFE/backward selection.

### Main output

The selected datasets are stored in separate worksheets:

```text
DataRF
DataLS
DataBC
```

These feature-selected datasets are subsequently evaluated using Random Forest and Support Vector Machine classifiers.

### Next steps

`05_RF_panel_building.R`

and

`06_SVM_panel_building.R`

---

## 05_RF_panel_building.R

### Purpose

Develops and evaluates CpG biomarker panels using Random Forest classification.

### Input

The script uses the three feature-selected datasets generated during the feature-selection stage:

```text
DataRF
DataLS
DataBC
```

### Main analytical workflow

Random Forest panel development is performed separately for each of the three feature-selection approaches.

The workflow includes:

- repeated analysis iterations;
- stratified 67:33 training/test partitioning;
- forward addition of candidate CpG predictors;
- repeated cross-validation within the training data;
- Random Forest model fitting;
- evaluation of candidate CpG panels;
- selection of the best-performing panel within each iteration; and
- resampling-based assessment of model performance.

Repeated cross-validation is performed using 10 folds and 5 repeats.

### Performance measures

The analysis evaluates multiple model-performance measures, including:

- Area Under the ROC Curve (AUC)
- Sensitivity
- Specificity
- F-score
- Accuracy
- Brier score
- Logarithmic score
- Area Under the Precision-Recall Curve (AUPRC)
- Divergence
- Kolmogorov-Smirnov statistic

### Output

Results are generated separately for CpG subsets originating from:

- Random Forest feature selection;
- LASSO feature selection; and
- RF-RFE/backward selection.

---

## 06_SVM_panel_building.R

### Purpose

Develops and evaluates CpG biomarker panels using Support Vector Machine classification.

### Input

The same three feature-selected datasets are evaluated:

```text
DataRF
DataLS
DataBC
```

### Main analytical workflow

The SVM analysis includes:

- 10 analysis iterations;
- stratified 67:33 training/test partitioning;
- forward addition of candidate CpGs;
- repeated cross-validation;
- radial-kernel Support Vector Machine modelling;
- evaluation of candidate panels;
- identification of the best-performing panel within each iteration; and
- repeated resampling for assessment of model performance.

Repeated cross-validation is performed using 10 folds and 5 repeats.

### Performance measures

Performance measures include:

- AUC
- Sensitivity
- Specificity
- F-score
- Accuracy
- Brier score
- Logarithmic score
- AUPRC
- Divergence
- Kolmogorov-Smirnov statistic

The results from the Random Forest and SVM analyses were used in the biomarker panel-development workflow.

---

# External Validation

The final 11-CpG biomarker panel was subsequently evaluated in three independent cohorts.

The external-validation analyses are provided as separate scripts because variable names and clinical covariates differ between the individual cohort datasets.

The final panel contains:

```text
cg10139015    PSMC3
cg11001216    A1BG
cg04539775    CLSTN3
cg04985016    TRIP12
cg24745753    CTNNA2
cg15980170    ZNF778
cg09518293    TTN
cg16311883    ZNF664
cg17967426    PANK3
cg09217522    RECQL
cg00078968    WNT3A
```

---

## 07_external_validation_GARBH_INi.R

### Purpose

Evaluates the final 11-CpG biomarker panel in the GARBH-INi external validation cohort.

### Models

Two logistic regression models are evaluated:

**Model 1a**

```text
11 CpGs
```

**Model 1b**

```text
11 CpGs + Age + BMI + family history of diabetes
```

### Evaluation

Predicted probabilities are obtained from the fitted models.

ROC analysis is performed using the `pROC` package, including:

- AUC;
- 95% confidence interval for AUC; and
- sensitivity and specificity at a probability threshold of 0.5.

ROC curves are generated for the CpG-only and CpG-plus-clinical models.

---

## 08_external_validation_PRiDE.R

### Purpose

Evaluates the final 11-CpG biomarker panel in the PRiDE external validation cohort.

### Outcome

The PRiDE-specific GDM group variable used in the analysis is:

```text
GDMgroup_Code
```

### Models

Two logistic regression models are evaluated:

**Model 1a**

```text
11 CpGs
```

**Model 1b**

```text
11 CpGs + age + BMI + family history of diabetes
```

The cohort-specific variable names present in the PRiDE dataset are retained in the analysis script.

### Evaluation

The script calculates:

- predicted probabilities;
- ROC curves;
- AUC;
- 95% confidence intervals for AUC; and
- sensitivity and specificity at a probability threshold of 0.5.

---

## 09_external_validation_REVAMP.R

### Purpose

Evaluates the final 11-CpG biomarker panel in the REVAMP external validation cohort.

### Models

Two logistic regression models are evaluated:

**Model 1a**

```text
11 CpGs
```

**Model 1b**

```text
11 CpGs + Age + BMI + family history of diabetes
```

### Evaluation

The external-validation analysis includes:

- logistic regression;
- calculation of predicted probabilities;
- ROC analysis;
- AUC;
- 95% confidence intervals for AUC; and
- sensitivity and specificity at a probability threshold of 0.5.

ROC curves are generated for both models.

---

## Data and Large Analysis Objects

Some input datasets and intermediate analysis objects are too large to store directly in this GitHub repository.

The corresponding data repository, file names, download links, and instructions for placing the files into the required directory structure will be documented separately in:

```text
data/README.md
```

Large analysis outputs will similarly be documented in:

```text
results/README.md
```

These files will be updated once the corresponding datasets and results have been deposited.

---

## Reproducibility Notes

The scripts retain the analytical workflow used for the study. Repository-relative file paths are used where possible to facilitate execution on systems other than the original analysis computer.

Users reproducing the complete workflow should execute the scripts in numerical order and ensure that the required input files are available before running each stage.

Software and package-version information will be provided separately to further support computational reproducibility.
