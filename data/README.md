# Data

## 1. Overview

This directory provides information about the datasets used in the DNA methylation-based biomarker discovery and validation workflow for gestational diabetes mellitus (GDM) prediction.

The reproducible analytical workflow begins with the preprocessed DNA methylation dataset (`01_DataProcessed.Rdata`). This repository does not include raw IDAT files. The preprocessing and quality-control procedures applied to the raw methylation data are documented separately in `R/01_Preprocessing.R`.

## 2. Discovery and Biomarker Development Datasets

### 2.1 `01_DataProcessed.Rdata`

This is the preprocessed DNA methylation dataset and represents the starting point for the reproducible analytical workflow.

**Used as input for:** `R/02_missingness_imputation.R`

The script performs CpG-level missingness filtering followed by imputation of missing methylation beta values.

**Next dataset generated:** `02_FDATAV1.Rdata`

### 2.2 `02_FDATAV1.Rdata`

This dataset is generated following missingness filtering and methylation-value imputation. It contains imputed CpG methylation data and the binary GDM outcome variable (`Gr`).

**Used as input for:** `R/03_EWAS.R`

The EWAS script performs the epigenome-wide association analysis and identifies CpGs meeting the study-defined significance criteria.

**Next dataset generated:** `03_d527.Rdata`

### 2.3 `03_d527.Rdata`

This dataset contains the 527 CpGs retained following EWAS filtering, together with the GDM outcome variable.

**EWAS filtering criterion:** Adjusted P-value < 1 × 10⁻⁶, with available gene annotation.

**Used as input for:** `R/04_feature_selection.R`

Three feature-selection approaches are subsequently applied: Random Forest Forward Selection (RFFS)-based feature ranking, LASSO regression, and Random Forest Recursive Feature Elimination (RF-RFE/backward selection).

**Next dataset generated:** `04_Data_FeatureSelection.xlsx`

### 2.4 `04_Data_FeatureSelection.xlsx`

This workbook contains the datasets generated from the three feature-selection approaches.

| Worksheet | Feature-selection approach | Number of selected CpGs |
|---|---|---:|
| `DataRF` | RFFS | 40 |
| `DataLS` | LASSO | 21 |
| `DataBC` | RF-RFE/backward selection | 40 |

These datasets are used as inputs for both classifier-based biomarker panel-development scripts:

- `R/05_RF_panel_building.R`
- `R/06_SVM_panel_building.R`

## 3. Data Flow

The discovery datasets and analysis scripts are summarized below:

```text
01_DataProcessed.Rdata
        ↓
02_missingness_imputation.R
        ↓
02_FDATAV1.Rdata
        ↓
03_EWAS.R
        ↓
03_d527.Rdata
        ↓
04_feature_selection.R
        ↓
04_Data_FeatureSelection.xlsx
   ├── DataRF
   ├── DataLS
   └── DataBC
        ↓
05_RF_panel_building.R and 06_SVM_panel_building.R
```

## 4. External Validation Datasets

The following datasets are used for external validation of the final CpG biomarker panel.

### 4.1 GARBH-INi 11CpG and clinical Data

This file contains methylation values for 11 CpGs, maternal risk factors, and group codes for ROC analysis.

### 4.2 REVAMP_11CpGs_Clinical data

This file contains methylation values for 11 CpGs, maternal risk factors, and group codes for ROC analysis.

### 4.3 PRiDE_GDM_Methylation validation_Clinical_DATA

This file contains methylation values for 11 CpGs, maternal risk factors, and group codes for ROC analysis.

## 5. External Validation Scripts

| Dataset | Analysis script |
|---|---|
| `GARBH_INi_11CpG_clinical_data.xlsx` | `R/07_external_validation_GARBH_INi.R` |
| `PRiDE_11CpG_clinical_data.xlsx` | `R/08_external_validation_PRiDE.R` |
| `REVAMP_11CpG_clinical_data.xlsx` | `R/09_external_validation_REVAMP.R` |
