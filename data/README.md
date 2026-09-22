# Data

## Overview

This directory provides information about the datasets used in the DNA methylation-based biomarker 
discovery and validation workflow for gestational diabetes mellitus (GDM).

The reproducible analytical workflow begins with the **preprocessed DNA methylation dataset (`01_DataProcessed.Rdata`)**. 
Raw IDAT files are not distributed through this repository. 
The preprocessing and quality-control procedures applied to the raw methylation data are documented separately in `R/01_Preprocessing.R`.

Because several datasets are too large for standard GitHub storage, the corresponding data files will be deposited separately. 
Download links will be provided below.

---

## Discovery and Biomarker Development Datasets

### 1. `01_DataProcessed.Rdata`

**Download:** [Link to be added]

This is the preprocessed DNA methylation dataset and represents the starting point for the reproducible analytical 
workflow provided in this repository.

It is used as the input for:

`R/02_missingness_imputation.R`

The script performs CpG-level missingness filtering followed by imputation of missing methylation beta values.

**Next dataset generated:** `02_FDATAV1.Rdata`

---

### 2. `02_FDATAV1.Rdata`

**Download:** [Link to be added]

This dataset is generated following missingness filtering and methylation-value imputation.

It contains the imputed CpG methylation data together with the binary GDM outcome variable (`Gr`).

It is used as the input for:

`R/03_EWAS.R`

The EWAS script performs the epigenome-wide association analysis and identifies CpGs meeting the study-defined 
significance criteria.

**Next dataset generated:** `03_d527.Rdata`

---

### 3. `03_d527.Rdata`

**Download:** [Link to be added]

This dataset contains the **527 CpGs retained following EWAS filtering**, together with the GDM outcome variable.

The EWAS filtering criterion used in the analysis was:

**Adjusted P-value < 1 × 10⁻⁶**

with available gene annotation.

This dataset is used as the input for:

`R/04_feature_selection.R`

Three feature-selection approaches are subsequently applied:

- Random Forest-based feature ranking
- LASSO regression
- Random Forest Recursive Feature Elimination (RF-RFE/backward selection)

**Next dataset generated:** `04_Data.xlsx`

---

### 4. `04_Data.xlsx`

**Download:** [Link to be added]

This workbook contains the datasets generated from the three feature-selection approaches.

The workbook contains three worksheets:

| Worksheet | Feature-selection approach | Number of selected CpGs |
|---|---|---:|
| `DataRF` | Random Forest-based feature selection | 40 |
| `DataLS` | LASSO | 21 |
| `DataBC` | RF-RFE/backward selection | 40 |

These datasets are used as inputs for both classifier-based biomarker panel-development scripts:

- `R/05_RF_panel_building.R`
- `R/06_SVM_panel_building.R`

---

## Data Flow

The relationship between the discovery datasets and analysis scripts is:

```text
01_DataProcessed.Rdata
        │
        ▼
02_missingness_imputation.R
        │
        ▼
02_FDATAV1.Rdata
        │
        ▼
03_EWAS.R
        │
        ▼
03_d527.Rdata
        │
        ▼
04_feature_selection.R
        │
        ▼
04_Data.xlsx
   ├── DataRF
   ├── DataLS
   └── DataBC
        │
        ├───────────────┐
        ▼               ▼
05_RF_panel_building.R  06_SVM_panel_building.R
```

---

# External Validation Datasets

The following datasets are used for external validation of the final CpG biomarker panel.

### GARBH-INi

**Dataset:** `GARBH_INi_11CpG_clinical_data.xlsx`

**Download:** [Link to be added]

---

### PRiDE

**Dataset:** `PRiDE_11CpG_clinical_data.xlsx`

**Download:** [Link to be added]

---

### REVAMP

**Dataset:** `REVAMP_11CpG_clinical_data.xlsx`

**Download:** [Link to be added]

---

## External Validation Scripts

The external validation datasets correspond to the following analysis scripts:

| Dataset | Analysis script |
|---|---|
| GARBH-INi | `R/07_external_validation_GARBH_INi.R` |
| PRiDE | `R/08_external_validation_PRiDE.R` |
| REVAMP | `R/09_external_validation_REVAMP.R` |

Detailed descriptions of the external validation datasets will be added once the final data files 
and associated data-sharing information are finalized.

---

## Data Availability

Large datasets required to reproduce the analyses are not stored directly in the GitHub repository. 
Links to the corresponding deposited datasets will be added to this README following finalization of the data repository.

Users should download the required dataset and place it in the appropriate local `data/` directory before executing the corresponding R script.
