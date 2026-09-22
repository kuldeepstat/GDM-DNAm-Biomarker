# Results

## Overview

This directory provides information about the main analysis outputs generated during the DNA methylation-based biomarker 
discovery and validation workflow for gestational diabetes mellitus (GDM).

Because some result files are relatively large or are maintained separately from the source-code repository, the corresponding 
files will be deposited separately. Download links will be provided below.

The results are organized according to the major stages of the analytical workflow: epigenome-wide association analysis (EWAS), 
identification of significant CpGs, Random Forest-based biomarker panel development, Support Vector Machine-based biomarker 
panel development, and external validation.

---

## 1. Epigenome-Wide Association Analysis (EWAS)

### `EWASFINAL.Rdata`

**Download:** [Link to be added]

This file contains the complete results from the epigenome-wide association analysis performed using `R/03_EWAS.R`.

The EWAS results include the statistical results and associated CpG annotation generated during differential methylation analysis.

The study-defined significance criterion used for downstream analysis was:

**Adjusted P-value < 1 × 10⁻⁶**

with available gene annotation.

---

## 2. Significant CpGs Identified from EWAS

### `Sgn_Cpg.csv`

**Download:** [Link to be added]

This file contains the list of **527 CpGs** retained following EWAS filtering.

The CpGs included in this file satisfy the study-defined significance criterion:

**Adjusted P-value < 1 × 10⁻⁶**

These 527 CpGs formed the candidate feature set for the subsequent feature-selection analysis.

---

## 3. Random Forest Panel-Building Results

Random Forest-based biomarker panel development was performed using the CpG subsets obtained from the 
three feature-selection approaches.

The corresponding analysis is implemented in:

`R/05_RF_panel_building.R`

Separate result files are available for each feature-selection approach.

### 3.1 `ResultRF_RF up.xlsx`

**Download:** [Link to be added]

This file contains the Random Forest panel-building results obtained using the CpG subset selected through 
Random Forest-based feature selection.

---

### 3.2 `ResultRF_LS up.xlsx`

**Download:** [Link to be added]

This file contains the Random Forest panel-building results obtained using the CpG subset selected through LASSO.

---

### 3.3 `ResultRF_BC up.xlsx`

**Download:** [Link to be added]

This file contains the Random Forest panel-building results obtained using the CpG subset selected through Random Forest 
Recursive Feature Elimination (RF-RFE/backward selection).

---

## 4. Support Vector Machine Panel-Building Results

Support Vector Machine-based biomarker panel development was performed using the CpG subsets obtained from the 
same three feature-selection approaches.

The corresponding analysis is implemented in:

`R/06_SVM_panel_building.R`

Separate result files are available for each feature-selection approach.

### 4.1 `SVMRFF.xlsx`

**Download:** [Link to be added]

This file contains the SVM panel-building results obtained using the CpG subset selected through Random 
Forest-based feature selection.

---

### 4.2 `SVMLASSO.xlsx`

**Download:** [Link to be added]

This file contains the SVM panel-building results obtained using the CpG subset selected through LASSO.

---

### 4.3 `SVMBCEL.xlsx`

**Download:** [Link to be added]

This file contains the SVM panel-building results obtained using the CpG subset selected through Random Forest 
Recursive Feature Elimination (RF-RFE/backward selection).

---

## 5. External Validation Results

The final CpG biomarker panel was evaluated in three external validation cohorts.

The corresponding validation scripts are:

- `R/07_external_validation_GARBH_INi.R`
- `R/08_external_validation_PRiDE.R`
- `R/09_external_validation_REVAMP.R`



### 5.1 `GARBH_INi_validation_results.xlsx`

**Download:** [Link to be added]

---

### 5.2 `PRiDE_validation_results.xlsx`

**Download:** [Link to be added]

---

### 5.3 `REVAMP_validation_results.xlsx`

**Download:** [Link to be added]

---

## 6. Summary of Result Files

| Analysis stage | Result file | Corresponding R script |
|---|---|---|
| EWAS | `EWASFINAL.Rdata` | `03_EWAS.R` |
| Significant CpGs | `Sgn_Cpg.csv` | `03_EWAS.R` |
| RF classifier – RF-selected CpGs | `ResultRF_RF up.xlsx` | `05_RF_panel_building.R` |
| RF classifier – LASSO-selected CpGs | `ResultRF_LS up.xlsx` | `05_RF_panel_building.R` |
| RF classifier – RF-RFE-selected CpGs | `ResultRF_BC up.xlsx` | `05_RF_panel_building.R` |
| SVM classifier – RF-selected CpGs | `SVMRFF.xlsx` | `06_SVM_panel_building.R` |
| SVM classifier – LASSO-selected CpGs | `SVMLASSO.xlsx` | `06_SVM_panel_building.R` |
| SVM classifier – RF-RFE-selected CpGs | `SVMBCEL.xlsx` | `06_SVM_panel_building.R` |
| GARBH-INi external validation | `GARBH_INi_validation_results.xlsx` | `07_external_validation_GARBH_INi.R` |
| PRiDE external validation | `PRiDE_validation_results.xlsx` | `08_external_validation_PRiDE.R` |
| REVAMP external validation | `REVAMP_validation_results.xlsx` | `09_external_validation_REVAMP.R` |

---

## 7. Relationship Between Data, Scripts, and Results

The major result-generation workflow is:

```text
02_FDATAV1.Rdata
        │
        ▼
03_EWAS.R
        │
        ├── EWASFINAL.Rdata
        │
        └── Sgn_Cpg.csv
                 │
                 ▼
          03_d527.Rdata
                 │
                 ▼
       04_feature_selection.R
                 │
                 ▼
           04_Data.xlsx
       ┌─────────┴─────────┐
       ▼                   ▼
05_RF_panel_building.R   06_SVM_panel_building.R
       │                   │
       ▼                   ▼
 RF result files       SVM result files
```

The final CpG biomarker panel derived from the biomarker-development workflow was subsequently evaluated in the three external validation cohorts.

---

## 8. Results Availability

Large or supplementary analysis outputs are not stored directly in the GitHub repository. The corresponding result files will be deposited separately, and permanent download links will be added to this README.

The GitHub repository contains the R scripts required to generate the corresponding results from the datasets described in `data/README.md`.
