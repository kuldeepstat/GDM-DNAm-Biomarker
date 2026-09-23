# Results

## 1. Overview

This directory provides information about the main analysis outputs generated during the DNA methylation-based biomarker discovery and validation workflow for gestational diabetes mellitus (GDM) prediction.

The results are organized according to the major stages of the analytical workflow: epigenome-wide association analysis (EWAS), identification of significant CpGs, feature selection, Random Forest-based biomarker panel development, Support Vector Machine-based biomarker panel development, and external validation.

## 2. Epigenome-Wide Association Analysis (EWAS)

### `EWAS_gdm_first_tri.xlsx`

This file contains the complete results from the epigenome-wide association analysis performed using `R/03_EWAS.R`.

The EWAS results include the statistical results and associated CpG annotation generated during differential methylation analysis.

**Study-defined significance criterion:** Adjusted P-value < 1 × 10⁻⁶, with available gene annotation.

## 3. Significant CpGs Identified from EWAS

### `Sgn_Cpg.csv`

This file contains the list of 527 CpGs retained following EWAS filtering.

**Filtering criterion:** Adjusted P-value < 1 × 10⁻⁶.

These 527 CpGs formed the candidate feature set for the subsequent feature-selection analysis.

## 4. Feature-Selected CpGs

### `Selected_CpGs_Feature_Selection.xlsx`

This file consists of three individual sheets with the list of CpGs selected via RFFS, LASSO, and RF-RFE/backward selection approaches.

## 5. Random Forest Panel-Building Results

Random Forest-based biomarker panel development was performed using the CpG subsets obtained from the three feature-selection approaches. The corresponding analysis is implemented in `R/05_RF_panel_building.R`.

### 5.1 `ResultNewRF`

This file contains the Random Forest panel-building results obtained using the CpG subset selected through Random Forest Forward Selection (RFFS)-based feature selection.

### 5.2 `ResultNewLS`

This file contains the Random Forest panel-building results obtained using the CpG subset selected through LASSO.

### 5.3 `ResultNewBC`

This file contains the Random Forest panel-building results obtained using the CpG subset selected through Random Forest Recursive Feature Elimination (RF-RFE/backward selection).

## 6. Support Vector Machine Panel-Building Results

Support Vector Machine-based biomarker panel development was performed using the CpG subsets obtained from the same three feature-selection approaches. The corresponding analysis is implemented in `R/06_SVM_panel_building.R`.

### 6.1 `SVMRFF.xlsx`

This file contains the SVM panel-building results obtained using the CpG subset selected through Random Forest Forward Selection (RFFS)-based feature selection.

### 6.2 `SVMLASSO.xlsx`

This file contains the SVM panel-building results obtained using the CpG subset selected through LASSO.

### 6.3 `SVMBCEL.xlsx`

This file contains the SVM panel-building results obtained using the CpG subset selected through Random Forest Recursive Feature Elimination (RF-RFE/backward selection).

## 7. Summary of Result Files

| Analysis stage | Result file | Corresponding R script |
|---|---|---|
| EWAS | `EWAS_gdm_first_tri.xlsx` | `03_EWAS.R` |
| Significant CpGs | `Sgn_Cpg.csv` | `03_EWAS.R` |
| Feature-selected CpGs | `Selected_CpGs_Feature_Selection.xlsx` | `04_feature_selection.R` |
| RF classifier - RF-selected CpGs | `ResultNewRF` | `05_RF_panel_building.R` |
| RF classifier - LASSO-selected CpGs | `ResultNewLS` | `05_RF_panel_building.R` |
| RF classifier - RF-RFE-selected CpGs | `ResultNewBC` | `05_RF_panel_building.R` |
| SVM classifier - RF-selected CpGs | `SVMRFF.xlsx` | `06_SVM_panel_building.R` |
| SVM classifier - LASSO-selected CpGs | `SVMLASSO.xlsx` | `06_SVM_panel_building.R` |
| SVM classifier - RF-RFE-selected CpGs | `SVMBCEL.xlsx` | `06_SVM_panel_building.R` |

## 8. Relationship Between Data, Scripts, and Results

The major result-generation workflow is summarized below:

```text
02_FDATAV1.Rdata
        ↓
03_EWAS.R
        ↓
├── EWAS_gdm_first_tri.xlsx
└── Sgn_Cpg.csv
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
┌────────────────────────┬────────────────────────┐
↓                                                 ↓
05_RF_panel_building.R                   06_SVM_panel_building.R
↓                                                 ↓
RF result files                            SVM result files
```

The final CpG biomarker panel derived from the biomarker-development workflow was subsequently evaluated in the three external validation cohorts.

## 9. Results Availability

Descriptions of the result files generated at different stages of the analytical workflow are provided in this README.

The GitHub repository contains the R scripts required to generate the corresponding results from the datasets described in `data/README.md`.
