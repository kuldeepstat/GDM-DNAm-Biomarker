# ==============================================================================
# 03_EWAS.R
#
# Purpose:
# Perform an epigenome-wide association study (EWAS) to identify CpG sites
# associated with gestational diabetes mellitus (GDM).
#
# Input:
# FDATAV1.Rdata
#   - 258 samples
#   - 597,670 CpG beta-value variables
#   - Binary outcome variable Gr (0 = Control, 1 = GDM)
#
# Procedure:
#   1. Separate CpG beta values and GDM outcome
#   2. Apply quantile normalization
#   3. Perform EWAS using ChAMP::champ.DMP
#   4. Apply Benjamini-Hochberg multiple-testing adjustment
#   5. Retain CpGs with adjusted P-value < 1e-6 and gene annotation
#   6. Construct the 527-CpG dataset for downstream feature selection
#
# Output:
# EWASFINAL.Rdata
#   - Complete EWAS results
#
# d527.Rdata
#   - 527 EWAS-selected CpGs
#   - Binary outcome variable Gr
#
# ==============================================================================


# ==============================================================================
# 1. Load required packages
# ==============================================================================

library(ChAMP)
library(dplyr)
library(preprocessCore)


# ==============================================================================
# 2. Load imputed methylation dataset
# ==============================================================================

load("data/FDATAV1.Rdata")


# ==============================================================================
# 3. Separate methylation data and outcome
# ==============================================================================

# Extract CpG beta values
data <- FDATAV1[, colnames(FDATAV1) != "Gr"]

# Extract GDM status
Group <- FDATAV1$Gr

# Check group distribution
table(Group)


# ==============================================================================
# 4. Quantile normalization
# ==============================================================================

nd <- normalize.quantiles(as.matrix(data))

nd <- as.data.frame(nd)

colnames(nd) <- colnames(data)
rownames(nd) <- rownames(data)


# ==============================================================================
# 5. Epigenome-wide association analysis
# ==============================================================================

DMP <- champ.DMP(
  beta = t(nd),
  pheno = Group,
  compare.group = NULL,
  adjPVal = 1,
  adjust.method = "BH",
  arraytype = "EPIC"
)


# ==============================================================================
# 6. Extract EWAS results
# ==============================================================================

scpg <- DMP[[1]]

# Add CpG identifier as an explicit variable
scpg$id <- rownames(scpg)


# ==============================================================================
# 7. Filter significant CpGs
# ==============================================================================

# Retain CpGs with:
#   BH-adjusted P-value < 1e-6
#   non-missing and non-empty gene annotation

scpg_filtered <- filter(
  scpg,
  adj.P.Val < 1e-6 &
    !is.na(gene) &
    gene != ""
)

# Number of EWAS-selected CpGs
nrow(scpg_filtered)


# ==============================================================================
# 8. Save complete EWAS results
# ==============================================================================

save(
  scpg,
  file = "results/EWASFINAL.Rdata"
)


# ==============================================================================
# 9. Construct dataset containing the 527 EWAS-selected CpGs
# ==============================================================================

# Verify that all selected CpGs are present in FDATAV1
all(rownames(scpg_filtered) %in% colnames(FDATAV1))

# Extract the selected CpGs
d527 <- FDATAV1[
  ,
  rownames(scpg_filtered),
  drop = FALSE
]

# Add GDM outcome
d527$Gr <- FDATAV1$Gr

# Check dimensions
dim(d527)

# Check outcome distribution
table(d527$Gr)


# ==============================================================================
# 10. Save dataset for downstream feature selection
# ==============================================================================

save(
  d527,
  file = "data/d527.Rdata"
)
