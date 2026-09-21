# ==============================================================================
# 04_feature_selection.R
#
# Purpose:
# Perform feature selection on the 527 CpGs identified through EWAS using
# three approaches:
#
#   1. Random Forest (RF)
#   2. LASSO
#   3. Random Forest Recursive Feature Elimination (RF-RFE)
#
# Input:
# d527.Rdata
#   - 258 samples
#   - 527 EWAS-selected CpGs
#   - Binary outcome variable Gr (0 = Control, 1 = GDM)
#
# Output:
# Data.xlsx containing three worksheets:
#
#   DataRF : 40 CpGs selected using Random Forest + Gr
#   DataLS : 21 CpGs selected using LASSO + Gr
#   DataBC : 40 CpGs selected using RF-RFE + Gr
#
# ==============================================================================


# ==============================================================================
# 1. Load required packages
# ==============================================================================

library(caret)
library(randomForest)
library(glmnet)
library(dplyr)
library(writexl)


# ==============================================================================
# 2. Load EWAS-selected dataset
# ==============================================================================

load("data/d527.Rdata")

# Convert outcome to factor
d527$Gr <- as.factor(d527$Gr)

# Check dataset
dim(d527)
table(d527$Gr)


# ==============================================================================
# 3. Prepare predictors and outcome
# ==============================================================================

# CpG predictors
x <- d527[, colnames(d527) != "Gr"]

# GDM outcome
y <- d527$Gr


# ==============================================================================
# 4. RANDOM FOREST FEATURE SELECTION
# ==============================================================================

# Define repeated cross-validation
cv_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5
)

# Train Random Forest model
set.seed(1623)

rf_first <- train(
  x = x,
  y = y,
  method = "rf",
  trControl = cv_control,
  importance = TRUE
)

# Extract variable importance
impo_df <- varImp(rf_first)

impo_scores <- impo_df$importance

impo_scores <- data.frame(
  CpG = rownames(impo_scores),
  Importance = impo_scores[, 1]
)

# Rank CpGs according to variable importance
impo_scores <- impo_scores[
  order(impo_scores$Importance, decreasing = TRUE),
]

# Select top 40 CpGs
rf_selected_cpgs <- impo_scores$CpG[1:40]

# Display selected CpGs
rf_selected_cpgs

# Check number of selected CpGs
length(rf_selected_cpgs)


# ==============================================================================
# 5. LASSO FEATURE SELECTION
# ==============================================================================

# Convert predictors to matrix for glmnet
x_lasso <- as.matrix(x)

# Binary outcome
y_lasso <- as.factor(y)

# Fit cross-validated LASSO model
set.seed(453)

lasso_first <- cv.glmnet(
  x_lasso,
  y_lasso,
  alpha = 1,
  family = "binomial"
)

# Extract variables with non-zero coefficients at lambda.min
selected_cpg_lasso <- rownames(
  coef(lasso_first, s = "lambda.min")
)[
  which(
    coef(lasso_first, s = "lambda.min") != 0
  )
]

# Remove intercept from the selected variables
lasso_selected_cpgs <- setdiff(
  selected_cpg_lasso,
  "(Intercept)"
)

# Display selected CpGs
lasso_selected_cpgs

# Check number of selected CpGs
length(lasso_selected_cpgs)


# ==============================================================================
# 6. RF-RFE / BACKWARD FEATURE SELECTION
# ==============================================================================

# Define repeated cross-validation for recursive feature elimination
rfe_control <- rfeControl(
  functions = rfFuncs,
  method = "repeatedcv",
  number = 10,
  repeats = 5
)

# Candidate subset sizes evaluated
sizes <- seq(10, 520, 30)

# Perform recursive feature elimination
set.seed(183)

rfe_first <- rfe(
  x = x,
  y = y,
  sizes = sizes,
  rfeControl = rfe_control
)

# Extract variable importance
impo_rfe <- varImp(rfe_first)

impo_rfe_scores <- data.frame(
  CpG = rownames(impo_rfe),
  Importance = impo_rfe[, 1]
)

# Rank CpGs according to variable importance
impo_rfe_scores <- impo_rfe_scores[
  order(impo_rfe_scores$Importance, decreasing = TRUE),
]

# Select top 40 CpGs
rfe_selected_cpgs <- impo_rfe_scores$CpG[1:40]

# Display selected CpGs
rfe_selected_cpgs

# Check number of selected CpGs
length(rfe_selected_cpgs)


# ==============================================================================
# 7. Verify selected CpGs
# ==============================================================================

# Confirm that all selected CpGs are present in the original
# 527-CpG EWAS dataset

all(rf_selected_cpgs %in% colnames(d527))

all(lasso_selected_cpgs %in% colnames(d527))

all(rfe_selected_cpgs %in% colnames(d527))


# ==============================================================================
# 8. Create the three feature-selected datasets
# ==============================================================================

# Random Forest:
# 40 selected CpGs + outcome
DataRF <- d527[
  ,
  c(rf_selected_cpgs, "Gr"),
  drop = FALSE
]


# LASSO:
# 21 selected CpGs + outcome
DataLS <- d527[
  ,
  c(lasso_selected_cpgs, "Gr"),
  drop = FALSE
]


# RF-RFE / backward selection:
# 40 selected CpGs + outcome
DataBC <- d527[
  ,
  c(rfe_selected_cpgs, "Gr"),
  drop = FALSE
]


# ==============================================================================
# 9. Check final datasets
# ==============================================================================

# Dimensions
dim(DataRF)
dim(DataLS)
dim(DataBC)

# Number of selected CpGs
ncol(DataRF) - 1
ncol(DataLS) - 1
ncol(DataBC) - 1

# Check outcome distributions
table(DataRF$Gr)
table(DataLS$Gr)
table(DataBC$Gr)


# ==============================================================================
# 10. Export feature-selected datasets
# ==============================================================================

# Create a single Excel workbook containing the three feature-selected
# datasets as separate worksheets.
#
# Original analysis file:
# E:/FilteredVariableAnalysis/Data.xlsx
#
# A relative path is used here for repository reproducibility.

write_xlsx(
  list(
    DataRF = DataRF,
    DataLS = DataLS,
    DataBC = DataBC
  ),
  path = "data/Data.xlsx"
)


# ==============================================================================
# 11. Final checks
# ==============================================================================

# Expected:
# DataRF = 258 samples x 41 columns (40 CpGs + Gr)
# DataLS = 258 samples x 22 columns (21 CpGs + Gr)
# DataBC = 258 samples x 41 columns (40 CpGs + Gr)

dim(DataRF)
dim(DataLS)
dim(DataBC)

# Selected CpG lists
rf_selected_cpgs
lasso_selected_cpgs
rfe_selected_cpgs
