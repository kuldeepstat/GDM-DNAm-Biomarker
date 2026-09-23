################################################################################
# 04_feature_selection.R
#
# Feature selection of EWAS-identified CpGs
#
# Methods:
#   1. Random Forest feature importance
#   2. LASSO regression
#   3. Random Forest Recursive Feature Elimination (RF-RFE)
#
# Input:
#   527 CpGs identified from the EWAS
#
# Output:
#   RF-selected CpGs    : 40
#   LASSO-selected CpGs : 21
#   RF-RFE-selected CpGs: 40
#
# This script retains the analytical workflow used in the study.
################################################################################


# ==============================================================================
# Packages
# ==============================================================================

library(ggplot2)
library(rio)
library(tidyverse)

options(digits = 4)
options(dplyr.print_min = Inf)
options(scipen = 999)  # turn off scientific notations

library(gtsummary)
library(gt)
library(multtest)
library(caret)
library(randomForest)
library(glmnet) # for lasso
library(gridExtra)


# ==============================================================================
# Data set
# ==============================================================================

# First trimester:
# Shortlisted CpGs (527 in Early Trimester; n = 258;
# outcome variable: GDM status, 0: No GDM; 1: GDM)

# importing the first trimester data

first_cpg <- read.table(
  'DATASET_FIRSTTRIMESTER527.txt',
  header = T,
  sep = ' ',
  row.names = 1
)


# ==============================================================================
# Modify the dataset
# ==============================================================================

# change the default colname
first <- first_cpg |> rename(SampleID = 'X.1')
# glimpse(first[, c(1, 525:528)])

# convert the response variable to factor
first <- first |> mutate_at(c('Gr'), as.factor)


# ==============================================================================
# Select only the numbers from sample ID
# ==============================================================================

first <- first |>
  mutate(Sample_ID = str_extract(SampleID, "\\d{5}")) |>
  mutate(Sample_ID = as.numeric(Sample_ID)) |>
  select(SampleID, Sample_ID, everything())


# ==============================================================================
# Prepare data matrix of variables and response
# ==============================================================================

# Prepare datasets for glmnet
# get the matrix of predictors and a response variable in both datasets

x_first <- first |> select(-c(SampleID, Samp, Gr))
y_first <- first |> select(Gr)


# ==============================================================================
# Embedded Feature Selection Methods: RF and LASSO
# ==============================================================================

# Embedded methods perform feature selection during the model training process.
# Regularized models like Lasso and Random Forest can be used for this purpose.


# ==============================================================================
# 1. Random Forest
# ==============================================================================

set.seed(453) # ??

# Define a control function for cross-validation

cv_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5
)


# Train Random Forest model using caret on the first trimester before filtering

set.seed(1623)

rf_first <- train(
  Gr ~ .,
  data = first,
  method = "rf",
  trControl = cv_control,
  importance = TRUE
)


# Extract cpgs importance

impo_df <- varImp(rf_first)
impo_scores <- impo_df$importance

impo_scores <- data.frame(
  CpG = rownames(impo_scores),
  Importance = impo_scores[, 1]
)


# Sort by importance

impo_scores <- impo_scores[
  order(
    impo_scores$Importance,
    decreasing = TRUE
  ),
]


# Select top 40 cpgs

top_cpgs <- impo_scores$CpG[1:40]


# top 40 cpgs selected using RF before filtering

rf_selected_cpgs <- c(
  "cg17967426", "cg24745753", "cg25425078",
  "cg07661849", "cg11001216", "cg11037466",
  "cg10139015", "cg22529952", "cg24489237",
  "cg07919162", "cg04539775", "cg09217522",
  "cg11021810", "cg05725489", "cg13653316",
  "cg03940688", "cg14730811", "cg19502700",
  "cg10685380", "cg26333513", "cg27104437",
  "cg10509965", "cg22064129", "cg10591475",
  "cg23507676", "cg08173730", "cg04391685",
  "cg25117809", "cg00078968", "cg08872579",
  "cg04985016", "cg09597829", "cg18153061",
  "cg09518293", "cg10413550", "cg25770783",
  "cg02489379", "cg04900672", "cg16311883",
  "cg07939646"
)


# ==============================================================================
# 2. LASSO
# ==============================================================================

# Applying Lasso (L1 regularization) using the glmnet package to select CpGs
# where the coefficients are non-zero.
#
# The Lasso method tends to select one variable from a group of highly
# correlated variables and shrink the others to exactly zero, effectively
# performing variable selection. This helps in reducing multicollinearity
# because it results in a sparse model where only one of the correlated
# variables remains.


# Prepare data for glmnet before filtering

x <- as.matrix(
  x_first |> select(-c(Sample_ID, SampleID))
)

y <- as.factor(unlist(y_first))


# Fit LASSO model

set.seed(453)

lasso_first <- cv.glmnet(
  x,
  y,
  alpha = 1,
  family = "binomial"
)


# Select non-zero coefficients

selected_cpg_lasso <-
  rownames(
    coef(
      lasso_first,
      s = "lambda.min"
    )
  )[
    which(
      coef(
        lasso_first,
        s = "lambda.min"
      ) != 0
    )
  ]

selected_cpg_lasso <- selected_cpg_lasso[-1]


# cpgs selected before removing correlated cpgs using lasso

lasso_selected_cpgs <- c(
  "cg16311883", "cg11000420", "cg10509965",
  "cg19749898", "cg07919162", "cg15738154",
  "cg11037466", "cg15980170", "cg16098718",
  "cg22454440", "cg11001216", "cg25425078",
  "cg16269199", "cg01400712", "cg17967426",
  "cg06821582", "cg03258272", "cg00542351",
  "cg14336654", "cg00078968", "cg09387867"
)


# ==============================================================================
# 3. Recursive Feature Elimination (RFE) (aka Backwards Selection)
# ==============================================================================

set.seed(8798)

# Define the control using random forest and cross-validation

rfe_control <- rfeControl(
  functions = rfFuncs,
  method = "repeatedcv",
  number = 10,
  repeats = 5
)


# Perform RFE

set.seed(183)

# sizes: an integer vector for the specific subset sizes that should be tested
# (which need not to include ncol(x))

sizes = seq(10, 520, 30)

rfe_first <- rfe(
  x,
  y,
  sizes = sizes,
  rfeControl = rfe_control
)


# Print the results

print(rfe_first)
print(predictors(rfe_first))
plot(rfe_first, type = c("g", "o"))

rfe_first$optsize  # number of CpG with max accuracy


# Selected CpGs from RFE

rfe_selected_cpgs <- c(
  "cg24745753", "cg25425078", "cg17967426",
  "cg07661849", "cg10139015", "cg09217522",
  "cg07919162", "cg26333513", "cg24489237",
  "cg05725489", "cg03940688", "cg10591475",
  "cg11021810", "cg13653316", "cg08872579",
  "cg04539775", "cg00078968", "cg09518293",
  "cg11037466", "cg14730811", "cg22529952",
  "cg11001216", "cg07939646", "cg10509965",
  "cg04985016", "cg15980170", "cg22797968",
  "cg10302336", "cg27104437", "cg26534993",
  "cg16311883", "cg13426133", "cg23507676",
  "cg16052686", "cg15149117", "cg06515734",
  "cg24652415", "cg13208845", "cg07159286",
  "cg04227758"
)


################################################################################
# Feature Selection Workflow
#
# Random Forest:
# We used 10-fold cross-validation repeated 5 times. CpGs were ranked by
# importance and the top 40 CpGs were selected.
#
# LASSO:
# LASSO was used for feature selection and selected 21 CpGs with non-zero
# coefficients.
#
# RFE:
# Recursive Feature Elimination (RFE) was performed using 10-fold
# cross-validation repeated 5 times. The highest cross-validated accuracy was
# achieved with 40 CpGs. Therefore, the top 40 CpGs were selected.
################################################################################
