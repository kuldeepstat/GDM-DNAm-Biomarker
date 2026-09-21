# ==============================================================================
# 05_RF_panel_building.R
# Random Forest biomarker panel building after three feature-selection methods
#
# Input:  data/Data.xlsx
#   - DataRF: Random Forest-selected CpGs
#   - DataLS: LASSO-selected CpGs
#   - DataBC: RF-RFE/backward-selection CpGs
#
# Outputs:
#   - results/ResultRF_UPDATED.xlsx
#   - results/ResultLS.xlsx
#   - results/ResultBC.xlsx
#
# The analysis logic follows the original study scripts. Repository-relative
# paths are used so the analysis can be run from the repository root.
# ==============================================================================

library(e1071)
library(pROC)
library(caret)
library(PRROC)
library(foreach)
library(doParallel)
library(readxl)
library(ggplot2)
library(tidyverse)
library(plyr)
library(writexl)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. RANDOM FOREST-SELECTED FEATURES (DataRF)
# Source: updated RF analysis supplied by the authors
# ==============================================================================

library(readxl)
Data <- read_excel("data/Data.xlsx", 
    sheet = "DataRF")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)


a<-Sys.time()
a
# Load your dataset (replace 'MAIN' with your actual dataset)
train_data <- MAIN

# Define the outcome variable
outcome_var <- "Gr"

# Convert the outcome variable levels to valid R variable names
train_data[[outcome_var]] <- factor(train_data[[outcome_var]], levels = unique(train_data[[outcome_var]]))
levels(train_data[[outcome_var]]) <- make.names(levels(train_data[[outcome_var]]))

# Define the predictors
predictors <- setdiff(names(train_data), outcome_var)

# Number of predictors to select
num_predictors_to_select <- 10

# Number of iterations
num_iterations <- 10

# Initialize variables to store overall best results
overall_best_auc <- 0
overall_best_model <- NULL
overall_best_predictors <- NULL
best_test_set <- NULL

# Initialize a data frame to store the results for each iteration
results_df <- data.frame(Iteration = integer(),
                         Round = integer(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         F_Score = numeric(),
                         Accuracy = numeric(),
                         Brier_Inaccuracy = numeric(),
                         Logarithmic_Score = numeric(),
                         AUPRC = numeric(),
                         Divergence = numeric(),
                         Kruskal_Statistic = numeric(),
                         Selected_Predictors = character(),
                         stringsAsFactors = FALSE)

# Function to calculate additional metrics
calc_metrics <- function(predictions, true_labels, prob_predictions) {
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate F-measure
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  f_measure <- 2 * (precision * recall) / (precision + recall)
  
  # Calculate Brier score
  true_labels_binary <- as.numeric(true_labels) - 1
  brier_score <- mean((prob_predictions - true_labels_binary) ^ 2)
  
  # Calculate logarithmic score (cross-entropy loss)
  cross_entropy <- -mean(true_labels_binary * log(prob_predictions) + (1 - true_labels_binary) * log(1 - prob_predictions))
  
  # Calculate area under the precision-recall curve (AUPRC)
  pr_curve <- pr.curve(scores.class0 = prob_predictions, weights.class0 = true_labels_binary, curve = TRUE)
  auprc <- pr_curve$auc.integral
  
  # Calculate divergence
  p_c1 <- prob_predictions[true_labels_binary == 1]
  p_c0 <- prob_predictions[true_labels_binary == 0]
  mean_c1 <- mean(p_c1)
  mean_c0 <- mean(p_c0)
  var_c1 <- var(p_c1)
  var_c0 <- var(p_c0)
  divergence <- ((mean_c1 - mean_c0)^2) / (0.5 * (var_c1 + var_c0))
  
  # Calculate Kolmogorov-Smirnov statistic
  ks_stat <- ks.test(p_c1, p_c0)$statistic
  
  # Calculate accuracy
  accuracy <- sum(predictions == true_labels) / length(true_labels)
  
  # Return metrics as a numeric vector
  return(c(f_measure, brier_score, cross_entropy, auprc, divergence, ks_stat, accuracy, recall, cm$byClass["Specificity"]))
}

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop for multiple iterations in parallel
results_list <- foreach(iter = 1:num_iterations, .packages = c("caret", "pROC", "PRROC")) %dopar% {
  # Set a different seed for each iteration
  set.seed(iter)
  
  # Create a data partition for cross-validation
  train_indices <- createDataPartition(train_data[[outcome_var]], p = 0.67, list = FALSE)
  
  # Split the data into training and testing sets
  train_set <- train_data[train_indices, ]
  test_set <- train_data[-train_indices, ]
  
  # Initialize variables to store results for the current iteration
  selected_predictors <- character(0)
  best_auc <- 0
  best_model <- NULL
  best_selected_predictors <- NULL
  best_round_metrics <- NULL
  round_results <- data.frame()
  
  # Stepwise selection for predictors using Random Forest with cross-validation
  for (i in 1:num_predictors_to_select) {
    # Initialize variables for this round
    best_round_auc <- 0
    best_round_predictor <- NULL
    best_round_model <- NULL
    
    # Loop over predictors to find the best one to add
    for (predictor in setdiff(predictors, selected_predictors)) {
      current_predictors <- c(selected_predictors, predictor)
      
      # Set up the control parameters for cross-validation
      ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary)  # 5-fold cross-validation with class probabilities
      
      # Create formula for Random Forest
      formula_rf <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Set up the grid for parameter tuning
      grid <- expand.grid(mtry = seq(2, 40, 2))  # Tuning mtry parameters
      
      # Train the Random Forest model with cross-validation and parameter tuning
      rf_model <- train(formula_rf, data = train_set, method = "rf", trControl = ctrl, tuneGrid = grid)
      
      # Make probability predictions on the test set
      prob_predictions <- predict(rf_model, newdata = test_set, type = "prob")[, 2]
      
      # Make class predictions on the test set
      predictions <- factor(ifelse(prob_predictions > 0.5, "X1", "X0"), levels = levels(test_set[[outcome_var]]))
      
      # Calculate AUC
      roc_curve <- roc(test_set[[outcome_var]], prob_predictions)
      round_auc <- auc(roc_curve)
      
      # Calculate additional metrics
      round_metrics <- calc_metrics(predictions, test_set[[outcome_var]], prob_predictions)
      
      # Update best round results
      if (round_auc > best_round_auc) {
        best_round_auc <- round_auc
        best_round_predictor <- predictor
        best_round_model <- rf_model
        best_round_metrics <- round_metrics
      }
    }
    
    # Update selected predictors
    selected_predictors <- c(selected_predictors, best_round_predictor)
    
    # Update best results
    if (best_round_auc > best_auc) {
      best_auc <- best_round_auc
      best_model <- best_round_model
      best_selected_predictors <- selected_predictors
      best_round_metrics <- best_round_metrics
    }
    
    # Store the results for the current round
    round_results <- rbind(round_results, data.frame(
      Iteration = iter,
      Round = i,
      AUC = best_round_auc,
      Sensitivity = best_round_metrics[8],  # Index for sensitivity
      Specificity = best_round_metrics[9],  # Index for specificity
      F_Score = best_round_metrics[1],      # Index for F-measure
      Accuracy = best_round_metrics[7],     # Index for accuracy
      Brier_Inaccuracy = best_round_metrics[2],  # Index for Brier score
      Logarithmic_Score = best_round_metrics[3],  # Index for cross-entropy
      AUPRC = best_round_metrics[4],        # Index for AUPRC
      Divergence = best_round_metrics[5],   # Index for divergence
      Kruskal_Statistic = best_round_metrics[6],  # Index for KS statistic
      Selected_Predictors = paste(selected_predictors, collapse = ", ")
    ))
    
    # Print round results
    cat("Iteration:", iter,
        "- Round:", i,
        "- Selected predictors:", paste(selected_predictors, collapse = ", "),
        "- AUC:", best_round_auc, "\n")
  }
  
  # Return the results for the current iteration
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
}

# Stop the cluster
stopCluster(cl)

# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}
RESULTRF<-results_df


# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}


# Assign the overall best model to bm
bm <- overall_best_model

if (!is.null(bm)) {
  # Make predictions with the best model
  predictions <- predict(bm, newdata = best_test_set)
  prob_predictions <- predict(bm, newdata = best_test_set, type = "prob")[,2]
  
  # Calculate AUC
  roc_curve <- roc(best_test_set[[outcome_var]], prob_predictions)
  round_auc <- auc(roc_curve)
  
  # Calculate sensitivity and specificity
  confusion_matrix <- table(predictions, best_test_set[[outcome_var]])
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
  specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
  
  # Print final results
  print(t(confusion_matrix))
  print(paste("AUC:", round_auc))
  print(paste("Specificity:", sensitivity))
  print(paste("Sensitivity:", specificity))
} else {
  cat("Overall best model is NULL\n")
}
bm

# Print the overall best-selected predictors
cat("Overall best-selected predictors corresponding to the highest AUC:", paste(overall_best_predictors, collapse = ", "), "\n")


rf_prediction <- predict(bm, best_test_set, type = "prob")
# ROC curves
library(pROC)
ROC_rf <- roc(best_test_set$Gr, rf_prediction[,2])
ROC_data<-as.data.frame(cbind(ROC_rf$sensitivities, ROC_rf$specificities, ROC_rf$thresholds))
auc<-round(round_auc, 3)
ggplot(ROC_data, aes(1-V2, V1)) + 
  geom_path(size=1.5, col="red") +
  geom_abline(intercept = 0, slope = 1, size = 1.5)+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  labs(x="1-Specificity", y="Sensitivity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face="bold"))+
  annotate("text", x=0.1, y=1, label= paste("AUC:", auc), size=4, col="brown") 


# Perform bootstrap validation
if (!is.null(overall_best_model)) {
  # Function to train the Random Forest model and return metrics
  train_rf_model <- function(train_data, test_data, best_model, outcome_var) {
    # Train the model
    rf_model <- train(as.formula(paste(outcome_var, "~", paste(best_model$finalModel$xNames, collapse = "+"))),
                      data = train_data, method = "rf", trControl = trainControl(method = "none"),
                      tuneGrid = best_model$bestTune)
    
    # Make probability predictions on the test set
    prob_predictions <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    
    # Make class predictions on the test set
    predictions <- predict(rf_model, newdata = test_data)
    
    # Calculate AUC
    roc_curve <- roc(test_data[[outcome_var]], prob_predictions)
    round_auc <- auc(roc_curve)
    
    # Calculate additional metrics
    metrics <- calc_metrics(predictions, test_data[[outcome_var]], prob_predictions)
    
    # Return metrics as numeric vector
    return(c(auc = round_auc, metrics))
  }
  
  # Perform bootstrap validation
  set.seed(123) # Set seed for reproducibility
  num_iterations <- 1000
  boot_results <- vector("list", num_iterations)
  
  for (i in 1:num_iterations) {
    # Sample indices with replacement
    indices <- sample(nrow(MAIN), replace = TRUE)
    
    # Extract bootstrap sample
    bootstrap_sample <- MAIN[indices, ]
    
    # Out-of-bootstrap sample
    out_of_bootstrap <- MAIN[-indices, ]
    
    # Ensure factor levels are consistent between training and testing sets
    out_of_bootstrap[[outcome_var]] <- factor(out_of_bootstrap[[outcome_var]], levels = levels(bootstrap_sample[[outcome_var]]))
    
    # Train the model on bootstrap sample and evaluate on out-of-bootstrap sample
    boot_results[[i]] <- train_rf_model(bootstrap_sample, out_of_bootstrap, overall_best_model, outcome_var)
  }
  
  # Combine bootstrap results into a dataframe and set column names
  boot_df <- do.call(rbind, boot_results)
  colnames(boot_df) <- c("AUC", "F_Score", "Brier_Inaccuracy", "Logarithmic_Score", "AUPRC", "Divergence", "Kruskal_Statistic", "Accuracy", "Sensitivity", "Specificity")
  
  # Summarize bootstrap results
  print(summary(boot_df))
} else {
  cat("Overall best model is NULL\n")
}
round<-seq(1, 1000, 1)
boot_dfRF<-as.data.frame(cbind(round, boot_df))


library(ggplot2)
library(tidyverse)
d<-sample_n(boot_dfRF[,c(1,2,9,10,11)], 500)
ld<-gather(d, Metric, Value, -round)

ggplot(ld, aes(round, Value))+
  geom_line(lwd=1)+
  ylim(0, 1)+
  theme_classic()+
  facet_wrap(Metric~.)


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ld[,-1], varname="Value", 
                    groupnames=c("Metric"))

library(ggplot2)
# Default bar plot
p<- ggplot(df2, aes(x=Metric, y=Value, fill = Metric)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme_classic()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=.2, size=0.8,
                position=position_dodge(.9)) 
print(p)


library(writexl)
list<-list("RoundsResult"= RESULTRF, "Bootstrap"=boot_dfRF, "RocOutput"=ROC_data)
write_xlsx(list, "results/ResultRF_UPDATED.xlsx")


# ==============================================================================
# 2. LASSO-SELECTED FEATURES (DataLS)
# Source: original RF analysis output
# ==============================================================================

library(readxl)
Data <- read_excel("data/Data.xlsx", 
    sheet = "DataLS")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)

a<-Sys.time()
a

# Load your dataset (replace 'MAIN' with your actual dataset)
train_data <- MAIN

# Define the outcome variable
outcome_var <- "Gr"

# Convert the outcome variable levels to valid R variable names
train_data[[outcome_var]] <- factor(train_data[[outcome_var]], levels = unique(train_data[[outcome_var]]))
levels(train_data[[outcome_var]]) <- make.names(levels(train_data[[outcome_var]]))

# Define the predictors
predictors <- setdiff(names(train_data), outcome_var)

# Number of predictors to select
num_predictors_to_select <- 10

# Number of iterations
num_iterations <- 10

# Initialize variables to store overall best results
overall_best_auc <- 0
overall_best_model <- NULL
overall_best_predictors <- NULL
best_test_set <- NULL

# Initialize a data frame to store the results for each iteration
results_df <- data.frame(Iteration = integer(),
                         Round = integer(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         F_Score = numeric(),
                         Accuracy = numeric(),
                         Brier_Inaccuracy = numeric(),
                         Logarithmic_Score = numeric(),
                         AUPRC = numeric(),
                         Divergence = numeric(),
                         Kruskal_Statistic = numeric(),
                         Selected_Predictors = character(),
                         stringsAsFactors = FALSE)

# Function to calculate additional metrics
calc_metrics <- function(predictions, true_labels, prob_predictions) {
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate F-measure
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  f_measure <- 2 * (precision * recall) / (precision + recall)
  
  # Calculate Brier score
  true_labels_binary <- as.numeric(true_labels) - 1
  brier_score <- mean((prob_predictions - true_labels_binary) ^ 2)
  
  # Calculate logarithmic score (cross-entropy loss)
  cross_entropy <- -mean(true_labels_binary * log(prob_predictions) + (1 - true_labels_binary) * log(1 - prob_predictions))
  
  # Calculate area under the precision-recall curve (AUPRC)
  pr_curve <- pr.curve(scores.class0 = prob_predictions, weights.class0 = true_labels_binary, curve = TRUE)
  auprc <- pr_curve$auc.integral
  
  # Calculate divergence
  p_c1 <- prob_predictions[true_labels_binary == 1]
  p_c0 <- prob_predictions[true_labels_binary == 0]
  mean_c1 <- mean(p_c1)
  mean_c0 <- mean(p_c0)
  var_c1 <- var(p_c1)
  var_c0 <- var(p_c0)
  divergence <- ((mean_c1 - mean_c0)^2) / (0.5 * (var_c1 + var_c0))
  
  # Calculate Kolmogorov-Smirnov statistic
  ks_stat <- ks.test(p_c1, p_c0)$statistic
  
  # Calculate accuracy
  accuracy <- sum(predictions == true_labels) / length(true_labels)
  
  # Return metrics as a numeric vector
  return(c(f_measure, brier_score, cross_entropy, auprc, divergence, ks_stat, accuracy, recall, cm$byClass["Specificity"]))
}

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop for multiple iterations in parallel
results_list <- foreach(iter = 1:num_iterations, .packages = c("caret", "pROC", "PRROC")) %dopar% {
  # Set a different seed for each iteration
  set.seed(iter)
  
  # Create a data partition for cross-validation
  train_indices <- createDataPartition(train_data[[outcome_var]], p = 0.67, list = FALSE)
  
  # Split the data into training and testing sets
  train_set <- train_data[train_indices, ]
  test_set <- train_data[-train_indices, ]
  
  # Initialize variables to store results for the current iteration
  selected_predictors <- character(0)
  best_auc <- 0
  best_model <- NULL
  best_selected_predictors <- NULL
  best_round_metrics <- NULL
  round_results <- data.frame()
  
  # Stepwise selection for predictors using Random Forest with cross-validation
  for (i in 1:num_predictors_to_select) {
    # Initialize variables for this round
    best_round_auc <- 0
    best_round_predictor <- NULL
    best_round_model <- NULL
    
    # Loop over predictors to find the best one to add
    for (predictor in setdiff(predictors, selected_predictors)) {
      current_predictors <- c(selected_predictors, predictor)
      
      # Set up the control parameters for cross-validation
      ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary)  # 5-fold cross-validation with class probabilities
      
      # Create formula for Random Forest
      formula_rf <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Set up the grid for parameter tuning
      grid <- expand.grid(mtry = seq(2, 27, 2))  # Tuning mtry parameters
      
      # Train the Random Forest model with cross-validation and parameter tuning
      rf_model <- train(formula_rf, data = train_set, method = "rf", trControl = ctrl, tuneGrid = grid)
      
      # Make probability predictions on the test set
      prob_predictions <- predict(rf_model, newdata = test_set, type = "prob")[, 2]
      
      # Make class predictions on the test set
      predictions <- factor(ifelse(prob_predictions > 0.5, "X1", "X0"), levels = levels(test_set[[outcome_var]]))
      
      # Calculate AUC
      roc_curve <- roc(test_set[[outcome_var]], prob_predictions)
      round_auc <- auc(roc_curve)
      
      # Calculate additional metrics
      round_metrics <- calc_metrics(predictions, test_set[[outcome_var]], prob_predictions)
      
      # Update best round results
      if (round_auc > best_round_auc) {
        best_round_auc <- round_auc
        best_round_predictor <- predictor
        best_round_model <- rf_model
        best_round_metrics <- round_metrics
      }
    }
    
    # Update selected predictors
    selected_predictors <- c(selected_predictors, best_round_predictor)
    
    # Update best results
    if (best_round_auc > best_auc) {
      best_auc <- best_round_auc
      best_model <- best_round_model
      best_selected_predictors <- selected_predictors
      best_round_metrics <- best_round_metrics
    }
    
    # Store the results for the current round
    round_results <- rbind(round_results, data.frame(
      Iteration = iter,
      Round = i,
      AUC = best_round_auc,
      Sensitivity = best_round_metrics[8],  # Index for sensitivity
      Specificity = best_round_metrics[9],  # Index for specificity
      F_Score = best_round_metrics[1],      # Index for F-measure
      Accuracy = best_round_metrics[7],     # Index for accuracy
      Brier_Inaccuracy = best_round_metrics[2],  # Index for Brier score
      Logarithmic_Score = best_round_metrics[3],  # Index for cross-entropy
      AUPRC = best_round_metrics[4],        # Index for AUPRC
      Divergence = best_round_metrics[5],   # Index for divergence
      Kruskal_Statistic = best_round_metrics[6],  # Index for KS statistic
      Selected_Predictors = paste(selected_predictors, collapse = ", ")
    ))
    
    # Print round results
    cat("Iteration:", iter,
        "- Round:", i,
        "- Selected predictors:", paste(selected_predictors, collapse = ", "),
        "- AUC:", best_round_auc, "\n")
  }
  
  # Return the results for the current iteration
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
}

# Stop the cluster
stopCluster(cl)

# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}
RESULTLASSO<-results_df

# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}


# Assign the overall best model to bm
bm <- overall_best_model

if (!is.null(bm)) {
  # Make predictions with the best model
  predictions <- predict(bm, newdata = best_test_set)
  prob_predictions <- predict(bm, newdata = best_test_set, type = "prob")[,2]
  
  # Calculate AUC
  roc_curve <- roc(best_test_set[[outcome_var]], prob_predictions)
  round_auc <- auc(roc_curve)
  
  # Calculate sensitivity and specificity
  confusion_matrix <- table(predictions, best_test_set[[outcome_var]])
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
  specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
  
  # Print final results
  print(t(confusion_matrix))
  print(paste("AUC:", round_auc))
  print(paste("Specificity:", sensitivity))
  print(paste("Sensitivity:", specificity))
} else {
  cat("Overall best model is NULL\n")
}

bm

# Print the overall best-selected predictors
cat("Overall best-selected predictors corresponding to the highest AUC:", paste(overall_best_predictors, collapse = ", "), "\n")

rf_prediction <- predict(bm, best_test_set, type = "prob")
# ROC curves
library(pROC)
ROC_rf <- roc(best_test_set$Gr, rf_prediction[,2])

ROC_data<-as.data.frame(cbind(ROC_rf$sensitivities, ROC_rf$specificities))
auc<-round(round_auc, 3)
ggplot(ROC_data, aes(1-V2, V1)) + 
  geom_path(size=1.5, col="red") +
  geom_abline(intercept = 0, slope = 1, size = 1.5)+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  labs(x="1-Specificity", y="Sensitivity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face="bold"))+
  annotate("text", x=0.1, y=1, label= paste("AUC:", auc), size=5.8, col="brown") 

# Perform bootstrap validation
if (!is.null(overall_best_model)) {
  # Function to train the Random Forest model and return metrics
  train_rf_model <- function(train_data, test_data, best_model, outcome_var) {
    # Train the model
    rf_model <- train(as.formula(paste(outcome_var, "~", paste(best_model$finalModel$xNames, collapse = "+"))),
                      data = train_data, method = "rf", trControl = trainControl(method = "none"),
                      tuneGrid = best_model$bestTune)
    
    # Make probability predictions on the test set
    prob_predictions <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    
    # Make class predictions on the test set
    predictions <- predict(rf_model, newdata = test_data)
    
    # Calculate AUC
    roc_curve <- roc(test_data[[outcome_var]], prob_predictions)
    round_auc <- auc(roc_curve)
    
    # Calculate additional metrics
    metrics <- calc_metrics(predictions, test_data[[outcome_var]], prob_predictions)
    
    # Return metrics as numeric vector
    return(c(auc = round_auc, metrics))
  }
  
  # Perform bootstrap validation
  set.seed(123) # Set seed for reproducibility
  num_iterations <- 1000
  boot_results <- vector("list", num_iterations)
  
  for (i in 1:num_iterations) {
    # Sample indices with replacement
    indices <- sample(nrow(MAIN), replace = TRUE)
    
    # Extract bootstrap sample
    bootstrap_sample <- MAIN[indices, ]
    
    # Out-of-bootstrap sample
    out_of_bootstrap <- MAIN[-indices, ]
    
    # Ensure factor levels are consistent between training and testing sets
    out_of_bootstrap[[outcome_var]] <- factor(out_of_bootstrap[[outcome_var]], levels = levels(bootstrap_sample[[outcome_var]]))
    
    # Train the model on bootstrap sample and evaluate on out-of-bootstrap sample
    boot_results[[i]] <- train_rf_model(bootstrap_sample, out_of_bootstrap, overall_best_model, outcome_var)
  }
  
  # Combine bootstrap results into a dataframe and set column names
  boot_df <- do.call(rbind, boot_results)
  colnames(boot_df) <- c("AUC", "F_Score", "Brier_Inaccuracy", "Logarithmic_Score", "AUPRC", "Divergence", "Kruskal_Statistic", "Accuracy", "Sensitivity", "Specificity")
  
  # Summarize bootstrap results
  print(summary(boot_df))
} else {
  cat("Overall best model is NULL\n")
}

round<-seq(1, 1000, 1)
boot_dfLS<-as.data.frame(cbind(round, boot_df))

library(ggplot2)
library(tidyverse)
d<-sample_n(boot_dfLS[,c(1,2,9,10,11)], 500)
ld<-gather(d, Metric, Value, -round)


ggplot(ld, aes(round, Value))+
  geom_line(lwd=1)+
  ylim(0, 1)+
  theme_classic()+
  facet_wrap(Metric~.)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ld[,-1], varname="Value", 
                    groupnames=c("Metric"))

library(ggplot2)
# Default bar plot
p<- ggplot(df2, aes(x=Metric, y=Value, fill = Metric)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme_classic()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=.2, size=0.8,
                position=position_dodge(.9)) 
print(p)

library(writexl)
list<-list("RoundsResult"= RESULTLASSO, "Bootstrap"=boot_dfLS, "ROC OUTPUT"=ROC_data)
write_xlsx(list, "results/ResultLS.xlsx")

# ==============================================================================
# 3. RF-RFE / BACKWARD-SELECTION FEATURES (DataBC)
# Source: original RF analysis output
# ==============================================================================

library(readxl)
Data <- read_excel("data/Data.xlsx", 
    sheet = "DataBC")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)

a<-Sys.time()
a

# Load your dataset (replace 'MAIN' with your actual dataset)
train_data <- MAIN

# Define the outcome variable
outcome_var <- "Gr"

# Convert the outcome variable levels to valid R variable names
train_data[[outcome_var]] <- factor(train_data[[outcome_var]], levels = unique(train_data[[outcome_var]]))
levels(train_data[[outcome_var]]) <- make.names(levels(train_data[[outcome_var]]))

# Define the predictors
predictors <- setdiff(names(train_data), outcome_var)

# Number of predictors to select
num_predictors_to_select <- 10

# Number of iterations
num_iterations <- 10

# Initialize variables to store overall best results
overall_best_auc <- 0
overall_best_model <- NULL
overall_best_predictors <- NULL
best_test_set <- NULL

# Initialize a data frame to store the results for each iteration
results_df <- data.frame(Iteration = integer(),
                         Round = integer(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         F_Score = numeric(),
                         Accuracy = numeric(),
                         Brier_Inaccuracy = numeric(),
                         Logarithmic_Score = numeric(),
                         AUPRC = numeric(),
                         Divergence = numeric(),
                         Kruskal_Statistic = numeric(),
                         Selected_Predictors = character(),
                         stringsAsFactors = FALSE)

# Function to calculate additional metrics
calc_metrics <- function(predictions, true_labels, prob_predictions) {
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate F-measure
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  f_measure <- 2 * (precision * recall) / (precision + recall)
  
  # Calculate Brier score
  true_labels_binary <- as.numeric(true_labels) - 1
  brier_score <- mean((prob_predictions - true_labels_binary) ^ 2)
  
  # Calculate logarithmic score (cross-entropy loss)
  cross_entropy <- -mean(true_labels_binary * log(prob_predictions) + (1 - true_labels_binary) * log(1 - prob_predictions))
  
  # Calculate area under the precision-recall curve (AUPRC)
  pr_curve <- pr.curve(scores.class0 = prob_predictions, weights.class0 = true_labels_binary, curve = TRUE)
  auprc <- pr_curve$auc.integral
  
  # Calculate divergence
  p_c1 <- prob_predictions[true_labels_binary == 1]
  p_c0 <- prob_predictions[true_labels_binary == 0]
  mean_c1 <- mean(p_c1)
  mean_c0 <- mean(p_c0)
  var_c1 <- var(p_c1)
  var_c0 <- var(p_c0)
  divergence <- ((mean_c1 - mean_c0)^2) / (0.5 * (var_c1 + var_c0))
  
  # Calculate Kolmogorov-Smirnov statistic
  ks_stat <- ks.test(p_c1, p_c0)$statistic
  
  # Calculate accuracy
  accuracy <- sum(predictions == true_labels) / length(true_labels)
  
  # Return metrics as a numeric vector
  return(c(f_measure, brier_score, cross_entropy, auprc, divergence, ks_stat, accuracy, recall, cm$byClass["Specificity"]))
}

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop for multiple iterations in parallel
results_list <- foreach(iter = 1:num_iterations, .packages = c("caret", "pROC", "PRROC")) %dopar% {
  # Set a different seed for each iteration
  set.seed(iter)
  
  # Create a data partition for cross-validation
  train_indices <- createDataPartition(train_data[[outcome_var]], p = 0.67, list = FALSE)
  
  # Split the data into training and testing sets
  train_set <- train_data[train_indices, ]
  test_set <- train_data[-train_indices, ]
  
  # Initialize variables to store results for the current iteration
  selected_predictors <- character(0)
  best_auc <- 0
  best_model <- NULL
  best_selected_predictors <- NULL
  best_round_metrics <- NULL
  round_results <- data.frame()
  
  # Stepwise selection for predictors using Random Forest with cross-validation
  for (i in 1:num_predictors_to_select) {
    # Initialize variables for this round
    best_round_auc <- 0
    best_round_predictor <- NULL
    best_round_model <- NULL
    
    # Loop over predictors to find the best one to add
    for (predictor in setdiff(predictors, selected_predictors)) {
      current_predictors <- c(selected_predictors, predictor)
      
      # Set up the control parameters for cross-validation
      ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary)  # 5-fold cross-validation with class probabilities
      
      # Create formula for Random Forest
      formula_rf <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Set up the grid for parameter tuning
      grid <- expand.grid(mtry = seq(2, 40, 2))  # Tuning mtry parameters
      
      # Train the Random Forest model with cross-validation and parameter tuning
      rf_model <- train(formula_rf, data = train_set, method = "rf", trControl = ctrl, tuneGrid = grid)
      
      # Make probability predictions on the test set
      prob_predictions <- predict(rf_model, newdata = test_set, type = "prob")[, 2]
      
      # Make class predictions on the test set
      predictions <- factor(ifelse(prob_predictions > 0.5, "X1", "X0"), levels = levels(test_set[[outcome_var]]))
      
      # Calculate AUC
      roc_curve <- roc(test_set[[outcome_var]], prob_predictions)
      round_auc <- auc(roc_curve)
      
      # Calculate additional metrics
      round_metrics <- calc_metrics(predictions, test_set[[outcome_var]], prob_predictions)
      
      # Update best round results
      if (round_auc > best_round_auc) {
        best_round_auc <- round_auc
        best_round_predictor <- predictor
        best_round_model <- rf_model
        best_round_metrics <- round_metrics
      }
    }
    
    # Update selected predictors
    selected_predictors <- c(selected_predictors, best_round_predictor)
    
    # Update best results
    if (best_round_auc > best_auc) {
      best_auc <- best_round_auc
      best_model <- best_round_model
      best_selected_predictors <- selected_predictors
      best_round_metrics <- best_round_metrics
    }
    
    # Store the results for the current round
    round_results <- rbind(round_results, data.frame(
      Iteration = iter,
      Round = i,
      AUC = best_round_auc,
      Sensitivity = best_round_metrics[8],  # Index for sensitivity
      Specificity = best_round_metrics[9],  # Index for specificity
      F_Score = best_round_metrics[1],      # Index for F-measure
      Accuracy = best_round_metrics[7],     # Index for accuracy
      Brier_Inaccuracy = best_round_metrics[2],  # Index for Brier score
      Logarithmic_Score = best_round_metrics[3],  # Index for cross-entropy
      AUPRC = best_round_metrics[4],        # Index for AUPRC
      Divergence = best_round_metrics[5],   # Index for divergence
      Kruskal_Statistic = best_round_metrics[6],  # Index for KS statistic
      Selected_Predictors = paste(selected_predictors, collapse = ", ")
    ))
    
    # Print round results
    cat("Iteration:", iter,
        "- Round:", i,
        "- Selected predictors:", paste(selected_predictors, collapse = ", "),
        "- AUC:", best_round_auc, "\n")
  }
  
  # Return the results for the current iteration
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
}

# Stop the cluster
stopCluster(cl)

# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}
RESULTBC<-results_df

# Combine results from all iterations
for (result in results_list) {
  results_df <- rbind(results_df, result$round_results)
  
  # Update overall best results if the current iteration's best model is better
  if (!is.null(result$best_auc) && result$best_auc > overall_best_auc) {
    overall_best_auc <- result$best_auc
    overall_best_model <- result$best_model
    overall_best_predictors <- result$best_selected_predictors
    best_test_set <- result$test_set
  }
}


# Assign the overall best model to bm
bm <- overall_best_model

if (!is.null(bm)) {
  # Make predictions with the best model
  predictions <- predict(bm, newdata = best_test_set)
  prob_predictions <- predict(bm, newdata = best_test_set, type = "prob")[,2]
  
  # Calculate AUC
  roc_curve <- roc(best_test_set[[outcome_var]], prob_predictions)
  round_auc <- auc(roc_curve)
  
  # Calculate sensitivity and specificity
  confusion_matrix <- table(predictions, best_test_set[[outcome_var]])
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
  specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
  
  # Print final results
  print(t(confusion_matrix))
  print(paste("AUC:", round_auc))
  print(paste("Specificity:", sensitivity))
  print(paste("Sensitivity:", specificity))
} else {
  cat("Overall best model is NULL\n")
}

bm

# Print the overall best-selected predictors
cat("Overall best-selected predictors corresponding to the highest AUC:", paste(overall_best_predictors, collapse = ", "), "\n")

rf_prediction <- predict(bm, best_test_set, type = "prob")
# ROC curves
library(pROC)
ROC_rf <- roc(best_test_set$Gr, rf_prediction[,2])

ROC_data<-as.data.frame(cbind(ROC_rf$sensitivities, ROC_rf$specificities))
auc<-round(round_auc, 3)
ggplot(ROC_data, aes(1-V2, V1)) + 
  geom_path(size=1.5, col="red") +
  geom_abline(intercept = 0, slope = 1, size = 1.5)+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  labs(x="1-Specificity", y="Sensitivity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face="bold"))+
  annotate("text", x=0.1, y=1, label= paste("AUC:", auc), size=5.8, col="brown") 

# Perform bootstrap validation
if (!is.null(overall_best_model)) {
  # Function to train the Random Forest model and return metrics
  train_rf_model <- function(train_data, test_data, best_model, outcome_var) {
    # Train the model
    rf_model <- train(as.formula(paste(outcome_var, "~", paste(best_model$finalModel$xNames, collapse = "+"))),
                      data = train_data, method = "rf", trControl = trainControl(method = "none"),
                      tuneGrid = best_model$bestTune)
    
    # Make probability predictions on the test set
    prob_predictions <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    
    # Make class predictions on the test set
    predictions <- predict(rf_model, newdata = test_data)
    
    # Calculate AUC
    roc_curve <- roc(test_data[[outcome_var]], prob_predictions)
    round_auc <- auc(roc_curve)
    
    # Calculate additional metrics
    metrics <- calc_metrics(predictions, test_data[[outcome_var]], prob_predictions)
    
    # Return metrics as numeric vector
    return(c(auc = round_auc, metrics))
  }
  
  # Perform bootstrap validation
  set.seed(123) # Set seed for reproducibility
  num_iterations <- 1000
  boot_results <- vector("list", num_iterations)
  
  for (i in 1:num_iterations) {
    # Sample indices with replacement
    indices <- sample(nrow(MAIN), replace = TRUE)
    
    # Extract bootstrap sample
    bootstrap_sample <- MAIN[indices, ]
    
    # Out-of-bootstrap sample
    out_of_bootstrap <- MAIN[-indices, ]
    
    # Ensure factor levels are consistent between training and testing sets
    out_of_bootstrap[[outcome_var]] <- factor(out_of_bootstrap[[outcome_var]], levels = levels(bootstrap_sample[[outcome_var]]))
    
    # Train the model on bootstrap sample and evaluate on out-of-bootstrap sample
    boot_results[[i]] <- train_rf_model(bootstrap_sample, out_of_bootstrap, overall_best_model, outcome_var)
  }
  
  # Combine bootstrap results into a dataframe and set column names
  boot_df <- do.call(rbind, boot_results)
  colnames(boot_df) <- c("AUC", "F_Score", "Brier_Inaccuracy", "Logarithmic_Score", "AUPRC", "Divergence", "Kruskal_Statistic", "Accuracy", "Sensitivity", "Specificity")
  
  # Summarize bootstrap results
  print(summary(boot_df))
} else {
  cat("Overall best model is NULL\n")
}

round<-seq(1, 1000, 1)
boot_dfBC<-as.data.frame(cbind(round, boot_df))

library(ggplot2)
library(tidyverse)
d<-sample_n(boot_dfBC[,c(1,2,9,10,11)], 500)
ld<-gather(d, Metric, Value, -round)


ggplot(ld, aes(round, Value))+
  geom_line(lwd=1)+
  ylim(0, 1)+
  theme_classic()+
  facet_wrap(Metric~.)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ld[,-1], varname="Value", 
                    groupnames=c("Metric"))

library(ggplot2)
# Default bar plot
p<- ggplot(df2, aes(x=Metric, y=Value, fill = Metric)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme_classic()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=.2, size=0.8,
                position=position_dodge(.9)) 
print(p)

library(writexl)
list<-list("RoundsResult"= RESULTBC, "Bootstrap"=boot_dfBC, "ROC OUTPUT"=ROC_data)
write_xlsx(list, "results/ResultBC.xlsx")
