
#--------------------------------------------- RF -----------------------------------------------

```{r, message=FALSE, warning=FALSE}
# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(PRROC)
library(readxl)
library(foreach)
library(doParallel)
```

```{r, message=FALSE, warning=FALSE}
Data <- read_excel("E:/FilteredVariableAnalysis/Data.xlsx", 
                   sheet = "DataRF")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)
```

```{r, warning=FALSE, message=FALSE}
a<-Sys.time()
a

set.seed(120)

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
      
      # Create formula for SVM
      formula_svm <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Train the SVM model with cross-validation
      svm_model <- train(formula_svm, data = train_set, method = "svmRadial", trControl = ctrl, metric = "ROC")
      # Make probability predictions on the test set
      prob_predictions <- predict(svm_model, newdata = test_set, type = "prob")[, 2]
      
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
        best_round_model <- svm_model
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
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_model_params = list(tuneGrid = best_model$bestTune), best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
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
Sys.time()-a
RESULT<-results_df
```


#---------- LIST OF BEST PANELS WITHIN EACH ITERATIONS ---------

```{r, warning=FALSE, message=FALSE}
library(dplyr)
best_auc <- RESULT %>%
  group_by(Iteration) %>%
  filter(AUC == max(AUC)) %>%
  ungroup()
BESTROUNDRESULTS<-best_auc
best_auc
```

#--------- BOOTSTRAP VALIDATION FOR EACH OF THE BEST PANEL -------------

```{r, warning=FALSE, message=FALSE}
# Number of bootstrap samples
num_bootstrap_samples <- 1000

# Initialize a data frame to store the bootstrap results
bootstrap_results_df <- data.frame(Bootstrap_Sample = integer(),
                                   Iteration = integer(),
                                   AUC = numeric(),
                                   Sensitivity = numeric(),
                                   Specificity = numeric(),
                                   F1_Score = numeric(),
                                   Accuracy = numeric(),
                                   stringsAsFactors = FALSE)

# Function to calculate metrics
calc_metrics <- function(predictions, true_labels) {
  # Ensure predictions and true_labels have the same levels
  predictions <- factor(predictions, levels = levels(true_labels))
  
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate Sensitivity, Specificity, and F1-Score
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Pos Pred Value"]
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # Calculate Accuracy
  accuracy <- cm$overall["Accuracy"]
  
  return(c(sensitivity, specificity, f1_score, accuracy))
}

# Function to perform bootstrap validation
# Function to perform bootstrap validation
bootstrap_validation <- function(iter_best_predictors,
                                 iter_best_model_params,
                                 data,
                                 outcome_var) {
  
  set.seed(123)  # For reproducibility
  
  iteration_bootstrap_results <- data.frame(
    Bootstrap_Sample = integer(),
    AUC = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    F1_Score = numeric(),
    Accuracy = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (b in 1:num_bootstrap_samples) {
    
    # Split the data into training (67%) and test (33%) sets
    train_indices <- sample(
      1:nrow(data),
      size = 0.67 * nrow(data)
    )
    
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    
    # Recode outcome in training data
    train_data$Gr <- as.numeric(train_data$Gr)
    train_data$Gr[train_data$Gr == 1] <- "X0"
    train_data$Gr[train_data$Gr == 2] <- "X1"
    train_data$Gr <- as.factor(train_data$Gr)
    
    
    # Recode outcome in test data
    test_data$Gr <- as.numeric(test_data$Gr)
    test_data$Gr[test_data$Gr == 1] <- "X0"
    test_data$Gr[test_data$Gr == 2] <- "X1"
    test_data$Gr <- as.factor(test_data$Gr)
    
    
    # Fit SVM using the best predictors and tuning parameters
    train_control <- trainControl(
      method = "none",
      classProbs = TRUE
    )
    
    iter_best_model <- train(
      as.formula(
        paste(
          outcome_var,
          "~",
          paste(iter_best_predictors, collapse = "+")
        )
      ),
      data = train_data,
      method = "svmRadial",
      trControl = train_control,
      tuneGrid = iter_best_model_params$tuneGrid
    )
    
    
    # Probability predictions
    prob_predictions <- predict(
      iter_best_model,
      newdata = test_data,
      type = "prob"
    )[, 2]
    
    
    # Class predictions
    predictions <- factor(
      ifelse(prob_predictions > 0.5, "X1", "X0"),
      levels = levels(test_data[[outcome_var]])
    )
    
    
    # Calculate AUC
    roc_curve <- roc(
      test_data[[outcome_var]],
      prob_predictions
    )
    
    auc_value <- auc(roc_curve)
    
    
    # Calculate classification metrics
    metrics <- calc_metrics(
      predictions,
      test_data[[outcome_var]]
    )
    
    
    # Store bootstrap results
    iteration_bootstrap_results <- rbind(
      iteration_bootstrap_results,
      data.frame(
        Bootstrap_Sample = b,
        AUC = as.numeric(auc_value),
        Sensitivity = as.numeric(metrics[1]),
        Specificity = as.numeric(metrics[2]),
        F1_Score = as.numeric(metrics[3]),
        Accuracy = as.numeric(metrics[4])
      )
    )
  }
  
  return(iteration_bootstrap_results)
}

# Loop for multiple iterations
for (iter in seq_along(results_list)) {
  iter_best_predictors <- results_list[[iter]]$best_selected_predictors
  iter_best_model_params <- results_list[[iter]]$best_model_params
  data <- MAIN
  
  # Perform bootstrap validation for the best model in this iteration
  iter_bootstrap_results <- bootstrap_validation(iter_best_predictors, iter_best_model_params, data, outcome_var = "Gr")
  iter_bootstrap_results$Iteration <- iter
  bootstrap_results_df <- rbind(bootstrap_results_df, iter_bootstrap_results)
}

BOOTSTRAP<-bootstrap_results_df
```

```{r, warning=FALSE, message=FALSE}
# Load necessary libraries
library(dplyr)

# Function to calculate summary statistics
calculate_summary_stats <- function(df, group_var, measure_vars) {
  df %>%
    group_by(!!sym(group_var)) %>%
    summarise(across(all_of(measure_vars), 
                     list(mean = ~mean(.), 
                          sd = ~sd(.), 
                          min = ~min(.), 
                          max = ~max(.), 
                          Q1 = ~quantile(., 0.25), 
                          median = ~median(.), 
                          Q3 = ~quantile(., 0.75)), 
                     .names = "{col}_{fn}"))
}

# Define the group variable and measure variables
group_var <- "Iteration"
measure_vars <- c("AUC", "Sensitivity", "Specificity", "F1_Score", "Accuracy")

# Calculate the summary statistics
summary_stats_df <- calculate_summary_stats(bootstrap_results_df, group_var, measure_vars)
# Load necessary libraries
library(tidyr)
library(dplyr)
# Reshape the data, keeping the Iteration column
reshaped_df <- summary_stats_df %>%
  pivot_longer(cols = -Iteration, 
               names_to = c("Metric", "Summary_Measure"), 
               names_sep = "_", 
               values_to = "Value") %>%
  select(Iteration, Metric, Summary_Measure, Value)
BOOTSTRAPSUMMARY<-reshaped_df
# View the result
summary_stats_df[,c(1, 2, 9, 16, 23, 30)]
```

```{r}
library(writexl)
list<-list("RoundsResult"= RESULT, "BestRoundResult"=BESTROUNDRESULTS, "Bootstrap"=BOOTSTRAP, "BootstrapSummary"=BOOTSTRAPSUMMARY)
write_xlsx(list, "E:\\FilteredVariableAnalysis\\SVM OUTPUT\\SVM 30082024\\SVMRFF.xlsx")
```


#--------------------------------------------- LASSO -----------------------------------------------

```{r, message=FALSE, warning=FALSE}
# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(PRROC)
library(readxl)
library(foreach)
library(doParallel)
```

```{r, message=FALSE, warning=FALSE}
Data <- read_excel("E:/FilteredVariableAnalysis/Data.xlsx", 
                   sheet = "DataLS")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)
```

```{r, warning=FALSE, message=FALSE}
a<-Sys.time()
a

set.seed(120)

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
      
      # Create formula for SVM
      formula_svm <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Train the SVM model with cross-validation
      svm_model <- train(formula_svm, data = train_set, method = "svmRadial", trControl = ctrl, metric = "ROC")
      # Make probability predictions on the test set
      prob_predictions <- predict(svm_model, newdata = test_set, type = "prob")[, 2]
      
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
        best_round_model <- svm_model
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
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_model_params = list(tuneGrid = best_model$bestTune), best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
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
Sys.time()-a
RESULT<-results_df
```


#---------- LIST OF BEST PANELS WITHIN EACH ITERATIONS ---------

```{r, warning=FALSE, message=FALSE}
library(dplyr)
best_auc <- RESULT %>%
  group_by(Iteration) %>%
  filter(AUC == max(AUC)) %>%
  ungroup()
BESTROUNDRESULTS<-best_auc
best_auc
```

#--------- BOOTSTRAP VALIDATION FOR EACH OF THE BEST PANEL -------------

```{r, warning=FALSE, message=FALSE}
# Number of bootstrap samples
num_bootstrap_samples <- 1000

# Initialize a data frame to store the bootstrap results
bootstrap_results_df <- data.frame(Bootstrap_Sample = integer(),
                                   Iteration = integer(),
                                   AUC = numeric(),
                                   Sensitivity = numeric(),
                                   Specificity = numeric(),
                                   F1_Score = numeric(),
                                   Accuracy = numeric(),
                                   stringsAsFactors = FALSE)

# Function to calculate metrics
calc_metrics <- function(predictions, true_labels) {
  # Ensure predictions and true_labels have the same levels
  predictions <- factor(predictions, levels = levels(true_labels))
  
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate Sensitivity, Specificity, and F1-Score
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Pos Pred Value"]
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # Calculate Accuracy
  accuracy <- cm$overall["Accuracy"]
  
  return(c(sensitivity, specificity, f1_score, accuracy))
}

# Function to perform bootstrap validation
# Function to perform bootstrap validation
bootstrap_validation <- function(iter_best_predictors,
                                 iter_best_model_params,
                                 data,
                                 outcome_var) {
  
  set.seed(123)  # For reproducibility
  
  iteration_bootstrap_results <- data.frame(
    Bootstrap_Sample = integer(),
    AUC = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    F1_Score = numeric(),
    Accuracy = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (b in 1:num_bootstrap_samples) {
    
    # Split the data into training (67%) and test (33%) sets
    train_indices <- sample(
      1:nrow(data),
      size = 0.67 * nrow(data)
    )
    
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    
    # Recode outcome in training data
    train_data$Gr <- as.numeric(train_data$Gr)
    train_data$Gr[train_data$Gr == 1] <- "X0"
    train_data$Gr[train_data$Gr == 2] <- "X1"
    train_data$Gr <- as.factor(train_data$Gr)
    
    
    # Recode outcome in test data
    test_data$Gr <- as.numeric(test_data$Gr)
    test_data$Gr[test_data$Gr == 1] <- "X0"
    test_data$Gr[test_data$Gr == 2] <- "X1"
    test_data$Gr <- as.factor(test_data$Gr)
    
    
    # Fit SVM using the best predictors and tuning parameters
    train_control <- trainControl(
      method = "none",
      classProbs = TRUE
    )
    
    iter_best_model <- train(
      as.formula(
        paste(
          outcome_var,
          "~",
          paste(iter_best_predictors, collapse = "+")
        )
      ),
      data = train_data,
      method = "svmRadial",
      trControl = train_control,
      tuneGrid = iter_best_model_params$tuneGrid
    )
    
    
    # Probability predictions
    prob_predictions <- predict(
      iter_best_model,
      newdata = test_data,
      type = "prob"
    )[, 2]
    
    
    # Class predictions
    predictions <- factor(
      ifelse(prob_predictions > 0.5, "X1", "X0"),
      levels = levels(test_data[[outcome_var]])
    )
    
    
    # Calculate AUC
    roc_curve <- roc(
      test_data[[outcome_var]],
      prob_predictions
    )
    
    auc_value <- auc(roc_curve)
    
    
    # Calculate classification metrics
    metrics <- calc_metrics(
      predictions,
      test_data[[outcome_var]]
    )
    
    
    # Store bootstrap results
    iteration_bootstrap_results <- rbind(
      iteration_bootstrap_results,
      data.frame(
        Bootstrap_Sample = b,
        AUC = as.numeric(auc_value),
        Sensitivity = as.numeric(metrics[1]),
        Specificity = as.numeric(metrics[2]),
        F1_Score = as.numeric(metrics[3]),
        Accuracy = as.numeric(metrics[4])
      )
    )
  }
  
  return(iteration_bootstrap_results)
}

# Loop for multiple iterations
for (iter in seq_along(results_list)) {
  iter_best_predictors <- results_list[[iter]]$best_selected_predictors
  iter_best_model_params <- results_list[[iter]]$best_model_params
  data <- MAIN
  
  # Perform bootstrap validation for the best model in this iteration
  iter_bootstrap_results <- bootstrap_validation(iter_best_predictors, iter_best_model_params, data, outcome_var = "Gr")
  iter_bootstrap_results$Iteration <- iter
  bootstrap_results_df <- rbind(bootstrap_results_df, iter_bootstrap_results)
}

BOOTSTRAP<-bootstrap_results_df
```

```{r, warning=FALSE, message=FALSE}
# Load necessary libraries
library(dplyr)

# Function to calculate summary statistics
calculate_summary_stats <- function(df, group_var, measure_vars) {
  df %>%
    group_by(!!sym(group_var)) %>%
    summarise(across(all_of(measure_vars), 
                     list(mean = ~mean(.), 
                          sd = ~sd(.), 
                          min = ~min(.), 
                          max = ~max(.), 
                          Q1 = ~quantile(., 0.25), 
                          median = ~median(.), 
                          Q3 = ~quantile(., 0.75)), 
                     .names = "{col}_{fn}"))
}

# Define the group variable and measure variables
group_var <- "Iteration"
measure_vars <- c("AUC", "Sensitivity", "Specificity", "F1_Score", "Accuracy")

# Calculate the summary statistics
summary_stats_df <- calculate_summary_stats(bootstrap_results_df, group_var, measure_vars)
# Load necessary libraries
library(tidyr)
library(dplyr)
# Reshape the data, keeping the Iteration column
reshaped_df <- summary_stats_df %>%
  pivot_longer(cols = -Iteration, 
               names_to = c("Metric", "Summary_Measure"), 
               names_sep = "_", 
               values_to = "Value") %>%
  select(Iteration, Metric, Summary_Measure, Value)
BOOTSTRAPSUMMARY<-reshaped_df
# View the result
summary_stats_df[,c(1, 2, 9, 16, 23, 30)]
```

```{r}
library(writexl)
list<-list("RoundsResult"= RESULT, "BestRoundResult"=BESTROUNDRESULTS, "Bootstrap"=BOOTSTRAP, "BootstrapSummary"=BOOTSTRAPSUMMARY)
write_xlsx(list, "E:\\FilteredVariableAnalysis\\SVM OUTPUT\\SVM 30082024\\SVMLASSO.xlsx")
```


#--------------------------------------------- BACKWARD ELIMINATION -----------------------------------------------

```{r, message=FALSE, warning=FALSE}
# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(PRROC)
library(readxl)
library(foreach)
library(doParallel)
```

```{r, message=FALSE, warning=FALSE}
Data <- read_excel("E:/FilteredVariableAnalysis/Data.xlsx", 
                   sheet = "DataBC")
MAIN<-Data
MAIN$Gr<-as.factor(MAIN$Gr)
```

```{r, warning=FALSE, message=FALSE}
a<-Sys.time()
a

set.seed(120)

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
      
      # Create formula for SVM
      formula_svm <- as.formula(paste(outcome_var, "~", paste(current_predictors, collapse = "+")))
      
      # Train the SVM model with cross-validation
      svm_model <- train(formula_svm, data = train_set, method = "svmRadial", trControl = ctrl, metric = "ROC")
      # Make probability predictions on the test set
      prob_predictions <- predict(svm_model, newdata = test_set, type = "prob")[, 2]
      
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
        best_round_model <- svm_model
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
  list(round_results = round_results, best_auc = best_auc, best_model = best_model, best_model_params = list(tuneGrid = best_model$bestTune), best_selected_predictors = best_selected_predictors, best_metrics = best_round_metrics, test_set = test_set)
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
Sys.time()-a
RESULT<-results_df
```


#---------- LIST OF BEST PANELS WITHIN EACH ITERATIONS ---------

```{r, warning=FALSE, message=FALSE}
library(dplyr)
best_auc <- RESULT %>%
  group_by(Iteration) %>%
  filter(AUC == max(AUC)) %>%
  ungroup()
BESTROUNDRESULTS<-best_auc
best_auc
```

#--------- BOOTSTRAP VALIDATION FOR EACH OF THE BEST PANEL -------------

```{r, warning=FALSE, message=FALSE}
# Number of bootstrap samples
num_bootstrap_samples <- 1000

# Initialize a data frame to store the bootstrap results
bootstrap_results_df <- data.frame(Bootstrap_Sample = integer(),
                                   Iteration = integer(),
                                   AUC = numeric(),
                                   Sensitivity = numeric(),
                                   Specificity = numeric(),
                                   F1_Score = numeric(),
                                   Accuracy = numeric(),
                                   stringsAsFactors = FALSE)

# Function to calculate metrics
calc_metrics <- function(predictions, true_labels) {
  # Ensure predictions and true_labels have the same levels
  predictions <- factor(predictions, levels = levels(true_labels))
  
  # Calculate confusion matrix
  cm <- confusionMatrix(predictions, true_labels)
  
  # Calculate Sensitivity, Specificity, and F1-Score
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Pos Pred Value"]
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  
  # Calculate Accuracy
  accuracy <- cm$overall["Accuracy"]
  
  return(c(sensitivity, specificity, f1_score, accuracy))
}

# Function to perform bootstrap validation
# Function to perform bootstrap validation
bootstrap_validation <- function(iter_best_predictors,
                                 iter_best_model_params,
                                 data,
                                 outcome_var) {
  
  set.seed(123)  # For reproducibility
  
  iteration_bootstrap_results <- data.frame(
    Bootstrap_Sample = integer(),
    AUC = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    F1_Score = numeric(),
    Accuracy = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (b in 1:num_bootstrap_samples) {
    
    # Split the data into training (67%) and test (33%) sets
    train_indices <- sample(
      1:nrow(data),
      size = 0.67 * nrow(data)
    )
    
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    
    # Recode outcome in training data
    train_data$Gr <- as.numeric(train_data$Gr)
    train_data$Gr[train_data$Gr == 1] <- "X0"
    train_data$Gr[train_data$Gr == 2] <- "X1"
    train_data$Gr <- as.factor(train_data$Gr)
    
    
    # Recode outcome in test data
    test_data$Gr <- as.numeric(test_data$Gr)
    test_data$Gr[test_data$Gr == 1] <- "X0"
    test_data$Gr[test_data$Gr == 2] <- "X1"
    test_data$Gr <- as.factor(test_data$Gr)
    
    
    # Fit SVM using the best predictors and tuning parameters
    train_control <- trainControl(
      method = "none",
      classProbs = TRUE
    )
    
    iter_best_model <- train(
      as.formula(
        paste(
          outcome_var,
          "~",
          paste(iter_best_predictors, collapse = "+")
        )
      ),
      data = train_data,
      method = "svmRadial",
      trControl = train_control,
      tuneGrid = iter_best_model_params$tuneGrid
    )
    
    
    # Probability predictions
    prob_predictions <- predict(
      iter_best_model,
      newdata = test_data,
      type = "prob"
    )[, 2]
    
    
    # Class predictions
    predictions <- factor(
      ifelse(prob_predictions > 0.5, "X1", "X0"),
      levels = levels(test_data[[outcome_var]])
    )
    
    
    # Calculate AUC
    roc_curve <- roc(
      test_data[[outcome_var]],
      prob_predictions
    )
    
    auc_value <- auc(roc_curve)
    
    
    # Calculate classification metrics
    metrics <- calc_metrics(
      predictions,
      test_data[[outcome_var]]
    )
    
    
    # Store bootstrap results
    iteration_bootstrap_results <- rbind(
      iteration_bootstrap_results,
      data.frame(
        Bootstrap_Sample = b,
        AUC = as.numeric(auc_value),
        Sensitivity = as.numeric(metrics[1]),
        Specificity = as.numeric(metrics[2]),
        F1_Score = as.numeric(metrics[3]),
        Accuracy = as.numeric(metrics[4])
      )
    )
  }
  
  return(iteration_bootstrap_results)
}

# Loop for multiple iterations
for (iter in seq_along(results_list)) {
  iter_best_predictors <- results_list[[iter]]$best_selected_predictors
  iter_best_model_params <- results_list[[iter]]$best_model_params
  data <- MAIN
  
  # Perform bootstrap validation for the best model in this iteration
  iter_bootstrap_results <- bootstrap_validation(iter_best_predictors, iter_best_model_params, data, outcome_var = "Gr")
  iter_bootstrap_results$Iteration <- iter
  bootstrap_results_df <- rbind(bootstrap_results_df, iter_bootstrap_results)
}

BOOTSTRAP<-bootstrap_results_df
```

```{r, warning=FALSE, message=FALSE}
# Load necessary libraries
library(dplyr)

# Function to calculate summary statistics
calculate_summary_stats <- function(df, group_var, measure_vars) {
  df %>%
    group_by(!!sym(group_var)) %>%
    summarise(across(all_of(measure_vars), 
                     list(mean = ~mean(.), 
                          sd = ~sd(.), 
                          min = ~min(.), 
                          max = ~max(.), 
                          Q1 = ~quantile(., 0.25), 
                          median = ~median(.), 
                          Q3 = ~quantile(., 0.75)), 
                     .names = "{col}_{fn}"))
}

# Define the group variable and measure variables
group_var <- "Iteration"
measure_vars <- c("AUC", "Sensitivity", "Specificity", "F1_Score", "Accuracy")

# Calculate the summary statistics
summary_stats_df <- calculate_summary_stats(bootstrap_results_df, group_var, measure_vars)
# Load necessary libraries
library(tidyr)
library(dplyr)
# Reshape the data, keeping the Iteration column
reshaped_df <- summary_stats_df %>%
  pivot_longer(cols = -Iteration, 
               names_to = c("Metric", "Summary_Measure"), 
               names_sep = "_", 
               values_to = "Value") %>%
  select(Iteration, Metric, Summary_Measure, Value)
BOOTSTRAPSUMMARY<-reshaped_df
# View the result
summary_stats_df[,c(1, 2, 9, 16, 23, 30)]
```

```{r}
library(writexl)
list<-list("RoundsResult"= RESULT, "BestRoundResult"=BESTROUNDRESULTS, "Bootstrap"=BOOTSTRAP, "BootstrapSummary"=BOOTSTRAPSUMMARY)
write_xlsx(list, "E:\\FilteredVariableAnalysis\\SVM OUTPUT\\SVM 30082024\\SVMBCEL.xlsx")
```

