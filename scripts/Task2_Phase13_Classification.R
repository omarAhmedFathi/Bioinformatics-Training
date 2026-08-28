# ==========================================
# Task 2 - Phase 13: Diabetes Classification Using Metabolites
# ==========================================

cat("=== Phase 13: Diabetes Classification Using Metabolites ===\n\n")

if (!require("caret", quietly=TRUE)) install.packages("caret")
if (!require("randomForest", quietly=TRUE)) install.packages("randomForest")
if (!require("pROC", quietly=TRUE)) install.packages("pROC")
if (!require("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!require("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!require("glmnet", quietly=TRUE)) install.packages("glmnet")

library(caret)
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)
library(glmnet)

# 1. Load Data
cat("Loading datasets...\n")
qc_data <- read.csv("../intermediate/Phase11_Metabolites_QC_Final.csv", stringsAsFactors = FALSE)

# We will build a Metabolite-Only classifier to specifically evaluate the predictive power of the metabolites.
cat("Model Type: Metabolite-Only Classifier\n")

# Format for ML
ml_data <- qc_data %>% select(-main_id, -mapped_id)
# Ensure Diabetes is a factor with valid R variable names for caret
ml_data$Diabetes <- factor(ml_data$Diabetes, levels = c(0, 1), labels = c("Control", "Diabetic"))

# 2. Stratified Train/Test Split (80/20)
cat("Splitting data into 80% train and 20% test sets...\n")
set.seed(123)
train_idx <- createDataPartition(ml_data$Diabetes, p = 0.8, list = FALSE)
train_set <- ml_data[train_idx, ]
test_set <- ml_data[-train_idx, ]

cat(sprintf("Training samples: %d (Control: %d, Diabetic: %d)\n", 
            nrow(train_set), sum(train_set$Diabetes == "Control"), sum(train_set$Diabetes == "Diabetic")))
cat(sprintf("Testing samples: %d (Control: %d, Diabetic: %d)\n", 
            nrow(test_set), sum(test_set$Diabetes == "Control"), sum(test_set$Diabetes == "Diabetic")))

# 3. Train Classification Models
cat("\nTraining Logistic Regression (Elastic Net) model...\n")
# Using glmnet to avoid convergence issues with many predictors
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

set.seed(123)
lr_model <- train(Diabetes ~ ., data = train_set, 
                  method = "glmnet", 
                  trControl = fitControl, 
                  metric = "ROC",
                  tuneLength = 5)

cat("Training Random Forest model...\n")
set.seed(123)
rf_model <- train(Diabetes ~ ., data = train_set, 
                  method = "rf", 
                  trControl = fitControl, 
                  metric = "ROC",
                  tuneLength = 5,
                  importance = TRUE)

# 4. Evaluate Performance on Test Set
cat("\nEvaluating models on the test set...\n")

evaluate_model <- function(model, test_data, model_name) {
  preds <- predict(model, test_data)
  probs <- predict(model, test_data, type = "prob")
  
  cm <- confusionMatrix(preds, test_data$Diabetes, positive = "Diabetic")
  
  roc_curve <- roc(test_data$Diabetes, probs$Diabetic, levels = c("Control", "Diabetic"))
  auc_val <- as.numeric(auc(roc_curve))
  
  metrics <- data.frame(
    Model = model_name,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"],
    AUC = auc_val,
    stringsAsFactors = FALSE
  )
  
  return(list(metrics = metrics, preds = data.frame(Sample_Index = rownames(test_data), True_Class = test_data$Diabetes, Predicted_Class = preds, Prob_Diabetic = probs$Diabetic)))
}

lr_eval <- evaluate_model(lr_model, test_set, "Logistic Regression (Elastic Net)")
rf_eval <- evaluate_model(rf_model, test_set, "Random Forest")

# Combine Metrics
metrics_table <- rbind(lr_eval$metrics, rf_eval$metrics)
rownames(metrics_table) <- NULL
write.csv(metrics_table, "../outputs/Phase13_Classification/Phase13_Classification_Metrics.csv", row.names = FALSE)

# Save Test Predictions
test_predictions <- rbind(
  cbind(Model = "Logistic Regression", lr_eval$preds),
  cbind(Model = "Random Forest", rf_eval$preds)
)
write.csv(test_predictions, "../intermediate/Phase13_TestSet_Predictions.csv", row.names = FALSE)

# 5. Extract Feature Importance (from Random Forest)
cat("\nExtracting feature importance...\n")
var_imp <- varImp(rf_model, scale = TRUE)
imp_matrix <- var_imp$importance

if ("Overall" %in% colnames(imp_matrix)) {
  imp_vals <- imp_matrix$Overall
} else if ("Diabetic" %in% colnames(imp_matrix)) {
  imp_vals <- imp_matrix$Diabetic
} else {
  imp_vals <- imp_matrix[,1]
}

imp_df <- data.frame(Metabolite = rownames(imp_matrix), Importance = imp_vals)
imp_df <- imp_df %>% arrange(desc(Importance))

write.csv(imp_df, "../outputs/Phase13_Classification/Phase13_Variable_Importance.csv", row.names = FALSE)

# Plot top 20 important metabolites
top_20_imp <- head(imp_df, 20)

pdf("../outputs/Phase13_Classification/Phase13_Variable_Importance_Plot.pdf", width=10, height=8)
print(
  ggplot(top_20_imp, aes(x = reorder(Metabolite, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Top 20 Predictive Metabolites (Random Forest)",
      x = "Metabolite",
      y = "Scaled Importance"
    )
)
invisible(dev.off())

# Save list of all metabolites used
write.table(imp_df$Metabolite, "../intermediate/Phase13_Predictor_Metabolites.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\nPhase 13 complete. Outputs saved to data/:\n")
cat("  - Phase13_Classification_Metrics.csv\n")
cat("  - Phase13_Variable_Importance.csv\n")
cat("  - Phase13_Variable_Importance_Plot.pdf\n")
cat("  - intermediate_files/Phase13_TestSet_Predictions.csv\n")
cat("  - intermediate_files/Phase13_Predictor_Metabolites.txt\n")
