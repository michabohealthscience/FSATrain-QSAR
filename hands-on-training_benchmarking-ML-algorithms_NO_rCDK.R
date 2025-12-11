# =============================
# Step 0: Load required packages
# =============================
packages <- c("recipes", "caret", "randomForest", "gbm", "FNN", "ggplot2",
              "dplyr", "reshape2", "e1071", "iml")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =============================
# Step 1: Load and preprocess data
# =============================
# df <- read.csv("m1-pki-qsar-ready.csv")  # Replace with your CSV file
# 
# # Ensure SMILES and activity are not missing
# df <- df[!is.na(df$SMILES) & !is.na(df$Activity), ]
# 
# # Parse SMILES and compute fingerprints
# mols <- parse.smiles(as.character(df$SMILES))
# valid_idx <- sapply(mols, Negate(is.null))
# df <- df[valid_idx, ]
# mols <- mols[valid_idx]
# 
# fps <- lapply(mols, get.fingerprint, type = "extended")
# fp_mat <- as.data.frame(fp.to.matrix(fps))
# fp_mat$Activity <- df$Activity  # Add activity column
# 
# saveRDS(fp_mat, 'fp_mat.RDS')

fp_mat <- readRDS('pre_generated/fp_mat.RDS')

# =============================
# Step 2: Create Cross-Validation Folds
# =============================
# Remove NAs from activity and check sample count
if (!"Activity" %in% colnames(fp_mat)) stop("Activity column missing")
fp_mat <- fp_mat[!is.na(fp_mat$Activity), ]
if (nrow(fp_mat) < 10) stop("Not enough rows for 10-fold CV")

set.seed(123)
k <- min(5, nrow(fp_mat))
folds <- createFolds(fp_mat$Activity, k = k, list = TRUE, returnTrain = FALSE)

# =============================
# Step 3: Train and evaluate models
# =============================
models <- c("rf", "gbm", "knn")
results <- list()

for (model in models) {
  cat("Training model:", model, "\n")
  fold_results <- data.frame()
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_data <- fp_mat[-test_idx, ]
    test_data  <- fp_mat[test_idx, ]
    
    x_train <- train_data[, -ncol(train_data)]
    y_train <- train_data$Activity
    x_test <- test_data[, -ncol(test_data)]
    y_test <- test_data$Activity
    
    # Train model
    if (model == "rf") {
      fit <- randomForest(x = x_train, y = y_train)
      pred <- predict(fit, x_test)
    } else if (model == "gbm") {
      fit <- gbm(Activity ~ ., data = train_data, distribution = "gaussian", n.trees = 100)
      pred <- predict(fit, x_test, n.trees = 100)
    } else if (model == "knn") {
      pred <- knn.reg(train = x_train, test = x_test, y = y_train, k = 5)$pred
    }
    
    # Metrics
    rmse <- RMSE(pred, y_test)
    mae <- MAE(pred, y_test)
    rsq <- cor(pred, y_test)^2
    
    fold_results <- rbind(fold_results, data.frame(Model = model, Fold = i, RMSE = rmse, MAE = mae, Rsquared = rsq))
    
    # Store residuals and model from first fold only
    if (i == 1) {
      results[[model]] <- list(
        Residuals = data.frame(Observed = y_test, Predicted = pred, Residual = y_test - pred),
        Model = fit
      )
    }
  }
  results[[model]]$Metrics <- fold_results
}

# =============================
# Step 4: Plot model metrics
# =============================
all_metrics <- do.call(rbind, lapply(results, function(x) x$Metrics))
melted_metrics <- melt(all_metrics, id.vars = c("Model", "Fold"))

ggplot(melted_metrics, aes(x = Model, y = value, fill = Model)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Model Performance (Cross-Validation)", y = "Metric Value", x = "") +
  theme_minimal()

# =============================
# Step 5: Residual Plot of Best Model
# =============================
avg_rmse <- sapply(results, function(x) mean(x$Metrics$RMSE))
best_model <- names(which.min(avg_rmse))
cat("Best model selected:", best_model, "\n")

residual_df <- results[[best_model]]$Residuals
ggplot(residual_df, aes(x = Predicted, y = Residual)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = paste("Residual Plot -", best_model), x = "Predicted", y = "Residual")

# =============================
# Step 6: Feature Importance Plot (Top 15)
# =============================
if (best_model == "rf") {
  imp <- importance(results[[best_model]]$Model)
  importance_df <- data.frame(Feature = rownames(imp), Importance = imp[, 1])
} else if (best_model == "gbm") {
  imp <- summary(results[[best_model]]$Model, plotit = FALSE)
  importance_df <- data.frame(Feature = imp$var, Importance = imp$rel.inf)
} else {
  importance_df <- data.frame(Feature = colnames(fp_mat)[-ncol(fp_mat)], Importance = rep(NA, ncol(fp_mat)-1))
}

# Top N features for clear display
top_n <- 15
importance_df <- importance_df[order(-importance_df$Importance), ][1:min(top_n, nrow(importance_df)), ]

ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = paste("Top", top_n, "Important Features -", best_model),
       x = "Feature", y = "Importance")

