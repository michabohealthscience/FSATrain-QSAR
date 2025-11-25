# ===============================================================
# QSAR Prediction Script for Unseen Molecules
# ===============================================================


# =============================
# 0. Packages (install & attach)
# =============================
required_pkgs <- c("rcdk", "fingerprint", "caret", "randomForest", "doParallel", "data.table", "dplyr")
for(pkg in required_pkgs){
  if(!pkg %in% installed.packages()[, "Package"]){
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# =============================
# 1. User settings
# =============================
train_csv <- "m1-pki-qsar-ready.csv"       # Training data: must contain SMILES + Activity
unseen_csv <- "list-of-tas_efsa-report.csv"       # Unseen data: must contain SMILES + Name
n_folds <- 5                                # CV folds
save_model_file <- "rf_maccs_model.rds"    # Save trained model
n_trees <- 500                              # Random Forest trees

# =============================
# 2. Load training data
# =============================
df <- data.table::fread(train_csv, data.table = FALSE)
if(!("SMILES" %in% names(df))) stop("Training data must contain 'SMILES'")
if(!("Activity" %in% names(df))) stop("Training data must contain 'Activity'")

df <- df[!is.na(df$SMILES) & !is.na(df$Activity), ]

# =============================
# 3. Parse SMILES -> molecules
# =============================
mols <- rcdk::parse.smiles(as.character(df$SMILES))
valid_idx <- sapply(mols, Negate(is.null))
if(any(!valid_idx)){
  cat("Dropping", sum(!valid_idx), "invalid SMILES rows\n")
  df <- df[valid_idx, , drop = FALSE]
  mols <- mols[valid_idx]
}
cat("Valid training molecules:", length(mols), "\n")

# =============================
# 4. Compute MACCS fingerprints
# =============================
fps <- lapply(mols, function(m) rcdk::get.fingerprint(m, type = "maccs"))
fp_mat <- fingerprint::fp.to.matrix(fps)
fp_df <- as.data.frame(fp_mat)
colnames(fp_df) <- paste0("MACCS_", seq_len(ncol(fp_df)))
rownames(fp_df) <- seq_len(nrow(fp_df))

# Combine features + target
train_data <- cbind(fp_df, Activity = df$Activity)

# =============================
# 5. Train Random Forest (caret)
# =============================
# Setup parallel
cores <- parallel::detectCores(logical = FALSE)
cl <- parallel::makePSOCKcluster(max(1, cores - 1))
doParallel::registerDoParallel(cl)

# Train control with 5-fold CV
fitControl <- caret::trainControl(method = "cv", number = n_folds,
                                  savePredictions = "final", allowParallel = TRUE)

# RF tuning: try a few mtry values
p <- ncol(fp_df)
tuneGrid <- expand.grid(mtry = unique(floor(c(sqrt(p), p/4, p/2))))

set.seed(123)
rf_cv <- caret::train(
  x = fp_df,
  y = df$Activity,
  method = "rf",
  trControl = fitControl,
  tuneGrid = tuneGrid,
  importance = TRUE,
  ntree = n_trees
)

# Stop parallel
parallel::stopCluster(cl)

cat("RF training completed. Best mtry:", rf_cv$bestTune$mtry, "\n")

# =============================
# 6. Save the model
# =============================
saveRDS(rf_cv, file = save_model_file)
cat("Saved RF model to", save_model_file, "\n")

# =============================
# 7. Load unseen data
# =============================
unseen <- data.table::fread(unseen_csv, data.table = FALSE)
if(!("SMILES" %in% names(unseen))) stop("Unseen data must contain 'SMILES'")
if(!("Name" %in% names(unseen))) stop("Unseen data must contain 'Name'")

# Parse molecules
mols_new <- rcdk::parse.smiles(as.character(unseen$SMILES))
valid_idx <- sapply(mols_new, Negate(is.null))
if(any(!valid_idx)){
  cat("Dropping", sum(!valid_idx), "invalid SMILES rows in unseen data\n")
}
unseen_valid <- unseen[valid_idx, , drop = FALSE]
mols_new <- mols_new[valid_idx]

# =============================
# 8. Compute MACCS fingerprints for unseen
# =============================
if(length(mols_new) == 0){
  cat("No valid molecules to predict.\n")
  pred_results <- data.frame()
} else {
  fps_new <- lapply(mols_new, function(m) rcdk::get.fingerprint(m, type = "maccs"))
  fp_mat_new <- fingerprint::fp.to.matrix(fps_new)
  fp_df_new <- as.data.frame(fp_mat_new)
  colnames(fp_df_new) <- paste0("MACCS_", seq_len(ncol(fp_df_new)))
  
  # Align columns with training set
  missing_cols <- setdiff(colnames(fp_df), colnames(fp_df_new))
  if(length(missing_cols) > 0){
    cat("Adding", length(missing_cols), "missing columns to unseen data with 0\n")
    for(col in missing_cols) fp_df_new[[col]] <- 0
  }
  fp_df_new <- fp_df_new[, colnames(fp_df), drop = FALSE]
  
  # Predict activity
  pred_values <- predict(rf_cv, newdata = fp_df_new)
  
  # Combine results
  pred_results <- data.frame(
    Name = unseen_valid$Name,
    SMILES = unseen_valid$SMILES,
    Predicted_Activity = pred_values
  )
}

# =============================
# 9. Save predictions
# =============================
if(nrow(pred_results) > 0){
  write.csv(pred_results, "predicted_activity_unseen.csv", row.names = FALSE)
  cat("Predictions saved to 'predicted_activity_unseen.csv'\n")
} else {
  cat("No predictions generated.\n")
}
