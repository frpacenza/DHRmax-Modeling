# ---------------------------------------------------------------------------- #
# FUNCTIONS
# ---------------------------------------------------------------------------- #
# Feature-selection methods

filtration <- function(order, cors, n, th = 0.85) {
  # ============================================================================
  # GREEDY FILTRATION ALGORITHM FOR ORDERED VECTORS
  #
  # Selects a subset of elements from an ordered named vector based on a 
  # maximum correlation threshold relative to already selected elements.
  #
  # @param order : Numeric named vector, sorted in descending order.
  # @param cors  : Symmetric correlation matrix with matching row/column names.
  # @param n     : Positive integer, max number of elements to select.
  # @param th    : Max allowed correlation threshold. Default is 0.85.
  #
  # @return      : A filtered subset (named vector) of 'order'.
  # ============================================================================
  
  # Initialization: take the name of the first element in 'order'
  filtered_names <- names(order)[1]
  
  # Loop through the remaining elements of 'order'
  for (i in 2:length(order)) {
    
    # Early exit check: stop if we have already reached 'n' elements
    if (length(filtered_names) == n) {
      break
    }
    
    current_candidate <- names(order)[i]
    
    # Extract correlations between the current candidate and ALL previously selected elements
    current_cors <- cors[current_candidate, filtered_names]
    
    # Condition: candidate must have a correlation < th with ALL selected elements
    if (all(current_cors < th)) {
      filtered_names <- c(filtered_names, current_candidate)
    }
  }
  
  # Extract the final named vector using the selected names
  filtered <- order[filtered_names]
  
  return(filtered)
}

corr_fs <- function(predictors, target, n, th = 0.85) {
  # ============================================================================
  # CORRELATION-BASED FEATURE SELECTION WITH REDUNDANCY FILTRATION
  #
  # Selects the top 'n' predictors by first ranking them based on their absolute
  # correlation with a target variable, and then filtering out redundant features
  # using a greedy filtration algorithm.
  # It assumes no NA values are present in predictors and target.
  #
  # @param predictors : Dataframe or matrix containing the predictor columns.
  # @param target     : Numeric vector representing the target/response variable.
  # @param n          : Positive integer, max number of features to select.
  # @param th         : Max allowed correlation between features. Default is 0.85.
  #
  # @return           : A named vector of the selected features and their 
  #                     absolute correlation with the target.
  # ============================================================================
  
  # --- STEP 1: Rank predictors by absolute correlation with target ---
  
  # Calculate correlation of each column in 'predictors' with 'target'
  # cor() returns a matrix/vector; as.vector() converts it to a standard vector
  raw_correlations <- as.vector(cor(predictors, target))
  
  # Assign column names to the correlation vector
  names(raw_correlations) <- colnames(predictors)
  
  # Take the absolute value (we care about the strength of the relationship, not the direction)
  abs_correlations <- abs(raw_correlations)
  
  # Sort in descending order to obtain the named vector "order"
  order_vector <- sort(abs_correlations, decreasing = TRUE)
  
  
  # --- STEP 2: Compute feature-to-feature correlation matrix ---
  
  # Calculate the symmetric correlation matrix among all predictors
  # Take the absolute value since the filtration is based on magnitude
  cors_matrix <- abs(cor(predictors))
  
  
  # --- STEP 3: Apply the filtration function ---
  
  # Call the previously defined function to get the filtered subset
  selected_features <- filtration(order = order_vector, cors = cors_matrix, n = n, th = th)
  
  return(selected_features)
}

enet_fs <- function(predictors, target, n, th = 0.85) {
  # ============================================================================
  # GLMNET-BASED FEATURE SELECTION WITH REDUNDANCY FILTRATION
  #
  # Selects the top 'n' predictors by ranking them based on the absolute beta 
  # coefficients of a GLMNET Elastic Net model tuned via caret.
  # It assumes no NA values are present in predictors and target.
  # TUNING DETAILS:
  # Method: Repeated 5-fold Cross-Validation (50 repeats) to minimize RMSE.
  # Grid Search: Optimizes '.alpha' (Lasso/Ridge mix: 0.1 to 1.0, by 0.1) and 
  # '.lambda' (Total penalty: 0 to 1, log-spaced standard values).
  #
  # @param predictors : Dataframe or matrix containing the predictor columns.
  # @param target     : Numeric vector representing the target/response variable.
  # @param n          : Positive integer, max number of features to select.
  # @param th         : Max allowed correlation between features. Default is 0.85.
  #
  # @return           : A named vector of the selected features and their 
  #                     absolute beta coefficients.
  # ============================================================================
  
  # Ensure required packages are available
  if (!requireNamespace("caret", "quietly" = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  if (!requireNamespace("glmnet", "quietly" = TRUE)) {
    stop("Package 'glmnet' is required but not installed.")
  }
  
  # --- STEP 0: Hyperparameter Tuning via caret ---
  
  # Define the parameter grid for glmnet
  # alpha = 0.1 (mostly Ridge) to 1.0 (pure Lasso)
  # lambda = controls the overall strength of the regularization penalty
  enet_grid <- expand.grid(
    .alpha = seq(0.1, 1.0, by = 0.1),
    .lambda = c(0, 0.0001, 0.001, 0.01, 0.1, 1)
  )
  
  # Set up the repeated cross-validation control
  set.seed(42)
  fitControl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50
  )
  
  # Train the Elastic Net models using "glmnet"
  # Note: glmnet strictly requires a matrix input for x, hence as.matrix()
  enet_model <- caret::train(
    x = as.matrix(predictors), 
    y = target,
    method = "glmnet",
    tuneGrid = enet_grid,
    trControl = fitControl
  )
  
  
  # --- STEP 1: Extract Beta Coefficients ---
  
  # Extract coefficients using glmnet's coef() function applied to the final model.
  # We must supply the optimal lambda value selected during the tuning process.
  best_lambda <- enet_model$bestTune$lambda
  raw_coefs <- coef(enet_model$finalModel, s = best_lambda)
  
  # raw_coefs is returned as a sparse matrix (dgCMatrix class).
  # Convert it to a standard numeric vector and map the feature names.
  beta_coefs <- as.vector(raw_coefs)
  names(beta_coefs) <- rownames(raw_coefs)
  
  # Remove the Intercept term, which glmnet explicitly labels as "(Intercept)"
  beta_coefs <- beta_coefs[names(beta_coefs) != "(Intercept)"]
  
  # Take the absolute value for ranking
  abs_betas <- abs(beta_coefs)
  
  # Sort in descending order to obtain the named vector "order"
  order_vector <- sort(abs_betas, decreasing = TRUE)
  
  
  # --- STEP 2: Compute feature-to-feature correlation matrix ---
  
  # Calculate the symmetric correlation matrix among all predictors
  # Take the absolute value since the filtration is based on magnitude
  cors_matrix <- abs(cor(predictors))
  
  
  # --- STEP 3: Apply the filtration function ---
  
  # Call the previously defined function to get the filtered subset
  selected_features <- filtration(order = order_vector, cors = cors_matrix, n = n, th = th)
  
  return(selected_features)
}

lasso_fs <- function(predictors, target, n, th = 0.85) {
  # ============================================================================
  # LASSO-BASED FEATURE SELECTION WITH REDUNDANCY FILTRATION
  #
  # Selects the top 'n' predictors by ranking them based on the absolute beta 
  # coefficients of a pure LASSO model (alpha = 1) tuned via caret.
  # It assumes no NA values are present in predictors and target.
  # TUNING DETAILS:
  # Method: Repeated 5-fold Cross-Validation (50 repeats) to minimize RMSE.
  # Grid Search: Forces '.alpha' = 1 (Pure Lasso) and optimizes '.lambda' 
  # over a sequence of regularization penalty values.
  #
  # @param predictors : Dataframe or matrix containing the predictor columns.
  # @param target     : Numeric vector representing the target/response variable.
  # @param n          : Positive integer, max number of features to select.
  # @param th         : Max allowed correlation between features. Default is 0.85.
  #
  # @return           : A named vector of the selected features and their 
  #                     absolute beta coefficients.
  # ============================================================================
  
  # Ensure required packages are available
  if (!requireNamespace("caret", "quietly" = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  if (!requireNamespace("glmnet", "quietly" = TRUE)) {
    stop("Package 'glmnet' is required but not installed.")
  }
  
  # --- STEP 0: Hyperparameter Tuning via caret ---
  
  # Define the parameter grid for glmnet
  # alpha = 0.1 (mostly Ridge) to 1.0 (pure Lasso)
  # lambda = controls the overall strength of the regularization penalty
  lasso_grid <- expand.grid(
    .alpha = 1,
    .lambda = c(0, 0.0001, 0.001, 0.01, 0.1, 1)
  )
  
  # Set up the repeated cross-validation control
  set.seed(42)
  fitControl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50
  )
  
  # Train the Lasso models using "glmnet"
  # Note: glmnet strictly requires a matrix input for x, hence as.matrix()
  lasso_model <- caret::train(
    x = as.matrix(predictors), 
    y = target,
    method = "glmnet",
    tuneGrid = lasso_grid,
    trControl = fitControl
  )
  
  
  # --- STEP 1: Extract Beta Coefficients ---
  
  # Extract coefficients using glmnet's coef() function applied to the final model.
  # We supply the optimal lambda value selected during the tuning process.
  best_lambda <- lasso_model$bestTune$lambda
  raw_coefs <- coef(lasso_model$finalModel, s = best_lambda)
  
  # raw_coefs is returned as a sparse matrix (dgCMatrix class).
  # Convert it to a standard numeric vector and map the feature names.
  beta_coefs <- as.vector(raw_coefs)
  names(beta_coefs) <- rownames(raw_coefs)
  
  # Remove the Intercept term, which glmnet explicitly labels as "(Intercept)"
  beta_coefs <- beta_coefs[names(beta_coefs) != "(Intercept)"]
  
  # Take the absolute value for ranking
  abs_betas <- abs(beta_coefs)
  
  # Sort in descending order to obtain the named vector "order"
  order_vector <- sort(abs_betas, decreasing = TRUE)
  
  
  # --- STEP 2: Compute feature-to-feature correlation matrix ---
  
  # Calculate the symmetric correlation matrix among all predictors
  # Take the absolute value since the filtration is based on magnitude
  cors_matrix <- abs(cor(predictors))
  
  
  # --- STEP 3: Apply the filtration function ---
  
  # Call the previously defined function to get the filtered subset
  selected_features <- filtration(order = order_vector, cors = cors_matrix, n = n, th = th)
  
  return(selected_features)
}

rforest_fs <- function(predictors, target, n, th = 0.85) {
  # ============================================================================
  # RANDOM FOREST-BASED FEATURE SELECTION WITH REDUNDANCY FILTRATION
  #
  # Selects the top 'n' predictors by ranking them based on the Mean Decrease 
  # in Impurity (Variable Importance) of a Random Forest model tuned via caret 
  # (using 'rf' method). It assumes no NA values are present in predictors and target.
  #
  # TUNING DETAILS:
  # - Method: Repeated 5-fold Cross-Validation (50 repeats) to minimize RMSE.
  # - Grid Search: Optimizes 'mtry' (number of variables randomly sampled as 
  #   candidates at each split) ranging from 1 to the total number of predictors.
  #
  # @param predictors : Dataframe or matrix containing the predictor columns.
  # @param target     : Numeric vector representing the target/response variable.
  # @param n          : Positive integer, max number of features to select.
  # @param th         : Max allowed correlation between features. Default is 0.85.
  #
  # @return           : A named vector of the selected features and their 
  #                     variable importance scores.
  # ============================================================================
  
  # Ensure required packages are available
  if (!requireNamespace("caret", "quietly" = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  if (!requireNamespace("randomForest", "quietly" = TRUE)) {
    stop("Package 'randomForest' is required but not installed.")
  }
  
  # --- STEP 0: Hyperparameter Tuning via caret ---
  
  # Define the parameter grid for mtry (from 1 up to the total number of predictors)
  rf_grid <- expand.grid(
    .mtry = seq(1, ncol(predictors), by = 1)
  )
  
  # Set up the repeated cross-validation control
  set.seed(42)
  fitControl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10
  )
  
  # Train the Random Forest models to tune the hyperparameter mtry
  rf_model <- caret::train(
    x = predictors, 
    y = target,
    method = "rf",
    tuneGrid = rf_grid,
    trControl = fitControl,
    importance = TRUE # Required to compute variable importance scores
  )
  
  # --- STEP 1: Extract Variable Importance ---
  
  # Extract variable importance from the final optimal model using caret's helper
  raw_importance <- caret::varImp(rf_model, scale = FALSE)
  
  # Extract the importance vector from the caret object
  importance_scores <- raw_importance$importance[, 1]
  
  # Ensure variable names are correctly preserved
  names(importance_scores) <- rownames(raw_importance$importance)
  
  # Take the absolute value (importance is naturally positive, kept for safety/consistency)
  abs_importance <- abs(importance_scores)
  
  # Sort in descending order to obtain the named vector "order"
  order_vector <- sort(abs_importance, decreasing = TRUE)
  
  
  # --- STEP 2: Compute feature-to-feature correlation matrix ---
  
  # Calculate the symmetric correlation matrix among all predictors
  # Take the absolute value since the filtration is based on magnitude
  cors_matrix <- abs(cor(predictors))
  
  
  # --- STEP 3: Apply the filtration function ---
  
  # Call the previously defined function to get the filtered subset
  selected_features <- filtration(order = order_vector, cors = cors_matrix, n = n, th = th)
  
  return(selected_features)
}

# ---------------------------------------------------------------------------- #
# Linear models assessment and construction

data_lm <- function(predictors, target, n, th = 0.85) {
  # ============================================================================
  # BENCHMARKING LINEAR MODELS BASED ON DIFFERENT FEATURE SELECTION METHODS
  #
  # 1. Applies 'corr_fs', 'enet_fs', and 'rforest_fs' to select features.
  # 2. Fits a final Linear Model using 100% of the data to extract Beta & p-values.
  # 3. Evaluates the feature subsets via Repeated CV, tracking out-of-fold RMSE and R2.
  #
  # CV DETAILS:
  # - Method: Repeated 5-fold Cross-Validation (50 repeats = 250 total splits).
  # - Metrics: Computed exclusively on out-of-fold validation test sets.
  #
  # @param predictors : Dataframe or matrix containing the predictor columns.
  # @param target     : Numeric vector representing the target/response variable.
  # @param n          : Positive integer, max number of features to select.
  # @param th         : Max allowed correlation between features. Default is 0.85.
  #
  # @return           : A data.frame summarizing the features, Final Betas (with p-values),
  #                     CV-RMSE, and CV-R2 with grouped rows for readability.
  # ============================================================================
  
  # Ensure required packages are available
  if (!requireNamespace("caret", "quietly" = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  
  # Map internal function names to their structured short names for the model ID
  method_mapping <- list(
    "corr_fs" = "corr",
    "enet_fs"        = "enet",
    "rforest_fs"     = "rforest"
  )
  
  # --- STEP 1: Apply Feature Selection Methods ---
  fs_results <- list(
    "corr_fs" = names(corr_fs(predictors, target, n = n, th = th)),
    "enet_fs"        = names(enet_fs(predictors, target, n = n, th = th)),
    "rforest_fs"     = names(rforest_fs(predictors, target, n = n, th = th))
  )
  
  # Set up CV Control (5-fold, 50 repeats = 250 total out-of-fold test evaluations)
  set.seed(42)
  cv_control <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50
  )
  
  # Initialize an empty list to collect data rows for the final table
  output_rows <- list()
  
  # --- STEP 2 & 3: Train LM, Extract Final Coefficients and CV Metrics ---
  
  for (method_name in names(fs_results)) {
    selected_vars <- fs_results[[method_name]]
    short_name <- method_mapping[[method_name]]
    model_id <- paste0(short_name, "_lm_", n)
    
    if (length(selected_vars) == 0) {
      output_rows[[method_name]] <- data.frame(
        `model name` = model_id,
        n = n,
        `feature-selection method` = method_name,
        `selected variables` = "None Selected",
        `final model beta (p-val)` = "N/A",
        `cv-RMSE (Mean pm SD)` = "N/A",
        `cv-R2 (Mean pm SD)` = "N/A",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Subset the predictors to include only the selected features
    x_subset <- predictors[, selected_vars, drop = FALSE]
    
    # Fit Linear Model via caret
    set.seed(42)
    lm_cv <- caret::train(
      x = x_subset,
      y = target,
      method = "lm",
      trControl = cv_control
    )
    
    # Extract CV Metrics
    rmse_mean <- mean(lm_cv$resample$RMSE)
    rmse_sd   <- sd(lm_cv$resample$RMSE)
    rmse_string <- paste0(round(rmse_mean, 4), " +/- ", round(rmse_sd, 4))
    
    r2_mean   <- mean(lm_cv$resample$Rsquared)
    r2_sd     <- sd(lm_cv$resample$Rsquared)
    r2_string <- paste0(round(r2_mean, 4), " +/- ", round(r2_sd, 4))
    
    # Extract Final Model Coefficients and P-values
    final_summary <- summary(lm_cv$finalModel)
    coef_matrix   <- final_summary$coefficients
    
    # Flag to identify the first row of the current model block
    is_first_row <- TRUE
    
    # Map over each selected variable to extract its beta and build the significance string
    for (var in selected_vars) {
      beta_val <- coef_matrix[var, "Estimate"]
      p_val    <- coef_matrix[var, "Pr(>|t|)"]
      p_str <- if(p_val < 0.001) "$<$ 0.001" else sprintf("%.3f", p_val)
      
      # Determine significance stars
      stars <- ""
      if (p_val < 0.001) { stars <- " ***" }
      else if (p_val < 0.01)  { stars <- " **" }
      else if (p_val < 0.05)  { stars <- " *" }
      else if (p_val < 0.1)   { stars <- " ." }
      
      # Build the final string: Beta (p-value *stars*)
      beta_string <- paste0(round(beta_val, 4), " (", p_str, stars, ")")
      
      # For rows after the first one within the same model, we blank out the global metrics
      row_model_id   <- if (is_first_row) model_id else ""
      row_n          <- if (is_first_row) as.character(n) else ""
      row_method     <- if (is_first_row) method_name else ""
      row_rmse       <- if (is_first_row) rmse_string else ""
      row_r2         <- if (is_first_row) r2_string else ""
      
      # Append row data to the list
      output_rows[[length(output_rows) + 1]] <- data.frame(
        `model name` = row_model_id,
        n = row_n,
        `feature-selection method` = row_method,
        `selected variables` = var,
        `final model beta (p-val)` = beta_string,
        `cv-RMSE (Mean pm SD)` = row_rmse,
        `cv-R2 (Mean pm SD)` = row_r2,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      # After the first iteration, switch the flag off for this model's remaining variables
      is_first_row <- FALSE
    }
  }
  
  # Combine all collected rows into a single structured data.frame
  final_table <- do.call(rbind, output_rows)
  rownames(final_table) <- NULL
  
  return(final_table)
}

physio_lm <- function(physiology_predictors, target, n) {
  # ============================================================================
  # BENCHMARKING PHYSIOLOGY-DRIVEN LINEAR MODELS VIA ALL POSSIBLE COMBINATIONS
  #
  # 1. Generates all possible combinations of 'n' predictors from the input set.
  # 2. Fits a final Linear Model using 100% of the data to extract Beta & p-values.
  # 3. Evaluates each model via Repeated CV (5-fold, 50 repeats).
  #
  # CV DETAILS:
  # - Method: Repeated 5-fold Cross-Validation (50 repeats = 250 total splits).
  # - Metrics: Computed exclusively on out-of-fold validation test sets.
  #
  # @param physiology_predictors : Dataframe or matrix containing ONLY the N 
  #                                physiological candidate columns.
  # @param target                : Numeric vector representing the target variable.
  # @param n                     : Positive integer, exact number of features per model.
  #
  # @return                      : A data.frame summarizing the features, Final Betas 
  #                                (with p-values), CV-RMSE, and CV-R2.
  # ============================================================================
  
  # Ensure required packages are available
  if (!requireNamespace("caret", "quietly" = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  
  # --- STEP 1: Generate All Combinations of size 'n' ---
  
  all_var_names <- colnames(physiology_predictors)
  N_total <- length(all_var_names)
  
  if (n > N_total) {
    stop("Requested 'n' cannot be greater than the total number of physiological predictors available.")
  }
  
  # Generate a matrix where each column is a combination of n variables
  comb_matrix <- combn(all_var_names, n)
  num_combinations <- ncol(comb_matrix)
  
  # Set up CV Control (5-fold, 50 repeats = 250 total out-of-fold test evaluations)
  set.seed(42)
  cv_control <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50
  )
  
  # Initialize an empty list to collect data rows for the final table
  output_rows <- list()
  
  # --- STEP 2 & 3: Train LM, Extract Final Coefficients and CV Metrics ---
  
  for (j in 1:num_combinations) {
    selected_vars <- comb_matrix[, j]
    
    # Unique model ID for this combination
    model_id <- paste0("physio_lm_", n, "_c", j)
    
    # Subset the predictors to include only the current combination
    x_subset <- physiology_predictors[, selected_vars, drop = FALSE]
    
    # Fit Linear Model via caret
    lm_cv <- caret::train(
      x = x_subset,
      y = target,
      method = "lm",
      trControl = cv_control
    )
    
    # Extract CV Metrics
    rmse_mean <- mean(lm_cv$resample$RMSE)
    rmse_sd   <- sd(lm_cv$resample$RMSE)
    rmse_string <- paste0(round(rmse_mean, 4), " +/- ", round(rmse_sd, 4))
    
    r2_mean   <- mean(lm_cv$resample$Rsquared)
    r2_sd     <- sd(lm_cv$resample$Rsquared)
    r2_string <- paste0(round(r2_mean, 4), " +/- ", round(r2_sd, 4))
    
    # Extract Final Model Coefficients and P-values
    final_summary <- summary(lm_cv$finalModel)
    coef_matrix   <- final_summary$coefficients
    
    # Flag to identify the first row of the current model block
    is_first_row <- TRUE
    
    # Map over each selected variable to extract its beta and build the significance string
    for (var in selected_vars) {
      beta_val <- coef_matrix[var, "Estimate"]
      p_val    <- coef_matrix[var, "Pr(>|t|)"]
      p_str <- if(p_val < 0.001) "$<$ 0.001" else sprintf("%.3f", p_val)
      
      # Determine significance stars
      stars <- ""
      if (p_val < 0.001) { stars <- " ***" }
      else if (p_val < 0.01)  { stars <- " **" }
      else if (p_val < 0.05)  { stars <- " *" }
      else if (p_val < 0.1)   { stars <- " ." }
      
      # Build the final string: Beta (p-value *stars*)
      beta_string <- paste0(round(beta_val, 4), " (", p_str, stars, ")")
      
      # For rows after the first one within the same model, we blank out the global metrics
      row_model_id   <- if (is_first_row) model_id else ""
      row_n          <- if (is_first_row) as.character(n) else ""
      row_method     <- if (is_first_row) "physiology_driven" else ""
      row_rmse       <- if (is_first_row) rmse_string else ""
      row_r2         <- if (is_first_row) r2_string else ""
      
      # Append row data to the list
      output_rows[[length(output_rows) + 1]] <- data.frame(
        `model name` = row_model_id,
        n = row_n,
        `feature-selection method` = row_method,
        `selected variables` = var,
        `final model beta (p-val)` = beta_string,
        `cv-RMSE (Mean pm SD)` = row_rmse,
        `cv-R2 (Mean pm SD)` = row_r2,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      # After the first iteration, switch the flag off for this model's remaining variables
      is_first_row <- FALSE
    }
  }
  
  # Combine all collected rows into a single structured data.frame
  final_table <- do.call(rbind, output_rows)
  rownames(final_table) <- NULL
  
  return(final_table)
}

# ---------------------------------------------------------------------------- #
# Latex code

to_latex_table <- function(data_df, physio_df, n) {
  # ============================================================================
  # AUTOMATED LATEX TABLE GENERATOR - WITH AUTO-COPY TO CLIPBOARD
  #
  # - Removes '_n' suffix from data-driven models.
  # - Generates sequential names for physiology models (e.g., physio_lm1).
  # - Formats metrics to 2 decimal places and wraps pm in math mode ($\pm$).
  # - Adds a title ABOVE the table and keeps the caption BELOW.
  # - AUTOMATICALLY COPIES THE OUTPUT TO YOUR CLIPBOARD.
  # ============================================================================
  
  escape_latex <- function(text_vector) {
    gsub("(?<!\\\\)_", "\\\\_", text_vector, perl = TRUE)
  }
  
  format_metrics <- function(text_vector) {
    converted <- gsub("\\+/-", "$\\\\pm$", text_vector)
    rounded <- gsub("(\\d+\\.\\d{2})\\d*", "\\1", converted)
    return(rounded)
  }
  
  latex_lines <- c()
  
  # --- Write Table Header ---
  latex_lines <- c(latex_lines, "\\begin{table}[htbp]")
  latex_lines <- c(latex_lines, "\\centering")
  
  # --- Added Heading ABOVE the table ---
  latex_lines <- c(latex_lines, paste0("\\textbf{Comparison between linear models with n = ", n, " predictors} \\\\[2mm]"))
  
  latex_lines <- c(latex_lines, "\\begin{tabular}{l l ccc}")
  latex_lines <- c(latex_lines, "\\toprule")
  latex_lines <- c(latex_lines, "\\textbf{Model Name} & \\textbf{Predictors} & \\makecell{$\\bm{\\beta}$\\\\ Estimate ($p$-val)} & \\makecell{\\textbf{cv-RMSE}\\\\ Mean $\\pm$ SD} & \\makecell{\\textbf{cv-$\\mathbf{R^2}$}\\\\ Mean $\\pm$ SD} \\\\")
  latex_lines <- c(latex_lines, "\\midrule")
  
  # --- STEP 1: Process Data-Driven Models ---
  latex_lines <- c(latex_lines, "% === DATA-DRIVEN APPROACHES ===")
  
  if (!is.null(data_df) && nrow(data_df) > 0) {
    for (i in 1:nrow(data_df)) {
      raw_model_name <- data_df$`model name`[i]
      if (raw_model_name != "") {
        raw_model_name <- gsub(paste0("_", n, "$"), "", raw_model_name)
      }
      
      model_name <- escape_latex(raw_model_name)
      variables  <- escape_latex(paste0("\\texttt{", data_df$`selected variables`[i], "}"))
      beta_pval_string <- sub("(", "\\hfill (", data_df$`final model beta (p-val)`[i], fixed = T)
      beta_pval  <- escape_latex(beta_pval_string)
      cv_rmse    <- format_metrics(data_df$`cv-RMSE (Mean pm SD)`[i])
      cv_r2      <- format_metrics(data_df$`cv-R2 (Mean pm SD)`[i])
      
      row_string <- paste(model_name, variables, beta_pval, cv_rmse, cv_r2, sep = " & ")
      row_string <- paste0(row_string, " \\\\")
      latex_lines <- c(latex_lines, row_string)
      
      if (i < nrow(data_df)) {
        if (data_df$`model name`[i + 1] != "") {
          latex_lines <- c(latex_lines, "\\addlinespace")
        }
      }
    }
  }
  
  # --- SECTION BREAK ---
  latex_lines <- c(latex_lines, "\\midrule")
  latex_lines <- c(latex_lines, "% === PHYSIOLOGY-DRIVEN APPROACHES ===")
  
  # --- STEP 2: Process Physiology-Driven Models ---
  if (!is.null(physio_df) && nrow(physio_df) > 0) {
    
    temp_model_ids <- physio_df$`model name`
    for(k in 1:length(temp_model_ids)) {
      if(temp_model_ids[k] == "") temp_model_ids[k] <- temp_model_ids[k-1]
    }
    
    unique_models <- unique(temp_model_ids)
    model_counter <- 1
    
    for (m_id in unique_models) {
      model_rows <- physio_df[temp_model_ids == m_id, ]
      
      sequential_model_name <- paste0("physio_lm", model_counter)
      model_counter <- model_counter + 1
      
      for (i in 1:nrow(model_rows)) {
        display_name <- if(i == 1) escape_latex(sequential_model_name) else ""
        
        variables  <- escape_latex(paste0("\\texttt{", model_rows$`selected variables`[i], "}"))
        beta_pval_string <- sub("(", "\\hfill (", model_rows$`final model beta (p-val)`[i], fixed = T)
        beta_pval  <- escape_latex(beta_pval_string)
        cv_rmse    <- format_metrics(model_rows$`cv-RMSE (Mean pm SD)`[i])
        cv_r2      <- format_metrics(model_rows$`cv-R2 (Mean pm SD)`[i])
        
        row_string <- paste(display_name, variables, beta_pval, cv_rmse, cv_r2, sep = " & ")
        row_string <- paste0(row_string, " \\\\")
        latex_lines <- c(latex_lines, row_string)
      }
      
      if (m_id != unique_models[length(unique_models)]) {
        latex_lines <- c(latex_lines, "\\addlinespace")
      }
    }
  }
  
  # --- Write Table Footer ---
  latex_lines <- c(latex_lines, "\\bottomrule")
  latex_lines <- c(latex_lines, "\\end{tabular}")
  
  # --- Caption BELOW the tabular block ---
  latex_lines <- c(latex_lines, paste0("\\caption{Performance comparison between data-driven and physiology-driven linear models with $n = ", n, "$ predictors, for the \\textbf{male cohort}. The $\\beta$ coefficients and corresponding $p$-values refer to the final model estimated on the complete dataset (of the male cohort). Cross-Validation metrics (cv-RMSE and cv-$R^2$) are computed exclusively on out-of-fold validation sets (Mean $\\pm$ SD).\\label{tab:performance_models_n", n, "}}"))
  
  latex_lines <- c(latex_lines, "\\end{table}")
  
  final_latex_code <- paste(latex_lines, collapse = "\n")
  
  # Print to console
  cat(final_latex_code, "\n\n")
  
  # --- CROSS-PLATFORM AUTO-COPY SYSTEM ---
  tryCatch({
    os <- Sys.info()[["sysname"]]
    if (os == "Windows") {
      con <- file("clipboard", open = "w")
      writeLines(final_latex_code, con)
      close(con)
      message(">>> LaTeX code automatically copied to clipboard! <<<")
    } else if (os == "Darwin") { # macOS
      con <- pipe("pbcopy", "w")
      writeLines(final_latex_code, con)
      close(con)
      message(">>> LaTeX code automatically copied to clipboard! <<<")
    } else if (os == "Linux") {
      if (system("which xclip", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        con <- pipe("xclip -selection clipboard", "w")
        writeLines(final_latex_code, con)
        close(con)
        message(">>> LaTeX code automatically copied to clipboard! <<<")
      } else {
        message("Note: Please install 'xclip' on Linux to enable auto-copy feature.")
      }
    }
  }, error = function(e) {
    message("Warning: Unable to access the system clipboard automatically.")
  })
  
  invisible(final_latex_code)
}

plot_rmse <- function(data_df, physio_df, n, filename = "plot_rmse.pdf") {
  # ============================================================================
  # R GGPLOT2 SYSTEM + AUTOMATED LATEX CODE GENERATOR & COPY (Multi-Asterisk)
  #
  # - Saves the vector PDF directly into the specified project folder.
  # - Fixed: Now highlights ALL models sharing the minimum mean RMSE with an asterisk.
  # - AUTOMATICALLY COPIES THE COMPLETE LATEX FIGURE CODE TO YOUR CLIPBOARD.
  # ============================================================================
  
  library(ggplot2)
  
  # Target directory configuration
  target_dir = "/Users/admin/Desktop/Project/Figures/"
  
  # If the user provided just a filename, prepend the full directory path
  if (!grepl("^/", filename)) {
    full_path <- file.path(target_dir, filename)
    latex_path <- file.path("Figures", filename) # Clean path for LaTeX document
  } else {
    full_path <- filename
    latex_path <- basename(filename) # Fallback
  }
  
  # Ensure the directory exists to avoid errors
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  # Helper to parse "Mean pm SD" or "Mean +/- SD" into a numeric data.frame
  parse_metrics <- function(df, type_label, n_val) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    valid_rows <- df$`model name` != ""
    sub_df <- df[valid_rows, ]
    
    means <- c()
    sds <- c()
    names <- c()
    
    physio_counter <- 1
    for(i in 1:nrow(sub_df)) {
      m_name <- sub_df$`model name`[i]
      if (type_label == "data") {
        m_name <- gsub(paste0("_", n_val, "$"), "", m_name)
      } else {
        m_name <- paste0("physio_lm", physio_counter)
        physio_counter <- physio_counter + 1
      }
      names <- c(names, m_name)
      
      clean_str <- gsub(" ", "", sub_df$`cv-RMSE (Mean pm SD)`[i])
      parts <- unlist(strsplit(clean_str, "pm|\\+/-"))
      
      means <- c(means, as.numeric(parts[1]))
      sds <- c(sds, as.numeric(parts[2]))
    }
    
    data.frame(model = names, mean = means, sd = sds, type = type_label, stringsAsFactors = FALSE)
  }
  
  # --- Process and combine data ---
  d_clean <- parse_metrics(data_df, "data", n)
  p_clean <- parse_metrics(physio_df, "physio", n)
  plot_data <- rbind(d_clean, p_clean)
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    stop("No valid data found to plot.")
  }
  
  plot_data$model <- factor(plot_data$model, levels = plot_data$model)
  
  # FIXED: Identify all models that match the absolute minimum (handles ties)
  min_rmse_val <- min(plot_data$mean)
  plot_data$is_best <- plot_data$mean == min_rmse_val
  
  # Calculate asterisk y-position
  plot_data$asterisk_y <- plot_data$mean + plot_data$sd + (max(plot_data$mean) * 0.03)
  
  # --- Build the Minimalist ggplot ---
  g <- ggplot(plot_data, aes(x = model, y = mean, color = type, shape = type)) +
    # Error bars
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.12, linewidth = 0.7) +
    # Mean points
    geom_point(size = 3) +
    # Best model asterisk marker (Plots on ALL rows where is_best == TRUE)
    geom_text(data = subset(plot_data, is_best == TRUE), 
              aes(y = asterisk_y, label = "*"), 
              size = 7, fontface = "bold", vjust = 0.3, show.legend = FALSE) +
    # Colors and Shapes (SteelBlue and Crimson)
    scale_color_manual(values = c("data" = "#4682B4", "physio" = "#C41E3A")) +
    scale_shape_manual(values = c("data" = 16, "physio" = 16)) +
    # Labels
    labs(x = NULL, y = "cv-RMSE (Mean \u00b1 SD)", title = NULL) +
    # Ultra-minimal Theme
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black", face = "plain", margin = margin(r = 8)),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.5),
      legend.position = "none"
    ) +
    ylim(0, max(plot_data$mean + plot_data$sd) * 1.1)
  
  # Print plot to R IDE
  print(g)
  
  # Save vector PDF to the absolute path
  ggsave(filename = full_path, plot = g, width = 6, height = 3.8, device = "pdf")
  message(paste(">>> Plot saved as vector PDF to:", full_path))
  
  # --- Generate LaTeX Code ---
  latex_lines <- c()
  latex_lines <- c(latex_lines, "\\begin{figure}[htbp]")
  latex_lines <- c(latex_lines, "  \\centering")
  latex_lines <- c(latex_lines, paste0("  \\textbf{Comparison between linear models with n = ", n, " predictors} \\\\[2mm]"))
  latex_lines <- c(latex_lines, paste0("  \\includegraphics[width=0.90\\linewidth]{", latex_path, "}"))
  # Updated caption text to reflect that multiple asterisks are possible in case of ex aequo
  latex_lines <- c(latex_lines, paste0("  \\caption{Cross-validated RMSE (Mean $\\pm$ SD) across different linear modeling approaches with $n = ", n, "$ predictors. Blue circles represent data-driven models, while red circles indicate physiology-driven models. The asterisk (*) highlights the best performing model(s) achieving the lowest average validation error. All evaluation metrics are computed out-of-fold.}"))
  latex_lines <- c(latex_lines, paste0("  \\label{fig:rmse_comparison_n", n, "}"))
  latex_lines <- c(latex_lines, "\\end{figure}")
  
  latex_code <- paste(latex_lines, collapse = "\n")
  
  # --- CROSS-PLATFORM AUTO-COPY SYSTEM ---
  tryCatch({
    os <- Sys.info()[["sysname"]]
    if (os == "Windows") { con <- file("clipboard", open = "w"); writeLines(latex_code, con); close(con) }
    else if (os == "Darwin") { con <- pipe("pbcopy", "w"); writeLines(latex_code, con); close(con) }
    else if (os == "Linux") {
      if (system("which xclip", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        con <- pipe("xclip -selection clipboard", "w"); writeLines(latex_code, con); close(con)
      }
    }
    message(">>> Complete LaTeX figure code automatically copied to clipboard! <<<")
  }, error = function(e) {})
  
  invisible(g)
}

plot_rmse_2 <- function(data_df, physio_df, n, filename = "plot_rmse.pdf") {
  # ============================================================================
  # SYSTEMA PLOT BASE R + AUTOMATED LATEX CODE GENERATOR & COPY
  # ============================================================================
  
  # Configuration directory di destinazione
  target_dir <- "/Users/admin/Desktop/Project/Figures/"
  
  if (!grepl("^/", filename)) {
    full_path <- file.path(target_dir, filename)
    latex_path <- file.path("Figures", filename)
  } else {
    full_path <- filename
    latex_path <- basename(filename)
  }
  
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  # --- 1. FUNZIONE INTERNA DI PARSING DEI DATI ---
  parse_metrics <- function(df, type_label, n_val) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    valid_rows <- df$`model name` != ""
    sub_df <- df[valid_rows, ]
    
    means <- c()
    sds <- c()
    names <- c()
    
    physio_counter <- 1
    for(i in 1:nrow(sub_df)) {
      m_name <- sub_df$`model name`[i]
      if (type_label == "data") {
        m_name <- gsub(paste0("_", n_val, "$"), "", m_name)
      } else {
        m_name <- paste0("physio_lm", physio_counter)
        physio_counter <- physio_counter + 1
      }
      names <- c(names, m_name)
      
      clean_str <- gsub(" ", "", sub_df$`cv-RMSE (Mean pm SD)`[i])
      parts <- unlist(strsplit(clean_str, "pm|\\+/-"))
      
      means <- c(means, as.numeric(parts[1]))
      sds <- c(sds, as.numeric(parts[2]))
    }
    
    data.frame(model = names, mean = means, sd = sds, type = type_label, stringsAsFactors = FALSE)
  }
  
  # Uniamo i dataset
  d_clean <- parse_metrics(data_df, "data", n)
  p_clean <- parse_metrics(physio_df, "physio", n)
  plot_data <- rbind(d_clean, p_clean)
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    stop("Nessun dato valido trovato per il plot.")
  }
  
  # Identifichiamo i modelli migliori (gestisce i pareggi)
  min_rmse_val <- min(plot_data$mean)
  plot_data$is_best <- plot_data$mean == min_rmse_val
  
  # Configurazione estetica (Grigio Antracite vs Grigio Medio)
  # MODIFICA: pch = 19 (cerchio pieno) per data, pch = 1 (cerchio vuoto) per physio
  plot_data$color <- ifelse(plot_data$type == "data", "#404040", "#808080")
  plot_data$pch   <- ifelse(plot_data$type == "data", 19, 1) 
  
  # --- 2. CODICE GRAFICO (FUNZIONE INTERNA) ---
  run_base_plot <- function() {
    par(mar = c(6, 5, 4, 2) + 0.1)
    
    y_min <- min(plot_data$mean - plot_data$sd, na.rm = TRUE) * 0.9
    y_max <- max(plot_data$mean + plot_data$sd, na.rm = TRUE) * 1.1
    x_coords <- 1:nrow(plot_data)
    
    # Creazione plot vuoto
    plot(x_coords, plot_data$mean, 
         type = "n", 
         xaxt = "n", 
         xlim = c(0.5, length(x_coords) + 0.5),
         ylim = c(y_min, y_max),
         xlab = "", 
         ylab = "cv-RMSE (Mean \u00b1 SD)",
         font.lab = 1, 
         las = 1)
    
    # Tacchette asse X
    axis(side = 1, at = x_coords, labels = FALSE, lwd = 1, lwd.ticks = 1)
    
    # 1. Barre di Errore (SD)
    arrows(x0 = x_coords, y0 = plot_data$mean - plot_data$sd, 
           x1 = x_coords, y1 = plot_data$mean + plot_data$sd, 
           code = 3, 
           angle = 90, 
           length = 0.05,
           col = plot_data$color,
           lwd = 1.2)
    
    # 2. Punti delle Medie (cerchio pieno vs cerchio vuoto)
    points(x_coords, plot_data$mean, 
           pch = plot_data$pch, 
           col = plot_data$color, 
           lwd = 1.5, # Un po' di spessore extra rende il cerchio vuoto più nitido
           cex = 1)
    
    # 3. Aggiunta dell'asterisco (*) per il/i modello/i migliore/i
    best_indices <- which(plot_data$is_best)
    if(length(best_indices) > 0) {
      asterisk_y <- plot_data$mean[best_indices] + plot_data$sd[best_indices] + (y_max - y_min) * 0.1
      
      # Ciclo per supportare l'assegnazione corretta del colore anche in caso di ex aequo
      for(idx in best_indices) {
        text(x = idx, 
             y = plot_data$mean[idx] + plot_data$sd[idx] + (y_max - y_min) * 0.1, 
             labels = "*", 
             col = plot_data$color[idx], # Cattura il colore specifico di quel modello
             cex = 1.5, 
             font = 1)
      }
    }
    # 4. Nomi dei modelli inclinati a 45 gradi
    text(x = x_coords, 
         y = par("usr")[3] - (y_max - y_min) * 0.1, 
         labels = plot_data$model, 
         srt = 45,             
         adj = 1,              
         xpd = TRUE,           
         font = 1,             
         cex = 1)
  }
  
  # --- 3. SALVATAGGIO IN APERTO E SU FILE ---
  run_base_plot()
  
  pdf(file = full_path, width = 6, height = 3.8)
  run_base_plot()
  dev.off()
  
  message(paste(">>> Plot salvato come PDF vettoriale in:", full_path))
  
  # --- 4. GENERAZIONE CODICE LATEX ---
  latex_lines <- c()
  latex_lines <- c(latex_lines, "\\begin{figure}[htbp]")
  latex_lines <- c(latex_lines, "  \\centering")
  latex_lines <- c(latex_lines, paste0("  \\textbf{Comparison between linear models with n = ", n, " predictors} \\\\[2mm]"))
  latex_lines <- c(latex_lines, paste0("  \\includegraphics[width=0.90\\linewidth]{", latex_path, "}"))
  # MODIFICA: Aggiornata la didascalia con "filled circles" e "open circles"
  latex_lines <- c(latex_lines, paste0("  \\caption{Cross-validated RMSE (Mean $\\pm$ SD) across different linear modeling approaches with $n = ", n, "$ predictors. Dark gray filled circles represent data-driven models, while light gray open circles indicate physiology-driven models. The asterisk (*) highlights the best performing model(s) achieving the lowest average validation error. All evaluation metrics are computed out-of-fold.}"))
  latex_lines <- c(latex_lines, paste0("  \\label{fig:rmse_comparison_n", n, "}"))
  latex_lines <- c(latex_lines, "\\end{figure}")
  
  latex_code <- paste(latex_lines, collapse = "\n")
  
  # --- 5. COPIA AUTOMATICA NEGLI APPUNTI ---
  tryCatch({
    os <- Sys.info()[["sysname"]]
    if (os == "Windows") { con <- file("clipboard", open = "w"); writeLines(latex_code, con); close(con) }
    else if (os == "Darwin") { con <- pipe("pbcopy", "w"); writeLines(latex_code, con); close(con) }
    else if (os == "Linux") {
      if (system("which xclip", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
        con <- pipe("xclip -selection clipboard", "w"); writeLines(latex_code, con); close(con)
      }
    }
    message(">>> Codice LaTeX completo copiato automaticamente negli appunti! <<<")
  }, error = function(e) {})
  
  invisible(plot_data)
}
