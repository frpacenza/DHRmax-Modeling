load("~/Desktop/Project/Data/multistudy_df.RData")
data <- multistudy_df %>%
  pivot_wider(
    names_from = condition,
    values_from = -c("ID", "condition", "method", "sex", "age", "height", "weight", "BMI")
  ) %>%
  mutate(
    DHRmax = HRmax_N - HRmax_H,
    altitude_H = factor(altitude_H, levels = c("2500", "3500", "4500"))
  ) %>%
  filter(
    sex == "M",
    altitude_H == 3500
  ) %>%
  select(
    age, height, weight, BMI,
    ends_with("_N"),
    SpO2rest_H,
    DHRmax
  ) %>%
  rename_with(~ sub("_N$", "", .x), ends_with("_N")) %>%
  rename(hSpO2rest = SpO2rest_H )
# md pattern per 3500 cohort: tanti missing values. togliendo le variabili con
# missing values restano solo 12 variabili, limitando l'analisi. Quindi per
# ora ci concentriamo sulla 4500 cohort.
md.pattern(data, rotate.names = T)

# ---------------------------------------------------------------------------- #
# 4500 COHORT OF THE MULTI-STUDY DATASET
# ---------------------------------------------------------------------------- #
data <- multistudy_df %>%
  pivot_wider(
    names_from = condition,
    values_from = -c("ID", "condition", "method", "sex", "age", "height", "weight", "BMI")
  ) %>%
  mutate(
    DHRmax = HRmax_N - HRmax_H,
    altitude_H = factor(altitude_H, levels = c("2500", "3500", "4500"))
  ) %>%
  filter(
    sex == "M",
    altitude_H == 4500
  ) %>%
  select(
    age, height, weight, BMI,
    ends_with("_N"), -altitude_N,
    SpO2rest_H,
    DHRmax
  ) %>%
  rename_with(~ sub("_N$", "", .x), ends_with("_N")) %>%
  rename(hSpO2rest = SpO2rest_H )
# md pattern per 4500 cohort: molto meglio. Ci sono solo quattro variabili con 
# NA (SpO2vt1, SpO2vt1, LArest, rPmax), che possiamo quindi rimuovere.
md.pattern(data, rotate.names = T)
data <- data %>% t() %>% as.data.frame() %>% na.omit() %>% t() %>% as.data.frame()

# ---------------------------------------------------------------------------- #
source("~/Desktop/Project/RFiles/functions.R")
get_sequential_cv_rmse <- function(selected, predictors, target) {
  # ============================================================================
  # SEQUENTIAL CROSS-VALIDATED RMSE EVALUATOR (via caret)
  #
  # - selected: Character vector with the ordered names of N predictors.
  # - predictors: Dataframe containing the predictor columns.
  # - target: Numeric vector containing the target variable (Y).
  # - Control: 5-fold Repeated CV with 50 repeats.
  #
  # OUTPUT: A dataframe with N+1 rows containing:
  #         Model_Size, Mean, SD, and determinants (the specific added variable).
  # ============================================================================
  
  library(caret)
  
  # Ensure the big dataframe is treated as a standard data.frame
  predictors <- as.data.frame(predictors)
  
  N <- length(selected)
  
  # Pre-allocate vectors for the final dataframe columns
  model_size_vec <- 0:N
  mean_vec       <- numeric(N + 1)
  sd_vec         <- numeric(N + 1)
  det_vec        <- character(N + 1) # Pre-allocate character vector for determinants
  
  # Define the resampling scheme: 5-fold CV repeated 50 times
  cv_control <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50
  )
  
  # --- Step 1: Model 0 (Intercept-only, i = 0) ---
  X_intercept <- data.frame(intercept = rep(1, nrow(predictors)))
  
  set.seed(42)
  model_0 <- caret::train(
    x = X_intercept,
    y = target,
    method = "lm",
    trControl = cv_control
  )
  
  # Extract Mean and SD across all 250 resamples
  mean_vec[1] <- round(mean(model_0$resample$RMSE), 4)
  sd_vec[1]   <- round(sd(model_0$resample$RMSE), 4)
  det_vec[1]  <- NA_character_ # Model 0 has no predictors
  
  # --- Step 2: Sequential Models (i = 1 to N) ---
  for (i in 1:N) {
    current_predictors <- selected[1:i]
    
    # Subset predictors safely for the mathematical model
    X_data <- predictors[, current_predictors, drop = FALSE]
    
    set.seed(42)
    model_i <- caret::train(
      x = X_data,
      y = target,
      method = "lm",
      trControl = cv_control
    )
    
    # Store metrics and ONLY the specific i-th variable added at this step
    mean_vec[i + 1] <- round(mean(model_i$resample$RMSE), 4)
    sd_vec[i + 1]   <- round(sd(model_i$resample$RMSE), 4)
    det_vec[i + 1]  <- selected[i] # Fixed: tracks only the specific new variable
  }
  
  # --- Step 3: Combine into a structured data.frame ---
  results_df <- data.frame(
    Model_Size   = model_size_vec,
    Mean         = mean_vec,
    SD           = sd_vec,
    determinants = det_vec,
    stringsAsFactors = FALSE
  )
  
  return(results_df)
}
plot_sequential_rmse <- function(res_df) {
  # ============================================================================
  # R BASE GRAPHICS FUNCTION - SEQUENTIAL cv-RMSE PLOT
  #
  # - Input: res_df (Dataframe with columns: Model_Size, Mean, SD)
  # - X-axis: Model_Size, Y-axis: Mean (dots) and SD (error bars).
  # - Minimalist: No main title, no X-axis label, no legend.
  # - Custom Y-axis label with the PM (\u00b1) symbol.
  # ============================================================================
  
  # 1. Calculate dynamic Y-axis limits to prevent error bars from clipping
  y_min <- min(res_df$Mean - res_df$SD) * 0.95
  y_max <- max(res_df$Mean + res_df$SD) * 1.05
  
  # 2. Generate the core plot for the Mean points
  plot(
    x = res_df$Model_Size, 
    y = res_df$Mean,
    pch = 16,                         # Solid circle
    col = "black",                  # Steel Blue
    cex = 1.3,                        # Point size
    ylim = c(y_min, y_max),           # Dynamic Y limits
    xaxt = "n",                       # Suppress automatic X-axis formatting
    xlab = "",                        # Remove X-axis label
    ylab = "cv-RMSE (Mean \u00b1 SD)", # Y-axis label with unicode ± symbol
    main = ""                         # Remove main title
  )
  
  # 3. Force exact integer ticks for the Model Size on the X-axis
  axis(side = 1, at = res_df$Model_Size, labels = res_df$Model_Size)
  
  # 4. Add the standard deviation (SD) error bars via base R arrows
  arrows(
    x0 = res_df$Model_Size, 
    y0 = res_df$Mean - res_df$SD, # Lower bound
    x1 = res_df$Model_Size, 
    y1 = res_df$Mean + res_df$SD, # Upper bound
    length = 0.05,                # Width of the horizontal cap
    angle = 90,                   # Perpendicular cap (90 degrees)
    code = 3,                     # Draw caps on both ends
    col = "black",              # Matching color
    lwd = 1.5                     # Line width
  )
}
generate_latex_table <- function(df_corr, df_enet, df_rforest) {
  # ============================================================================
  # LATEX TABLE GENERATOR & AUTO-COPY SYSTEM (PDF-Landscape - Syntax Fixed)
  #
  # - Inputs: 3 dataframes from 'get_sequential_cv_rmse'
  # - Formatting: Centered first column, \texttt{} for determinants.
  # - Fixed: Replaced ALL structural '\\//' with standard '\\\\' row endings.
  # - Output: Standard \begin{table} wrapped inside \begin{landscape}
  # - AUTOMATICALLY COPIES THE COMPLETE LATEX TABLE CODE TO YOUR CLIPBOARD.
  # ============================================================================
  
  if (!all(df_corr$Model_Size == df_enet$Model_Size) || !all(df_corr$Model_Size == df_rforest$Model_Size)) {
    stop("Error: The three dataframes must have the same Model_Size sequence.")
  }
  
  sizes <- df_corr$Model_Size
  N_rows <- length(sizes)
  
  build_cell_string <- function(mean_val, sd_val) {
    mean_fmt <- format(round(mean_val, 2), nsmall = 2)
    sd_fmt   <- format(round(sd_val, 2), nsmall = 2)
    paste0("$", mean_fmt, " \\pm ", sd_fmt, "$")
  }
  
  corr_cv    <- build_cell_string(df_corr$Mean, df_corr$SD)
  enet_cv    <- build_cell_string(df_enet$Mean, df_enet$SD)
  rforest_cv <- build_cell_string(df_rforest$Mean, df_rforest$SD)
  
  # Format determinants with \texttt{} or return --- if NA
  clean_det <- function(det_vector) {
    ifelse(is.na(det_vector), "---", paste0("\\texttt{", as.character(det_vector), "}"))
  }
  
  corr_det    <- clean_det(df_corr$determinants)
  enet_det    <- clean_det(df_enet$determinants)
  rforest_det <- clean_det(df_rforest$determinants)
  
  # Generate the LaTeX Code using pdflscape environment
  latex_lines <- c()
  latex_lines <- c(latex_lines, "\\begin{landscape}")
  latex_lines <- c(latex_lines, "\\vspace*{\\fill}")
  latex_lines <- c(latex_lines, "\\begin{table}[htbp]")
  latex_lines <- c(latex_lines, "  \\centering")
  latex_lines <- c(latex_lines, "  \\textbf{\\Large Sequential cv-RMSE performance and selected determinants across methods} \\\\[3mm]")
  latex_lines <- c(latex_lines, "  \\begin{tabular}{ccccccc}")
  latex_lines <- c(latex_lines, "    \\toprule")
  
  # Row 1: Macro-headers
  # FIXED: Corrected syntax to standard '\\\\'
  macro_header <- paste0(
    "    \\multirow{2}{*}{\\textbf{Model size}} & ",
    "\\multicolumn{2}{c}{\\textbf{Correlation}} & ",
    "\\multicolumn{2}{c}{\\textbf{Elastic Net}} & ",
    "\\multicolumn{2}{c}{\\textbf{Random Forest}} \\\\"
  )
  latex_lines <- c(latex_lines, macro_header)
  
  latex_lines <- c(latex_lines, "    \\cmidrule(r){2-3} \\cmidrule(lr){4-5} \\cmidrule(l){6-7}")
  
  # Row 2: Sub-headers
  # FIXED: Corrected syntax to standard '\\\\'
  sub_header <- paste0(
    "    & \\textbf{Determinant} & \\shortstack{\\textbf{cv-RMSE} \\\\ (Mean $\\pm$ SD)} ",
    "& \\textbf{Determinant} & \\shortstack{\\textbf{cv-RMSE} \\\\ (Mean $\\pm$ SD)} ",
    "& \\textbf{Determinant} & \\shortstack{\\textbf{cv-RMSE} \\\\ (Mean $\\pm$ SD)} \\\\"
  )
  latex_lines <- c(latex_lines, sub_header)
  latex_lines <- c(latex_lines, "    \\midrule")
  
  # Populate rows
  for (i in 1:N_rows) {
    size_label <- as.character(sizes[i])
    
    # FIXED: Corrected syntax to standard '\\\\'
    row_line <- paste0(
      "    ", size_label, " & ",
      corr_det[i], " & ", corr_cv[i], " & ",
      enet_det[i], " & ", enet_cv[i], " & ",
      rforest_det[i], " & ", rforest_cv[i], " \\\\"
    )
    latex_lines <- c(latex_lines, row_line)
  }
  
  latex_lines <- c(latex_lines, "    \\bottomrule")
  latex_lines <- c(latex_lines, "  \\end{tabular}")
  latex_lines <- c(latex_lines, "  \\caption{Comparison of cv-RMSE obtained via repeated 5-fold cross-validation (50 repeats) as a function of model size $n$. In the compared linear models, features are entered sequentially based on the order established by Correlation, Elastic Net, and Random Forest based feature-selection methods, respectively. $\\mathbf{n = 0}$ \\textbf{corresponds to the Null Model.}}")
  latex_lines <- c(latex_lines, "  \\label{tab:sequential_rmse}")
  latex_lines <- c(latex_lines, "\\end{table}")
  latex_lines <- c(latex_lines, "\\vspace*{\\fill}")
  latex_lines <- c(latex_lines, "\\end{landscape}")
  
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
    message(">>> Complete LaTeX table code automatically copied to clipboard! <<<")
  }, error = function(e) {})
  
  invisible(latex_code)
}


target_4500M <- data$DHRmax
predictors_4500M <- data %>% select(-DHRmax)

corr_predictors_4500M <- corr_fs(predictors = predictors_4500M, target = target_4500M, n = 10)
corr_selected_4500M <- names(corr_predictors_4500M)
corr_res_4500M <- get_sequential_cv_rmse(selected = corr_selected_4500M, target = target_4500M, predictors = predictors_4500M)

enet_predictors_4500M <- enet_fs(predictors = predictors_4500M, target = target_4500M, n = 10)
enet_selected_4500M <- names(enet_predictors_4500M)
enet_res_4500M <- get_sequential_cv_rmse(selected = enet_selected_4500M, target = target_4500M, predictors = predictors_4500M)

rforest_predictors_4500M <- rforest_fs(predictors = predictors_4500M, target = target_4500M, n = 10)
rforest_selected_4500M <- names(rforest_predictors_4500M)
rforest_res_4500M <- get_sequential_cv_rmse(selected = rforest_selected_4500M, target = target_4500M, predictors = predictors_4500M)

plot_sequential_rmse(corr_res_4500M)
plot_sequential_rmse(enet_res_4500M)
plot_sequential_rmse(rforest_res_4500M)

generate_latex_table(corr_res_4500M, enet_res_4500M, rforest_res_4500M)

# ---------------------------------------------------------------------------- #
# 4500-UNIFIED MODEL
source("~/Desktop/Project/RFiles/altitude_interactions.R")
raw_multistudy_dataset_v2 <- read_excel("Desktop/Project/Data/multistudy_dataset_v2.xlsx")
multistudy_VErest_v2 <- raw_multistudy_dataset_v2$`VE REST`
multistudy_SpO2rest_v2 <- raw_multistudy_dataset_v2$`SPO2 REST`
multistudy_df_v2 <- multistudy_df %>%
  mutate(
    VErest = multistudy_VErest_v2,
    SpO2rest = multistudy_SpO2rest_v2
  )

data <- multistudy_df_v2 %>%
  filter(
    sex == "M"
  ) %>%
  select(
    c("altitude", "Pmax", "HRrvt1", "VErest", "SpO2rest", "SpO2end", "LApeak", "HRmax", "condition", "ID")
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = c("altitude", "Pmax", "HRrvt1", "VErest", "SpO2rest", "SpO2end", "LApeak", "HRmax")
  ) %>%
  mutate(
    DHRmax = HRmax_H - HRmax_N,
    altitude_H = as.factor(altitude_H)
  ) %>%
  rename(
    altitude = altitude_H,
    Pmax = Pmax_N,
    HRrvt1 = HRrvt1_N,
    VErest = VErest_N,
    SpO2rest = SpO2rest_N,
    SpO2end = SpO2end_N,
    LApeak = LApeak_N,
    HRmax = HRmax_N
  ) %>%
  select(
    c("altitude", "Pmax", "HRrvt1", "VErest", "SpO2rest", "SpO2end", "LApeak", "HRmax", "DHRmax")
  ) %>% na.omit()

model_id <- "4500-Driven Model"
covariates <- c("altitude", "Pmax", "HRrvt1", "VErest", "SpO2rest", "SpO2end", "LApeak", "HRmax") 
uni_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

altitude_tab <- rbind(alt_tab, add_tab, uni_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)
