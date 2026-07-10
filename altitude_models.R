# ---------------------------------------------------------------------------- #
# MODELS WITH ALTITUDE
# ---------------------------------------------------------------------------- #
# Data organization

load("~/Desktop/Project/Data/big_df.RData")
# sviluppo e seleziono (e DHRmax) le covariate per i modelli con l'altitudine
# per la male cohort del big dataset
altitude_data <- big_df %>%
  filter(
    sex == "M"
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = c(altitude, HRmax, Pmax, VEmax, Pvt2, VEvt2, SpO2end)
  ) %>%
  mutate(
    DHRmax = HRmax_H - HRmax_N,
    VErvt2 = VEvt2_N / VEmax_N,
    Prvt2 = Pvt2_N / Pmax_N
  ) %>%
  rename(
    VEmax = VEmax_N,
    altitude = altitude_H,
    hSpO2end = SpO2end_H
  ) %>%
  select(
    ID, altitude, age, BMI, Prvt2, VErvt2, VEmax, DHRmax, hSpO2end
  )

# ---------------------------------------------------------------------------- #
# Funzioni

summary_lmcv <- function(lm_cv, model_id = "cv_lm") {
  output_rows <- list()
  
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
  coefficients <- names(lm_cv$finalModel$coefficients)
  for (coef in coefficients) {
    beta_val <- coef_matrix[coef, "Estimate"]
    p_val    <- coef_matrix[coef, "Pr(>|t|)"]
    p_str <- if(p_val < 0.001) "< 0.001" else sprintf("%.3f", p_val)
    
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
    row_rmse       <- if (is_first_row) rmse_string else ""
    row_r2         <- if (is_first_row) r2_string else ""
    
    # Append row data to the list
    output_rows[[length(output_rows) + 1]] <- data.frame(
      `model name` = row_model_id,
      `predictors` = coef,
      `final model beta (p-val)` = beta_string,
      `cv-RMSE (Mean pm SD)` = row_rmse,
      `cv-R2 (Mean pm SD)` = row_r2,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    # After the first iteration, switch the flag off for this model's remaining variables
    is_first_row <- FALSE
  }
  
  final_table <- do.call(rbind, output_rows)
  rownames(final_table) <- NULL
  
  return(final_table)
}

lm_table <- function(data, model_id, covariates) {
  
  # fit Linear Model via caret
  x_subset <- data[, covariates, drop = FALSE]
  target <- data$DHRmax
  set.seed(42)
  # fisso gli split, per coerenza tra le cv sui diversi modelli
  folds <- createMultiFolds(
    y = target,
    k = 5,
    times = 50
  )
  cv_control <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 50,
    index = folds
  )
  set.seed(42)
  lm_cv <- caret::train(
    x = x_subset,
    y = target,
    method = "lm",
    trControl = cv_control
  )
  
  results_tab <- summary_lmcv(lm_cv = lm_cv, model_id = model_id)
  
  return(results_tab)
}

mourot_table <- function(data) {
  
  beta0 = 0.7296
  beta1 = -0.0024
  
  altitude_gain <- data[, "altitude", drop = FALSE] - 200
  target <- data$DHRmax
  set.seed(42)
  # fisso gli split, per coerenza tra le cv sui diversi modelli
  folds <- createMultiFolds(
    y = target,
    k = 5,
    times = 50
  )
  
  # inizializzo il df dove salvare le metriche di ogni split
  cv_metrics <- data.frame(
    RMSE = numeric(),
    Rsquared = numeric()
  )
  
  # calcolo e salvo le metriche di ogni split
  for(i in seq_along(folds)) {
    # stime di mourot sugli stessi split (test) dei nuovi modelli con altitudine
    train_idx <- folds[[i]]
    test_idx <- setdiff(
      seq_len(length(target)),
      train_idx
    )
    x_test <- altitude_gain[test_idx, ]
    y_test <- target[test_idx]
    mourot_prediction <- beta0 + beta1 * x_test

    # calcolo e salvo le metriche di cv
    fold_res <- postResample(
      pred = mourot_prediction,
      obs = y_test
    )
    cv_metrics <- rbind(
      cv_metrics,
      data.frame(
        RMSE = fold_res["RMSE"],
        Rsquared = fold_res["Rsquared"]
      )
    )
  }
  
  # output table (simile a summary_lmcv)
  
  rmse_mean <- mean(cv_metrics$RMSE)
  rmse_sd   <- sd(cv_metrics$RMSE)
  rmse_string <- paste0(round(rmse_mean, 4), " +/- ", round(rmse_sd, 4))
  
  r2_mean   <- mean(cv_metrics$Rsquared)
  r2_sd     <- sd(cv_metrics$Rsquared)
  r2_string <- paste0(round(r2_mean, 4), " +/- ", round(r2_sd, 4))
    
  first_row <- data.frame(
    `model name` = "mourot equation",
    `predictors` = "(Intercept)",
    `final model beta (p-val)` = beta0,
    `cv-RMSE (Mean pm SD)` = rmse_string,
    `cv-R2 (Mean pm SD)` = r2_string,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  second_row <- data.frame(
    `model name` = "",
    `predictors` = "altitude gain",
    `final model beta (p-val)` = beta1,
    `cv-RMSE (Mean pm SD)` = "",
    `cv-R2 (Mean pm SD)` = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  final_table <- rbind(first_row, second_row)
  rownames(final_table) <- NULL
  
  return(final_table)
}

to_latex <- function(altitude_tab) {
  
  print(altitude_tab)
  
  latex_lines <- c()
  
  # Intestazione del blocco
  latex_lines <- c(latex_lines, "\\midrule")
  
  for(i in 1:nrow(altitude_tab)) {
    row <- altitude_tab[i, ]
    model_name <- gsub("_", "\\\\_", as.character(row$`model name`))
    
    # GESTIONE DEI SEPARATORI ORIZZONTALI (\midrule o \addlinespace)
    if (i > 1 && model_name != "") {
      if (model_name == "mourot equation") {
        latex_lines <- c(latex_lines, "      \\midrule") # \midrule specifico per mourot
      } else {
        latex_lines <- c(latex_lines, "      \\addlinespace") # \addlinespace per gli altri modelli
      }
    }
    
    predictor <- sprintf("\\texttt{%s}", as.character(row$predictors))
    
    # GESTIONE DEI COEFFICIENTI BETA E P-VALUE
    beta_pval <- as.character(row$`final model beta (p-val)`)
    
    if (grepl("\\(", beta_pval)) {
      # Per i modelli data/physiology-driven e lm standard
      parts <- unlist(strsplit(beta_pval, " \\("))
      beta_val <- parts[1]
      pval_stars <- paste0("(", parts[2])
      pval_stars <- gsub("<", "$<$", pval_stars) # Fix per p-val < 0.001
      beta_col <- sprintf("%s \\hfill %s", beta_val, pval_stars)
    } else {
      # Per la mourot Equation: restituisce solo il numero, senza (-) e centrato nella cella
      beta_col <- as.character(beta_pval)
    }
    
    # GESTIONE E ARROTONDAMENTO A 2 CIFRE PER LE METRICHE CV
    format_metric <- function(val_string) {
      if (val_string == "" || is.na(val_string)) return("")
      
      # Estrae i numeri (sia interi che decimali, positivi o negativi)
      numeri <- as.numeric(unlist(regmatches(val_string, gregexpr("-?[0-9]+\\.[0-9]+|-?[0-9]+", val_string))))
      
      if (length(numeri) == 2) {
        # Arrotonda a 2 cifre decimali e formatta con la virgola se necessario
        mean_val <- sprintf("%.2f", numeri[1])
        sd_val   <- sprintf("%.2f", numeri[2])
        return(sprintf("%s $\\pm$ %s", mean_val, sd_val))
      }
      return(val_string)
    }
    
    rmse_val <- format_metric(as.character(row$`cv-RMSE (Mean pm SD)`))
    r2_val   <- format_metric(as.character(row$`cv-R2 (Mean pm SD)`))
    
    # Costruzione della riga
    latex_row <- sprintf("      %s & %s & %s & %s & %s \\\\", 
                         model_name, predictor, beta_col, rmse_val, r2_val)
    latex_lines <- c(latex_lines, latex_row)
  }
  
  latex_lines <- c(latex_lines, "      \\bottomrule")
  
  # Unisce tutto il testo con i ritorni a capo
  full_output <- paste(latex_lines, collapse = "\n")
  
  # COPIA AUTOMATICA NEGLI APPUNTI (CLIPBOARD)
  os <- Sys.info()["sysname"]
  if (os == "Windows") {
    writeClipboard(full_output)
    cat("✔️ Codice LaTeX copiato negli appunti di Windows!\n")
  } else if (os == "Darwin") { # macOS
    con <- pipe("pbcopy", "w")
    writeLines(full_output, con)
    close(con)
    cat("✔️ Codice LaTeX copiato negli appunti di macOS!\n")
  } else {
    cat("⚠️ Sistema operativo non supportato per la copia automatica. Copia manualmente l'output sopra.\n")
  }
}

plot_rmse <- function(altitude_tab) {
  # --- 1. ESTRAZIONE E PULIZIA DEI DATI ---
  
  # Filtriamo via le righe vuote dei predittori secondari
  plot_data <- altitude_tab[altitude_tab$`model name` != "", ]
  
  # Separiamo l'equazione di mourot dagli altri modelli
  mourot_row <- plot_data[plot_data$`model name` == "mourot equation", ]
  models_row <- plot_data[plot_data$`model name` != "mourot equation", ]
  
  # Funzione helper interna per estrarre Media e SD numeriche
  extract_rmse <- function(string_vector) {
    matches <- gregexpr("[0-9]+\\.[0-9]+|[0-9]+", string_vector)
    reg_list <- regmatches(string_vector, matches)
    
    mean_val <- sapply(reg_list, function(x) as.numeric(x[1]))
    sd_val   <- sapply(reg_list, function(x) as.numeric(x[2]))
    
    return(list(mean = mean_val, sd = sd_val))
  }
  
  # Estraiamo i valori per i modelli standard
  rmse_models <- extract_rmse(models_row$`cv-RMSE (Mean pm SD)`)
  model_names <- models_row$`model name`
  model_names <- gsub(" Model", "\n Model", model_names)
  means       <- rmse_models$mean
  sds         <- rmse_models$sd
  
  # Estraiamo il valore medio di mourot
  mourot_rmse_mean <- extract_rmse(mourot_row$`cv-RMSE (Mean pm SD)`)$mean
  
  # --- 2. CONFIGURAZIONE SPAZI E MARGINI ---
  
  # Aumentiamo il margine inferiore (las = 2 richiede spazio per il testo verticale/obliquo)
  old_par <- par(no.readonly = TRUE) # Salviamo i parametri grafici attuali
  par(mar = c(6, 5, 4, 2) + 0.1)
  
  # Calcoliamo i limiti dell'asse Y
  y_min <- min(c(means - sds, mourot_rmse_mean), na.rm = TRUE) * 0.9
  y_max <- max(c(means + sds, mourot_rmse_mean), na.rm = TRUE) * 1.1
  
  x_coords <- 1:length(model_names)
  
  # --- 3. COSTRUZIONE DEL GRAFICO ---
  
  plot(x_coords, means, 
       type = "n", 
       xaxt = "n", 
       xlim = c(0.5, length(model_names) + 0.5),
       ylim = c(y_min, y_max),
       xlab = "", # Rimosso per non sovrapporsi ai nomi dei modelli
       ylab = "cv-RMSE (Mean \u00b1 SD)",
       font.lab = 1, 
       las = 1)
  
  axis(side = 1, at = x_coords, labels = FALSE, lwd = 1, lwd.ticks = 1)
  
  # 1. Linea orizzontale per mourot Equation
  abline(h = mourot_rmse_mean, lty = "dotted", lwd = 1)
  text(x = length(model_names) + 0.5, 
       y = mourot_rmse_mean + (y_max - y_min) * 0.02, 
       labels = "Mourot", 
       adj = c(1, 0), 
       cex = 0.75, 
       font = 3)
  
  # 2. Barre di errore (SD)
  arrows(x0 = x_coords, y0 = means - sds, 
         x1 = x_coords, y1 = means + sds, 
         code = 3, 
         angle = 90, 
         length = 0.05,
         lwd = 1.2)
  
  # 3. Punti delle medie
  points(x_coords, means, pch = 19, cex = 1)
  
  # 4. Asse X con testo inclinato a 45 gradi (per evitare sovrapposizioni se i nomi sono lunghi)
  text(x = x_coords, 
       y = par("usr")[3] - (y_max - y_min) * 0.2, # Posiziona il testo appena sotto l'asse Y minimo
       labels = model_names, 
       # srt = 45,             # Inclinazione a 45 gradi
       # adj = 0.5,              # Allineamento a destra del punto di ancoraggio
       xpd = TRUE,           # Permette di disegnare fuori dall'area del plot
       font = 1,             # Grassetto
       cex = 1)            # Dimensione del testo leggermente ridotta
  
  # Ripristiniamo i parametri grafici originali di R
  par(old_par)
}

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: altitude itself
# - data: raw altitude dataset (con NA per VErvt2 e Prvt2, quindi non costruiamo il Complex Model)

data <- altitude_data

model_id <- "Altitude Only Model"
covariates <- c("altitude")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("altitude", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = data)

altitude_tab <- rbind(alt_tab, simple_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: altitude itself
# - data: complete altitude dataset (senza NA per VErvt2 e Prvt2, possiamo costruire il Complex Model)

data <- na.omit(altitude_data)

model_id <- "Altitude Only Model"
covariates <- c("altitude")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("altitude", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Complex Model"
covariates <- c("altitude", "VEmax", "age", "VErvt2", "Prvt2", "BMI")
complex_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = data)

altitude_tab <- rbind(alt_tab, simple_tab, complex_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: altitude itself (as a factor!)
# - data: raw altitude dataset (con NA per VErvt2 e Prvt2, quindi non costruiamo il Complex Model)

num_data <- altitude_data # Morour equation vuole comunque altitude as numeric
data <- num_data %>% mutate(altitude = factor(altitude))

model_id <- "Altitude Only Model"
covariates <- c("altitude")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("altitude", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = num_data)

altitude_tab <- rbind(alt_tab, simple_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: altitude itself (as a factor!)
# - data: complete altitude dataset (senza NA per VErvt2 e Prvt2, possiamo costruire il Complex Model)

num_data <- na.omit(altitude_data) # Morour equation vuole comunque altitude as numeric
data <- num_data %>% mutate(altitude = factor(altitude))

model_id <- "Altitude Only Model"
covariates <- c("altitude")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("altitude", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Complex Model"
covariates <- c("altitude", "VEmax", "age", "VErvt2", "Prvt2", "BMI")
complex_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = num_data)

altitude_tab <- rbind(alt_tab, simple_tab, complex_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: hSpO2end
# - data: raw altitude dataset (con NA per VErvt2 e Prvt2, quindi non costruiamo il Complex Model)

num_data <- altitude_data # Morour equation vuole comunque altitude as numeric
data <- num_data %>% mutate(altitude = factor(altitude))

model_id <- "Saturation Only Model"
covariates <- c("hSpO2end")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("hSpO2end", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = num_data)

altitude_tab <- rbind(alt_tab, simple_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

# ---------------------------------------------------------------------------- #
# Risultati e codice LaTeX
# - Effect of altitude: hSpO2end
# - data: complete altitude dataset (senza NA per VErvt2 e Prvt2, possiamo costruire il Complex Model)

num_data <- na.omit(altitude_data) # Morour equation vuole comunque altitude as numeric
data <- num_data %>% mutate(altitude = factor(altitude))

model_id <- "Altitude-Only Model"
covariates <- c("hSpO2end")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("hSpO2end", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Complex Model"
covariates <- c("hSpO2end", "VEmax", "age", "VErvt2", "Prvt2", "BMI")
complex_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = num_data)

altitude_tab <- rbind(alt_tab, simple_tab, complex_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)
