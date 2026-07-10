# ---------------------------------------------------------------------------- #
# ALTITUDE EDA
# ---------------------------------------------------------------------------- #
# Libraries and data loading

library(dplyr)
library(tidyr)
library(pcaMethods)
load("~/Desktop/Project/Data/multistudy_df.RData")

# ---------------------------------------------------------------------------- #
# 1. Summary table of included studies

multistudy_tab <- function(df) {
  if (!requireNamespace("clipr", quietly = TRUE)) {
    stop("Il pacchetto 'clipr' è richiesto per copiare negli appunti. Installalo con install.packages('clipr')")
  }
  
  latex_rows <- df %>%
    # 1. Isolo l'acronimo finale dall'ID
    mutate(
      acronym = gsub("^S0[0-9]+", "", ID),
      study = case_when(
        acronym == "Cal2025" ~ "Callovini et al. 2025",
        acronym == "Cal2026" ~ "Callovini et al. 2026",
        acronym == "Dor"     ~ "Dorelli et al. (unpublished data)",
        acronym == "For2018" ~ "Fornasiero et al. 2018",
        acronym %in% c("Sky2024a", "Sky2024b", "Sky2025a", "Sky2025b", "Sky2026a") ~ "Skyrunning Study",
        TRUE ~ acronym
      )
    ) %>%
    # 2. Trovo l'altitudine massima ipossica per ogni soggetto
    group_by(ID) %>%
    mutate(hypoxic_altitude = max(altitude)) %>%
    ungroup() %>%
    # 3. Filtro in NORMOSSIA (1 riga per soggetto)
    filter(condition == "N") %>%
    # 4. Raggruppo per studio e calcolo i riassunti
    group_by(study) %>%
    summarise(
      n_subjects = n(),
      pct_female = mean(sex == "F", na.rm = TRUE) * 100,
      methods = paste(unique(method), collapse = ", "),
      altitudes = paste(sort(unique(hypoxic_altitude)), collapse = ", "),
      
      age_mean = mean(age, na.rm = TRUE),
      age_sd   = sd(age, na.rm = TRUE),
      bmi_mean = mean(BMI, na.rm = TRUE),
      bmi_sd   = sd(BMI, na.rm = TRUE),
      
      rvo2_mean = mean(rVO2max, na.rm = TRUE),
      rvo2_sd   = sd(rVO2max, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # 5. Generazione stringhe LaTeX con l'escape corretto per la barra rovesciata (\\\\\\%)
    mutate(
      latex_line = sprintf(
        "  %s & %d & %.1f\\%% & %s & %s & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f & %.1f $\\pm$ %.1f \\\\",
        study, n_subjects, pct_female, methods, altitudes, age_mean, age_sd, bmi_mean, bmi_sd, rvo2_mean, rvo2_sd
      )
    ) %>%
    # Sostituzione di sicurezza per garantire al 100% che rimanga \\% negli appunti
    mutate(latex_line = gsub("%", "\\\\%", latex_line, fixed = TRUE)) %>%
    # Ripristiniamo i dollari della matematica che il gsub ha raddoppiato per errore
    mutate(latex_line = gsub("$\\\\\\pm$", "$\\pm$", latex_line, fixed = TRUE)) %>%
    pull(latex_line) %>%
    paste(collapse = "\n")
  
  clipr::write_clip(latex_rows)
  
  cat("Successo! Righe con \\% corretto copiate negli appunti.\n\n")
  cat(latex_rows, "\n")
}
multistudy_tab(multistudy_df)

# ---------------------------------------------------------------------------- #
# 2. Visual Insight: Boxplot and distribution estimate of DHRmax for each study

# Data organization
data <- multistudy_df %>%
  pivot_wider(
    names_from = condition,
    values_from = -c("ID", "condition", "method", "sex", "age", "height", "weight", "BMI")
  ) %>%
  mutate(
    DHRmax = HRmax_H - HRmax_N,
    altitude_H = factor(altitude_H, levels = c("2500", "3500", "4500"))
  )

colors <- c("#d8b4fe", "#a855f7", "#5b127d")
names(colors) <- levels(data$altitude_H)

# Plot
par(mfrow = c(1, 2))
# Grafico A: Boxplot per Studio
boxplot(DHRmax ~ altitude_H, 
        data = data,
        main = "Boxplot",
        col = colors,
        ylab = "DHRmax (bpm)", 
        xlab = "Altitude (m)")#,
        #las = 2) # Ruota le etichette degli studi di 90 gradi per leggerle
# Grafico B: Densità sovrapposte
density_list <- lapply(split(data$DHRmax, data$altitude), density)
xlim_global <- range(sapply(density_list, function(d) d$x))
ylim_global <- range(sapply(density_list, function(d) d$y))
plot(density_list[[1]], 
     main = "Probability distribution", 
     xlab = "DHRmax (bpm)", 
     ylab = "probability",
     xlim = xlim_global, 
     ylim = ylim_global, 
     col = colors[1], 
     lwd = 2)
for(i in 2:length(density_list)) {
  lines(density_list[[i]], col = colors[i], lwd = 2)
}
# Ripristiniamo il layout singolo
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------- #
# 3. Visual Insight: PCA plot

data <- multistudy_df %>%
  pivot_wider(
    names_from = condition,
    values_from = -c("ID", "condition", "method", "sex", "age", "height", "weight", "BMI")
  ) %>%
  mutate(
    altitude_H = factor(altitude_H, levels = c("2500", "3500", "4500"))
  ) %>%
  select(-c("ID", "method", "sex", "altitude_N"))

# Risultati PPCA (con pcaMethods, poiché ci sono molti missing data)
ppca_res <- pca(data %>% select(-altitude_H), method = "ppca", nPcs = 3, scale = "vector", center = TRUE)

# Dati per il plot:
var_cum <- ppca_res@R2cum * 100
var_marginal <- c(var_cum[1], diff(var_cum))
lbl_pc1 <- paste0("PC1 (", round(var_marginal[1], 1), "%)")
lbl_pc2 <- paste0("PC2 (", round(var_marginal[2], 1), "%)")
coordinate_pc <- scores(ppca_res)
pca_colors <- colors[as.character(data$altitude_H)]

# Plot PPCA 
par(mar = c(4.1, 4.1, 3.1, 2.1))
plot(-coordinate_pc[, 1], -coordinate_pc[, 2],
     col = pca_colors,
     pch = 16,
     cex = 1.3,
     main = "PPCA on anthropometric and exercise variables",
     xlab = lbl_pc1,
     ylab = lbl_pc2)
par(mar = c(5.1, 4.1, 4.1, 2.1))

# ---------------------------------------------------------------------------- #
# Effect on the predictors dataset

# 1. Table for anthropometric variables

# Prepariamo una matrice vuota per i risultati
anthro_vars <- c("age", "height", "weight", "BMI")
anthro_tab <- matrix(NA, nrow = length(anthro_vars), ncol = 6)
rownames(anthro_tab) <- anthro_vars
colnames(anthro_tab) <- c("2500 m (Mean $\\pm$ SD)", "3500 m (Mean $\\pm$ SD)", "4500 m (Mean $\\pm$ SD)", "t", "p-value", " ")

# riempiamo la tabella con i risultati
for(var in anthro_vars) {
  aov <- aov(data[[var]] ~ data$altitude_H)
  
  m_2500 <- mean(data[[var]][data$altitude_H == 2500])
  sd_2500 <- sd(data[[var]][data$altitude_H == 2500])
  m_3500 <- mean(data[[var]][data$altitude_H == 3500])
  sd_3500 <- sd(data[[var]][data$altitude_H == 3500])
  m_4500 <- mean(data[[var]][data$altitude_H == 4500])
  sd_4500 <- sd(data[[var]][data$altitude_H == 4500])
  
  str_2500 <- paste0(round(m_2500, 1), " $\\pm$ ", round(sd_2500, 1))
  str_3500 <- paste0(round(m_3500, 1), " $\\pm$ ", round(sd_3500, 1))
  str_4500 <- paste0(round(m_4500, 1), " $\\pm$ ", round(sd_4500, 1))
  F_val <- summary(aov)[[1]][["F value"]][1]
  p_val <- summary(aov)[[1]][["Pr(>F)"]][1]
  p_str <- if(p_val < 0.001) "$<$ 0.001" else sprintf("%.3f", p_val)
  simbolo <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else if(p_val < 0.1) "." else ""
  
  anthro_tab[var, ] <- c(str_2500, str_3500, str_4500, round(F_val, 2), p_str, simbolo)
}

testo_latex <- ""
for(i in 1:nrow(anthro_tab)) {
  testo_latex <- paste0(testo_latex, "\\texttt{", rownames(anthro_tab)[i], "} & ", 
                        paste(anthro_tab[i, ], collapse = " & "), " \\\\\n")
}
# Copia diretta negli appunti
writeLines(testo_latex, pipe("pbcopy"))
