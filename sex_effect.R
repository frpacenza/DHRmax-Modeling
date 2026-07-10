# ---------------------------------------------------------------------------- #
# EFFECT OF SEX
# ---------------------------------------------------------------------------- #
# Data organization

load("~/Desktop/Project/Data/predictors_df_with_target.RData")

data <- predictors_df_with_target %>% select(
  -ID, -subject, -subjectID, -session, -hHRmax, -Lapeak
)
complete_idx <- which(complete.cases(data))
complete_data <- data[complete_idx,]

target_M <- complete_data$DHRmax[complete_data$sex == "M"]
target_F <- complete_data$DHRmax[complete_data$sex == "F"]

# ---------------------------------------------------------------------------- #
# Effect on target variable DHRmax

# 1. Visual insight: boxplot and probability distribution estimate

# layout doppio
par(mfrow = c(1, 2))

# Boxplot affiancati (sx)
boxplot(DHRmax ~ sex, data = data,
        main = "Boxplot",
        col = c("#ffb7b2", "#a0c4ff"),
        ylab = "DHRmax (bpm)",
        xlab = "sex",
        frame.plot = T)

# Distribuzioni stimate (dx)
dens_M <- density(target_M)
dens_F <- density(target_F)

# Base plot with male density
plot(dens_M, 
     main = "Probability distribution",
     col = "#a0c4ff", 
     lwd = 2,
     ylab = "probability",
     xlab = "DHRmax (bpm)",
     xlim = range(c(dens_M$x, dens_F$x)), # Assicura che l'asse X copra entrambi
     ylim = range(c(dens_M$y, dens_F$y))) # Assicura che l'asse Y copra entrambi
# Add female density
lines(dens_F, col = "#ffb7b2", lwd = 2)

# Ripristina il layout a schermo singolo (buona pratica)
par(mfrow = c(1, 1))

# 2. t-test
print(t.test(DHRmax ~ sex, data = data))

# ---------------------------------------------------------------------------- #
# Effect on the predictors dataset

# 1. Table for anthropometric variables

# Prepariamo una matrice vuota per i risultati
anthro_vars <- c("age", "height", "weight", "BMI")
anthro_tab <- matrix(NA, nrow = length(anthro_vars), ncol = 5)
rownames(anthro_tab) <- anthro_vars
colnames(anthro_tab) <- c("Male (Mean $\\pm$ SD)", "Female (Mean $\\pm$ SD)", "t", "p-value", " ")

# riempiamo la tabella con i risultati
for(var in anthro_vars) {
  ttest <- t.test(complete_data[[var]] ~ complete_data$sex)
  
  m_M <- mean(complete_data[[var]][complete_data$sex == "M"])
  sd_M <- sd(complete_data[[var]][complete_data$sex == "M"])
  m_F <- mean(complete_data[[var]][complete_data$sex == "F"])
  sd_F <- sd(complete_data[[var]][complete_data$sex == "F"])
  
  str_M <- paste0(round(m_M, 1), " $\\pm$ ", round(sd_M, 1))
  str_F <- paste0(round(m_F, 1), " $\\pm$ ", round(sd_F, 1))
  p_val <- ttest$p.value
  p_str <- if(p_val < 0.001) "$<$ 0.001" else sprintf("%.3f", p_val)
  simbolo <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else if(p_val < 0.1) "." else ""

  anthro_tab[var, ] <- c(str_M, str_F, round(ttest$statistic, 2), p_str, simbolo)
}

# trascriviamo la tabella in codice LaTeX
testo_latex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\begin{tabular}{lccccc}\n",
  "\\toprule\n"
)
testo_latex <- paste0(testo_latex, "\\textbf{Variable} & \\makecell{\\textbf{Male}\\\\ (N = ", length(target_M), ")} & \\makecell{\\textbf{Female}\\\\ (N = ", length(target_F), ")} & $t$ & $p$-val & \\\\ \n")
testo_latex <- paste0(testo_latex, "\\midrule\n")
for(i in 1:nrow(anthro_tab)) {
  testo_latex <- paste0(testo_latex, "\\texttt{", rownames(anthro_tab)[i], "} & ", 
                        paste(anthro_tab[i, ], collapse = " & "), " \\\\\n")
}
testo_latex <- paste0(testo_latex, "\\bottomrule\n\\end{tabular}\n")
testo_latex <- paste0(testo_latex, "\\caption{Comparative analysis of anthropometric variables by biological sex. Results are presented as Mean $\\pm$ SD, along with Welch $t$-test statistics and associated $p$-values.}\n", "\\end{table}\n")

# Copia diretta negli appunti
writeLines(testo_latex, pipe("pbcopy"))

# 2. PCA for test metrics

# Eseguiamo la PCA sulle sole metriche di esercizio (i.e., non quelle antropometriche)
pca_vars <- setdiff(names(complete_data), c("sex", anthro_vars))
pca_data  <- complete_data[, pca_vars]
pca_res <- prcomp(pca_data, scale. = TRUE)

# Calcoliamo la varianza spiegata dalle prime due componenti per inserirla negli assi del grafico
var_spiegata <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
lbl_pc1 <- paste0("PC1 (", round(var_spiegata[1], 1), "%)")
lbl_pc2 <- paste0("PC2 (", round(var_spiegata[2], 1), "%)")

# Impostiamo i colori
col_M <- "#a0c4ff"
col_F <- "#ffb7b2"
colors <- ifelse(complete_data$sex == "M", col_M, col_F)

# Impostiamo i margini
par(mar = c(4.1, 4.1, 3.1, 2.1))

# PCA plot
plot(pca_res$x[, 1], pca_res$x[, 2],
     col = colors,
     pch = 16,
     cex = 1.2,
     main = "PCA on exercise metrics",
     xlab = lbl_pc1,
     ylab = lbl_pc2)