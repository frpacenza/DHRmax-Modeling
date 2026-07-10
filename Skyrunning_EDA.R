# ---------------------------------------------------------------------------- #
# Skyrunning Dataset EDA
# ---------------------------------------------------------------------------- #
# Data organization:

load("~/Desktop/Project/Data/predictors_df_with_target.RData")
data <- predictors_df_with_target %>%
  select(-Lapeak) %>%
  na.omit() %>%
  mutate(
    subjectID = factor(subjectID)
  )

colors <- c(
  "#2B5C8F", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
  "#E6AB02", "#A6761D", "#666666", "#1B9E77", "#A6CEE3", 
  "#FB9A99", "#FDBF6F", "#CAB2D6", "#B2DF8A"
)
names(colors) <- levels(data$subjectID)
id_colors <- colors[data$subjectID]

# ---------------------------------------------------------------------------- #
# 1. PCA for exercise metrics

# Eseguiamo la PCA sulle sole metriche di esercizio (i.e., non quelle antropometriche)
# same data as for the sex EDA:
pca_data <- predictors_df_with_target %>% 
  select(-c(ID, subject, subjectID, session, hHRmax, Lapeak, sex, age, height, weight, BMI)) %>% 
  na.omit()
pca_res <- prcomp(pca_data, scale. = TRUE)

# Calcoliamo la varianza spiegata dalle prime due componenti per inserirla negli assi del grafico
var_spiegata <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
lbl_pc1 <- paste0("PC1 (", round(var_spiegata[1], 1), "%)")
lbl_pc2 <- paste0("PC2 (", round(var_spiegata[2], 1), "%)")


# Impostiamo i margini
par(mar = c(4.1, 4.1, 3.1, 2.1))

# PCA plot
plot(pca_res$x[, 1], pca_res$x[, 2],
     col = id_colors,
     pch = 16,
     cex = 1.2,
     main = "PCA on exercise metrics",
     xlab = lbl_pc1,
     ylab = lbl_pc2)

for (soggetto in levels(data$subjectID)) {
  # Trova gli indici delle righe corrispondenti a quel determinato atleta
  indices <- which(data$subjectID == soggetto)
  
  # Se l'atleta ha fatto almeno due test session, unisci i punti nell'ordine in cui compaiono
  if (length(indices) > 1) {
    lines(pca_res$x[indices, 1], pca_res$x[indices, 2],
          col = id_colors[soggetto], # Usa lo stesso colore assegnato al soggetto
          lty = 1,                              # 2 significa linea tratteggiata ("dashed")
          lwd = 1)                            # Spessore della linea
  }
}

# ---------------------------------------------------------------------------- #
# 2.Spaghetti plot:

plot(NULL, xlim = c(0.7, 2.3), 
     ylim = c(min(c(data$HRmax, data$hHRmax), na.rm = TRUE) - 5, 
              max(c(data$HRmax, data$hHRmax), na.rm = TRUE) + 5),
     xaxt = "n",
     xlab = "",
     ylab = "HRmax (bpm)",
     main = "Hypoxia-induced HRmax decline")

# Aggiungi le etichette personalizzate sull'asse X
axis(1, at = c(1, 2), labels = c("Normoxia", "Hypoxia"))
grid(nx = NA, ny = NULL, col = "gray90", lty = "solid") # griglia orizzontale

# --- CICLO PER DISEGNARE LE LINEE E I PUNTI ---
for (i in 1:nrow(data)) {
  soggetto <- data$subjectID[i]
  colore_atleta <- id_colors[soggetto]
  
  # Coordinate: x=1 è Normossia (HRmax), x=2 è Ipossia (hHRmax)
  x_coords <- c(1, 2)
  y_coords <- c(data$HRmax[i], data$hHRmax[i])
  
  # Disegna il segmento che unisce i due test della sessione i-esima
  lines(x_coords, y_coords, col = colore_atleta, lwd = 1.5, lty = 1)
  
  # Disegna i punti nei rispettivi estremi
  points(x_coords, y_coords, col = colore_atleta, pch = 16, cex = 1.2)
}

mean_normoxia <- mean(data$HRmax, na.rm = TRUE)
mean_hypoxia  <- mean(data$hHRmax, na.rm = TRUE)

# Disegnamo la linea nera marcata sopra tutti gli altri spaghetti
lines(c(1, 2), c(mean_normoxia, mean_hypoxia), col = "black", lwd = 4.0, lty = 2)

# Aggiungiamo due punti quadrati o cerchiati grandi per la media
points(c(1, 2), c(mean_normoxia, mean_hypoxia), col = "black", pch = 16, bg = "black", cex = 1.8)

