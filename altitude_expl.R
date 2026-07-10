# ---------------------------------------------------------------------------- #
# ALTITUDE EXPLORATION
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

plot_altitude_levels <- function(df, cat_var, cont_var) {
  # Assicuriamoci che i nomi delle colonne siano gestiti correttamente
  # cat_var: nome della colonna categoriale (es. "altitudine")
  # cont_var: nome della colonna continua (es. "valore")
  
  # --- 1. AGGREGAZIONE E CALCOLO STATISTICHE ---
  # Calcoliamo media e SD per ogni livello
  means_df <- aggregate(df[[cont_var]] ~ df[[cat_var]], FUN = mean)
  sds_df   <- aggregate(df[[cont_var]] ~ df[[cat_var]], FUN = sd)
  
  # Estraiamo i vettori ordinati
  levels_x <- as.numeric(as.character(means_df[, 1])) # Convertiamo in numerico per l'asse X
  means    <- means_df[, 2]
  sds      <- sds_df[, 2]
  
  # --- 2. CONFIGURAZIONE SPAZI E MARGINI ---
  old_par <- par(no.readonly = TRUE) # Salviamo i parametri grafici attuali
  par(mar = c(5, 5, 4, 2) + 0.1)     # Margini coerenti con il tuo standard
  
  # Calcoliamo i limiti degli assi per non tagliare le barre di errore
  y_min <- min(means - sds, na.rm = TRUE) * 0.9
  y_max <- max(means + sds, na.rm = TRUE) * 1.1
  
  # Definiamo i limiti X lasciando un po' di spazio laterale (es. 500 unità a destra e sinistra)
  x_min <- min(levels_x) - 500
  x_max <- max(levels_x) + 500
  
  # --- 3. COSTRUZIONE DEL GRAFICO ---
  # Inizializziamo il grafico vuoto (type = "n")
  plot(levels_x, means, 
       type = "n", 
       xaxt = "n", # Disattiviamo l'asse X di default per personalizzarlo
       xlim = c(x_min, x_max),
       ylim = c(y_min, y_max),
       main = "hSpO2end at different levels of altitude",
       xlab = "Altitude (m)", 
       ylab = paste(cont_var, "(%)"),
       font.lab = 1, 
       las = 1)
  
  # Disegniamo l'asse X personalizzato esattamente sui tuoi tre livelli
  axis(side = 1, at = levels_x, labels = levels_x, lwd = 1, lwd.ticks = 1)
  
  # 1. Barre di errore (SD) usando arrows() con i tuoi stessi parametri
  arrows(x0 = levels_x, y0 = means - sds, 
         x1 = levels_x, y1 = means + sds, 
         code = 3, 
         angle = 90, 
         length = 0.05,
         lwd = 1.2)
  
  # 2. Punti delle medie (pch = 19)
  points(levels_x, means, pch = 19, cex = 1)
  
  # Ripristiniamo i parametri grafici originali di R
  par(old_par)
}

# ---------------------------------------------------------------------------- #
# Risultati

data <- altitude_data %>% mutate(altitude = as.factor(altitude))
plot_altitude_levels(data, "altitude", "hSpO2end")