# ---------------------------------------------------------------------------- #
# ALTITUDE: NON-ADDITIVE AND INTERACTION EFFECTS
# ---------------------------------------------------------------------------- #
# Libraries, functions and data loading:

load("~/Desktop/Project/Data/multistudy_df.RData")
source("~/Desktop/Project/RFiles/altitude_models.R")

library(dplyr)
library(tidyr)
library(caret)

# ---------------------------------------------------------------------------- #
# Data organization:

data <- multistudy_df %>%
  filter(sex == "M") %>%
  select(condition, ID, VEmax, altitude, HRmax, age, SpO2end) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(altitude, HRmax, VEmax, SpO2end)
  ) %>%
  rename(
    VEmax = VEmax_N,
    altitude = altitude_H,
    hSpO2end = SpO2end_H
  ) %>%
  mutate(
    DHRmax = HRmax_H - HRmax_N,
    altitude = factor(altitude)
  )  %>%
  select(
    ID, age, altitude, VEmax, DHRmax, hSpO2end
  ) %>%
  na.omit()

# ---------------------------------------------------------------------------- #
# Results
# ---------------------------------------------------------------------------- #
# Additive Model:

model_id <- "Additive Model"
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
  DHRmax ~ VEmax + altitude,
  data = data,
  method = "lm",
  trControl = cv_control
)

add_tab <- summary_lmcv(lm_cv = lm_cv, model_id = model_id)

colors <- c("#d8b4fe", "#a855f7", "#5b127d")
altitude_colors <- colors[factor(data$altitude, levels = c(2500, 3500, 4500))]
plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression lines for the Additive Model"
)
abline(a = lm_cv$finalModel$coefficients[1], 
       b = lm_cv$finalModel$coefficients[2],
       col = "#d8b4fe", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[3],
       b = lm_cv$finalModel$coefficients[2], 
       col = "#a855f7", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[4], 
       b = lm_cv$finalModel$coefficients[2],
       col = "#5b127d", lwd = 1)

# ---------------------------------------------------------------------------- #
# Interaction Model:

model_id <- "Interaction Model"
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
  DHRmax ~ VEmax * altitude,
  data = data,
  method = "lm",
  trControl = cv_control
)
int_tab <- summary_lmcv(lm_cv = lm_cv, model_id = model_id)

plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression lines for the Interaction Model"
)
abline(a = lm_cv$finalModel$coefficients[1], 
       b = lm_cv$finalModel$coefficients[2],
       col = "#d8b4fe", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[3],
       b = lm_cv$finalModel$coefficients[2] + lm_cv$finalModel$coefficients[5], 
       col = "#a855f7", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[4], 
       b = lm_cv$finalModel$coefficients[2] + lm_cv$finalModel$coefficients[6],
       col = "#5b127d", lwd = 1)

# ---------------------------------------------------------------------------- #
# Summary table

to_latex(rbind(alt_tab, add_tab, simple_tab, mourot_tab))
plot_rmse(rbind(alt_tab, add_tab, simple_tab, mourot_tab))

to_latex(rbind(add_tab, int_tab))
plot_rmse(rbind(alt_tab, add_tab, int_tab, simple_tab, mourot_tab))

# ---------------------------------------------------------------------------- #
# Specific Models:
# ---------------------------------------------------------------------------- #
# 2500-Specific Model:

model_id <- "2500-Specific Model"
model_id <- "2500-Sp. Model"
target <- data %>% filter(altitude == 2500) %>% pull(DHRmax)
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
lm2500_cv <- caret::train(
  DHRmax ~ VEmax,
  data = data %>% filter(altitude == 2500),
  method = "lm",
  trControl = cv_control
)
s2500_tab <- summary_lmcv(lm_cv = lm2500_cv, model_id = model_id)

plot(data %>% filter(altitude == 2500) %>% pull(VEmax),
     data %>% filter(altitude == 2500) %>% pull(DHRmax),
     pch = 19,
     col = "#d8b4fe",
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression line for the 2500-Specific Model"
)
abline(a = lm2500_cv$finalModel$coefficients[1], 
       b = lm2500_cv$finalModel$coefficients[2],
       col = "#d8b4fe", lwd = 1)

# ---------------------------------------------------------------------------- #
# 3500-Specific Model:

model_id <- "3500-Specific Model"
model_id <- "3500-Sp. Model"
target <- data %>% filter(altitude == 3500) %>% pull(DHRmax)
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
lm3500_cv <- caret::train(
  DHRmax ~ VEmax,
  data = data %>% filter(altitude == 3500),
  method = "lm",
  trControl = cv_control
)
s3500_tab <- summary_lmcv(lm_cv = lm3500_cv, model_id = model_id)

plot(data %>% filter(altitude == 3500) %>% pull(VEmax),
     data %>% filter(altitude == 3500) %>% pull(DHRmax),
     pch = 19,
     col = "#a855f7",
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression line for the 3500-Specific Model"
)
abline(a = lm3500_cv$finalModel$coefficients[1], 
       b = lm3500_cv$finalModel$coefficients[2],
       col = "#a855f7", lwd = 1)

# ---------------------------------------------------------------------------- #
# 4500-Specific Model:

model_id <- "4500-Specific Model"
model_id <- "4500-Sp. Model"
target <- data %>% filter(altitude == 4500) %>% pull(DHRmax)
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
lm4500_cv <- caret::train(
  DHRmax ~ VEmax,
  data = data %>% filter(altitude == 4500),
  method = "lm",
  trControl = cv_control
)
s4500_tab <- summary_lmcv(lm_cv = lm4500_cv, model_id = model_id)

plot(data %>% filter(altitude == 4500) %>% pull(VEmax),
     data %>% filter(altitude == 4500) %>% pull(DHRmax),
     pch = 19,
     col = "#5b127d",
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression line for the 4500-Specific Model"
)
abline(a = lm4500_cv$finalModel$coefficients[1], 
       b = lm4500_cv$finalModel$coefficients[2],
       col = "#5b127d", lwd = 1)

# ---------------------------------------------------------------------------- #
# All together now:

to_latex(rbind(s2500_tab, s3500_tab, s4500_tab, int_tab))
plot_rmse(rbind(s2500_tab, s3500_tab, s4500_tab, int_tab, mourot_tab))

plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Regression lines for the Altitude-Specific Models"
)
abline(a = lm2500_cv$finalModel$coefficients[1], 
       b = lm2500_cv$finalModel$coefficients[2],
       col = "#d8b4fe", lwd = 1)
abline(a = lm3500_cv$finalModel$coefficients[1],
       b = lm3500_cv$finalModel$coefficients[2],
       col = "#a855f7", lwd = 1)
abline(a = lm4500_cv$finalModel$coefficients[1], 
       b = lm4500_cv$finalModel$coefficients[2],
       col = "#5b127d", lwd = 1)

