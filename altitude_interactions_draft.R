load("~/Desktop/Project/Data/multistudy_df.RData")
source("~/Desktop/Project/RFiles/altitude_models.R")

library(dplyr)
library(tidyr)
library(caret)

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

model_id <- "Altitude Only Model"
covariates <- c("altitude")
alt_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Additive Model"
covariates <- c("altitude", "VEmax")
add_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

model_id <- "Simple Model"
covariates <- c("altitude", "VEmax", "age")
simple_tab <- lm_table(data = data, model_id = model_id, covariates = covariates)

mourot_tab <- mourot_table(data = data %>% mutate(altitude = as.numeric(as.character(altitude))))

altitude_tab <- rbind(alt_tab, add_tab, simple_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)


model_id <- "Specific-4500 Model"
covariates <- c("VEmax")
s4500_tab <- lm_table(data = data %>% filter(altitude == 4500), model_id = model_id, covariates = covariates)

model_id <- "Specific-3500 Model"
covariates <- c("VEmax")
s3500_tab <- lm_table(data = data %>% filter(altitude == 3500), model_id = model_id, covariates = covariates)

model_id <- "Specific-4500 Model"
covariates <- c("VEmax")
s4500_tab <- lm_table(data = data %>% filter(altitude == 4500), model_id = model_id, covariates = covariates)

altitude_tab <- rbind(alt_tab, s2500_tab, s3500_tab, s4500_tab, mourot_tab)
to_latex(altitude_tab)
plot_rmse(altitude_tab)

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

results_tab <- summary_lmcv(lm_cv = lm_cv, model_id = model_id)
to_latex(results_tab)

colors <- c("lightblue", "blue", "darkblue")
altitude_colors <- colors[factor(data$altitude, levels = c(2500, 3500, 4500))]
plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Interaction Model regression lines"
)
abline(a = lm_cv$finalModel$coefficients[1], 
       b = lm_cv$finalModel$coefficients[2],
       col = "lightblue", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[3],
       b = lm_cv$finalModel$coefficients[2], 
       col = "blue", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[4], 
       b = lm_cv$finalModel$coefficients[2],
       col = "darkblue", lwd = 1)

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

results_tab <- summary_lmcv(lm_cv = lm_cv, model_id = model_id)
to_latex(results_tab)

plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Interaction Model regression lines"
)
abline(a = lm_cv$finalModel$coefficients[1], 
       b = lm_cv$finalModel$coefficients[2],
       col = "lightblue", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[3],
       b = lm_cv$finalModel$coefficients[2] + lm_cv$finalModel$coefficients[5], 
       col = "blue", lwd = 1)
abline(a = lm_cv$finalModel$coefficients[1] + lm_cv$finalModel$coefficients[4], 
       b = lm_cv$finalModel$coefficients[2] + lm_cv$finalModel$coefficients[6],
       col = "darkblue", lwd = 1)

# ---------------------------------------------------------------------------- #
# Specific Models:

model_id <- "Specific-2500 Model"
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
  data = target <- data %>% filter(altitude == 2500),
  method = "lm",
  trControl = cv_control
)
s2500_tab <- summary_lmcv(lm_cv = lm2500_cv, model_id = model_id)

model_id <- "Specific-3500 Model"
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
  data = target <- data %>% filter(altitude == 3500),
  method = "lm",
  trControl = cv_control
)
s3500_tab <- summary_lmcv(lm_cv = lm3500_cv, model_id = model_id)

model_id <- "Specific-4500 Model"
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
  data = target <- data %>% filter(altitude == 4500),
  method = "lm",
  trControl = cv_control
)
s4500_tab <- summary_lmcv(lm_cv = lm4500_cv, model_id = model_id)

plot(data$VEmax, data$DHRmax, 
     pch = 19,
     col = altitude_colors,
     xlab = "VEmax (L/min)",
     ylim = c(min(data$DHRmax, na.rm = T), 
              max(data$DHRmax, na.rm = T) * 1.15),
     ylab = "DHRmax (bpm)",
     main = "Interaction Model regression lines"
)
abline(a = lm2500_cv$finalModel$coefficients[1], 
       b = lm2500_cv$finalModel$coefficients[2],
       col = "lightblue", lwd = 1)
abline(a = lm3500_cv$finalModel$coefficients[1],
       b = lm3500_cv$finalModel$coefficients[2],
       col = "blue", lwd = 1)
abline(a = lm4500_cv$finalModel$coefficients[1], 
       b = lm4500_cv$finalModel$coefficients[2],
       col = "darkblue", lwd = 1)



