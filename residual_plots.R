# ---------------------------------------------------------------------------- #
# RESIDUAL PLOTS
# ---------------------------------------------------------------------------- #
# Import data and functions

load("~/Desktop/Project/Data/predictors_df_with_target.RData")
# ---------------------------------------------------------------------------- #
# Simple and Complex Models

data <- predictors_df_with_target %>%
  select(-Lapeak) %>%
  na.omit()
original_predictors <- predictors_df_with_target %>% select(
  -ID, -subject, -subjectID, -session, -sex, -hHRmax, -DHRmax, -Lapeak
)
complete_idx <- which(complete.cases(original_predictors))
predictors <- original_predictors[complete_idx,]

original_target <- predictors_df_with_target$DHRmax
target <- original_target[complete_idx]

simple <- lm(DHRmax ~ VEmax + age, data = data)
plot(simple, which = 1, sub = "")
title(main = "Residual plot for the Simple Model")

complex <- lm(DHRmax ~ VEmax + VErvt2 + age + BMI + Prvt2, data = data)
plot(complex, which = 1, sub = "")
title(main = "Residual plot for the Complex Model")

# ---------------------------------------------------------------------------- #
# Additive and Interaction Models
# ---------------------------------------------------------------------------- #
# Libraries, functions and data loading:

load("~/Desktop/Project/Data/multistudy_df.RData")

library(dplyr)
library(tidyr)
library(caret)

# ---------------------------------------------------------------------------- #
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

addittive <- lm(DHRmax ~ VEmax + altitude, data = data)
plot(addittive, which = 1, sub = "")
title(main = "Residual plot for the Additive Model")

interaction <- lm(DHRmax ~ VEmax * altitude, data = data)
plot(interaction, which = 1, sub = "")
title(main = "Residual plot for the Interaction Model")
