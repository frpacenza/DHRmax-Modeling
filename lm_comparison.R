# ---------------------------------------------------------------------------- #
# COMPARISON OF LINEAR MODELS WITH THE SAME NUMBER OF PREDICTORS
# ---------------------------------------------------------------------------- #
# Import data and functions

load("~/Desktop/Project/predictors_df_with_target.RData")
source("~/Desktop/Project/RFiles/functions.R")

# ---------------------------------------------------------------------------- #
# Data organization (full cohort)

original_predictors <- predictors_df_with_target %>% select(
  -ID, -subject, -subjectID, -session, -sex, -hHRmax, -DHRmax, -Lapeak
)
complete_idx <- which(complete.cases(original_predictors))
predictors <- original_predictors[complete_idx,]

original_target <- predictors_df_with_target$DHRmax
target <- original_target[complete_idx]

physiology_predictors <- as.data.frame(predictors[, c("rVO2max", "hSpO2rest", "SpO2max")])

# ---------------------------------------------------------------------------- #
# Results and latex code (full cohort)

data_lm_3 <- data_lm(predictors, target, n = 3)
physio_lm_3 <- physio_lm(physiology_predictors, target, n = 3)
save(data_lm_3, file = "~/Desktop/Project/Results/data_lm_3.RData")
save(physio_lm_3, file = "~/Desktop/Project/Results/physio_lm_3.RData")
to_latex_table(data_lm_3, physio_lm_3, 3)
to_latex_plot(data_lm_3, physio_lm_3, 3)
plot_rmse_2(data_lm_3, physio_lm_3, 3, filename = "plot_rmse_n3.pdf")

data_lm_2 <- data_lm(predictors, target, n = 2)
physio_lm_2 <- physio_lm(physiology_predictors, target, n = 2)
save(data_lm_2, file = "~/Desktop/Project/Results/data_lm_2.RData")
save(physio_lm_2, file = "~/Desktop/Project/Results/physio_lm_2.RData")
to_latex_table(data_lm_2, physio_lm_2, 2)
plot_rmse_2(data_lm_2, physio_lm_2, 2, filename = "plot_rmse_n2.pdf")

data_lm_1 <- data_lm(predictors, target, n = 1)
physio_lm_1 <- physio_lm(physiology_predictors, target, n = 1)
save(data_lm_1, file = "~/Desktop/Project/Results/data_lm_1.RData")
save(physio_lm_1, file = "~/Desktop/Project/Results/physio_lm_1.RData")
to_latex_table(data_lm_1, physio_lm_1, 1)
plot_rmse_2(data_lm_1, physio_lm_1, 1, filename = "plot_rmse_n1.pdf")

# ---------------------------------------------------------------------------- #
# Data organization (male cohort)

original_predictors_M <- predictors_df_with_target %>% filter(
  sex == "M"
) %>% select(
  -ID, -subject, -subjectID, -session, -sex, -hHRmax, -DHRmax, -Lapeak
)
complete_idx_M <- which(complete.cases(original_predictors_M))
predictors_M <- original_predictors_M[complete_idx_M,]

original_target_M <- predictors_df_with_target %>% filter(
  sex == "M"
) %>% pull(DHRmax)
target_M <- original_target_M[complete_idx_M]

physiology_predictors_M <- as.data.frame(predictors_M[, c("rVO2max", "hSpO2rest", "SpO2max")])

# ---------------------------------------------------------------------------- #
# Results and latex code (male cohort)

data_lm_3_M <- data_lm(predictors_M, target_M, n = 3)
physio_lm_3_M <- physio_lm(physiology_predictors_M, target_M, n = 3)
save(data_lm_3_M, file = "~/Desktop/Project/Results/data_lm_3_M.RData")
save(physio_lm_3_M, file = "~/Desktop/Project/Results/physio_lm_3_M.RData")
to_latex_table(data_lm_3_M, physio_lm_3_M, 3)
plot_rmse_2(data_lm_3_M, physio_lm_3_M, 3, filename = "plot_rmse_n3_M.pdf")

data_lm_2_M <- data_lm(predictors_M, target_M, n = 2)
physio_lm_2_M <- physio_lm(physiology_predictors_M, target_M, n = 2)
save(data_lm_2_M, file = "~/Desktop/Project/Results/data_lm_2_M.RData")
save(physio_lm_2_M, file = "~/Desktop/Project/Results/physio_lm_2_M.RData")
to_latex_table(data_lm_2_M, physio_lm_2_M, 2)
plot_rmse_2(data_lm_2_M, physio_lm_2_M, 2, filename = "plot_rmse_n2_M.pdf")

data_lm_1_M <- data_lm(predictors_M, target_M, n = 1)
physio_lm_1_M <- physio_lm(physiology_predictors_M, target_M, n = 1)
save(data_lm_1_M, file = "~/Desktop/Project/Results/data_lm_1_M.RData")
save(physio_lm_1_M, file = "~/Desktop/Project/Results/physio_lm_1_M.RData")
to_latex_table(data_lm_1_M, physio_lm_1_M, 1)
plot_rmse_2(data_lm_1_M, physio_lm_1_M, 1, filename = "plot_rmse_n1_M.pdf")
