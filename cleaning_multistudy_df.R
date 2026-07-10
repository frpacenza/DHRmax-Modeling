# ---------------------------------------------------------------------------- #
# RAW MULTI-STUDY DATASET CLEANING
# ---------------------------------------------------------------------------- #
# Pulizia e riorganizzazione di Raw Multi-Study Dataset.
# Carico in R il file spreadsheet .xlsx e lo salvo in formato .RData.
# Definisco un nuovo ID alfanumerico che indica sia soggetto che studio. 
# Elimino le colonne ridontanti e/o non di interesse. 
# Per coerenza su tutto il progetto, rinomino le colonne restanti e aggiusto 
# (arrotondo) le variabili antropometriche.

# ---------------------------------------------------------------------------- #
# Libraries:

library(readxl)
library(dplyr)

# ---------------------------------------------------------------------------- #
# Pulizia Raw Multi-Study Dataset

raw_multistudy_df <- read_excel("Desktop/Project/Data/raw_multistudy_dataset.xlsx")
load("~/Desktop/Project/Data/multistudy_prj_map.RData")

# ID: S + numero unico per soggetto + abbreviazione project name
newID <- raw_multistudy_df %>% 
  rename(project = PROJ) %>%
  left_join(multistudy_prj_map, by = "project") %>%
  mutate(ID = paste0("S0", ID, projectID)) %>%
  pull(ID)

# aggiungo il nuovo ID, rinomino e seleziono le variabili che mi interessano
multistudy_df <- raw_multistudy_df %>%
  mutate(
    ID = newID
  ) %>%
  rename(
    method = "MODALITà",
    condition = COND,
    altitude = QUOTA,
    sex = SEX,
    LArest = "LA REST",
    HRrest = "HR REST",
    SpO2rest = "SPO2 REST",
    RFrest = "RF REST",
    VErest = "VE REST",
    HRmax = "HRMAX...22",
    VO2max = "VO2max...23",
    rVO2max = "VO2max/kg...24",
    Pmax = "PERF...25",
    rPmax = "PERF/kg...26",
    LApeak = "Lapeak...27",
    RFmax = "RF max...28",
    VEmax = "VE max...29",
    SpO2end = "SPO2 END...30",
    Pvt2 = "Perf VT2...31",
    SpO2vt2 = "SPO2 VT2...32",
    VEvt2 = "VEVT2...33",
    RFvt2 = "RfVT2...34",
    VO2vt2 = "VO2VT2...35",
    rVO2vt2 = "VO2VT2/kg...36",
    VO2rvt2 = "VO2VT2%MAX...37",
    HRvt2 = "HRVT2...38",
    HRrvt2 = "HRVT2%MAX...39",
    Pvt1 = "PERF VT1...40",
    SpO2vt1 = "SPO2 VT1...41",
    VEvt1 = "VEVT1...42",
    RFvt1 = "RfVT1...43",
    VO2vt1 = "VO2VT1...44",
    rVO2vt1 = "VO2VT1/kg...45", 
    VO2rvt1 = "VO2VT1%MAX...46",
    HRvt1 = "HRVT1...47",
    HRrvt1 = "HRVT1%MAX...48"
  ) %>%
  # per coerenza con il dataset predictors_df_with_target aggiusto age, height, weight e BMI
  mutate( 
    age = floor(age),
    height = floor(height + 0.5), # trucco per arrodntare il .5 per eccesso (stranamente round() non lo fa...)
    weight = floor(weight + 0.5),
    HRmax = floor(HRmax + 0.5)
  ) %>%
  mutate(
    BMI = weight / (height/100)^2
  ) %>%
  select(
    ID, condition, altitude, method, # protocol variables
    sex, age, height, weight, BMI,   # anthropometric variables
    LArest, LApeak,                  # exercise variables
    HRrest, HRvt1, HRrvt1, HRvt2, HRrvt2, HRmax, 
    SpO2rest, SpO2vt1, SpO2vt2, SpO2end,
    RFrest, RFvt1, RFvt2, RFmax,
    VErest, VEvt1, VEvt2, VEmax,
    VO2vt1, rVO2vt1, VO2rvt1, VO2vt2, rVO2vt2, VO2rvt2, VO2max, rVO2max,
    Pvt1, Pvt2, Pmax, rPmax
  ) %>%
  # trasformo variabili in factor
  mutate(
    condition = factor(condition),
    method = factor(method),
    sex = factor(sex)
  )
save(multistudy_df, file = "~/Desktop/Project/Data/multistudy_df.RData")