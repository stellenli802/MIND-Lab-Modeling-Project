setwd("~/Desktop/MindLabWork")
library(tidyverse)
library(survival)

toydata <- readRDS("toy_data.rds")
# Full equation consists of interaction term: Age0 +CAG+(Age0×CAG)
aft_full <- survreg(Surv(W, delta) ~ age_0 * allele2l, data = toydata, dist = "loglogistic") 
summary(aft_full)

# I am also fitting this model based ont the CAP score and CAP-scaled score.
toydata$CAP <- toydata$age_0 * (toydata$allele2l - 33.66)

# Fit AFT model using CAP
aft_cap <- survreg(Surv(W, delta) ~ CAP, data = toydata, dist = "loglogistic")
summary(aft_cap)

# Fit AFT model using CAP-scaled
toydata$CAPs <- toydata$CAP/432.3326 #empirically derived CAP value from PREDICT-HD 
aft_caps <- survreg(Surv(W, delta) ~ CAPs, data = toydata, dist = "loglogistic")
summary(aft_caps)