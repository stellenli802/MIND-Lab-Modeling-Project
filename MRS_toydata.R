setwd("~/Desktop/MIND Lab")
library(tidyverse)
library(survival)

toydata <- readRDS("toydata.rds")
#column names for the data (X, C, W, delta, age_0, sdmt1, allele2l, motscore, swrt1, scnt1, sit1, diagconf):
#colnames(toydata)

#Factoring DCL into 3 levels
toydata$diagconf <- factor(toydata$diagconf, levels = c("DCL1", "DCL2", "DCL3"), ordered = TRUE)

#Fitting the Cox regression model
model_MRS <- coxph(Surv(W, delta) ~ diagconf + motscore + I(motscore^2) + scnt1 + swrt1
                   + sit1 + sdmt1 + allele2l + age_0 + allele2l:motscore + allele2l:age_0, data=toydata)
summary(model_MRS)