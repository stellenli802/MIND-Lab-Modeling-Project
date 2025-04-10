setwd("~/Desktop/MIND Lab")
library(tidyverse)
library(survival)

toydata <- readRDS("toydata.rds")
#column names for the data (X, C, W, delta, age_0, sdmt1, allele2l, motscore, swrt1, scnt1, sit1, diagconf):
#colnames(toydata)

#Creating the CAP covariate
toydata$CAP <- toydata$age_0 * (toydata$allele2l - 34)

#Fitting Cox regression model on TMS+SDMT+CAP
fit_PI_HD = coxph(formula = Surv(W, delta) ~ motscore + sdmt1 + CAP, data=toydata)
summary(fit_PI_HD)

#Computing PI_HD and normalizing for each individual
#PI_HD = 51*TMS + (-34)*SDMT + 7*Age*(CAG-34)
PI_HD = 51*toydata$motscore + (-34)*toydata$sdmt1 + 7*toydata$age_0*(toydata$allele2l-34)

PIN_HD <- (PI_HD - 883)/1044
PIN_HD