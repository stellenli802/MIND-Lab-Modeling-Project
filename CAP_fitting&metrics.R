setwd("~/Desktop/MindLabWork")
library(tidyverse)
library(survival)
library(pec)
library(rms)
library(timeROC)

toydata <- readRDS("toydata.rds")
# Full equation consists of interaction term: Age0 + (Age0:CAG)
toydata$age_allele2l <- toydata$age_0 * toydata$allele2l
aft_full <- survreg(Surv(W, delta) ~ age_0 + toydata$age_allele2l, data = toydata, dist = "loglogistic") 
summary(aft_full)

# I am also fitting this model based on the CAP score and CAP-scaled score.
toydata$CAP <- toydata$age_0 * (toydata$allele2l - 33.66)

# Fit AFT model using CAP
aft_cap <- survreg(Surv(W, delta) ~ CAP, data = toydata, dist = "loglogistic")
summary(aft_cap)

# Fit AFT model using CAP-scaled
toydata$CAPs <- toydata$CAP/432.3326 #empirically derived CAP value from PREDICT-HD 
aft_caps <- survreg(Surv(W, delta) ~ CAPs, data = toydata, dist = "loglogistic")
summary(aft_caps)

# Harrell's C (seeing how many pairs are concordant: the subject predicted to fail later actually fails later)
lp <- predict(aft_full, type = "lp")
conc <- survConcordance(Surv(W, delta) ~ lp, data = toydata)
c_index <- conc$concordance
c_index

# Compute the time-dependent ROC object at time t = 1 (How well can my model predict who will experience the event by time = 1?)
risk_score <- predict(aft_full, type = "lp")
roc1 <- timeROC(
  T      = toydata$W,        
  delta  = toydata$delta,    
  marker = risk_score,       
  cause  = 1,                
  times  = c(1,2,3),               
  iid    = FALSE             
)
auc_at_1 <- roc1$AUC[1]
auc_at_1

# Calibration Plots (compare the predicted survival probability with actual observed survival rate)
dd <- datadist(toydata)
options(datadist = "dd")
aft_full_1 <- psm(Surv(W, delta) ~ age_0 + age_allele2l, data = toydata, dist = "loglogistic") 

t_pred <- 1

calib_obj <- pec(
  object = aft_full_1,
  formula = Surv(W, delta) ~ 1,     
  data = toydata,
  times = t_pred,
  exact = TRUE,
  cens.model = "marginal",          
)

plot(calib_obj,
     xlab = "Predicted survival probability at time t=1",
     ylab = "Observed survival probability at time t=1",
     legend = FALSE
     )
title("Calibration Plot at Time t = 1")

legend("topright",
       legend = c("Reference", "AFT Model (Full)"),
       col = c("black", "red"),
       lwd = 2)

# Brier score (deviation of predicted survival probabilities at a specific time)
bs <- pec(
  object = list("AFT Model (Full)" = aft_full_1),
  data = toydata,
  formula = Surv(W,delta) ~ 1,
  times = c(1,2,3),
  exact = TRUE,
  cens.model = "marginal",
)
bs
idx <- which.min( abs(bs$time - 1) )   # index of the closest time to 1
bs$time[idx] 
brier_t1 <- bs$AppErr[["AFT Model (Full)"]][ idx ]
brier_t1
