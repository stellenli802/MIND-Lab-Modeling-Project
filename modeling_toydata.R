# set working directory
setwd("../Desktop")

# load packages
library(tidyverse)
library(survival)

# load the toy data set
toydata = readRDS("toydata.RDS")

# what are the column names of the toy data set?
colnames(toydata)

# what does the survival response look like?
Surv(toydata$W, toydata$delta)

# Cox proportional hazards model
coxph(formula = Surv(W, delta) ~ age_0 + sdmt1,
      data = toydata)
?coxph

# AFT model
survreg(formula = Surv(W, delta) ~ age_0 + sdmt1,
        data = toydata, dist = "weibull")
?survreg