library(survival)
library(tidyverse)

# simulate data from a proportional hazards model
# using a Weibull distribution
# n: sample size
# beta: regression coefficients on covariates
# shape: shape parameter of Weibull distribution
# rate_C: rate parameter for censoring variable C
simulate_weibull_data <- function(n, beta, shape, rate_C) {
  # Generate covariates
  Age <- rnorm(n)  # Age
  SDMT <- rnorm(n) # Symbol Digit Modalities Test score
  CAG <- rnorm(n) # CAG repeat length
  TMS <- rnorm(n)  # Total Motor Score
  Stroop_Word <- rnorm(n)   # Stroop Word Test score
  Stroop_Color <- rnorm(n)  # Stroop Color Test score
  Stroop_Interf <- rnorm(n) # Stroop Interference Test score
  DCL <- sample(c("DCL1", "DCL2", "DCL3"), n, replace = TRUE)  # diagnostic confidence level
  
  # Quadratic term for Total Motor Score
  TMS_sq <- TMS^2
  
  # Interaction terms
  CAG_TMS <- CAG * TMS
  CAG_Age <- CAG * Age
  
  # Convert levels of diagnostic confidence levels to binary variables
  DCL_2 <- as.numeric(DCL == "DCL2")
  DCL_3 <- as.numeric(DCL == "DCL3")
  
  # Design matrix (including categorical variables and higher order terms)
  design_matrix <- cbind(1, Age, SDMT, CAG, TMS, TMS_sq, CAG_TMS, CAG_Age, 
                         Stroop_Word, Stroop_Color, Stroop_Interf, DCL_2, DCL_3)  
  linear_predictor <- design_matrix %*% beta
  
  # Generate uniform random variable for inverse transform
  U <- runif(n)
  
  # Generate event times X using the proportional hazards model
  X <- (-log(U) / exp(linear_predictor))^(1 / shape)
  
  # Generate censoring times
  C <- rexp(n = n, rate = rate_C)
  
  # Observed time and censoring indicator
  W <- pmin(X, C)
  delta <- as.numeric(X <= C) 
  
  # Create a data frame
  df <- data.frame(X, C, W, delta, Age, SDMT, CAG, TMS, 
                   Stroop_Word, Stroop_Color, Stroop_Interf, DCL)
  
  return(df)
}

# Example usage
n <- 10000
beta <- rep(0.05, 13)  # Regression coefficients
shape <- 1.5  # Weibull shape parameter

# Simulate proportional hazards data
sim_data <- simulate_weibull_data(n, beta, shape, rate_C = 2.75)
head(sim_data)

# Check censoring rate
mean(sim_data$delta)

# Estimate parameters of Cox PH model
fit_results = coxph(formula = Surv(W, delta) ~ Age + SDMT + CAG + TMS + TMS^2 + CAG*TMS + CAG*Age + Stroop_Word + Stroop_Color + Stroop_Interf + factor(DCL),
                      data = sim_data)

# Check parameter estimates
fit_results$coefficients

# Rename variables to match ENROLL column names
sim_data = sim_data %>%
  rename(age_0 = Age,
         sdmt1 = SDMT,
         motscore = TMS,
         scnt1 = Stroop_Color,
         swrt1 = Stroop_Word,
         sit1 = Stroop_Interf,
         diagconf = DCL,
         allele2l = CAG)

head(sim_data)