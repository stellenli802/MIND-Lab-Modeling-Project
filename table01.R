## R script to create "Table 1" for ENROLL-HD data for Modeling TTD paper

library(tidyverse)
library(tableone)
library(knitr)


## data and directory definitions
ENROLL <- readRDS("G:/projects/tpgarcia/enroll-hd/02-Kyle-Comparing-Time-To-Diagnosis-Models/enrollhd_baseline_survival_data.rds") %>%
  ungroup()

enrollhd_data <- ENROLL %>%
  filter(W != 0) %>%
  filter(complete.cases(select(., W, delta, CAG, age_0,
                               diagconf, motscore, scnt1, swrt1,
                               sit1, sdmt1))) %>%
  filter(CAG > 40 & CAG < 57) %>%
  mutate(
    W_age = age_0 + W,
    motscore2 = motscore^2,
    CAP = age_0 * (CAG - 34),
    diagconf1 = ifelse(diagconf == 1, 1, 0),
    diagconf2 = ifelse(diagconf == 2, 1, 0),
    diagconf3 = ifelse(diagconf == 3, 1, 0)
  )



## summary table via tableone with stratification by censoring status
tab1_strat_delta <- CreateTableOne(vars = c("age_0", "motscore", "sdmt1", "diagconf",
                                "swrt1", "scnt1", "sit1", "sex", "race", "CAG"),
                       strata = "delta",
                       data = enrollhd_data,
                       factorVars = c("diagconf", "sex", "race"))
tab1_strat_delta_print <- print(tab1_strat_delta, printToggle = FALSE, showAllLevels = TRUE)
colnames(tab1_strat_delta_print) <- c("level", "censored", "observed", "p", "test")
rownames(tab1_strat_delta_print) <- c("n", "Age at study entry", "Motor score", "SDMT", "Diagnostic Confidence Level", "", "", "", 
                                      "Stroop Word Reading Test", "Stroop Color Naming Test", "Stroop Interference Test",
                                      "Sex", "", "Race", "", "CAG")
kab_t1_str <- kable(tab1_strat_delta_print[, 1:3], format = "latex", booktabs = TRUE, caption = "Table 1: Baseline Characteristics Stratified by Censoring")

## summary table via tableone with no stratification
tab1 <- CreateTableOne(vars = c("age_0", "motscore", "sdmt1", "diagconf",
                                "swrt1", "scnt1", "sit1", "sex", "race", "CAG"),
                       data = enrollhd_data,
                       factorVars = c("diagconf", "sex", "race"))
tab1_print <- print(tab1, printToggle = FALSE, showAllLevels = TRUE)
colnames(tab1_print) <- c("level", "overall")
rownames(tab1_print) <- c("n", "Age at study entry", "Motor score", "SDMT", "Diagnostic Confidence Level", "", "", "", 
                                      "Stroop Word Reading Test", "Stroop Color Naming Test", "Stroop Interference Test",
                                      "Sex", "", "Race", "", "CAG")
kab_t1 <- kable(tab1_print, format = "latex", booktabs = TRUE, caption = "Table 1: Baseline Characteristics")

