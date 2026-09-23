## STAT 7900 Final Project / SEASUG 2026 Paper #103
## Causal effect of maternal smoking during pregnancy on preterm birth
## 2023 NVSS/NCHS U.S. natality data (singleton births)
## Data source: https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(MatchIt); library(tidyr)
})
options(width = 120)
set.seed(79001)

dat0 <- read.csv("nvss2023_extract.csv", stringsAsFactors = FALSE,
                  colClasses = c(dmar = "character"))
cat("Raw extract rows (singleton, known gestation, known smoking):", nrow(dat0), "\n")

## ---- recode ----
dat0 <- dat0 %>%
  mutate(
    smoker    = ifelse(cig_rec == "Y", 1L, 0L),
    preterm   = ifelse(oegest_comb < 37, 1L, 0L),
    ga_weeks  = oegest_comb,
    mager     = as.numeric(mager),
    bmi       = ifelse(bmi >= 99.9, NA, as.numeric(bmi)),
    meduc     = ifelse(meduc == 9, NA, meduc),
    mrace6    = factor(mrace6, levels = 1:6,
                        labels = c("White","Black","AIAN","Asian","NHOPI","MultiRace")),
    meduc_f   = factor(meduc, levels = 1:8,
                        labels = c("<9th","9-12,no dip","HS grad/GED","Some college",
                                   "Associate","Bachelor's","Master's","Doct/Prof")),
    precare5  = ifelse(precare5 == 5, NA, precare5),
    precare_f = factor(precare5, levels = 1:4,
                        labels = c("1st tri","2nd tri","3rd tri","No care")),
    wic       = ifelse(wic == "U", NA, wic),
    wic01     = ifelse(wic == "Y", 1L, 0L),
    parity0   = ifelse(priorlive == 0 & priorterm == 0, 1L, 0L),  # nulliparous indicator
    rf_ppterm = ifelse(rf_ppterm == "U", NA, rf_ppterm),
    ppterm01  = ifelse(rf_ppterm == "Y", 1L, 0L),
    rf_pdiab  = ifelse(rf_pdiab == "U", NA, rf_pdiab),
    pdiab01   = ifelse(rf_pdiab == "Y", 1L, 0L),
    rf_phype  = ifelse(rf_phype == "U", NA, rf_phype),
    phype01   = ifelse(rf_phype == "Y", 1L, 0L),
    pay_f     = factor(pay_rec, levels = c(1,2,3,4,9),
                        labels = c("Medicaid","Private","Self-pay","Other","Unknown")),
    pay_f     = ifelse(pay_f == "Unknown", NA, as.character(pay_f))
  )

covs <- c("mager","mrace6","meduc_f","precare_f","wic01","parity0",
          "ppterm01","pdiab01","phype01","bmi","pay_f")

cc <- dat0 %>% select(smoker, preterm, ga_weeks, all_of(covs)) %>% na.omit()
cat("Complete cases on outcome + DAG-adjustment covariates:", nrow(cc),
    " (", round(100*nrow(cc)/nrow(dat0),1), "% of extract )\n")
cat("Dropped for missing covariate data:", nrow(dat0) - nrow(cc), "\n")

cc$mrace6  <- factor(cc$mrace6)
cc$meduc_f <- factor(cc$meduc_f)
cc$precare_f <- factor(cc$precare_f)
cc$pay_f   <- factor(cc$pay_f)

cat("\n=== Table 1 drivers: N, treated(smoker), preterm counts ===\n")
cat("N total:", nrow(cc), "\n")
cat("Smokers (T=1):", sum(cc$smoker==1), sprintf("(%.2f%%)", 100*mean(cc$smoker==1)), "\n")
cat("Non-smokers (T=0):", sum(cc$smoker==0), sprintf("(%.2f%%)", 100*mean(cc$smoker==0)), "\n")
cat("Preterm overall:", sum(cc$preterm==1), sprintf("(%.2f%%)", 100*mean(cc$preterm==1)), "\n")
tab <- cc %>% group_by(smoker) %>% summarise(n=n(), preterm_n=sum(preterm), rate=mean(preterm))
print(as.data.frame(tab))

saveRDS(cc, "cc.rds")
write.csv(tab, "table_crude_rates.csv", row.names = FALSE)
cat("\nStep 1 (load/clean/restrict) complete.\n")
