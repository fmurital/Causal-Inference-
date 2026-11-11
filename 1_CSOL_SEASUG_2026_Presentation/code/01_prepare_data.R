###############################################################################
## 01_prepare_data.R
## STEP 1: read the analysis extract, clean it, and build the analytic cohort.
##
## ---------------------------------------------------------------------------
## THE DATA  (public-use; no registration; no personal identifiers)
## ---------------------------------------------------------------------------
## U.S. National Center for Health Statistics (NCHS), National Vital Statistics
## System, 2023 Natality Public Use File (every U.S. birth certificate of 2023).
##   Official page  : https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm
##   Official file  : https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip
##   User guide     : https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/UserGuide2023.pdf
##   CSV version    : https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip   (NBER mirror, 198 MB)
##
## The full file is 3.6 million births and 200+ columns (2.2 GB unzipped), too
## large to keep in a code repository.  This project therefore uses a small
## extract of 24 columns and 3,476,180 births (singletons with known gestational
## age and known smoking status).  It ships with this code as
##       data/nvss2023_extract.csv.gz          (36 MB, nothing to download)
## so you can run everything from this script onward immediately.
##
## To REBUILD the extract yourself from the original download, follow the
## instructions at the top of 00_extract_data.R (it checks that the rebuilt file
## matches this one).
##
## Input  : data/nvss2023_extract.csv.gz  (or the unzipped data/nvss2023_extract.csv)
## Output : output/analytic_cohort.rds  (the cleaned cohort used by every later step)
##          output/cohort_flow.csv      (how many births remain after each rule)
##          output/missingness.csv      (how much data is missing, by variable)
##
## HOW TO RUN (set the working directory to the folder that contains this file):
##     Rscript 01_prepare_data.R          # or open in RStudio and click "Source"
## Needs only the data.table package:  install.packages("data.table")
## Run time: about 1 minute.
###############################################################################

set.seed(1234)                       # project-wide seed (nothing random happens here)
source("utils_causal.R")

DATA_PATH <- "data/nvss2023_extract.csv"
## If only the compressed file is present, unzip it once (base R, works on Windows, Mac, Linux).
if (!file.exists(DATA_PATH)) {
  if (!file.exists("data/nvss2023_extract.csv.gz")) stop("Data file not found. See the DATA notes at the top of this script.")
  cat("Decompressing data/nvss2023_extract.csv.gz ...\n")
  inp <- gzfile("data/nvss2023_extract.csv.gz", "rb"); out <- file(DATA_PATH, "wb")
  repeat { chunk <- readBin(inp, "raw", 5e7); if (length(chunk) == 0) break; writeBin(chunk, out) }
  close(inp); close(out)
}
dir.create("output", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Read the raw extract
## ---------------------------------------------------------------------------
## Several columns hold letters (Y / N / U) so we read them as text.
raw <- fread(DATA_PATH, colClasses = list(character =
  c("dmar", "wic", "rf_ppterm", "rf_pdiab", "rf_phype", "rf_gdiab", "rf_ghype", "cig_rec", "sex")))
flow <- data.table(step = "Singleton live births with known gestational age and known smoking status (extract)",
                   n = nrow(raw))

## ---------------------------------------------------------------------------
## 2. Residency rule
## ---------------------------------------------------------------------------
## NCHS convention: births to mothers who are NOT U.S. residents (restatus == 4)
## are left out of national statistics.  We follow that convention.
raw <- raw[restatus != 4]
flow <- rbind(flow, data.table(step = "Exclude births to non-U.S.-resident mothers (restatus = 4)", n = nrow(raw)))

## ---------------------------------------------------------------------------
## 3. Recode every variable.   "Unknown" codes become NA (missing).
## ---------------------------------------------------------------------------
yn <- function(x) factor(ifelse(x == "Y", "Yes", ifelse(x == "N", "No", NA)), levels = c("No", "Yes"))

d <- raw[, .(
  ## ---- treatment and outcome ----
  smoker  = as.integer(cig_rec == "Y"),        # T: smoked any cigarettes during pregnancy
  preterm = as.integer(oegest_comb < 37),      # Y: delivered before 37 completed weeks
  ga_weeks = oegest_comb,                      # obstetric estimate of gestational age

  ## ---- the ten confounders in the final DAG ----
  mager,                                       # maternal age (years), kept for plots/continuous spec
  age_grp = cut(mager, c(-Inf, 19, 24, 29, 34, 39, Inf),
                labels = c("<20", "20-24", "25-29", "30-34", "35-39", "40+")),
  race    = factor(mrace6, levels = 1:6,
                   labels = c("White", "Black", "AIAN", "Asian", "NHOPI", "Multiracial")),
  educ    = factor(ifelse(meduc == 9, NA, meduc), levels = 1:8,
                   labels = c("<9th grade", "9-12th, no diploma", "HS grad/GED", "Some college",
                              "Associate", "Bachelor's", "Master's", "Doctorate/Prof.")),
  precare = factor(ifelse(precare5 == 5, NA, precare5), levels = 1:4,
                   labels = c("1st trimester", "2nd trimester", "3rd trimester", "No care")),
  wic     = yn(ifelse(wic == "U", NA, wic)),
  ## parity = number of previous live births (born alive, whether or not still living)
  parity  = {pl <- ifelse(priorlive == 99 | priordead == 99, NA, priorlive + priordead)
             factor(pmin(pl, 3), levels = 0:3, labels = c("0", "1", "2", "3+"))},
  ppterm  = yn(ifelse(rf_ppterm == "U", NA, rf_ppterm)),   # previous preterm birth
  pdiab   = yn(ifelse(rf_pdiab  == "U", NA, rf_pdiab)),    # pre-pregnancy diabetes
  phype   = yn(ifelse(rf_phype  == "U", NA, rf_phype)),    # pre-pregnancy hypertension
  bmi     = ifelse(bmi >= 99.9, NA, bmi),                  # 99.9 = unknown
  bmi_grp = cut(ifelse(bmi >= 99.9, NA, bmi), c(-Inf, 18.5, 25, 30, 35, 40, Inf), right = FALSE,
                labels = c("<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+")),

  ## ---- insurance (NOT in the final DAG; used only in a sensitivity analysis) ----
  pay     = factor(pay_rec, levels = c(1, 2, 3, 4, 9),
                   labels = c("Medicaid", "Private", "Self-pay", "Other", "Unknown"))
)]

## ---------------------------------------------------------------------------
## 4. Missing data summary, then complete-case restriction
## ---------------------------------------------------------------------------
adj_vars <- c("age_grp", "race", "educ", "precare", "wic", "parity", "ppterm", "pdiab", "phype", "bmi_grp")
miss <- data.table(variable = adj_vars,
                   n_missing = sapply(adj_vars, function(v) sum(is.na(d[[v]]))))
miss[, pct_missing := round(100 * n_missing / nrow(d), 2)]
fwrite(miss, "output/missingness.csv")
print(miss)

d <- d[complete.cases(d[, ..adj_vars])]
flow <- rbind(flow, data.table(step = "Keep births with complete data on all ten confounders (analytic cohort)", n = nrow(d)))
flow[, pct_of_extract := round(100 * n / flow$n[1], 2)]
fwrite(flow, "output/cohort_flow.csv")
print(flow)

## ---------------------------------------------------------------------------
## 5. Save the analytic cohort
## ---------------------------------------------------------------------------
saveRDS(d, "output/analytic_cohort.rds")
cat("\nAnalytic cohort:", format(nrow(d), big.mark = ","), "births;",
    format(sum(d$smoker), big.mark = ","), "smokers (",
    round(100 * mean(d$smoker), 2), "%);",
    format(sum(d$preterm), big.mark = ","), "preterm (", round(100 * mean(d$preterm), 2), "%)\n")
cat("Crude preterm risk: smokers", round(100 * mean(d$preterm[d$smoker == 1]), 2), "%  vs  non-smokers",
    round(100 * mean(d$preterm[d$smoker == 0]), 2), "%\n")
