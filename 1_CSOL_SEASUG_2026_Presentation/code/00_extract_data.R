###############################################################################
## 00_extract_data.R
## STEP 0 (OPTIONAL): build the small analysis file from the full public data.
##
## YOU DO NOT NEED THIS STEP to reproduce the paper.  The repository already
## contains the small file this script creates (data/nvss2023_extract.csv.gz).
## Run this script only if you want to rebuild that file from the original
## government download, to check that nothing was altered.
##
## ---------------------------------------------------------------------------
## WHERE THE DATA COME FROM  (public-use, no registration, no personal identifiers)
## ---------------------------------------------------------------------------
## Source : U.S. National Center for Health Statistics (NCHS), National Vital
##          Statistics System, 2023 Natality Public Use File (all U.S. births).
## Official page (links to every year):
##          https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm
## Official file   (fixed-width text, 224 MB zip):
##          https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip
## Official user guide (variable definitions and codes):
##          https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/UserGuide2023.pdf
##
## This project used the COMMA-SEPARATED (CSV) version of the same file, with
## lower-case variable names, distributed by the National Bureau of Economic
## Research (NBER) at
##          https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip
## (about 198 MB zipped, 2.2 GB unzipped, 3.6 million rows, 200+ columns).
## Both files contain the same births.  NBER names are used throughout this code
## (for example mager = mother's age, oegest_comb = obstetric gestation in weeks).
##
## ---------------------------------------------------------------------------
## HOW TO USE
## ---------------------------------------------------------------------------
## 1. Download natality2023us.zip from the NBER link above and put it in the
##    folder  data/  (so the path is  data/natality2023us.zip).
##    Do NOT unzip it by hand; this script does it.  You need about 2.5 GB of
##    free disk space while it runs; the large temporary file is deleted at the end.
## 2. In R, set the working directory to the folder containing this file, then:
##        source("00_extract_data.R")        # or:  Rscript 00_extract_data.R
##    Needs only the data.table package.  Run time: 3 to 8 minutes.
## 3. Output: data/nvss2023_extract.csv  (24 columns, 3,476,180 births).
##
## If R cannot unzip the 2.2 GB file on your computer (a rare zip64 limitation),
## unzip it with 7-Zip or Windows Explorer instead, place natality2023us.csv in
## data/, and re-run this script; it will skip the unzip step.
###############################################################################

set.seed(1234)
suppressPackageStartupMessages(library(data.table))

ZIP_FILE <- "data/natality2023us.zip"
CSV_FILE <- "data/natality2023us.csv"
OUT_FILE <- "data/nvss2023_extract.csv"

## The 24 columns the analysis needs (meaning in comments)
keep <- c("dob_yy",       # year of birth (2023)
          "mager",        # mother's age in years
          "mrace6",       # mother's race: 6 groups
          "dmar",         # marital status
          "meduc",        # mother's education (8 levels, 9 = unknown)
          "precare5",     # month prenatal care began, in 4 groups (5 = unknown)
          "wic",          # participated in WIC (Y/N/U)
          "priorlive",    # previous births still living
          "priordead",    # previous births now deceased
          "priorterm",    # previous other terminations
          "rf_ppterm",    # prior preterm birth (Y/N/U)
          "rf_pdiab",     # pre-pregnancy diabetes (Y/N/U)
          "rf_phype",     # pre-pregnancy hypertension (Y/N/U)
          "rf_gdiab",     # gestational diabetes (not used as a confounder)
          "rf_ghype",     # gestational hypertension (not used as a confounder)
          "bmi",          # pre-pregnancy body mass index
          "cig_0",        # cigarettes per day before pregnancy
          "cig_rec",      # smoked at any time during pregnancy (Y/N)  = exposure T
          "dplural",      # plurality (1 = singleton)
          "sex",          # infant sex
          "oegest_comb",  # obstetric estimate of gestational age, weeks = outcome Y
          "dbwt",         # birth weight in grams
          "pay_rec",      # principal source of payment (insurance), sensitivity analysis only
          "restatus")     # residence status (4 = not a U.S. resident)

if (!file.exists(CSV_FILE)) {
  if (!file.exists(ZIP_FILE)) stop("Put natality2023us.zip in the data/ folder first (see the link above).")
  cat("Unzipping (this takes a few minutes)...\n")
  utils::unzip(ZIP_FILE, files = "natality2023us.csv", exdir = "data")
  made_csv <- TRUE
} else made_csv <- FALSE

cat("Reading the 24 needed columns...\n")
raw <- fread(CSV_FILE, select = keep, colClasses = list(character = c("wic", "rf_ppterm", "rf_pdiab",
             "rf_phype", "rf_gdiab", "rf_ghype", "cig_rec", "sex", "dmar")))
cat("Rows in the full 2023 file:", format(nrow(raw), big.mark = ","), "\n")

## Inclusion rules (the same ones described in the paper)
##   singleton births, known gestational age (oegest_comb < 99), known smoking status
raw <- raw[dplural == 1 & oegest_comb < 99 & cig_rec %in% c("Y", "N")]
cat("Rows after the rules:", format(nrow(raw), big.mark = ","), "(expected 3,476,180)\n")

fwrite(raw[, ..keep], OUT_FILE)
if (made_csv) file.remove(CSV_FILE)          # remove the 2.2 GB temporary file
cat("Wrote", OUT_FILE, "\nDone.\n")
