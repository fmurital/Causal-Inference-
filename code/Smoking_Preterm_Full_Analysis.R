###############################################################################
##  Smoking in pregnancy and preterm birth
##  Doubly robust (AIPTW) analysis of the 2023 U.S. natality file
##
##  TO RUN: select everything (Ctrl+A) and run (Ctrl+Enter), or click Source.
##  Both data files are downloaded from this GitHub repository on the first
##  run. Tables are written to results/tables and figures to results/figures.
##  The random seed is 1234, so every run gives the same numbers.
##
##  SECTIONS
##    1   Load the data, build the analytic cohort
##    2   Data dictionary
##    3   Exploratory analysis, missing data, association matrix
##    4   Which confounders matter most
##    5   Main estimate, comparison estimators, adjustment-set checks
##    6   Overlap, weights and covariate balance
##    7   Variance: influence function versus full-sample bootstrap
##    8   Flexible models with cross-fitting (LASSO, random forest)
##    9   Smoking history and cigarettes per day
##    10  Trimester of smoking
##    11  Unmeasured confounding: E-value and benchmarking
##    12  Public-health translation and final summary table
###############################################################################

## Work in the folder that holds this script, wherever it is run from.
local({
  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (!nzchar(p)) p <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
  if (!nzchar(p)) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) p <- sub("^--file=", "", a[1]) }
  if (nzchar(p) && file.exists(p)) setwd(dirname(p))
})

## ======= SECTION 0. SETTINGS AND HELPER FUNCTIONS ===========================

QUICK         <- FALSE        # TRUE = short demo (small B, one cross-fit repeat); FALSE = final numbers
RUN_TRIMESTER <- TRUE         # Section 10 (trimester of smoking); set FALSE to skip it
SEED          <- 1234         # project-wide random seed
CORES         <- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)   # every core except one, so the computer stays usable
B             <- if (QUICK) 40L else 1000L               # number of bootstrap resamples (Section 7)
K             <- 5L                                      # number of cross-fitting folds (Section 8)
REPS          <- if (QUICK) 1L else 2L                   # repeats of the random fold split (Section 8)

## Packages. Missing ones are installed once.
need <- c("data.table", "ggplot2", "glmnet", "ranger", "Matrix")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(glmnet); library(ranger); library(Matrix); library(parallel)
})
set.seed(SEED)
options(width = 160, datatable.print.nrows = 60)

dir.create("data",            showWarnings = FALSE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
save_tab <- function(x, name) fwrite(x, file.path("results/tables", name))   # one place to write tables
rd2 <- function(x, k = 2) x[, lapply(.SD, function(v) if (is.numeric(v)) round(v, k) else v)]   # round numeric columns for printing
cat("Cores used:", CORES, "| QUICK:", QUICK, "| bootstrap B:", B, "| cross-fit repeats:", REPS, "\n")

NAVY <- "#1F3864"; ORANGE <- "#C1611A"; GRAY <- "#6B7F99"             # colours used in every figure
theme_proj <- theme_minimal(base_size = 15) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

## The ten confounders of the causal diagram (names used throughout).
CONF     <- c("age_grp", "race", "educ", "precare", "wic", "parity", "ppterm", "pdiab", "phype", "bmi_grp")
RHS_MAIN <- paste(CONF, collapse = " + ")                            # model formula for the main specification

## ---------------------------------------------------------------------------
## The core idea used by almost every section.
## The cohort has 3.26 million births, but the ten confounders take only about
## 65,787 distinct combinations ("cells"). For every cell four counts are
## enough: smokers (n1), preterm births among smokers (y1), non-smokers (n0),
## preterm births among non-smokers (y0). A logistic regression depends on the
## data only through these counts, so fitting it on the cells gives exactly the
## same answer as fitting it on all 3.26 million births, in seconds.
## ---------------------------------------------------------------------------
make_cells <- function(d, vars = CONF) {
  d[, .(n1 = sum(smoker == 1L), n0 = sum(smoker == 0L),
        y1 = sum(preterm[smoker == 1L]), y0 = sum(preterm[smoker == 0L])), by = vars]
}

## Weighted quantile: value of x below which a share p of all births lies (each cell counts n times).
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); x[o][which(cw >= p)[1]] }

## Clamp the propensity score at its lo and hi percentiles (keeps everyone, limits extreme weights).
clamp_e <- function(e, n, lo = 0.01, hi = 0.99) {
  if (is.null(lo)) return(pmin(pmax(e, 1e-6), 1 - 1e-6))
  pmin(pmax(e, wq(e, n, lo), 1e-6), pmin(wq(e, n, hi), 1 - 1e-6))
}

## The three working models of AIPTW, fitted on cell counts:
##   e(X)   probability of smoking given the confounders (propensity score)
##   mu1(X) probability of preterm birth given X among smokers
##   mu0(X) probability of preterm birth given X among non-smokers
## `rhs` is the right-hand side of the model formula. Predictions are returned for every cell of `newdata`.
fit_glm_cells <- function(train, rhs = RHS_MAIN, newdata = train) {
  ps <- glm(as.formula(paste("cbind(n1, n0) ~", rhs)),       data = train,         family = binomial())
  o1 <- glm(as.formula(paste("cbind(y1, n1 - y1) ~", rhs)),  data = train[n1 > 0], family = binomial())
  o0 <- glm(as.formula(paste("cbind(y0, n0 - y0) ~", rhs)),  data = train[n0 > 0], family = binomial())
  list(e = predict(ps, newdata, type = "response"), mu1 = predict(o1, newdata, type = "response"),
       mu0 = predict(o0, newdata, type = "response"), ps_fit = ps, o1_fit = o1, o0_fit = o0)
}

## Every estimator is an average of one value per birth ("psi"), and its standard error is
## sd(psi - average) / sqrt(N) (the influence-function or sandwich standard error).
## Inside a cell there are only four kinds of birth, so psi is computed for the four kinds
## (A smoker preterm, B smoker term, C non-smoker preterm, D non-smoker term) and weighted by counts.
## p1 and p0 are lists with elements A, B, C, D: psi for "everyone smokes" and "no one smokes".
summ_kinds <- function(cells, p1, p0) {
  cnt <- list(A = cells$y1, B = cells$n1 - cells$y1, C = cells$y0, D = cells$n0 - cells$y0)
  N   <- sum(unlist(cnt))
  avg <- function(p) sum(cnt$A * p$A + cnt$B * p$B + cnt$C * p$C + cnt$D * p$D) / N
  m1 <- avg(p1); m0 <- avg(p0)
  ## influence values for the risk difference and the log risk ratio
  sq <- function(f) sum(cnt$A * f$A^2 + cnt$B * f$B^2 + cnt$C * f$C^2 + cnt$D * f$D^2) / N
  fd <- Map(function(a, b) (a - m1) - (b - m0), p1, p0)
  fl <- Map(function(a, b) (a - m1) / m1 - (b - m0) / m0, p1, p0)
  rd <- m1 - m0
  list(rd = 100 * rd, se_rd = 100 * sqrt(sq(fd) / N), rr = m1 / m0, se_lrr = sqrt(sq(fl) / N),
       risk1 = 100 * m1, risk0 = 100 * m0, N = N)
}

## AIPTW (augmented inverse-probability weighting, doubly robust):
##   psi1 = mu1 + T (Y - mu1) / e          psi0 = mu0 + (1 - T)(Y - mu0) / (1 - e)
aiptw_cells <- function(cells, e, mu1, mu0, clamp = c(0.01, 0.99)) {
  n <- cells$n1 + cells$n0
  e <- if (is.null(clamp)) clamp_e(e, n, NULL) else clamp_e(e, n, clamp[1], clamp[2])
  r <- summ_kinds(cells,
        p1 = list(A = mu1 + (1 - mu1) / e, B = mu1 - mu1 / e, C = mu1, D = mu1),
        p0 = list(A = mu0, B = mu0, C = mu0 + (1 - mu0) / (1 - e), D = mu0 - mu0 / (1 - e)))
  r$e_used <- e; r
}

## Weighted-mean estimators (Hajek form): w1 weights the smokers, w0 the non-smokers.
##   inverse-probability weights: w1 = 1/e, w0 = 1/(1-e)     overlap weights: w1 = 1-e, w0 = e
hajek_cells <- function(cells, w1, w0) {
  m1 <- sum(w1 * cells$y1) / sum(w1 * cells$n1); m0 <- sum(w0 * cells$y0) / sum(w0 * cells$n0)
  N  <- sum(cells$n1 + cells$n0)
  mw1 <- sum(w1 * cells$n1) / N; mw0 <- sum(w0 * cells$n0) / N
  summ_kinds(cells,
        p1 = list(A = m1 + w1 * (1 - m1) / mw1, B = m1 - w1 * m1 / mw1, C = m1, D = m1),
        p0 = list(A = m0, B = m0, C = m0 + w0 * (1 - m0) / mw0, D = m0 - w0 * m0 / mw0))
}

## ATT: effect among the mothers who actually smoked. Smokers' own outcomes are observed, so only
## their "no smoking" risk is built from non-smokers, with the weight e/(1-e).
att_cells <- function(cells, e, mu0, clamp = c(0.01, 0.99)) {
  n <- cells$n1 + cells$n0; N <- sum(n); N1 <- sum(cells$n1)
  e <- if (is.null(clamp)) clamp_e(e, n, NULL) else clamp_e(e, n, clamp[1], clamp[2])
  w <- e / (1 - e)
  aA <- 1 - mu0; aB <- -mu0; aC <- -w * (1 - mu0); aD <- w * mu0        # value of one birth of each kind
  cA <- cells$y1; cB <- cells$n1 - cells$y1; cC <- cells$y0; cD <- cells$n0 - cells$y0
  tau <- sum(cA * aA + cB * aB + cC * aC + cD * aD) / N1
  p   <- N1 / N
  v   <- sum(cA * (aA - tau)^2 + cB * (aB - tau)^2 + cC * aC^2 + cD * aD^2) / N / p^2
  list(att = 100 * tau, se = 100 * sqrt(v / N))
}

fmt_ci <- function(est, se) sprintf("%.2f (%.2f, %.2f)", est, est - 1.96 * se, est + 1.96 * se)
evalue_rr <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }   # VanderWeele and Ding (2017)


## ======= SECTION 1. DATA ====================================================
## Source: National Center for Health Statistics, National Vital Statistics System,
## 2023 Natality Public Use File (public, no personal identifiers).
##   Official file : https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip
##   User guide    : https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/UserGuide2023.pdf
## The full file has 3.6 million births and more than 200 columns. The analysis uses a
## 24-column extract (3,476,180 singleton births with known gestational age and known smoking
## status) stored in the project repository.
##
## WHERE THE DATA ARE READ FROM (first one found is used):
##   1. data/nvss2023_extract.csv        (already unzipped)
##   2. data/nvss2023_extract.csv.gz     (compressed copy)
##   3. EXTRACT_URL below                (the copy stored in the project repository on GitHub; downloaded once
##                                        into data/, so the script runs on any computer with no other setup).

EXTRACT_CSV <- "data/nvss2023_extract.csv"
EXTRACT_GZ  <- paste0(EXTRACT_CSV, ".gz")
EXTRACT_URL <- "https://github.com/fmurital/Causal-Inference-/raw/main/data/nvss2023_extract.csv.gz"   # file kept in the repository folder data/

if (!file.exists(EXTRACT_CSV)) {
  if (!file.exists(EXTRACT_GZ)) {
    cat("Downloading the extract (36 MB) ...\n")
    download.file(EXTRACT_URL, EXTRACT_GZ, mode = "wb")
  }
  cat("Decompressing ...\n")
  inp <- gzfile(EXTRACT_GZ, "rb"); out <- file(EXTRACT_CSV, "wb")       # base R, works on Windows, Mac and Linux
  repeat { chunk <- readBin(inp, "raw", 5e7); if (length(chunk) == 0) break; writeBin(chunk, out) }
  close(inp); close(out)
}

## Columns holding letters (Y / N / U) are read as text so that nothing is silently converted.
CHAR_COLS <- c("dmar", "wic", "rf_ppterm", "rf_pdiab", "rf_phype", "rf_gdiab", "rf_ghype", "cig_rec", "sex")
raw <- fread(EXTRACT_CSV, colClasses = list(character = CHAR_COLS))
flow <- data.table(step = "Singleton live births, known gestational age and smoking status (extract)", n = nrow(raw))

## Residency rule: births to mothers who are not U.S. residents (restatus = 4) are excluded, as in NCHS reports.
raw  <- raw[restatus != 4]
flow <- rbind(flow, data.table(step = "After excluding non-U.S.-resident mothers", n = nrow(raw)))

## build_cohort(): recode the raw columns into analysis variables. "Unknown" codes become NA.
## It is a function so the same recoding can be applied to the full natality file in Section 10.
yn <- function(x) factor(ifelse(x == "Y", "Yes", ifelse(x == "N", "No", NA)), levels = c("No", "Yes"))
build_cohort <- function(raw) {
  d <- raw[, .(
    smoker   = as.integer(cig_rec == "Y"),                 # exposure T: smoked at any time during pregnancy
    preterm  = as.integer(oegest_comb < 37),               # outcome Y: born before 37 completed weeks
    ga_weeks = oegest_comb,                                # gestational age in weeks (obstetric estimate)
    mager    = mager,                                      # maternal age in whole years (a number)
    age_grp  = cut(mager, c(-Inf, 19, 24, 29, 34, 39, Inf),            # the same age, in six groups
                   labels = c("<20", "20-24", "25-29", "30-34", "35-39", "40+")),
    race     = factor(mrace6, levels = 1:6, labels = c("White", "Black", "AIAN", "Asian", "NHOPI", "Multiracial")),
    educ     = factor(ifelse(meduc == 9, NA, meduc), levels = 1:8,
                      labels = c("<9th grade", "9-12th, no diploma", "HS grad/GED", "Some college",
                                 "Associate", "Bachelor's", "Master's", "Doctorate/Prof.")),
    precare  = factor(ifelse(precare5 == 5, NA, precare5), levels = 1:4,
                      labels = c("1st trimester", "2nd trimester", "3rd trimester", "No care")),
    wic      = yn(ifelse(wic == "U", NA, wic)),
    parity   = {pl <- ifelse(priorlive == 99 | priordead == 99, NA, priorlive + priordead)   # previous live births
                factor(pmin(pl, 3), levels = 0:3, labels = c("0", "1", "2", "3+"))},
    ppterm   = yn(ifelse(rf_ppterm == "U", NA, rf_ppterm)),            # previous preterm birth
    pdiab    = yn(ifelse(rf_pdiab  == "U", NA, rf_pdiab)),             # pre-pregnancy diabetes
    phype    = yn(ifelse(rf_phype  == "U", NA, rf_phype)),             # pre-pregnancy hypertension
    bmi      = ifelse(bmi >= 99.9, NA, bmi),                           # pre-pregnancy BMI (a number)
    bmi_grp  = cut(ifelse(bmi >= 99.9, NA, bmi), c(-Inf, 18.5, 25, 30, 35, 40, Inf), right = FALSE,
                   labels = c("<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+")),
    pay      = factor(pay_rec, levels = c(1, 2, 3, 4, 9), labels = c("Medicaid", "Private", "Self-pay", "Other", "Unknown")),
    cig_0    = cig_0                                                   # cigarettes per day before pregnancy
  )]
  for (v in intersect(c("cig_1", "cig_2", "cig_3"), names(raw))) d[[v]] <- raw[[v]]   # trimester columns, full file only
  d
}
d_all <- build_cohort(raw)                         # every resident birth, before dropping missing values
ok    <- complete.cases(d_all[, ..CONF])           # TRUE when all ten confounders are known
d     <- d_all[ok]                                 # analytic cohort: complete data on all ten confounders
flow  <- rbind(flow, data.table(step = "Complete data on all ten confounders (analytic cohort)", n = nrow(d)))
flow[, pct_of_extract := round(100 * n / flow$n[1], 2)]
print(flow); save_tab(flow, "cohort_flow.csv")
stopifnot(nrow(d) == 3259496L)                     # the cohort size reported in the paper
cat(sprintf("\nAnalytic cohort: %s births, %s smokers (%.2f%%), %s preterm (%.2f%%)\n",
            format(nrow(d), big.mark = ","), format(sum(d$smoker), big.mark = ","), 100 * mean(d$smoker),
            format(sum(d$preterm), big.mark = ","), 100 * mean(d$preterm)))
cat(sprintf("Crude preterm risk: smokers %.2f%%, non-smokers %.2f%%, difference %.2f points\n",
            100 * mean(d$preterm[d$smoker == 1]), 100 * mean(d$preterm[d$smoker == 0]),
            100 * (mean(d$preterm[d$smoker == 1]) - mean(d$preterm[d$smoker == 0]))))

cells <- make_cells(d)                             # the covariate cells used from here on
N     <- sum(cells$n1 + cells$n0)
cat("Distinct covariate cells:", format(nrow(cells), big.mark = ","), "\n")


## ======= SECTION 2. DATA DICTIONARY: WHAT EACH VARIABLE IS ==================
## Question answered here: is age a number or a category, and what are the levels of every input?

## 2.1 Raw columns in the extract, as they appear in the file
dict_raw <- data.table(
  column = c("mager", "mrace6", "meduc", "precare5", "wic", "priorlive", "priordead", "rf_ppterm", "rf_pdiab",
             "rf_phype", "bmi", "cig_0", "cig_rec", "oegest_comb", "pay_rec", "restatus"),
  stored_as = c("number (years)", "code 1-6", "code 1-9", "code 1-5", "letter Y/N/U", "number", "number",
                "letter Y/N/U", "letter Y/N/U", "letter Y/N/U", "number (kg/m2)", "number (cigarettes/day)",
                "letter Y/N", "number (weeks)", "code 1-4, 9", "code 1-4"),
  meaning = c("Mother's age in whole years", "Mother's race (6 groups)", "Mother's education (8 levels, 9 = unknown)",
              "Prenatal care began: 1 = first trimester ... 4 = none, 5 = unknown", "Received WIC food assistance",
              "Previous births still living", "Previous births now deceased", "Previous preterm birth",
              "Diabetes before pregnancy", "Hypertension before pregnancy", "Pre-pregnancy body mass index (99.9 = unknown)",
              "Cigarettes per day in the 3 months before pregnancy (99 = unknown)", "Smoked at any time during pregnancy (the exposure)",
              "Gestational age at birth in weeks (the outcome is < 37)", "Payer: 1 Medicaid 2 Private 3 Self-pay 4 Other",
              "Residence: 4 = not a U.S. resident (excluded)"))
cat("\n2.1 RAW COLUMNS\n"); print(dict_raw); save_tab(dict_raw, "dictionary_raw_columns.csv")

## 2.2 Analysis variables after recoding: type, number of levels, the levels themselves, and missing share
describe_var <- function(v) {
  x <- d_all[[v]]
  data.table(variable = v, type = class(x)[1],
             n_levels = if (is.factor(x)) nlevels(x) else NA_integer_,
             levels_or_range = if (is.factor(x)) paste(levels(x), collapse = " | ") else paste(round(range(x, na.rm = TRUE), 1), collapse = " to "),
             pct_missing = round(100 * mean(is.na(x)), 2))
}
dict <- rbindlist(lapply(c("smoker", "preterm", "ga_weeks", "mager", "age_grp", "race", "educ", "precare", "wic", "parity",
                           "ppterm", "pdiab", "phype", "bmi", "bmi_grp", "pay", "cig_0"), describe_var))
cat("\n2.2 ANALYSIS VARIABLES (before dropping missing values)\n"); print(dict); save_tab(dict, "dictionary_analysis_variables.csv")

## 2.3 Age: reported as a whole number of years (mager) and also grouped (age_grp)
cat("\n2.3 AGE: mager is a number; age_grp is the same value cut into six groups\n")
print(summary(d$mager)); print(head(d[, .(mager, age_grp)], 8))
print(d[, .(n = .N, min_age = min(mager), max_age = max(mager)), by = age_grp][order(min_age)])

## 2.4 Number of births in every level of every confounder, by smoking status
lev_counts <- rbindlist(lapply(CONF, function(v)
  d[, .(n = .N, n_smokers = sum(smoker), pct_smokers = round(100 * mean(smoker), 2),
        pct_preterm = round(100 * mean(preterm), 2)), by = .(level = get(v))][order(level)][, `:=`(level = as.character(level), variable = v)][]))
setcolorder(lev_counts, "variable")
cat("\n2.4 LEVELS OF EVERY CONFOUNDER (all rows saved to results/tables/levels_by_confounder.csv)\n")
print(lev_counts[1:20]); save_tab(lev_counts, "levels_by_confounder.csv")


## ======= SECTION 3. EXPLORATORY ANALYSIS (EDA) ==============================

## 3.1 Headline numbers
d[, group := factor(ifelse(smoker == 1, "Smoked in pregnancy", "Did not smoke"), levels = c("Did not smoke", "Smoked in pregnancy"))]
eda_numbers <- data.table(
  quantity = c("Births in cohort", "Smokers", "Smoking prevalence (%)", "Preterm, all (%)", "Preterm, smokers (%)",
               "Preterm, non-smokers (%)", "Mean gestational age, smokers (weeks)", "Mean gestational age, non-smokers (weeks)",
               "Mean maternal age, smokers", "Mean maternal age, non-smokers"),
  value = c(nrow(d), sum(d$smoker), 100 * mean(d$smoker), 100 * mean(d$preterm), 100 * mean(d$preterm[d$smoker == 1]),
            100 * mean(d$preterm[d$smoker == 0]), mean(d$ga_weeks[d$smoker == 1]), mean(d$ga_weeks[d$smoker == 0]),
            mean(d$mager[d$smoker == 1]), mean(d$mager[d$smoker == 0])))
cat("\n3.1 HEADLINE NUMBERS\n"); print(rd2(eda_numbers)); save_tab(eda_numbers, "eda_numbers.csv")

## 3.2 Table 1: each confounder level, share of smokers and of non-smokers, and preterm risk within each group
table1 <- rbindlist(lapply(CONF, function(v) {
  tab <- d[, .N, by = .(level = as.character(get(v)), smoker)]
  tab <- dcast(tab, level ~ smoker, value.var = "N", fill = 0)
  setnames(tab, c("0", "1"), c("n_nonsmoker", "n_smoker"))
  tab[, `:=`(variable = v, pct_of_nonsmokers = 100 * n_nonsmoker / sum(n_nonsmoker), pct_of_smokers = 100 * n_smoker / sum(n_smoker))]
  pt <- d[, .(r = 100 * mean(preterm)), by = .(level = as.character(get(v)), smoker)]
  tab[, preterm_pct_nonsmoker := pt[smoker == 0][match(tab$level, level), r]]
  tab[, preterm_pct_smoker    := pt[smoker == 1][match(tab$level, level), r]]
  tab
}))
setcolorder(table1, "variable")
cat("\n3.2 TABLE 1 (first rows)\n"); print(rd2(table1, 1)[1:12]); save_tab(table1, "table1_by_smoking.csv")

## 3.3 Is the complete-case cohort representative? Compare births kept with births dropped for missing data.
ok_rows <- complete.cases(d_all[, ..CONF]); d_all[, kept := ok_rows]
cc <- d_all[, .(births = .N, pct_smokers = round(100 * mean(smoker), 2), pct_preterm = round(100 * mean(preterm), 2)), by = .(kept)]
cat("\n3.3 BIRTHS KEPT (TRUE) VERSUS DROPPED FOR MISSING DATA (FALSE)\n"); print(cc); save_tab(cc, "complete_case_comparison.csv")
miss <- data.table(variable = CONF, n_missing = sapply(CONF, function(v) sum(is.na(d_all[[v]]))))
miss[, pct_missing := round(100 * n_missing / nrow(d_all), 2)]
cat("Missing values by confounder:\n"); print(miss); save_tab(miss, "missingness.csv")

## 3.4 Figures. EDA 2: gestational age by smoking. EDA 1: who smokes by education, age and race.
ga <- d[ga_weeks >= 24, .N, by = .(group, ga_weeks)]; ga[, share := 100 * N / sum(N), by = group]
p_ga <- ggplot(ga, aes(ga_weeks, share, fill = group)) + geom_col(position = "dodge", width = 0.85) +
  geom_vline(xintercept = 36.5, linetype = "dashed") +
  annotate("text", x = 36.3, y = max(ga$share) * 0.95, label = "Preterm\n(< 37 weeks)", hjust = 1, size = 4) +
  scale_fill_manual(values = c(GRAY, ORANGE), name = NULL) +
  labs(x = "Gestational age at delivery (completed weeks)", y = "% of births within group") + theme_proj
save_tab(ga, "eda_gestational_age.csv")
ggsave("results/figures/fig_eda_gestational_age.png", p_ga, width = 9, height = 4.8, dpi = 200)

prev <- function(v, lab) d[, .(prev = 100 * mean(smoker)), by = .(level = as.character(get(v)))][, var := lab][]
pv <- rbind(prev("educ", "Education"), prev("age_grp", "Maternal age"), prev("race", "Race/ethnicity"))
pv[, level := factor(level, levels = unlist(lapply(c("educ", "age_grp", "race"), function(v) levels(d[[v]]))))]
p_who <- ggplot(pv, aes(level, prev)) + geom_col(fill = NAVY, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", prev)), vjust = -0.4, size = 4) + facet_wrap(~var, scales = "free_x", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) + labs(x = NULL, y = "% who smoked in pregnancy") +
  theme_proj + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("results/figures/fig_eda_who_smokes.png", p_who, width = 13, height = 5.2, dpi = 200)

## 3.5 Association matrix (the "correlation matrix" for categorical data).
## Cramer's V runs from 0 (no association) to 1 (perfect association): V = sqrt(chi2 / (N (min(rows, cols) - 1))).
vars12 <- c("smoker", "preterm", "educ", "race", "age_grp", "precare", "wic", "parity", "ppterm", "pdiab", "phype", "bmi_grp")
cramers_v <- function(a, b) {
  tab  <- table(a, b)
  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  sqrt(as.numeric(chi2) / (sum(tab) * (min(dim(tab)) - 1)))
}
fac <- lapply(vars12, function(v) factor(d[[v]])); names(fac) <- vars12
M <- matrix(NA_real_, 12, 12, dimnames = list(vars12, vars12))
for (i in 1:12) for (j in 1:12) M[i, j] <- if (i == j) 1 else if (j < i) M[j, i] else cramers_v(fac[[i]], fac[[j]])
cat("\n3.5 ASSOCIATION MATRIX (Cramer's V)\n"); print(round(M, 2))
save_tab(as.data.table(M, keep.rownames = "var"), "cramers_v.csv")
rm(fac); invisible(gc())


## ======= SECTION 4. WHICH CONFOUNDERS MATTER MOST ===========================
## Method: drop one confounder at a time and refit. Two readings are reported:
##   (a) share of the model's explanatory power (deviance) lost, for the smoking model and for the
##       preterm model among non-smokers: this says how strongly the variable sorts smokers
##       or predicts preterm birth;
##   (b) the change in the AIPTW risk difference when the variable is left out: this says how much
##       the variable matters for the answer.
dev_expl <- function(m) m$null.deviance - m$deviance
fit_models <- function(vs) {
  f <- if (length(vs)) paste(vs, collapse = " + ") else "1"
  list(ps = glm(as.formula(paste("cbind(n1, n0) ~", f)), data = cells, family = binomial()),
       o0 = glm(as.formula(paste("cbind(y0, n0 - y0) ~", f)), data = cells[n0 > 0], family = binomial()),
       o1 = glm(as.formula(paste("cbind(y1, n1 - y1) ~", f)), data = cells[n1 > 0], family = binomial()))
}
rd_of <- function(m) {
  e <- predict(m$ps, cells, type = "response"); m1 <- predict(m$o1, cells, type = "response"); m0 <- predict(m$o0, cells, type = "response")
  aiptw_cells(cells, e, m1, m0)$rd
}
full <- fit_models(CONF); rd_full <- rd_of(full)
contrib <- rbindlist(lapply(CONF, function(v) {
  m <- fit_models(setdiff(CONF, v))
  data.table(variable = v,
             share_of_smoking_model_pct   = 100 * (dev_expl(full$ps) - dev_expl(m$ps)) / dev_expl(full$ps),
             share_of_preterm_model_pct   = 100 * (dev_expl(full$o0) - dev_expl(m$o0)) / dev_expl(full$o0),
             rd_without_variable          = rd_of(m),
             change_in_rd                 = rd_of(m) - rd_full)
}))
cat("\nAIPTW risk difference with all ten confounders:", round(rd_full, 3), "(reported 3.26)\n")
cat("\n4. CONFOUNDER RANKING\n"); print(rd2(contrib[order(-share_of_smoking_model_pct)], 2))
save_tab(contrib, "confounder_contributions.csv")


## ======= SECTION 5. MAIN ESTIMATE AND COMPARISON ESTIMATORS =================
## Estimand: average treatment effect (ATE) on the risk-difference scale = preterm risk if every
## mother smoked minus preterm risk if none did. The propensity score is clamped at its 1st and 99th
## percentiles (everyone stays in the analysis; only extreme weights are limited).

## 5.1 Fit the three working models and compute every estimator
fit <- fit_glm_cells(cells)                                  # e(X), mu1(X), mu0(X) for every cell
cells[, `:=`(e = fit$e, mu1 = fit$mu1, mu0 = fit$mu0)]       # keep the predictions with the cells
n_c <- cells$n1 + cells$n0
e_cl <- clamp_e(cells$e, n_c, 0.01, 0.99)                    # clamped propensity score

aip  <- aiptw_cells(cells, cells$e, cells$mu1, cells$mu0)
iptw <- hajek_cells(cells, 1 / e_cl, 1 / (1 - e_cl))
ato  <- hajek_cells(cells, 1 - cells$e, cells$e)           # overlap weights (need no clamping)
att  <- att_cells(cells, cells$e, cells$mu0)
crude_rd <- 100 * (sum(cells$y1) / sum(cells$n1) - sum(cells$y0) / sum(cells$n0))
crude_se <- 100 * sqrt(mean(d$preterm[d$smoker == 1]) * (1 - mean(d$preterm[d$smoker == 1])) / sum(cells$n1) +
                       mean(d$preterm[d$smoker == 0]) * (1 - mean(d$preterm[d$smoker == 0])) / sum(cells$n0))
gcomp_rd <- 100 * sum(n_c * (cells$mu1 - cells$mu0)) / N      # outcome models only (no weighting)

main_tab <- data.table(
  estimator = c("Crude (unadjusted)", "Outcome regression (g-computation)", "IPTW (stabilized)", "AIPTW (doubly robust)",
                "AIPTW, effect among smokers (ATT)", "Overlap weights (ATO)"),
  RD = c(crude_rd, gcomp_rd, iptw$rd, aip$rd, att$att, ato$rd),
  SE = c(crude_se, NA, iptw$se_rd, aip$se_rd, att$se, ato$se_rd))
main_tab[, `:=`(lo = RD - 1.96 * SE, hi = RD + 1.96 * SE)]
main_tab[, RR := c(sum(cells$y1) / sum(cells$n1) / (sum(cells$y0) / sum(cells$n0)), NA, iptw$rr, aip$rr, NA, ato$rr)]
cat("\n5.1 RISK DIFFERENCE IN PREGNANCY SMOKING, PERCENTAGE POINTS (95% CI)\n"); print(rd2(main_tab, 3)); save_tab(main_tab, "results_main.csv")
cat(sprintf("\nAIPTW risk difference %s | risk ratio %.2f (%.2f, %.2f)\n", fmt_ci(aip$rd, aip$se_rd), aip$rr,
            exp(log(aip$rr) - 1.96 * aip$se_lrr), exp(log(aip$rr) + 1.96 * aip$se_lrr)))
cat(sprintf("Risk if all smoked %.2f%% | risk if none smoked %.2f%%\n", aip$risk1, aip$risk0))

## 5.2 Which half of AIPTW carries the number? Outcome models alone and weighting alone should agree.
cat(sprintf("\n5.2 Outcome-model part (g-computation): %.2f | weighting part (IPTW): %.2f | combined AIPTW: %.2f\n",
            gcomp_rd, iptw$rd, aip$rd))

## 5.3 Odds ratio on the marginal scale (the scale most earlier studies report)
m1 <- aip$risk1 / 100; m0 <- aip$risk0 / 100
cat(sprintf("Marginal odds ratio from AIPTW risks: %.2f\n", (m1 / (1 - m1)) / (m0 / (1 - m0))))

## 5.4 CONCERN 1. THE ADJUSTMENT SET.
## Prenatal-care timing (precare) and WIC are recorded during pregnancy. If smoking changes them, they would be
## mediators or colliders and adjusting for them could bias the estimate. Two checks are run:
##   (i)  refit with pre-pregnancy variables only (age, race, education, parity, previous preterm birth, diabetes,
##        hypertension, BMI): precare and WIC removed;
##   (ii) add insurance type, and add two-way interactions, to see whether the answer depends on the specification.
run_spec <- function(rhs, vars, label, clamp = c(0.01, 0.99)) {
  cl <- make_cells(d, vars); f <- suppressWarnings(fit_glm_cells(cl, rhs))   # interaction terms for empty combinations trigger a harmless rank warning
  r  <- aiptw_cells(cl, f$e, f$mu1, f$mu0, clamp); a <- att_cells(cl, f$e, f$mu0, clamp)
  data.table(analysis = label, RD = r$rd, lo = r$rd - 1.96 * r$se_rd, hi = r$rd + 1.96 * r$se_rd, ATT = a$att)
}
PRE <- setdiff(CONF, c("precare", "wic"))
spec_tab <- rbind(
  run_spec(RHS_MAIN, CONF, "Main: ten confounders"),
  run_spec(paste(PRE, collapse = " + "), PRE, "Pre-pregnancy set only (no care timing, no WIC)"),
  run_spec(paste(RHS_MAIN, "+ pay"), c(CONF, "pay"), "Main plus insurance type"),
  run_spec(paste(RHS_MAIN, "+ race:educ + age_grp:educ + age_grp:parity + educ:wic + educ:precare"), CONF, "Main plus two-way interactions"))
cat("\n5.4 ADJUSTMENT-SET CHECKS\n"); print(rd2(spec_tab, 2))
cat(sprintf("Dropping care timing and WIC moves the risk difference by %+.2f points (%.2f to %.2f); the confidence intervals overlap.\n",
            spec_tab$RD[2] - spec_tab$RD[1], spec_tab$RD[1], spec_tab$RD[2]))

## 5.5 Propensity-score handling: no clamp, 0.1/99.9, 1/99 (main), 5/95, and deleting births outside [0.01, 0.99]
p_smk <- sum(cells$n1) / N                                    # overall share of smokers (stabilizing factor)
clamp_tab <- rbindlist(lapply(list(c(NA, NA), c(0.001, 0.999), c(0.01, 0.99), c(0.05, 0.95)), function(cl) {
  cc <- if (is.na(cl[1])) NULL else cl
  r  <- aiptw_cells(cells, cells$e, cells$mu1, cells$mu0, cc)
  ee <- r$e_used; w <- p_smk / ee[cells$n1 > 0]                   # stabilized weight of each smoker cell
  data.table(analysis = if (is.null(cc)) "No clamp" else sprintf("Clamp at %g / %g percentile", 100 * cl[1], 100 * cl[2]),
             RD = r$rd, lo = r$rd - 1.96 * r$se_rd, hi = r$rd + 1.96 * r$se_rd,
             max_weight_stabilized = max(w), effective_n_smokers = sum(cells$n1[cells$n1 > 0] * w)^2 / sum(cells$n1[cells$n1 > 0] * w^2))
}))
keep <- cells$e >= 0.01 & cells$e <= 0.99
r_trim <- aiptw_cells(cells[keep], cells$e[keep], cells$mu1[keep], cells$mu0[keep], clamp = NULL)
clamp_tab <- rbind(clamp_tab, data.table(analysis = sprintf("Delete births with score outside [0.01, 0.99] (keeps %.1f%%)", 100 * sum(n_c[keep]) / N),
                   RD = r_trim$rd, lo = r_trim$rd - 1.96 * r_trim$se_rd, hi = r_trim$rd + 1.96 * r_trim$se_rd,
                   max_weight_stabilized = NA_real_, effective_n_smokers = NA_real_))
cat("\n5.5 PROPENSITY-SCORE HANDLING\n"); print(rd2(clamp_tab, 2)); save_tab(clamp_tab, "sens_ps_handling.csv")

## 5.6 Everything in one table for the robustness figure
robust <- rbind(
  data.table(analysis = "Crude, unadjusted", RD = crude_rd, lo = crude_rd - 1.96 * crude_se, hi = crude_rd + 1.96 * crude_se),
  spec_tab[, .(analysis, RD, lo, hi)],
  clamp_tab[analysis == "Clamp at 5 / 95 percentile", .(analysis, RD, lo, hi)],
  data.table(analysis = "Overlap weights", RD = ato$rd, lo = ato$rd - 1.96 * ato$se_rd, hi = ato$rd + 1.96 * ato$se_rd),
  data.table(analysis = "Among smokers (ATT)", RD = att$att, lo = att$att - 1.96 * att$se, hi = att$att + 1.96 * att$se))
cat("\n5.6 ROBUSTNESS SUMMARY\n"); print(rd2(robust, 2)); save_tab(robust, "robustness_summary.csv")
cat(sprintf("Range of adjusted estimates across the analytic choices: %.2f to %.2f points\n",
            min(robust$RD[-1]), max(robust$RD[-1])))
p_rob <- ggplot(robust[, analysis := factor(analysis, levels = rev(analysis))], aes(RD, analysis)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = NAVY) + geom_point(size = 3, colour = ORANGE) +
  labs(x = "Risk difference in preterm birth (percentage points, 95% CI)", y = NULL) + theme_proj
ggsave("results/figures/fig_robustness.png", p_rob, width = 9, height = 4.8, dpi = 200)


## ======= SECTION 6. OVERLAP, WEIGHTS AND COVARIATE BALANCE ======
## Positivity needs every kind of mother to have some chance of smoking and of not smoking. It can
## fail in practice: some non-smokers have almost no smoker who resembles them. This section measures how
## often, who those mothers are, what the weights look like, and whether weighting balances the covariates.

## 6.1 Propensity-score distribution by group (unclamped score, counted per birth)
qs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
ps_q <- data.table(quantile = qs,
  smokers     = sapply(qs, function(p) wq(cells$e[cells$n1 > 0], cells$n1[cells$n1 > 0], p)),
  nonsmokers  = sapply(qs, function(p) wq(cells$e[cells$n0 > 0], cells$n0[cells$n0 > 0], p)))
cat("\n6.1 PROPENSITY SCORE QUANTILES\n"); print(rd2(ps_q, 4))
share_low_ns <- 100 * sum(cells$n0[cells$e < 0.01]) / sum(cells$n0)
share_low_s  <- 100 * sum(cells$n1[cells$e < 0.01]) / sum(cells$n1)
cat(sprintf("Non-smokers with score below 0.01: %.1f%% | smokers with score below 0.01: %.1f%%\n", share_low_ns, share_low_s))
cat("Reading: about", round(share_low_ns), "% of non-smokers have essentially no smoker who resembles them, which limits what can be said about them.\n")

## 6.2 Who are the non-smokers with no look-alike?
cells[, low := e < 0.01]
who <- rbindlist(lapply(c("educ", "age_grp", "race"), function(v) {
  tt <- cells[n0 > 0, .(low_n = sum(n0[low]), other_n = sum(n0[!low])), by = .(level = as.character(get(v)))]
  tt[, `:=`(variable = v, pct_of_low_group = 100 * low_n / sum(low_n), pct_of_rest = 100 * other_n / sum(other_n))]
  tt[, .(variable, level, pct_of_low_group, pct_of_rest)]
}))
cat("\n6.2 NON-SMOKERS WITH SCORE < 0.01 VERSUS THE REST (share within each group)\n"); print(rd2(who, 1)); save_tab(who, "who_has_no_lookalike.csv")

## 6.3 Weights after clamping: largest weight and effective sample size, (sum w)^2 / sum(w^2)
p_s   <- sum(cells$n1) / N
w1    <- p_s / e_cl; w0 <- (1 - p_s) / (1 - e_cl)               # stabilized ATE weights for smokers and non-smokers
ess   <- function(w, n) sum(w * n)^2 / sum(w^2 * n)
wt <- data.table(group = c("Smokers", "Non-smokers"), births = c(sum(cells$n1), sum(cells$n0)),
                 max_weight = c(max(w1[cells$n1 > 0]), max(w0[cells$n0 > 0])),
                 mean_weight = c(sum(w1 * cells$n1) / sum(cells$n1), sum(w0 * cells$n0) / sum(cells$n0)),
                 effective_n = c(ess(w1[cells$n1 > 0], cells$n1[cells$n1 > 0]), ess(w0[cells$n0 > 0], cells$n0[cells$n0 > 0])))
cat("\n6.3 STABILIZED WEIGHTS\n"); print(rd2(wt, 1)); save_tab(wt, "weight_summary.csv")

## 6.4 Covariate balance: standardized mean difference (SMD) of every covariate level, before and after weighting.
## SMD = (p1 - p0) / sqrt((p1 (1 - p1) + p0 (1 - p0)) / 2); below 0.10 is the usual criterion for good balance.
smd_w <- function(x, wa, wb) {                              # x: 0/1 per cell; wa, wb: weights for smokers, non-smokers
  p1 <- sum(wa * cells$n1 * x) / sum(wa * cells$n1); p0 <- sum(wb * cells$n0 * x) / sum(wb * cells$n0)
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
}
ov1 <- 1 - cells$e; ov0 <- cells$e                            # overlap weights
bal <- rbindlist(lapply(CONF, function(v) rbindlist(lapply(levels(cells[[v]]), function(lv) {
  x <- as.numeric(cells[[v]] == lv)
  data.table(covariate = paste0(v, ": ", lv), before = smd_w(x, 1, 1), after_iptw = smd_w(x, w1, w0), after_overlap = smd_w(x, ov1, ov0))
}))))
bal <- bal[!is.na(before) & !(grepl(": No$", covariate))]    # keep "Yes" levels of yes/no variables only
cat("\n6.4 COVARIATE BALANCE (absolute SMD)\n")
cat(sprintf("Largest before weighting: %.3f (%s)\n", max(abs(bal$before)), bal$covariate[which.max(abs(bal$before))]))
cat(sprintf("Largest after IPTW:      %.3f (%s)\n", max(abs(bal$after_iptw)), bal$covariate[which.max(abs(bal$after_iptw))]))
cat(sprintf("Largest after overlap weights: %.3f\n", max(abs(bal$after_overlap))))
cat(sprintf("Levels with |SMD| > 0.10: before %d of %d, after IPTW %d, after overlap weights %d\n",
            sum(abs(bal$before) > 0.1), nrow(bal), sum(abs(bal$after_iptw) > 0.1), sum(abs(bal$after_overlap) > 0.1)))
save_tab(bal, "balance_table.csv")

save_tab(cells[, .(e, n1, n0, n_births = n1 + n0)], "cell_scores.csv")   # one row per cell: used to draw the overlap figure elsewhere

## 6.5 Figure: propensity score overlap (log scale). Weights in the density are births per cell.
dens <- rbind(data.table(e = cells$e[cells$n1 > 0], w = cells$n1[cells$n1 > 0] / sum(cells$n1), group = "Smoked"),
              data.table(e = cells$e[cells$n0 > 0], w = cells$n0[cells$n0 > 0] / sum(cells$n0), group = "Did not smoke"))
p_ov <- ggplot(dens, aes(e, weight = w, fill = group, colour = group)) + geom_density(alpha = 0.45, bw = 0.15) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1), labels = c("0.001", "0.01", "0.1")) +
  geom_vline(xintercept = 0.01, linetype = "dashed", colour = GRAY) +
  scale_fill_manual(values = c(GRAY, ORANGE), name = NULL) + scale_colour_manual(values = c(GRAY, ORANGE), name = NULL) +
  labs(x = "Estimated probability of smoking (log scale)", y = "Density") + theme_proj
ggsave("results/figures/fig_ps_overlap.png", p_ov, width = 8, height = 4.8, dpi = 200)


## ======= SECTION 7. VARIANCE: INFLUENCE FUNCTION VERSUS BOOTSTRAP
## The standard error in Section 5 comes from the influence function (a sandwich formula) that treats the
## fitted propensity and outcome models as known. Two checks:
##   7.1 the fast cell-level calculation must equal a direct calculation on all 3.26 million births;
##   7.2 a nonparametric bootstrap on the FULL cohort: every resample redraws 3.26 million births, refits all
##       three models and re-clamps the propensity score, so model-fitting uncertainty is included.
##       Resampling people with replacement is the same as drawing the four counts of every cell from a
##       multinomial distribution, which takes a fraction of a second.

## 7.1 Direct person-level calculation, using the same fitted scores
cells[, e_used := aip$e_used]
dd <- cells[, c(CONF, "e_used", "mu1", "mu0"), with = FALSE][d[, c("smoker", "preterm", CONF), with = FALSE], on = CONF]
psi1 <- dd$mu1 + dd$smoker * (dd$preterm - dd$mu1) / dd$e_used
psi0 <- dd$mu0 + (1 - dd$smoker) * (dd$preterm - dd$mu0) / (1 - dd$e_used)
rd_person <- 100 * (mean(psi1) - mean(psi0)); se_person <- 100 * sd(psi1 - psi0) / sqrt(length(psi1))
cat(sprintf("\n7.1 PERSON-LEVEL CHECK  RD %.4f SE %.4f   |   CELL-LEVEL  RD %.4f SE %.4f\n", rd_person, se_person, aip$rd, aip$se_rd))
stopifnot(abs(rd_person - aip$rd) < 1e-6, abs(se_person - aip$se_rd) < 1e-6)
rm(dd, psi1, psi0); invisible(gc())

## 7.2 Full-sample bootstrap, run on several cores
cnt  <- c(cells$y1, cells$n1 - cells$y1, cells$y0, cells$n0 - cells$y0)      # the four counts of every cell
prob <- cnt / sum(cnt)                                                        # share of all births in each count
C    <- nrow(cells)
one_boot <- function(b) {
  k  <- as.vector(rmultinom(1, size = N, prob = prob))                        # one resampled cohort
  bc <- cells[, CONF, with = FALSE]                                           # same cells, new counts
  bc[, `:=`(y1 = k[1:C], n1 = k[1:C] + k[(C + 1):(2 * C)],
            y0 = k[(2 * C + 1):(3 * C)], n0 = k[(2 * C + 1):(3 * C)] + k[(3 * C + 1):(4 * C)])]
  f  <- fit_glm_cells(bc, RHS_MAIN)                                           # refit all three models
  r  <- aiptw_cells(bc, f$e, f$mu1, f$mu0)                                    # re-clamp and re-estimate
  c(rd = r$rd, rr = r$rr, se_if = r$se_rd)
}
safe_boot <- function(b) tryCatch(one_boot(b), error = function(e) c(rd = NA, rr = NA, se_if = NA))

## Run a function over 1..n on `cores` cores (works on Windows, Mac and Linux).
run_parallel <- function(n, fn, cores, export) {
  if (cores <= 1) return(lapply(seq_len(n), fn))
  cl <- makeCluster(cores); on.exit(stopCluster(cl))
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)) })
  clusterExport(cl, export, envir = globalenv())
  clusterSetRNGStream(cl, SEED)
  parLapplyLB(cl, seq_len(n), fn)
}
t0 <- Sys.time()
bs <- do.call(rbind, run_parallel(B, safe_boot, CORES,
        c("cells", "prob", "N", "C", "CONF", "RHS_MAIN", "one_boot", "fit_glm_cells", "aiptw_cells", "summ_kinds", "clamp_e", "wq")))
cat(sprintf("\n7.2 %d bootstrap resamples on %d cores took %.1f minutes; %d failed\n", B, CORES,
            as.numeric(Sys.time() - t0, units = "mins"), sum(is.na(bs[, "rd"]))))
bs <- bs[!is.na(bs[, "rd"]), , drop = FALSE]
if (nrow(bs) == 0) stop("Every bootstrap resample failed; run one_boot(1) on its own to see the error message.")

## 7.3 Comparison
var_tab <- data.table(
  quantity = c("AIPTW risk difference (points)", "Standard error, influence function", "Standard error, full-sample bootstrap",
               "Ratio: bootstrap SE / influence-function SE", "95% CI, influence function (lower)", "95% CI, influence function (upper)",
               "95% CI, bootstrap percentile (lower)", "95% CI, bootstrap percentile (upper)", "Bootstrap mean of the estimate",
               "Bootstrap resamples used"),
  value = c(aip$rd, aip$se_rd, sd(bs[, "rd"]), sd(bs[, "rd"]) / aip$se_rd, aip$rd - 1.96 * aip$se_rd, aip$rd + 1.96 * aip$se_rd,
            quantile(bs[, "rd"], 0.025), quantile(bs[, "rd"], 0.975), mean(bs[, "rd"]), nrow(bs)))
cat("\n7.3 INFLUENCE FUNCTION VERSUS BOOTSTRAP\n"); print(rd2(var_tab, 4)); save_tab(var_tab, "variance_check.csv")
saveRDS(bs, "results/tables/bootstrap_draws.rds")
cat("Reading: a ratio near 1 means the influence-function standard error is reliable; below 1 means it is slightly conservative.\n")


## ======= SECTION 8. FLEXIBLE MODELS WITH CROSS-FITTING ==========
## The main analysis uses logistic regressions with main effects only. Here the same AIPTW estimator is run with
## flexible learners that find curved and interacting patterns without being told where to look:
##   Logistic_main        main-effects logistic regression (the benchmark)
##   Lasso_interactions   logistic regression with all two-way interactions, shrunk by the LASSO penalty
##   Random_forest        probability forest
## Flexible learners can over-fit, so they are cross-fitted (double machine learning, Chernozhukov et al. 2018):
##   1. split births at random into K folds;
##   2. for fold k, train all three models on the other K - 1 folds;
##   3. predict for fold k only, so nobody is scored by a model that saw them;
##   4. pool the folds and compute AIPTW once.
cx  <- make_cells(d); cx[, cid := .I]                         # cells with an id
d[cx, on = CONF, cid := i.cid]                                # cell id of every birth
sm  <- d$smoker; pt <- d$preterm; cid <- d$cid; NN <- length(cid)
cxDF <- as.data.frame(cx[, ..CONF])
X2 <- sparse.model.matrix(~ .^2, data = cxDF)[, -1]           # all dummy variables and all two-way interactions
cat("\n8. LASSO design matrix:", nrow(X2), "cells x", ncol(X2), "columns\n")

make_folds <- function(K) {                                   # random fold of every birth, then four counts per cell per fold
  fold_id <- sample.int(K, NN, replace = TRUE)
  cnt1 <- function(mask, k) tabulate(cid[mask & fold_id == k], nbins = nrow(cx))
  s1 <- sm == 1L; p1 <- pt == 1L
  list(fold_id = fold_id, folds = lapply(1:K, function(k)
    list(y1 = cnt1(s1 & p1, k), n1 = cnt1(s1, k), y0 = cnt1(!s1 & p1, k), n0 = cnt1(!s1, k))))
}
clip <- function(p, lo = 1e-4) pmin(pmax(p, lo), 1 - lo)       # keep probabilities strictly inside (0, 1)

learner_logit <- function(tr, in_train) {                     # (a) main-effects logistic
  tt <- cbind(cx[, ..CONF], tr)
  fit_glm_cells(tt, RHS_MAIN, tt)[c("e", "mu1", "mu0")]
}
learner_lasso <- function(tr, in_train) {                     # (b) LASSO with all two-way interactions
  one <- function(y, n) {                                     # y events out of n births, per cell
    ok <- n > 0
    cv <- cv.glmnet(X2[ok, ], cbind(n[ok] - y[ok], y[ok]), family = "binomial", nfolds = 3, nlambda = 25)
    clip(as.vector(predict(cv, X2, s = "lambda.min", type = "response")))
  }
  list(e = one(tr$n1, tr$n1 + tr$n0), mu1 = one(tr$y1, tr$n1), mu0 = one(tr$y0, tr$n0))
}
learner_forest <- function(tr, in_train) {                    # (c) random forest, trained on individual births
  frame <- function(rows, y) { x <- cxDF[cid[rows], , drop = FALSE]; x$y <- factor(y[rows]); x }
  grow  <- function(dat, sf, mns) ranger(y ~ ., data = dat, probability = TRUE, num.trees = 100, mtry = 5,
             min.node.size = mns, sample.fraction = sf, respect.unordered.factors = "order", num.threads = CORES, seed = SEED)
  pr <- function(rf) clip(predict(rf, cxDF)$predictions[, "1"])
  ## Propensity score: all training smokers but only 10% of training non-smokers (saves memory); the odds are
  ## corrected afterwards by the factor 0.10, the standard correction for this kind of sampling.
  f    <- 0.10
  keep <- in_train & (sm == 1L | runif(NN) < f)
  e_s  <- pr(grow(frame(keep, sm), 0.25, 100))
  e    <- clip((e_s / (1 - e_s)) * f / (1 + (e_s / (1 - e_s)) * f))
  mu1  <- pr(grow(frame(in_train & sm == 1L, pt), 0.50, 50))                          # outcome among smokers
  mu0  <- pr(grow(frame(in_train & sm == 0L & runif(NN) < 0.25, pt), 0.25, 100))      # outcome among non-smokers (25% sample)
  list(e = e, mu1 = mu1, mu0 = mu0)
}
learners <- list(Logistic_main = learner_logit, Lasso_interactions = learner_lasso, Random_forest = learner_forest)

cf_rows <- list()
for (rep in seq_len(REPS)) {
  mf   <- make_folds(K)
  tot  <- list(y1 = cx$y1, n1 = cx$n1, y0 = cx$y0, n0 = cx$n0)
  pred <- lapply(learners, function(f) vector("list", K))
  for (k in 1:K) {
    tr <- lapply(names(tot), function(nm) tot[[nm]] - mf$folds[[k]][[nm]]); names(tr) <- names(tot)   # counts of the training folds
    for (ln in names(learners)) {
      t1 <- Sys.time(); pred[[ln]][[k]] <- learners[[ln]](as.data.table(tr), mf$fold_id != k)
      cat(sprintf("  repeat %d fold %d %-18s %.0f s\n", rep, k, ln, as.numeric(Sys.time() - t1, units = "secs")))
    }
  }
  pred$Average_of_3 <- lapply(1:K, function(k) lapply(c(e = "e", mu1 = "mu1", mu0 = "mu0"), function(q)
    Reduce(`+`, lapply(names(learners), function(ln) pred[[ln]][[k]][[q]])) / length(learners)))   # simple ensemble
  for (ln in c(names(learners), "Average_of_3")) {
    st <- do.call(rbind, lapply(1:K, function(k) data.table(                                       # stack the K held-out tables
      y1 = mf$folds[[k]]$y1, n1 = mf$folds[[k]]$n1, y0 = mf$folds[[k]]$y0, n0 = mf$folds[[k]]$n0,
      e = pred[[ln]][[k]]$e, mu1 = pred[[ln]][[k]]$mu1, mu0 = pred[[ln]][[k]]$mu0)))
    for (cl in list(c(0.01, 0.99), NULL)) {
      r  <- aiptw_cells(st, st$e, st$mu1, st$mu0, clamp = cl)
      nn <- st$n1 + st$n0; ee <- r$e_used; ws <- 1 / ee
      cf_rows[[length(cf_rows) + 1]] <- data.table(
        rep = rep, learner = ln, clamp = if (is.null(cl)) "none" else "1%/99%",
        RD = r$rd, SE = r$se_rd, lo = r$rd - 1.96 * r$se_rd, hi = r$rd + 1.96 * r$se_rd, RR = r$rr,
        max_weight_smokers = max(ws[st$n1 > 0]), ESS_smokers = sum(st$n1 * ws)^2 / sum(st$n1 * ws^2),
        RD_outcome_only   = 100 * sum(nn * (st$mu1 - st$mu0)) / sum(nn),
        RD_weighting_only = 100 * (sum(st$y1 / ee) / sum(st$n1 / ee) - sum(st$y0 / (1 - ee)) / sum(st$n0 / (1 - ee))),
        pred_risk_smokers = 100 * sum(st$n1 * st$mu1) / sum(st$n1), obs_risk_smokers = 100 * sum(st$y1) / sum(st$n1),
        pred_risk_nonsmokers = 100 * sum(st$n0 * st$mu0) / sum(st$n0), obs_risk_nonsmokers = 100 * sum(st$y0) / sum(st$n0),
        pred_smoking_prob = 100 * sum(nn * st$e) / sum(nn), obs_smoking_prob = 100 * sum(st$n1) / sum(nn))
    }
  }
}
cf <- rbindlist(cf_rows); save_tab(cf, "crossfit_by_repeat.csv")
cf_sum <- cf[, .(RD = mean(RD), SE = mean(SE), lo = mean(lo), hi = mean(hi), RR = mean(RR),
                 range_over_repeats = max(RD) - min(RD), max_weight = mean(max_weight_smokers), ESS = mean(ESS_smokers),
                 RD_outcome_only = mean(RD_outcome_only), RD_weighting_only = mean(RD_weighting_only),
                 pred_risk_smokers = mean(pred_risk_smokers), obs_risk_smokers = mean(obs_risk_smokers),
                 pred_risk_nonsmokers = mean(pred_risk_nonsmokers), obs_risk_nonsmokers = mean(obs_risk_nonsmokers),
                 pred_smoking_prob = mean(pred_smoking_prob), obs_smoking_prob = mean(obs_smoking_prob)), by = .(learner, clamp)]
cat("\n8. CROSS-FITTED AIPTW\n"); print(rd2(cf_sum, 3)); save_tab(cf_sum, "crossfit_summary.csv")
cat("Reading: similar risk differences across learners mean the answer does not depend on the main-effects logistic form.\n",
    "Predicted and observed risks should match (calibration). A learner whose weighting part departs from the others\n",
    "(see RD_weighting_only) is estimating the propensity score poorly and is better read as a stress test.\n", sep = "")
rm(X2); invisible(gc())


## ======= SECTION 9. SMOKING HISTORY AND CIGARETTES PER DAY ======
## The extract has two smoking variables:
##   cig_0   = cigarettes per day in the 3 months BEFORE pregnancy (0 to 98; 99 = unknown)
##   smoker  = smoked at any time DURING pregnancy
## They support two exploratory questions:
##   (A) Cessation: among women who smoked before pregnancy, does continuing (versus quitting) change preterm risk?
##   (B) Intensity: among women who kept smoking, does risk rise with the number of cigarettes per day before
##       pregnancy, compared with never-smokers? (Pre-pregnancy intensity, not intensity during pregnancy.)
## All results in Sections 9 and 10 are exploratory: the questions were added after the main analysis was fixed.

d[, pre_smk := fifelse(cig_0 == 99, NA_integer_, as.integer(cig_0 > 0))]        # smoked before pregnancy; 99 = unknown -> NA
d[, hist := fcase(pre_smk == 0 & smoker == 0, "Never",
                  pre_smk == 1 & smoker == 0, "Quit",
                  pre_smk == 1 & smoker == 1, "Continued",
                  pre_smk == 0 & smoker == 1, "Started in pregnancy")]
d[, int_grp := cut(cig_0, c(0, 5, 10, 20, 98), labels = c("1-5", "6-10", "11-20", "21+"))]   # cigarettes per day before pregnancy
crude_h <- d[!is.na(hist), .(births = .N, preterm_pct = round(100 * mean(preterm), 2)), by = hist][order(-preterm_pct)]
cat("\n9.1 CRUDE PRETERM RISK BY SMOKING HISTORY\n"); print(crude_h); save_tab(crude_h, "history_crude.csv")

## contrast(): AIPTW for "group tt = 1 versus tt = 0" inside any subset of births, adjusting for `vars`.
contrast <- function(dsub, vars, label) {
  x  <- copy(dsub)[, smoker := tt]                                 # make_cells() reads the column called smoker
  cl <- make_cells(x, vars); f <- paste(vars, collapse = " + ")
  fit <- fit_glm_cells(cl, f)
  r <- aiptw_cells(cl, fit$e, fit$mu1, fit$mu0); a <- att_cells(cl, fit$e, fit$mu0)
  data.table(contrast = label, n_group = sum(cl$n1), n_reference = sum(cl$n0),
             crude_RD = 100 * (sum(cl$y1) / sum(cl$n1) - sum(cl$y0) / sum(cl$n0)),
             RD = r$rd, lo = r$rd - 1.96 * r$se_rd, hi = r$rd + 1.96 * r$se_rd,
             ATT = a$att, ATT_lo = a$att - 1.96 * a$se, ATT_hi = a$att + 1.96 * a$se)
}
hist_out <- list()
pre <- d[pre_smk == 1 & !is.na(int_grp)]; pre[, tt := as.integer(smoker == 1)]              # (A) women who smoked before pregnancy
hist_out[[1]] <- contrast(pre, c(CONF, "int_grp"), "A. Continued vs quit (also adjusted for cigarettes/day before)")
hist_out[[2]] <- contrast(pre, setdiff(c(CONF, "int_grp"), c("precare", "wic")), "A2. Continued vs quit, pre-pregnancy confounders only")
never <- d[hist == "Never"]
for (g in levels(d$int_grp)) {                                                              # (B) each dose group versus never-smokers
  sub <- rbind(d[hist == "Continued" & int_grp == g], never, fill = TRUE)
  sub[, tt := as.integer(hist != "Never")]
  hist_out[[length(hist_out) + 1]] <- contrast(sub, CONF, paste0("B. Continued, ", g, " cigarettes/day before pregnancy, vs never"))
}
sub <- rbind(d[hist == "Quit"], never, fill = TRUE); sub[, tt := as.integer(hist != "Never")]
hist_out[[length(hist_out) + 1]] <- contrast(sub, CONF, "B2. Quit vs never-smokers")
hist_tab <- rbindlist(hist_out)
cat("\n9.2 CESSATION AND INTENSITY (risk difference in points)\n"); print(rd2(hist_tab, 2)); save_tab(hist_tab, "history_intensity_results.csv")
## How to read the table. RD is the average treatment effect over the combined group (smoking group plus reference group).
## When one group is very small compared with the other (for example 21+ cigarettes/day against 3 million never-smokers) the
## ATE leans on mothers who resemble nobody in the small group and becomes unstable (see the wide interval for 21+).
## ATT, the effect among the mothers in the smoking group itself, is stable and is the reading used for the dose and the
## quit-versus-never comparisons on the poster. Row A (continued versus quit) has two groups of similar size, so RD is fine there.


## ======= SECTION 10. TRIMESTER OF SMOKING ===============
## Trimester-specific smoking (cig_1, cig_2, cig_3 = cigarettes per day in each trimester) is in the FULL natality
## file but not in the 24-column extract. This section downloads the full file once (about 200 MB zipped, 2.2 GB
## unzipped), keeps only the columns needed, and rebuilds the same cohort. It can take 5 to 10 minutes.
## To skip it, set RUN_TRIMESTER <- FALSE in Section 0.
FULL_ZIP <- "data/natality2023us.zip"; FULL_CSV <- "data/natality2023us.csv"
FULL_URL <- "https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip"
FULL_COLS <- c("dob_yy", "mager", "mrace6", "dmar", "meduc", "precare5", "wic", "priorlive", "priordead", "priorterm",
               "rf_ppterm", "rf_pdiab", "rf_phype", "rf_gdiab", "rf_ghype", "bmi", "cig_0", "cig_1", "cig_2", "cig_3",
               "cig_rec", "dplural", "sex", "oegest_comb", "dbwt", "pay_rec", "restatus")
MANUAL_MSG <- paste0("To do it by hand: open https://data.nber.org/nvss/natality/csv/2023/ in your browser, click natality2023us.zip ",
                     "(about 190 MB), save it in the folder 'data' next to this script (or in your Downloads folder), and run the script again.")

## find_full_file(): use a copy that is already on the computer (data folder, working folder or Downloads) before downloading.
find_full_file <- function() {
  places <- c("data", ".", "~/Downloads", file.path(Sys.getenv("USERPROFILE"), "Downloads"))
  for (p in places) for (f in c("natality2023us.csv", "natality2023us.zip")) {
    fp <- file.path(p, f); if (file.exists(fp) && file.size(fp) > 1e8) return(fp)
  }
  NULL
}
step <- function(msg) cat(format(Sys.time(), "%H:%M:%S"), msg, "\n")

## The compact trimester file (the 27 columns needed, all 2023 births) is kept in the repository so that nobody has to download
## the 2 GB national file. It is built once from the full NCHS file (below) and uploaded to the repository folder data/.
TRI_GZ  <- "data/nvss2023_trimester_extract.csv.gz"
TRI_URL <- "https://github.com/fmurital/Causal-Inference-/raw/main/data/nvss2023_trimester_extract.csv.gz"

read_full_raw <- function() {
  if (!file.exists(TRI_GZ)) {                                              # 1. the compact file stored in the repository
    step("Looking for the compact trimester file in the repository ...")
    tryCatch(suppressWarnings(download.file(TRI_URL, TRI_GZ, mode = "wb", quiet = TRUE)), error = function(e) NULL)
    if (file.exists(TRI_GZ) && file.size(TRI_GZ) < 1e6) file.remove(TRI_GZ)     # an error page, not the data
  }
  if (file.exists(TRI_GZ)) { step("Reading the compact trimester file ..."); return(fread(TRI_GZ, colClasses = list(character = CHAR_COLS))) }
  loc <- find_full_file(); csv <- if (!is.null(loc) && grepl("\\.csv$", loc)) loc else NULL      # 2. the full NCHS file
  if (is.null(csv)) {
    zp <- if (!is.null(loc)) loc else FULL_ZIP
    if (is.null(loc)) {
      step("Downloading natality2023us.zip (about 190 MB) ..."); options(timeout = 7200)
      tryCatch(download.file(FULL_URL, FULL_ZIP, mode = "wb"), error = function(e) message("Download error: ", conditionMessage(e)))
      if (!file.exists(FULL_ZIP) || file.size(FULL_ZIP) < 1e8)
        stop("the download did not complete (file missing or too small). ", MANUAL_MSG)
    } else step(paste("Using the file already on this computer:", loc))
    inside <- utils::unzip(zp, list = TRUE)$Name; inside <- inside[grepl("\\.csv$", inside, ignore.case = TRUE)]
    if (!length(inside)) stop("the zip file holds no .csv file. ", MANUAL_MSG)
    step(paste("Unzipping", inside[1], "(a few minutes, about 2 GB) ..."))
    utils::unzip(zp, files = inside[1], exdir = "data", junkpaths = TRUE)
    csv <- file.path("data", basename(inside[1]))
    if (!file.exists(csv)) stop("unzipping failed; check free disk space (about 3 GB needed).")
  }
  hdr <- names(fread(csv, nrows = 0))
  lack <- setdiff(FULL_COLS, hdr)
  if (length(lack)) stop("columns not found in the full file: ", paste(lack, collapse = ", "),
                         ". Similar names: ", paste(grep("cig", hdr, value = TRUE), collapse = ", "))
  step("Reading the 27 needed columns ...")
  fr <- fread(csv, select = FULL_COLS, colClasses = list(character = CHAR_COLS))
  fwrite(fr, TRI_GZ, compress = "gzip")                                    # keep the compact copy for the repository
  step(sprintf("Saved %s (%.0f MB). Upload this file to the repository folder data/ so that others can skip the large download.",
               TRI_GZ, file.size(TRI_GZ) / 1e6))
  fr
}

load_full_file <- function() {
  fr <- read_full_raw()
  cat("Rows in the full file:", format(nrow(fr), big.mark = ","), "\n")
  fr <- fr[dplural == 1 & oegest_comb < 99 & cig_rec %in% c("Y", "N") & as.integer(restatus) != 4L]     # same rules as the extract
  dfull <- build_cohort(fr)
  ok <- complete.cases(dfull[, ..CONF])
  dfull[ok]
}
tri_ok <- RUN_TRIMESTER && tryCatch({ dT <- load_full_file(); TRUE },
  error = function(e) { cat("\n*** Section 10 skipped: ", conditionMessage(e), "\n", sep = ""); FALSE })

if (tri_ok) {
  ## safe_contrast(): one group that is too small or fails to fit is reported and skipped; it does not stop the run
  safe_contrast <- function(sub, vars, label) {
    if (sum(sub$tt == 1) < 500) { cat("  skipped (fewer than 500 births in the group):", label, "\n"); return(NULL) }
    tryCatch(contrast(sub, vars, label), error = function(e) { cat("  skipped (", conditionMessage(e), "):", label, "\n"); NULL })
  }
  cat("Cohort rebuilt from the full file:", format(nrow(dT), big.mark = ","), "births (extract cohort: 3,259,496)\n")
  dT <- dT[cig_0 < 99 & cig_1 < 99 & cig_2 < 99 & cig_3 < 99]                         # known cigarettes in every period
  dT[, pattern := fcase(
        cig_0 == 0 & cig_1 == 0 & cig_2 == 0 & cig_3 == 0,                "Never smoked",
        cig_1 > 0  & cig_2 == 0 & cig_3 == 0,                              "First trimester only",
        cig_1 > 0  & cig_2 > 0  & cig_3 == 0,                              "First and second trimesters only",
        cig_1 > 0  & cig_2 > 0  & cig_3 > 0,                               "All three trimesters",
        cig_1 == 0 & (cig_2 > 0 | cig_3 > 0),                              "Started after the first trimester")]
  cat("\n10.1 PRETERM RISK BY TRIMESTER PATTERN (crude)\n")
  tri_crude <- dT[!is.na(pattern), .(births = .N, preterm_pct = round(100 * mean(preterm), 2)), by = pattern][order(preterm_pct)]
  print(tri_crude); save_tab(tri_crude, "trimester_crude.csv")
  ref <- dT[pattern == "Never smoked"]
  tri_out <- list()
  for (pt_name in c("First trimester only", "First and second trimesters only", "All three trimesters", "Started after the first trimester")) {
    sub <- rbind(dT[pattern == pt_name], ref, fill = TRUE); sub[, tt := as.integer(pattern != "Never smoked")]
    tri_out[[length(tri_out) + 1]] <- safe_contrast(sub, CONF, paste0(pt_name, " vs never"))
  }
  dT[, cig3_grp := cut(cig_3, c(0, 5, 10, 20, 98), labels = c("1-5", "6-10", "11-20", "21+"))]   # cigarettes per day, third trimester
  for (g in levels(dT$cig3_grp)) {
    sub <- rbind(dT[pattern == "All three trimesters" & cig3_grp == g], ref, fill = TRUE); sub[, tt := as.integer(pattern != "Never smoked")]
    tri_out[[length(tri_out) + 1]] <- safe_contrast(sub, CONF, paste0("Smoked all trimesters, ", g, " cigarettes/day in 3rd trimester, vs never"))
  }
  tri_tab <- rbindlist(Filter(Negate(is.null), tri_out))
  cat("\n10.2 TRIMESTER TIMING AND THIRD-TRIMESTER DOSE (AIPTW, points)\n"); print(rd2(tri_tab, 2)); save_tab(tri_tab, "trimester_results.csv")
  rm(dT, ref); invisible(gc())
} else cat("\nSection 10 not run (RUN_TRIMESTER = FALSE or the full file was not available).\n")


## ======= SECTION 11. UNMEASURED CONFOUNDING (E-value and benchmarking) ======
## E-value: the smallest risk-ratio association that an unmeasured confounder would need with BOTH smoking and
## preterm birth, beyond the measured confounders, to fully explain away the estimate (VanderWeele and Ding 2017).
rr_ci_lo <- exp(log(aip$rr) - 1.96 * aip$se_lrr)
ev <- data.table(estimate = c("AIPTW risk ratio", "Lower 95% limit of the risk ratio"), risk_ratio = c(aip$rr, rr_ci_lo),
                 E_value = c(evalue_rr(aip$rr), evalue_rr(rr_ci_lo)))
cat("\n11.1 E-VALUES\n"); print(rd2(ev, 2)); save_tab(ev, "evalues.csv")

## Benchmark: how strong are the measured confounders on the same scale? For each, the risk ratio with smoking and with
## preterm birth (among non-smokers), and how much the adjusted estimate changes when it is left out (Section 4).
defs <- list("Education: high school or less" = d$educ %in% c("<9th grade", "9-12th, no diploma", "HS grad/GED"),
             "WIC participation" = d$wic == "Yes", "Previous preterm birth" = d$ppterm == "Yes",
             "Pre-pregnancy diabetes" = d$pdiab == "Yes", "Pre-pregnancy hypertension" = d$phype == "Yes",
             "BMI 40 or more" = d$bmi_grp == "40+", "Parity 3 or more" = d$parity == "3+",
             "Care: none or 3rd trimester" = d$precare %in% c("No care", "3rd trimester"), "Age under 20" = d$age_grp == "<20")
bench <- rbindlist(lapply(names(defs), function(nm) { x <- defs[[nm]]
  data.table(covariate = nm, pct_of_cohort = 100 * mean(x),
             RR_with_smoking = mean(d$smoker[x] == 1) / mean(d$smoker[!x] == 1),
             RR_with_preterm = mean(d$preterm[x & d$smoker == 0]) / mean(d$preterm[!x & d$smoker == 0]))
}))
bench[, exceeds_E_value_on_both := RR_with_smoking >= ev$E_value[1] & RR_with_preterm >= ev$E_value[1]]
bench[, var_key := c("educ", "wic", "ppterm", "pdiab", "phype", "bmi_grp", "parity", "precare", "age_grp")]
bench <- merge(bench, contrib[, .(var_key = variable, change_in_rd_when_left_out = change_in_rd)], by = "var_key", sort = FALSE)[, var_key := NULL]
cat("\n11.2 MEASURED CONFOUNDERS ON THE E-VALUE SCALE\n"); print(rd2(bench, 2)); save_tab(bench, "evalue_benchmark.csv")
cat(sprintf("Reading: left out one at a time, the measured confounders move the estimate by %.2f points at most (%s). Previous preterm birth is the only\n",
            max(abs(bench$change_in_rd_when_left_out)), bench$covariate[which.max(abs(bench$change_in_rd_when_left_out))]),
    sprintf("one above the E-value on both associations, and leaving it out moves the estimate by %.2f points.\n", bench[covariate == "Previous preterm birth"]$change_in_rd_when_left_out), sep = "")


## ======= SECTION 12. PUBLIC-HEALTH TRANSLATION AND FINAL SUMMARY ============
n_sm <- sum(d$smoker)
ph <- data.table(quantity = c("Excess preterm births per 1,000 smoking pregnancies (ATT)", "Smoking pregnancies in the cohort",
                              "Preterm births attributable to smoking in the cohort (ATT x smokers)", "One extra preterm birth for every ... smoking pregnancies",
                              "Excess preterm births per 1,000 pregnancies if all mothers smoked versus none (ATE)"),
                 value = c(10 * att$att, n_sm, att$att / 100 * n_sm, 100 / att$att, 10 * aip$rd))
cat("\n12.1 PUBLIC-HEALTH TRANSLATION\n"); print(rd2(ph, 1)); save_tab(ph, "public_health_translation.csv")

summary_numbers <- data.table(item = c("Births in the analytic cohort", "Smokers", "Crude risk difference (points)",
  "AIPTW risk difference (points)", "AIPTW 95% CI lower", "AIPTW 95% CI upper", "AIPTW risk ratio", "ATT (points)",
  "SE, influence function", "SE, bootstrap", "E-value (risk ratio)", "E-value (lower limit)",
  "Share of non-smokers with score below 0.01 (%)", "Estimates across analytic choices, minimum", "Estimates across analytic choices, maximum"),
  value = c(N, n_sm, crude_rd, aip$rd, aip$rd - 1.96 * aip$se_rd, aip$rd + 1.96 * aip$se_rd, aip$rr, att$att,
            aip$se_rd, sd(bs[, "rd"]), ev$E_value[1], ev$E_value[2], share_low_ns, min(robust$RD[-1]), max(robust$RD[-1])))
cat("\n12.2 NUMBERS USED ON THE POSTER AND IN THE MANUSCRIPT\n"); print(rd2(summary_numbers, 2)); save_tab(summary_numbers, "summary_numbers.csv")

concern_map <- data.table(
  concern = c("1 Adjustment set", "2 Overlap and imbalance", "3 Variance and bootstrap", "4 Flexible, cross-fitted models",
              "5 Novelty and public health", "6 Intensity, timing, cessation"),
  where   = c("Section 5.4, 5.6", "Section 6", "Section 7", "Section 8", "Section 12.1 (numbers); introduction text", "Sections 9 and 10"),
  result  = c(sprintf("Pre-pregnancy set: %.2f vs main %.2f", spec_tab$RD[2], spec_tab$RD[1]),
              sprintf("%.0f%% of non-smokers have score < 0.01; largest SMD after weighting %.3f", share_low_ns, max(abs(bal$after_iptw))),
              sprintf("IF SE %.3f, bootstrap SE %.3f (ratio %.2f)", aip$se_rd, sd(bs[, "rd"]), sd(bs[, "rd"]) / aip$se_rd),
              sprintf("Cross-fitted RD: %s", paste(sprintf("%s %.2f", cf_sum[clamp == "1%/99%"]$learner, cf_sum[clamp == "1%/99%"]$RD), collapse = "; ")),
              sprintf("%.0f excess preterm births per 1,000 smoking pregnancies", 10 * att$att),
              sprintf("Continued vs quit: %.2f; quit vs never-smokers: %.2f among quitters (whole-population version %.2f)", hist_tab$RD[1], hist_tab[contrast == "B2. Quit vs never-smokers"]$ATT, hist_tab[contrast == "B2. Quit vs never-smokers"]$RD)))
cat("\nCONCERN MAP\n"); print(concern_map, right = FALSE); save_tab(concern_map, "concern_map.csv")
cat("\nAll tables are in results/tables and all figures in results/figures.\n")
sessionInfo()
