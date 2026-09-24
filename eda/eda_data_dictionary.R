## =============================================================================
## STAT 7900 -- Maternal Smoking & Preterm Birth
## Standalone Data Dictionary + Exploratory Data Analysis script
## Author: Faruk Muritala (with Claude assistance; see AI-use disclosure in paper)
##
## WHAT THIS SCRIPT DOES
##   1. Loads the 2023 NVSS natality extract (fast, memory-efficient read).
##   2. Recodes NCHS "unknown / not stated" sentinel codes to real NA, and
##      derives the treatment (T), outcome (Y), and nulliparous indicator.
##   3. Builds a full data dictionary (variable, label, type, coding,
##      % missing) and writes it to a CSV you can open directly in Excel.
##   4. Produces missingness and univariate summaries for every variable.
##   5. Produces a correlation matrix + heatmap across T, Y, and the 10
##      adjustment-set covariates.
##   6. Saves everything (CSVs + PNG plots) into OUTPUT_DIR.
##
## HOW TO RUN THIS ON YOUR OWN LAPTOP
##   1. Open this file in RStudio.
##   2. Edit ONLY the two lines under "USER SETTINGS" below -- set DATA_PATH
##      to wherever nvss2023_extract.csv actually sits on YOUR computer.
##   3. Click "Source" (or Ctrl+Shift+Enter). That's it -- everything else
##      is self-contained: no other files, no working-directory assumptions,
##      no relative paths. This is exactly why earlier code broke for you:
##      it depended on files/paths that only existed in the assistant's own
##      session. This script depends on nothing but the one CSV you point it
##      at and CRAN (for package installation, first run only).
## =============================================================================


## ------------------------------- USER SETTINGS ------------------------------
## EDIT THESE TWO LINES ONLY. Everything below this section runs unattended.

# Full path to the NVSS extract on YOUR machine. Use forward slashes "/" even
# on Windows (R accepts this, and it avoids the backslash-escaping trap).
DATA_PATH  <- "C:/Users/fmurital/OneDrive - Kennesaw State University/MS KSU/Fall 2026/SpTp Causal Analysis XLS Group STAT 7900 Fall Semester 2026 CO - 9202026 - 711 AM/Project Folder/nvss2023_extract.csv"

# Folder where the data dictionary CSV, missingness table, correlation matrix,
# and plots will be saved. Defaults to a new "eda_output" subfolder next to
# the data file -- change it if you want the results somewhere else.
OUTPUT_DIR <- file.path(dirname(DATA_PATH), "eda_output")

## -----------------------------------------------------------------------------


## ------------------------------- 0. SETUP ------------------------------------
cat("=== STEP 0: checking/installing required packages (first run only) ===\n")

required_pkgs <- c("data.table", "ggplot2", "reshape2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  installing '%s' ...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
})

if (!file.exists(DATA_PATH)) {
  stop(
    "Cannot find the data file at:\n  ", DATA_PATH, "\n",
    "Fix: open this .R file, find the line that starts with DATA_PATH <-, and\n",
    "replace the path with the exact location of nvss2023_extract.csv on your\n",
    "own computer (right-click the file -> Properties/Get Info to copy the\n",
    "full path)."
  )
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat(sprintf("Output will be saved to: %s\n", OUTPUT_DIR))


## ------------------------------- 1. LOAD DATA --------------------------------
cat("\n=== STEP 1: reading data (fread is fast even for multi-million-row files) ===\n")
t0 <- Sys.time()
df <- fread(DATA_PATH, showProgress = TRUE)
cat(sprintf("Loaded %s rows x %s columns in %.1f seconds.\n",
            format(nrow(df), big.mark = ","), ncol(df),
            as.numeric(Sys.time() - t0, units = "secs")))


## ------------------------- 2. RECODE SENTINEL "UNKNOWN" CODES ----------------
## NCHS natality files use specific numeric/letter codes to mean "not stated"
## instead of a blank cell. Left un-recoded, these masquerade as real values
## (e.g. bmi = 99.9 would silently drag your mean BMI way up). This step turns
## every documented "unknown" code into a real NA, so summaries, correlations,
## and complete-case counts are all correct.
cat("\n=== STEP 2: recoding NCHS 'unknown/not stated' sentinel codes to NA ===\n")

df[, meduc_c    := fifelse(meduc == 9, NA_integer_, meduc)]
df[, bmi_c      := fifelse(bmi == 99.9, NA_real_, bmi)]
df[, precare5_c := fifelse(precare5 == 5, NA_integer_, precare5)]
df[, pay_rec_c  := fifelse(pay_rec == 9, NA_integer_, pay_rec)]
df[, priorlive_c := fifelse(priorlive == 99, NA_integer_, priorlive)]
df[, priordead_c := fifelse(priordead == 99, NA_integer_, priordead)]
df[, priorterm_c := fifelse(priorterm == 99, NA_integer_, priorterm)]

yn_to_bin <- function(x) fifelse(x == "Y", 1L, fifelse(x == "N", 0L, NA_integer_))
for (v in c("wic", "rf_ppterm", "rf_pdiab", "rf_phype", "rf_gdiab", "rf_ghype", "cig_rec")) {
  df[[paste0(v, "_c")]] <- yn_to_bin(df[[v]])
}

## Derived variables actually used in the analysis
df[, nulliparous := fifelse(
  priorlive_c == 0 & priordead_c == 0 & priorterm_c == 0, 1L,
  fifelse(is.na(priorlive_c) | is.na(priordead_c) | is.na(priorterm_c), NA_integer_, 0L)
)]
df[, preterm := as.integer(oegest_comb < 37)]   # outcome Y
df[, smoker  := cig_rec_c]                       # treatment T


## ------------------------------ 3. DATA DICTIONARY ----------------------------
cat("\n=== STEP 3: building the data dictionary ===\n")

## Coding reference is per the published NCHS Natality Public-Use File User
## Guide's standard category scheme (stable across recent data years). If you
## are citing this in a submitted document, cross-check the exact 2023 User
## Guide PDF for any year-specific changes -- flagged here for transparency,
## not verified against that PDF in this script.
dict <- data.table(
  variable = c("dob_yy","mager","mrace6","dmar","meduc","precare5","wic",
               "priorlive","priordead","priorterm","rf_ppterm","rf_pdiab","rf_phype",
               "rf_gdiab","rf_ghype","bmi","cig_0","cig_rec","dplural","sex",
               "oegest_comb","dbwt","pay_rec","restatus"),
  label = c(
    "Birth year",
    "Mother's age (single years)",
    "Mother's bridged race (6-category recode)",
    "Marital status",
    "Mother's education (2003 revision recode)",
    "Month prenatal care began (5-category recode)",
    "WIC participation during pregnancy",
    "Prior births now living",
    "Prior births now dead",
    "Prior other pregnancy terminations",
    "Risk factor: previous preterm birth",
    "Risk factor: pre-pregnancy diabetes",
    "Risk factor: pre-pregnancy hypertension",
    "Risk factor: gestational diabetes",
    "Risk factor: gestational hypertension",
    "Pre-pregnancy BMI",
    "Cigarettes/day before pregnancy",
    "Cigarette smoker recode (any smoking during pregnancy)",
    "Plurality (this extract is pre-restricted to singletons)",
    "Infant sex",
    "Obstetric estimate of gestation, combined (completed weeks)",
    "Birth weight (grams)",
    "Payment source for delivery",
    "Residence status relative to state of occurrence"
  ),
  role = c("context","confounder (adjustment set)","confounder (adjustment set)",
           "not used","confounder (adjustment set)","confounder (adjustment set)",
           "confounder (adjustment set)","used to derive parity","used to derive parity",
           "used to derive parity","confounder (adjustment set)","confounder (adjustment set)",
           "confounder (adjustment set)","excluded (gestational, possible mediator)",
           "excluded (gestational, possible mediator)","confounder (adjustment set)",
           "not used directly (see cig_rec)","TREATMENT (T)","restriction variable",
           "not used","used to derive OUTCOME (Y)","not used","considered & excluded -- see proposal",
           "not used"),
  coding = c(
    "2023 (single value)",
    "12-50 (years)",
    "1=NH White, 2=NH Black, 3=NH AIAN, 4=NH Asian, 5=NH NHOPI, 6=Hispanic",
    "1=Married, 2=Unmarried (blank in some states -- see missingness)",
    "1=8th grade or less ... 8=Doctorate/Professional, 9=Unknown",
    "1=1st trimester, 2=2nd, 3=3rd, 4=No care, 5=Unknown",
    "Y/N/U",
    "0-30 (count), 99=Not stated",
    "0-30 (count), 99=Not stated",
    "0-30 (count), 99=Not stated",
    "Y/N/U", "Y/N/U", "Y/N/U", "Y/N/U", "Y/N/U",
    "13.0-99.8 (kg/m^2), 99.9=Not stated",
    "0-98 (count), 99=Unknown",
    "Y = smoked at any point during pregnancy, N = did not",
    "1 in this extract (singletons only)",
    "M/F",
    "17-47 (weeks)",
    "227-9998 (grams), 9999=Not stated",
    "1=Medicaid, 2=Private, 3=Self-pay, 4=Other, 9=Unknown",
    "1=Resident, 2=Intrastate nonresident, 3=Interstate nonresident, 4=Foreign resident"
  ),
  source_note = "Coding per the published NCHS Natality Public-Use File User Guide's standard scheme; cross-check the exact-year PDF before citing precise category text in a submission."
)

## Attach REAL, computed missingness (post sentinel-recode where applicable)
miss_map <- c(
  dob_yy = "dob_yy", mager = "mager", mrace6 = "mrace6", dmar = "dmar",
  meduc = "meduc_c", precare5 = "precare5_c", wic = "wic_c",
  priorlive = "priorlive_c", priordead = "priordead_c", priorterm = "priorterm_c",
  rf_ppterm = "rf_ppterm_c", rf_pdiab = "rf_pdiab_c", rf_phype = "rf_phype_c",
  rf_gdiab = "rf_gdiab_c", rf_ghype = "rf_ghype_c", bmi = "bmi_c",
  cig_0 = "cig_0", cig_rec = "cig_rec_c", dplural = "dplural", sex = "sex",
  oegest_comb = "oegest_comb", dbwt = "dbwt", pay_rec = "pay_rec_c", restatus = "restatus"
)
dict[, pct_missing := sapply(miss_map[variable], function(col) round(100 * mean(is.na(df[[col]])), 3))]
dict[, n_missing := sapply(miss_map[variable], function(col) sum(is.na(df[[col]])))]

fwrite(dict, file.path(OUTPUT_DIR, "data_dictionary.csv"))
cat(sprintf("  wrote data_dictionary.csv (%d variables)\n", nrow(dict)))
print(dict[, .(variable, label, pct_missing)])


## ---------------------------- 4. UNIVARIATE SUMMARIES -------------------------
cat("\n=== STEP 4: univariate summaries (top categories / distribution) ===\n")
summary_lines <- c()
for (col in names(df)) {
  summary_lines <- c(summary_lines, sprintf("=== %s ===", col))
  if (is.numeric(df[[col]])) {
    s <- summary(df[[col]])
    summary_lines <- c(summary_lines, paste(capture.output(print(s)), collapse = "\n"))
  } else {
    tb <- sort(table(df[[col]], useNA = "ifany"), decreasing = TRUE)
    summary_lines <- c(summary_lines, paste(capture.output(print(head(tb, 15))), collapse = "\n"))
  }
  summary_lines <- c(summary_lines, "")
}
writeLines(summary_lines, file.path(OUTPUT_DIR, "univariate_summaries.txt"))
cat("  wrote univariate_summaries.txt\n")


## ---------------------------- 5. CORRELATION MATRIX ---------------------------
cat("\n=== STEP 5: correlation matrix across T, Y, and the adjustment-set covariates ===\n")

corr_vars <- data.table(
  label = c("smoker (T)", "preterm (Y)", "mager", "mrace6", "meduc", "precare5",
            "wic", "nulliparous", "prior_preterm", "pre_diabetes",
            "pre_hypertension", "bmi"),
  col   = c("smoker", "preterm", "mager", "mrace6", "meduc_c", "precare5_c",
            "wic_c", "nulliparous", "rf_ppterm_c", "rf_pdiab_c",
            "rf_phype_c", "bmi_c")
)
corr_df <- as.data.frame(df)[, corr_vars$col]
names(corr_df) <- corr_vars$label
corr_mat <- cor(corr_df, use = "pairwise.complete.obs")

write.csv(round(corr_mat, 4), file.path(OUTPUT_DIR, "correlation_matrix.csv"))
cat("  wrote correlation_matrix.csv\n")

## Heatmap
cm <- melt(corr_mat)
p <- ggplot(cm, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "#C1611A", mid = "white", high = "#1F3864",
                        midpoint = 0, limits = c(-1, 1), name = "r") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Correlation matrix: T, Y, and the 10 adjustment-set covariates",
       subtitle = "2023 NVSS extract, pairwise-complete observations")
ggsave(file.path(OUTPUT_DIR, "correlation_heatmap.png"), p, width = 9, height = 7.5, dpi = 200)
cat("  wrote correlation_heatmap.png\n")


## ------------------------------- 6. MISSINGNESS PLOT ---------------------------
cat("\n=== STEP 6: missingness plot ===\n")
miss_plot_df <- dict[order(-pct_missing)]
miss_plot_df$variable <- factor(miss_plot_df$variable, levels = rev(miss_plot_df$variable))
p2 <- ggplot(miss_plot_df, aes(x = variable, y = pct_missing)) +
  geom_col(fill = "#1F3864") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "Missingness by variable (2023 NVSS extract)",
       x = NULL, y = "% missing (after sentinel-code recode)")
ggsave(file.path(OUTPUT_DIR, "missingness_plot.png"), p2, width = 8, height = 7, dpi = 200)
cat("  wrote missingness_plot.png\n")


## ------------------------------- 7. COMPLETE-CASE COUNT ------------------------
cc_cols <- corr_vars$col[!(corr_vars$label %in% c("smoker (T)", "preterm (Y)"))]
cc_cols <- c("smoker", "preterm", cc_cols)
n_total <- nrow(df)
n_complete <- sum(complete.cases(as.data.frame(df)[, cc_cols]))
cat(sprintf(
  "\n=== SUMMARY ===\nTotal rows: %s\nComplete cases (T, Y + 10 adjustment-set covariates): %s (%.2f%% dropped)\n",
  format(n_total, big.mark=","), format(n_complete, big.mark=","),
  100 * (n_total - n_complete) / n_total
))

cat(sprintf("\nAll done. Open the folder below to see every output file:\n  %s\n", OUTPUT_DIR))
