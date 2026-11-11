# Maternal smoking in pregnancy and preterm birth: doubly robust (AIPTW) analysis

Code for the SEASUG 2026 paper by Faruk Muritala and Dhrubajyoti Ghosh, Kennesaw State University.
Data: 2023 U.S. Natality Public Use File (NCHS), 3,476,180 singleton births, 3,259,496 with complete data on all ten confounders.

Everything uses R (tested with R 4.3.3) and a fixed random seed of 1234.

## What you need

- R 4.3 or newer.
- R packages: `data.table` and `ggplot2` only (everything else is base R).
  Install once in R: `install.packages(c("data.table", "ggplot2"))`
- About 8 GB of RAM and 2 CPU cores. The main analysis takes a few minutes; the sensitivity script takes 45 to 60 minutes; the bootstrap check about 40 minutes.

## Fastest way to reproduce the paper

1. Download this folder and keep the structure unchanged. The data file must be in `data/` as either `nvss2023_extract.csv.gz` (35.8 MB, included in the GitHub copy) or the unzipped `nvss2023_extract.csv` (241 MB, 3,476,181 lines including the header). If your copy of the code has no data file, place the extract there, or build it with `00_extract_data.R`.
2. In a terminal, go into the folder and run: `bash run_all.sh`
   On Windows without bash, open R, set the working directory to this folder with `setwd("path/to/this/folder")`, then run the scripts below one at a time with `source("01_prepare_data.R")` and so on.
3. Results are written to `output/` (CSV tables), `figures/` (PNG graphics) and `logs/`.

## Scripts, in run order

| Script | What it does |
|---|---|
| `00_extract_data.R` | Optional. Rebuilds the 24-column extract from the original 2.2 GB government file. Not needed to reproduce results. |
| `01_prepare_data.R` | Cleans variables, drops non-U.S. residents and incomplete records, writes the cohort flow table. |
| `02_eda.R` | Table 1 and the exploratory figures. |
| `03_main_analysis.R` | Propensity model, five estimators (crude, g-computation, stabilized IPTW, PS stratification, AIPTW), ATT, overlap weights, balance, E-value, main figures. |
| `04_sensitivity.R` | Propensity score handling grid, four alternative covariate specifications, subgroups (slow). |
| `04b_sensitivity_figures.R` | Draws the robustness and subgroup figures from the output of script 04. |
| `06_draw_dag.R` | Draws the causal diagram. |
| `07_paper_numbers.R` | Collects every number quoted in the paper into `output/paper_numbers.csv`. |
| `05_variance_check.R` | Bootstrap check of the influence-function standard errors on a 200,000-birth subsample (slow; runs last). |
| `utils_causal.R`, `pipeline_functions.R` | Helper functions (chunked logistic regression, estimators, balance). |

## Data links

- Official file (fixed-width, 224 MB zip): https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip
- User guide with all variable codes: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/UserGuide2023.pdf
- NCHS data page: https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm
- CSV version used here (NBER): https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip
  (check that the link still works before relying on it; the header of `00_extract_data.R` explains each step)

The data are public-use and contain no personal identifiers.

## Notes on reading the results

- The primary estimate is the AIPTW risk difference, 3.26 percentage points (95% CI 2.47 to 4.04).
- Overlap between smokers and non-smokers is limited. The paper reports the effect under several propensity score handling rules and with overlap weights; read Section 6 before quoting the main number.
- Research code. Smoking is self-reported on the birth certificate; see the Limitations section of the paper.
