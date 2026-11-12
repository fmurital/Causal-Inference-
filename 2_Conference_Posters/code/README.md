# Code

`Smoking_Preterm_Full_Analysis.R` is the single script behind both posters. It runs from top to bottom, writes every table to `results/tables` and every figure to `results/figures`, and prints each result as it goes. The random seed is 1234.

## How to run

1. Install R 4.3 or newer. Packages (`data.table`, `ggplot2`, `glmnet`, `ranger`, `Matrix`) are installed on the first run if missing.
2. Save the script in a folder of its own. The script creates `data` and `results` inside that folder.
3. In R or RStudio, set the working directory to that folder (RStudio: Session, Set Working Directory, Choose Directory) and run the whole script. On the first run it downloads `data/nvss2023_extract.csv.gz` from this repository.
4. Settings are in Section 0: `K` (cross-fitting folds), `REPS` (random splits), `B` (bootstrap resamples) and `RUN_TRIMESTER`.

The run takes about one to two hours on a computer with several cores; most of the time is the cross-fitted random forest and the bootstrap.

## Contents of the script

| Section | Content |
|---|---|
| 0 | Settings and helper functions |
| 1 | Load the data and build the analytic cohort |
| 2 | Data dictionary |
| 3 | Exploratory analysis, missing data, association matrix |
| 4 | Which confounders matter most |
| 5 | Main estimate and comparison estimators, adjustment sets |
| 6 | Overlap, weights and covariate balance |
| 7 | Variance: influence function against full-sample bootstrap |
| 8 | Flexible models with cross-fitting |
| 9 | Smoking history and cigarettes per day |
| 10 | Trimester of smoking (reads the trimester file in `../../data`; see `../../data/README.md`) |
| 11 | Unmeasured confounding (E-value and benchmarking) |
| 12 | Public-health translation and final summary |
