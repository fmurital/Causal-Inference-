# Data

Both files come from the 2023 U.S. Natality Public Use File of the National Center for Health Statistics. They are public, contain no personal identifiers, and are the only inputs the R code needs.

| File | Rows | Used for |
|---|---|---|
| `nvss2023_extract.csv.gz` (36 MB) | 3,476,180 singleton births, 24 columns | Every main analysis, the posters and the paper |
| `nvss2023_trimester_extract.csv.gz` (42 MB) | 3,605,081 births (all 2023 deliveries), 27 columns | The trimester analysis (Section 10 of the script) |

## How the code finds the data

`Smoking_Preterm_Full_Analysis.R` looks in its own `data` folder first. If the file is not there it downloads it from this repository, so a fresh copy of the script runs with no manual download. The links it uses are:

- `https://github.com/fmurital/Causal-Inference-/raw/main/data/nvss2023_extract.csv.gz`
- `https://github.com/fmurital/Causal-Inference-/raw/main/data/nvss2023_trimester_extract.csv.gz`

## The trimester file

The trimester of smoking is recorded in the national file (`cig_1`, `cig_2`, `cig_3`: cigarettes per day in each trimester) but not in the 24-column extract. The national file is about 190 MB zipped and 2.2 GB unzipped, which is too large for GitHub. The script therefore reads the national file once, keeps the 27 columns it needs, and writes them to `data/nvss2023_trimester_extract.csv.gz`. That compact file is the one stored in this repository. It is stored in this repository, so Section 10 runs without any large download. Rebuilding the analysis cohort from it gives the same 3,259,496 births as the extract.

National file: https://data.nber.org/nvss/natality/csv/2023/natality2023us.zip (the same data are at https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip).
