# Smoking During Pregnancy and Preterm Birth

**A Doubly Robust Propensity Score Weighting Analysis of Millions of U.S. Births**

Faruk Muritala, Kennesaw State University. Advisor and co-author: Dhrubajyoti Ghosh, Ph.D.

This repository holds one research project at three stages: the SEASUG (CSOL) 2026 presentation, the two Kennesaw State conference posters that followed it, and the journal paper in preparation. Each stage has its own numbered folder with its own code and update notes. Every number comes from the R code in the repository, and the data the code reads are stored here, so anyone can run the analysis without a separate download.

## The question

Mothers who smoke during pregnancy deliver before 37 weeks more often than mothers who do not (13.9% against 8.3% in the 2023 U.S. natality file, a raw gap of 5.6 percentage points). Smokers also differ from non-smokers in education, age and health, which affect preterm birth too. The question is how much of the gap smoking itself causes. The estimate comes from Augmented Inverse Probability of Treatment Weighting (AIPTW), a doubly robust method that gives a valid estimate if either the smoking model or the birth-outcome model is correct.

**Main result (3,259,496 singleton births, 97,187 to mothers who smoked).** After adjustment the risk of preterm birth is 3.26 percentage points higher among mothers who smoked (bootstrap standard error 0.39; 95% interval 2.5 to 4.1), a gap two fifths smaller than the raw 5.6 points. Among mothers who smoked this is about 34 extra preterm births per 1,000 pregnancies.

## Where to start

| Folder | What it holds | Open it for |
|---|---|---|
| [`1_CSOL_SEASUG_2026_Presentation`](1_CSOL_SEASUG_2026_Presentation) | The SEASUG (CSOL) 2026 paper, slides, original R scripts and outputs, unchanged | the work as presented at the conference |
| [`2_Conference_Posters`](2_Conference_Posters) | The C Day 2026 and Analytic Day 2026 posters, the single analysis script, result tables, and a list of what changed after CSOL | the newest analyses and the posters |
| [`3_Journal_Paper`](3_Journal_Paper) | The journal manuscript (LaTeX and PDF), its code, tables and figures | the full write-up |
| [`data`](data) | The data files the code reads | running the code yourself |

## Run it yourself

1. Install R 4.3 or newer. The script installs the packages it needs (`data.table`, `ggplot2`, `glmnet`, `ranger`, `Matrix`) if they are missing.
2. Download `2_Conference_Posters/code/Smoking_Preterm_Full_Analysis.R`.
3. Open it in R or RStudio, set the working directory to the folder that holds the script, and run the whole script. On the first run it downloads `data/nvss2023_extract.csv.gz` from this repository, then writes every table to `results/tables` and every figure to `results/figures`. The run takes about one to two hours on a computer with several cores.

## Data

2023 U.S. Natality Public Use File, National Center for Health Statistics (National Vital Statistics System). The file is public and has no personal identifiers.

- `data/nvss2023_extract.csv.gz` is a 24-column extract (3,476,180 singleton births with known gestational age and smoking status) used for every main analysis.
- `data/nvss2023_trimester_extract.csv.gz` holds the cigarettes-per-day columns for each trimester, which the trimester analysis needs. See `data/README.md`.
- Official source: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2023us.zip. The user guide with every variable code is at https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/UserGuide2023.pdf.

## Contact

Faruk Muritala, Kennesaw State University. LinkedIn: https://www.linkedin.com/in/faruk-muritala-22b36b195
