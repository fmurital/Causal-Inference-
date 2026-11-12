# Smoking During Pregnancy and Preterm Birth

A doubly robust propensity score weighting analysis of 3.26 million U.S. births (2023 natality file).

## Contribution

Mothers who smoke in pregnancy deliver before 37 weeks more often (13.9% against 8.3%), but they also differ in education, age and health. This project estimates how much of the 5.6-point raw gap smoking itself causes, using Augmented Inverse Probability of Treatment Weighting (AIPTW), which is valid if either the smoking model or the birth-outcome model is correct.

**Result.** The adjusted risk of preterm birth is 3.26 percentage points higher among mothers who smoked (bootstrap SE 0.39; 95% interval 2.5 to 4.1), about 34 extra preterm births per 1,000 smoking pregnancies. The estimate holds across alternative adjustment sets, weighting rules and cross-fitted flexible models (LASSO, random forest), and an unmeasured confounder would need a risk ratio of about 2.1 with both smoking and preterm birth to explain it away.

## What has been developed

| Folder | Contents |
|---|---|
| [`1_CSOL_SEASUG_2026_Presentation`](1_CSOL_SEASUG_2026_Presentation) | SEASUG (CSOL) 2026 paper, slides and original code |
| [`2_Conference_Posters`](2_Conference_Posters) | C Day and Analytic Day posters (PowerPoint and PDF) |
| [`3_Journal_Paper`](3_Journal_Paper) | Journal manuscript (PDF) and its source |
| [`code`](code) | The single analysis script behind the posters and the paper |
| [`data`](data) | The two data files the script reads |

## Code

Open [`code/Smoking_Preterm_Full_Analysis.R`](code/Smoking_Preterm_Full_Analysis.R) in R or RStudio, select all and run. The data download from this repository on the first run, and every table and figure is written to `results/`.

## Data

2023 U.S. Natality Public Use File, National Center for Health Statistics. Public, no personal identifiers. See [`data/README.md`](data/README.md).
