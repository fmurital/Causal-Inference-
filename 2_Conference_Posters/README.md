# Conference Posters

Two posters were prepared after the SEASUG (CSOL) 2026 presentation. This folder holds both, the single analysis script behind them, and the result tables from the latest full run. The section below lists what changed since CSOL.

| Poster | Folder | Format |
|---|---|---|
| C Day 2026 | [`C_Day_2026`](C_Day_2026) | Native PowerPoint with an editable flow diagram, and PDF |
| Analytic Day 2026 | [`Analytic_Day_2026`](Analytic_Day_2026) | Analytic Day template, two layouts, PowerPoint and PDF |

Both posters carry the same text, figures and numbers. Only the layout differs.

## Headline numbers

| Quantity | Value |
|---|---|
| Births in the analysis | 3,259,496 (97,187 to mothers who smoked) |
| Raw gap in preterm birth | 13.9% against 8.3%, 5.6 points |
| Adjusted risk difference (AIPTW) | 3.26 points, bootstrap SE 0.39, 95% interval 2.5 to 4.1 |
| Extra preterm births among mothers who smoked | about 34 per 1,000 pregnancies |
| Range across analytic choices | 3.0 to 3.7 points |
| E-value | 2.12 |

## What changed since CSOL

At CSOL the analysis was a ten-script pipeline with a logistic-regression AIPTW estimate, an influence-function interval and a bootstrap check on a 200,000-birth subsample. After the conference the committee comments were worked through one at a time. The full list with the numbers is in [`UPDATES_SINCE_CSOL.md`](UPDATES_SINCE_CSOL.md). In short:

1. **Variance.** The standard error is now checked with a full-sample bootstrap of 1,000 resamples in which every model is refitted. Bootstrap SE 0.388 against 0.399 from the influence function.
2. **Flexible models.** AIPTW is repeated with 5-fold cross-fitting and three learners: logistic regression, LASSO with all two-way interactions, and random forest. Estimates: 3.22, 2.80 and 2.30 points.
3. **Overlap and balance.** The effect among smokers, overlap weights, weight diagnostics and a balance table for every covariate level. Levels out of balance fall from 21 to 3 after weighting and 0 with overlap weights; 41% of non-smokers have almost no smoker like them.
4. **Adjustment set.** The ten-variable model is compared with a pre-pregnancy-only set (3.65 points), because prenatal-care timing and WIC are measured during pregnancy.
5. **Quitting and dose.** Four smoking-history groups (never, quit before pregnancy, started in pregnancy, smoked throughout) and cigarettes per day before pregnancy. Quitters look like never-smokers after adjustment (-0.07 points); the extra risk is similar at every cigarette level.
6. **Who is compared.** An association matrix and a ranking of which confounders drive smoking and preterm birth, with the change in the estimate when each is left out.
7. **Public-health translation.** The effect among smokers expressed as births per 1,000 pregnancies.
8. **Code.** One script, `Smoking_Preterm_Full_Analysis.R`, now produces every table and figure, with seed 1234, and reads its data from this repository.

Trimester of smoking: run after the posters were finalized, in Section 10 of the script, from the full national file (see [`../data/README.md`](../data/README.md)). The tables are `results/tables/trimester_crude.csv` and `trimester_results.csv`, and the results are in the journal paper. They are not on either poster.

## Folder contents

| Path | Contents |
|---|---|
| `C_Day_2026/` | Poster, editable flow diagram, README with the poster's section-by-section map |
| `Analytic_Day_2026/` | Poster in two layouts, README |
| `code/` | `Smoking_Preterm_Full_Analysis.R` and instructions to run it |
| `results/tables/` | The 26 result tables from the full run (CSV) |
| `results/figures/` | Figures written by the script |
| `UPDATES_SINCE_CSOL.md` | Comment-by-comment record of the changes |
