# Updates Since CSOL

This note records what was added or changed after the SEASUG (CSOL) 2026 presentation, one item at a time, with the numbers and the place in the code where each one is produced. The CSOL materials are unchanged in [`../1_CSOL_SEASUG_2026_Presentation`](../1_CSOL_SEASUG_2026_Presentation).

All numbers below come from the full run of `code/Smoking_Preterm_Full_Analysis.R` (seed 1234). The tables are in `results/tables`.

## Summary

| # | Concern raised | What was done | Result | Script section |
|---|---|---|---|---|
| 1 | Adjustment set: prenatal-care timing and WIC are measured during pregnancy | Main ten-variable model compared with a pre-pregnancy-only set and other sets | 3.26 against 3.65 points; all sets between 3.0 and 3.7 | 5 |
| 2 | Overlap and residual imbalance | Effect among smokers, overlap weights, weight and effective-sample-size table, balance for every covariate level | 41% of non-smokers have a smoking score below 0.01; levels out of balance 21, 3, 0 | 6 |
| 3 | Check the influence-function variance and strengthen the bootstrap | Full-sample bootstrap of 1,000 resamples with every model refitted | Bootstrap SE 0.388 against 0.399; ratio 0.97; 95% interval 2.55 to 4.08 | 7 |
| 4 | Flexible models, possibly cross-fitted | 5-fold cross-fitted AIPTW with logistic, LASSO (all two-way interactions) and random forest, two random splits | 3.22, 2.80, 2.30; average of the three 2.84 | 8 |
| 5 | Novelty and public-health meaning | Effect among smokers expressed in births | about 34 extra preterm births per 1,000 smoking pregnancies | 12 |
| 6 | Intensity, timing and cessation | Quit and dose comparisons from cigarettes before pregnancy; trimester of smoking from the full natality file | Continued against quit 3.08 points; quitters against never-smokers -0.07; trimester results in the journal paper | 9, 10 |

## What changed in the analysis

**Variance (concern 3).** The 3.26 million births collapse into 65,787 cells of identical confounder values. A logistic regression on the cell counts gives the same coefficients as one on all births, and resampling births is the same as drawing the cell counts from a multinomial distribution. This makes 1,000 full-sample bootstrap resamples practical. At CSOL the bootstrap used a 200,000-birth subsample with about 6,000 smokers.

**Flexible models (concern 4).** Each birth is scored by models that did not see it. The births are split at random into five folds. The models are fitted on four folds and used to predict the fifth, once for each fold, and the five sets of predictions are stacked before the AIPTW estimate is computed. The whole procedure is repeated with a second random split. Cross-fitting itself changes little (logistic 3.22 against 3.26 without it). The random forest gives a lower estimate mostly because its weights are more concentrated.

**Overlap (concern 2).** The effect among smokers (3.39 points) is reported beside the whole-population effect, since it needs only non-smokers who resemble smokers. Overlap weights give 3.30. The largest standardized mean difference falls from 0.64 before weighting to 0.13 after weighting and 0.00 with overlap weights.

**Adjustment set (concern 1).** Dropping prenatal-care timing and WIC moves the estimate from 3.26 to 3.65 points. The ranking of confounders shows prenatal-care timing changes the estimate by +0.64 points when left out and WIC by -0.13.

**Who smokes and what drives the estimate.** A Cramér's V association matrix shows education and WIC are strongly linked (0.42). The ranking of confounders shows education best separates smokers from non-smokers (46% of the smoking model's fit) and a previous preterm birth best predicts preterm birth (38% of that model's fit). Leaving education out would move the estimate by +1.03 points.

**Quitting and dose (concern 6).** Preterm birth was 8.3% in never-smokers, 9.4% in mothers who quit before pregnancy, 12.5% in mothers who started during pregnancy and 13.9% in those who smoked throughout. After adjustment, mothers who quit look like never-smokers (-0.07 points, 95% CI -0.40 to 0.27). Continuing against quitting differs by 3.08 points. The extra risk is 3.55, 3.48, 3.27 and 3.39 points at 1 to 5, 6 to 10, 11 to 20 and 21 or more cigarettes per day before pregnancy. These results are exploratory because intensity during pregnancy and its timing are not in the extract.

**Sensitivity analysis.** An unmeasured factor would need a risk ratio of about 2.12 with both smoking and preterm birth to explain the estimate away (lower confidence limit 1.92).

## What changed in the code

| At CSOL | Now |
|---|---|
| Ten scripts (`00` to `07`) run in order with `run_all.sh` | One script, `Smoking_Preterm_Full_Analysis.R`, with a table of contents and a data dictionary |
| Data extract inside `code/data` | Data in the repository `data` folder, downloaded automatically on the first run |
| Standard errors from the influence function; bootstrap on a subsample | Both reported; bootstrap on the full sample |
| Logistic regression only | Logistic, LASSO and random forest with cross-fitting |
| Quitting and dose not analysed | Sections 9 and 10 |

## What changed in the presentation

The C Day and Analytic Day posters add, in this order: a problem statement, a methodology diagram that follows the data through the two models to the estimate and its checks, nine numbered figures with short notes inside each, a results section that quotes the numbers, and a limitations section. The poster contains no trimester results.

## Trimester of smoking (Section 10)

The full national file records cigarettes per day in each trimester. Rebuilding the cohort from it gives the same 3,259,496 births, and 3,218,258 have known cigarettes in every period. Preterm risk was 7.94% in mothers who never smoked, 10.73% for the first trimester only, 13.18% for the first and second trimesters, 13.93% for all three, and 10.75% for those who started after the first trimester. The effect among those who smoked was 1.38 points (0.86 to 1.90) for the first trimester only, 3.43 (2.55 to 4.31) for the first and second, and 3.68 (3.43 to 3.93) for all three. Among mothers who smoked throughout, it rose with third-trimester intensity: 3.26, 3.72, 4.44 and 5.03 points at 1 to 5, 6 to 10, 11 to 20 and 21 or more cigarettes per day. The 2,288 mothers who started after the first trimester give 0.87 points (-0.36 to 2.11), too imprecise to read. The time a mother can smoke is limited by when she delivers, so these are patterns of reported smoking and not the effect of stopping at a chosen time. The results are exploratory, are in the journal paper (Appendix Table A7), and are not on the posters.
