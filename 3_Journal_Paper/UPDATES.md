# Updates for the Journal Version

This note lists what the journal draft adds to the SEASUG (CSOL) 2026 paper and what remains.

## Added or changed

| Area | CSOL paper | Journal draft |
|---|---|---|
| Standard error | Influence function; bootstrap on a 200,000-birth subsample | Influence function and a full-sample bootstrap of 1,000 resamples (SE 0.388 against 0.399; percentile interval 2.55 to 4.08), reported as the primary interval |
| Models | Logistic regression | Logistic, LASSO with all two-way interactions and random forest, each with 5-fold cross-fitting and two random splits |
| Target populations | Whole population, effect among smokers | Adds overlap weights, a weight and effective-sample-size table, and balance for every covariate level |
| Adjustment set | Ten confounders | Ten confounders compared with a pre-pregnancy-only set (3.65 against 3.26 points) |
| Smoking history | Not analysed | Quit, started and continued groups; continued against quit 3.08 points; quitters against never-smokers -0.07 points; cigarettes per day before pregnancy |
| Introduction | Problem and method | Adds what is known, three limits of existing work, and what this study adds |
| Public-health reading | Effect in points | Effect among smokers in births per 1,000 pregnancies (about 34) |
| Code | Ten scripts | One script with a data dictionary and a table of contents |

## Still open

1. **Trimester of smoking.** Needs the national natality file; the script and the manuscript section are ready (see `../data/README.md`). The manuscript marks the section in red.
2. **Journal choice and formatting.** The draft follows a general article layout. It will be moved into the journal's template once the journal is chosen.
3. **Preterm subtypes** (extremely, very and late preterm) and a **STROBE/RECORD checklist** are not yet included.
4. **Smoking misclassification.** Self-reported smoking is discussed qualitatively; a quantitative bias analysis is not yet done.
