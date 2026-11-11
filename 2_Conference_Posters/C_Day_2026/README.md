# C Day 2026 Poster

**Smoking During Pregnancy and Preterm Birth: A Doubly Robust Propensity Score Weighting Analysis of Millions of U.S. Births**

Author: Faruk Muritala. Advisor: Dhrubajyoti Ghosh, Ph.D., Kennesaw State University.

| File | Description |
|---|---|
| `Poster_CDay2026_Muritala_version_7.pptx` | The poster in PowerPoint. The methodology diagram is built from native shapes and connectors, so every box, node and arrow can be edited. |
| `Poster_CDay2026_Muritala_version_7.pdf` | The same poster as PDF. |
| `Flow_Diagram_Editable_v7.pptx` | The methodology diagram on its own, for reuse. |

## How the poster reads

Left column: Problem Statement, Methods (two models, cross-fitting, uncertainty, checks, with the estimator equation), Results and Findings, Limitations.
Right panel: the headline result, the methodology diagram (data, ten confounders, two working models, estimator, output, checks) and nine figures.

| Figure | Shows |
|---|---|
| 1 | Causal diagram: each confounder affects both smoking and preterm birth |
| 2 | Percentage of mothers who smoked, by education |
| 3 | Births at each week of delivery, smokers and non-smokers |
| 4 | Association between each pair of variables |
| 5 | What drives confounding: share of each model's fit lost when a variable is dropped, and the change in the estimate |
| 6 | Overlap of estimated smoking probabilities and covariate balance |
| 7 | Risk difference under alternative analytic choices |
| 8 | Flexible cross-fitted models and variance |
| 9 | Preterm births by smoking history, and extra risk by cigarettes per day |

## Code and data

The numbers come from [`../code/Smoking_Preterm_Full_Analysis.R`](../code/Smoking_Preterm_Full_Analysis.R); the result tables are in [`../results/tables`](../results/tables). What changed since CSOL is in [`../UPDATES_SINCE_CSOL.md`](../UPDATES_SINCE_CSOL.md). The poster contains no trimester results.
