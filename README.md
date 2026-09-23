# Causal Inference — STAT 7900 / SEASUG 2026

Reproducible causal-inference work from Faruk Muritala's STAT 7900 (Kennesaw State University, Fall 2026), including a full propensity-score analysis of maternal smoking and preterm birth built on the 2023 NVSS natality data, submitted as SEASUG 2026 Paper #103.

## Contents

- **`paper/`** — the SEASUG 2026 conference paper / STAT 7900 final project (`SEASUG2026_Paper103_Muritala.docx`), including the full write-up, DAG, results, discussion, AI-use disclosure, and complete R + SAS code appendices.
- **`code/R/`** — the executed R analysis pipeline, run in this order:
  1. `analysis.R` — load, recode, and restrict the NVSS extract to the complete-case analytic cohort
  2. `analysis2_ps.R` — propensity-score model (logistic regression, exposure-stratified subsample) + overlap diagnostics
  3. `analysis3_iptw.R` — positivity trimming, stabilized IPTW weights, PS-quintile strata, covariate-balance diagnostics
  4. `analysis3_match.R` — the abandoned 1:1 nearest-neighbor matching attempt, kept for transparency (see paper §5.2 — it never completed; no numbers from this script appear anywhere in the results)
  5. `analysis4_outcomes.R` — the five outcome models (crude, DAG-adjusted, IPTW, doubly robust, PS-quintile Mantel–Haenszel)
  6. `fix_mh.R` — manual double-precision Mantel–Haenszel calculation (works around a base-R integer-overflow bug on large stratified tables)
  7. `make_dag.py`, `make_forest.py` — the DAG and forest-plot figures (Python/matplotlib)
- **`code/SAS/`** — `nvss_smoking_preterm.sas`, a SAS translation of the same pipeline (PROC LOGISTIC, PROC SURVEYLOGISTIC, PROC FREQ/CMH). **Not yet executed or verified** — no SAS environment was available when this was written; run it and check its output against `code/R` before relying on it.
- **`figures/`** — the DAG, propensity-score overlap plots, covariate-balance love plot, and results forest plot used in the paper and slides.
- **`funding/`** — draft SEASUG 2026 travel funding request (`SEASUG2026_Funding_Request_DRAFT.docx`), mirroring the FTC 2026 request template; several fields still need confirmation (see the orange callouts in the document).
- **`assignment1/`** — STAT 7900 Assignment 1 solutions (potential outcomes, DAGs, backdoor paths, colliders, simulation, Simpson's paradox), included for course-record completeness.

## Slide deck and poster

- **SEASUG 2026 slide deck** ("Applying Propensity Score Matching to Evaluate Interventions"): https://claude.ai/artifact/AyBa1KWEQAaeZRRyYbCUC5
- **KSU Analytics Day poster**: in progress — link to be added here once published.

## Data

This project uses the 2023 U.S. National Vital Statistics System (NVSS) natality public-use file (NCHS/CDC, distributed via NBER: `data.nber.org/nvss/natality/csv/2023/natality2023us.zip`). The raw file (~2GB unzipped) and the working extract are **not** included in this repository because of size; `code/R/analysis.R` documents the exact filtering and recoding used to go from the raw file to the analytic cohort (N = 3,251,950).

## Reproducibility

Every number in `paper/SEASUG2026_Paper103_Muritala.docx` was produced by executing the scripts in `code/R/` against the NVSS extract, in the order listed above. The SAS code in `code/SAS/` mirrors the same pipeline but has not itself been executed — see the note in the paper's SAS appendix. No result anywhere in this repository is estimated, assumed, or hand-calculated.

## AI-use disclosure

Generative AI assistance (Claude, Anthropic) was used for code drafting, debugging, document formatting, and literature-search support throughout this project, under the author's direction and review, consistent with the STAT 7900 course AI-use policy. See `paper/SEASUG2026_Paper103_Muritala.docx` §9 for the full disclosure statement.
