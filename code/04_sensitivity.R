###############################################################################
## 04_sensitivity.R
## STEP N: sensitivity and robustness analyses.
##
## Question answered: "Does the AIPTW result change if we change the analyst's
## choices?"  We vary (A) how extreme propensity scores are handled, (B) which
## covariates and which model forms are used, and (C) we look at subgroups.
##
## Input  : output/analytic_cohort.rds
## Output : output/sens_truncation.csv, sens_specifications.csv, subgroup_aiptw.csv
##          (the figures are drawn from these files by 04b_sensitivity_figures.R)
## HOW TO RUN:  Rscript 04_sensitivity.R      (slow: about 45-60 minutes on a laptop with 8 GB RAM,
##              because specification D has about 150 columns; the other parts take about 10 minutes)
###############################################################################

set.seed(1234)
source("pipeline_functions.R")
d <- readRDS("output/analytic_cohort.rds")
n <- nrow(d); y <- d$preterm; t <- d$smoker
fmt <- function(x) x[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
options(width = 200)

## ---------------------------------------------------------------------------
## A. Handling of extreme propensity scores (main specification)
## ---------------------------------------------------------------------------
cat("== A. Propensity-score truncation, main specification ==\n")
fit <- run_spec(d, specs$main$rhs)
rows <- list()
for (cut in list(c(0, 1), c(0.001, 0.999), c(0.01, 0.99), c(0.05, 0.95))) {
  lab <- if (cut[1] == 0) "No truncation" else sprintf("Clamp PS at %g / %g percentile", 100 * cut[1], 100 * cut[2])
  r <- estimate_all(d, fit, cut[1], cut[2])$table
  r <- r[estimator %in% c("IPTW (stabilized)", "AIPTW (doubly robust)")]
  r[, scenario := lab]
  e_used <- estimate_all(d, fit, cut[1], cut[2])$e_used
  w <- ifelse(t == 1, mean(t) / e_used, (1 - mean(t)) / (1 - e_used))
  r[, max_weight := max(w)][, smoker_ess := sum(w[t == 1])^2 / sum(w[t == 1]^2)]
  rows[[length(rows) + 1]] <- r
}
## Dropping (rather than clamping) births with PS outside [0.01, 0.99]:
## this CHANGES the target population, so it is shown only for comparison.
keep <- fit$e >= 0.01 & fit$e <= 0.99
ra <- est_aiptw(y[keep], t[keep], fit$e[keep], fit$mu1[keep], fit$mu0[keep])
ra[, `:=`(scenario = sprintf("Delete births with PS outside [0.01, 0.99] (keeps %.1f%% of cohort)", 100 * mean(keep)),
          max_weight = NA_real_, smoker_ess = NA_real_)]
rows[[length(rows) + 1]] <- ra
## Overlap weights (ATO): target population = mothers with the most comparable
## smokers and non-smokers; bounded weights, so no truncation is needed.
ato <- est_wmean(y, t, t * (1 - fit$e), (1 - t) * fit$e, "Overlap-weighted (ATO)")
ato[, `:=`(scenario = "Overlap weights (target = ATO)", max_weight = NA_real_, smoker_ess = NA_real_)]
rows[[length(rows) + 1]] <- ato
trunc_tab <- rbindlist(rows, fill = TRUE)
setcolorder(trunc_tab, c("scenario", "estimator"))
print(fmt(trunc_tab[, .(scenario, estimator, RD_pp, RD_lo, RD_hi, OR, OR_lo, OR_hi, max_weight, smoker_ess)]))
fwrite(trunc_tab, "output/sens_truncation.csv")

## ---------------------------------------------------------------------------
## B. Alternative covariate sets and model forms (primary truncation)
## ---------------------------------------------------------------------------
cat("\n== B. Alternative specifications ==\n")
spec_rows <- list()
for (nm in names(specs)) {
  cat(" ", specs[[nm]]$label, "\n")
  f <- if (nm == "main") fit else run_spec(d, specs[[nm]]$rhs)
  r <- estimate_all(d, f, 0.01, 0.99)
  tt <- r$table[estimator %in% c("Outcome regression (g-computation)", "IPTW (stabilized)", "AIPTW (doubly robust)")]
  tt <- cbind(specification = specs[[nm]]$label, tt)
  tt[, ATT_AIPTW_pp := r$att$ATT_RD_pp]
  spec_rows[[length(spec_rows) + 1]] <- tt
  if (nm != "main") rm(f); gc(FALSE)
}
spec_tab <- rbindlist(spec_rows)
print(fmt(spec_tab[, .(specification, estimator, RD_pp, RD_lo, RD_hi, OR, OR_lo, OR_hi)]))
fwrite(spec_tab, "output/sens_specifications.csv")

## ---------------------------------------------------------------------------
## C. Subgroup (effect-modification) analysis with the main-specification AIPTW
## ---------------------------------------------------------------------------
## AIPTW is an average of per-mother values (psi), so a subgroup estimate is just
## the average of those values over the subgroup.  Exploratory, not confirmatory.
cat("\n== C. AIPTW by subgroup (exploratory) ==\n")
e1 <- clamp_ps(fit$e, 0.01, 0.99)
psi1 <- fit$mu1 + t * (y - fit$mu1) / e1
psi0 <- fit$mu0 + (1 - t) * (y - fit$mu0) / (1 - e1)
sg <- list()
for (v in c("age_grp", "race", "educ", "bmi_grp", "parity", "ppterm")) {
  for (lv in levels(d[[v]])) {
    ii <- which(d[[v]] == lv)
    if (sum(t[ii] == 1) < 500) next                  # skip tiny cells
    s <- summarize_if(psi1[ii], psi0[ii], lv)
    sg[[length(sg) + 1]] <- cbind(variable = v, level = lv, n = length(ii), n_smokers = sum(t[ii]), s[, -1])
  }
}
sg <- rbindlist(sg)
print(fmt(sg[, .(variable, level, n, n_smokers, risk_nosmoke, RD_pp, RD_lo, RD_hi, OR)]))
fwrite(sg, "output/subgroup_aiptw.csv")

cat("\nDone. Figures are drawn by 04b_sensitivity_figures.R\n")
