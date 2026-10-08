###############################################################################
## pipeline_functions.R
## The model specifications and the one function (run_spec) that carries a
## specification from "model formula" to "all effect estimates".
## Sourced by 03_main_analysis.R and 04_sensitivity.R.  Not run on its own.
###############################################################################
source("utils_causal.R")

## ---------------------------------------------------------------------------
## Model specifications = the right-hand side of the regression formulas.
## ---------------------------------------------------------------------------
## MAIN: the ten confounders of the final DAG.  Age and BMI enter as the
## standard clinical categories so the model is flexible (no straight-line
## assumption) in those two variables.
rhs_main <- ~ age_grp + race + educ + precare + wic + parity + ppterm + pdiab + phype + bmi_grp

specs <- list(
  main       = list(label = "Main: ten DAG confounders", rhs = rhs_main),
  ## Sensitivity A: drop the two variables measured DURING pregnancy
  ## (prenatal-care timing and WIC), keeping only strictly pre-pregnancy factors.
  prepreg    = list(label = "A: pre-pregnancy confounders only (no prenatal-care timing, no WIC)",
                    rhs = ~ age_grp + race + educ + parity + ppterm + pdiab + phype + bmi_grp),
  ## Sensitivity B: add insurance (the covariate removed from the DAG).
  insurance  = list(label = "B: main set + insurance type",
                    rhs = update(rhs_main, ~ . + pay)),
  ## Sensitivity C: age and BMI as smooth quadratic curves instead of categories.
  continuous = list(label = "C: age and BMI as continuous (quadratic)",
                    rhs = ~ mager + I(mager^2) + race + educ + precare + wic + parity + ppterm +
                           pdiab + phype + bmi + I(bmi^2)),
  ## Sensitivity D: richer models with two-way interactions among the strongest
  ## predictors, to check that conclusions do not depend on main-effects-only models.
  rich       = list(label = "D: main set + two-way interactions",
                    rhs = update(rhs_main, ~ . + race:educ + age_grp:educ + age_grp:parity +
                                   educ:wic + educ:precare))
)

## ---------------------------------------------------------------------------
## run_spec(): fit the three nuisance models for one specification.
##   1. propensity-score model   P(smoker = 1 | X)
##   2. outcome model among smokers      P(preterm = 1 | X, smoker = 1)  -> mu1
##   3. outcome model among non-smokers  P(preterm = 1 | X, smoker = 0)  -> mu0
## and return the predictions for EVERY mother in the cohort.
## ---------------------------------------------------------------------------
run_spec <- function(d, rhs, verbose = TRUE) {
  tic <- Sys.time()
  smk <- which(d$smoker == 1); non <- which(d$smoker == 0)
  ps_fit <- fit_logit(rhs, d, "smoker")
  e      <- predict_logit(ps_fit, d)
  o1_fit <- fit_logit(rhs, d, "preterm", rows = smk)
  o0_fit <- fit_logit(rhs, d, "preterm", rows = non)
  mu1 <- predict_logit(o1_fit, d)
  mu0 <- predict_logit(o0_fit, d)
  if (verbose) cat(sprintf("   fitted 3 models in %.0f s (iterations: PS %d, outcome|smoke %d, outcome|no smoke %d)\n",
                           as.numeric(Sys.time() - tic, units = "secs"),
                           ps_fit$iterations, o1_fit$iterations, o0_fit$iterations))
  list(e = e, mu1 = mu1, mu0 = mu0, ps_fit = ps_fit, o1_fit = o1_fit, o0_fit = o0_fit,
       corr1 = gcomp_correction(o1_fit, d, mu1, "preterm", smk),
       corr0 = gcomp_correction(o0_fit, d, mu0, "preterm", non))
}

## ---------------------------------------------------------------------------
## estimate_all(): given fitted nuisance models, compute every effect estimator
## for one choice of propensity-score truncation (lo, hi quantiles).
## ---------------------------------------------------------------------------
estimate_all <- function(d, fit, lo = 0.01, hi = 0.99, label_suffix = "") {
  y <- d$preterm; t <- d$smoker
  e <- if (lo == 0 && hi == 1) fit$e else clamp_ps(fit$e, lo, hi)
  e <- pmin(pmax(e, 1e-6), 1 - 1e-6)           # numerical guard only
  st <- est_strata(y, t, fit$e)                # strata use the raw PS ranks
  list(
    table = rbind(
      est_iptw(y, t, rep(mean(t), length(y)), label = "Crude (unadjusted)"),
      est_gcomp(y, t, fit$mu1, fit$mu0, fit$corr1, fit$corr0),
      est_iptw(y, t, e),
      st$marginal,
      est_aiptw(y, t, e, fit$mu1, fit$mu0)
    ),
    mh = st$mantel_haenszel,
    strata = st$strata_table,
    att = est_aiptw_att(y, t, e, fit$mu0),
    e_used = e
  )
}
