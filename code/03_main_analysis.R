###############################################################################
## 03_main_analysis.R
## STEP N: the main causal analysis (AIPTW, IPTW, g-computation, PS strata).
##
## Input  : output/analytic_cohort.rds   (made by 01_prepare_data.R)
## Output : output/results_main.csv, ps_model_odds_ratios.csv, balance_table.csv,
##          weight_summary.csv, evalues.csv, figures/*.png
##
## HOW TO RUN:   Rscript 03_main_analysis.R      (about 5-10 minutes)
## Packages  :   data.table, ggplot2
###############################################################################

set.seed(1234)
suppressPackageStartupMessages(library(ggplot2))
source("pipeline_functions.R")
dir.create("figures", showWarnings = FALSE)

NAVY <- "#1F3864"; ORANGE <- "#C1611A"; GRAY <- "#6B7F99"
d <- readRDS("output/analytic_cohort.rds")
n <- nrow(d); y <- d$preterm; t <- d$smoker
cat("Analytic cohort:", format(n, big.mark = ","), "births\n\n")

## ---------------------------------------------------------------------------
## SELF-TEST: prove that our chunked logistic regression equals R's glm().
## We use a random 100,000-row sample (seed 1234) so glm() runs quickly.
## ---------------------------------------------------------------------------
cat("== Self-test: fit_logit() versus glm() on a 100,000-row sample ==\n")
idx   <- sample(n, 1e5)
dsub  <- d[idx]
g_ref <- glm(update(rhs_main, smoker ~ .), data = dsub, family = binomial())
g_new <- fit_logit(rhs_main, dsub, "smoker", chunk = 30000L)
cat("   largest coefficient difference between glm() and fit_logit():",
    format(max(abs(coef(g_ref) - g_new$beta)), digits = 3), "\n\n")
stopifnot(max(abs(coef(g_ref) - g_new$beta)) < 1e-5)
rm(dsub, g_ref, g_new)

## ---------------------------------------------------------------------------
## 1. Fit the three nuisance models for the MAIN specification
## ---------------------------------------------------------------------------
cat("== Fitting propensity-score and outcome models (main specification) ==\n")
fit <- run_spec(d, specs$main$rhs)

## Propensity-score model coefficients, shown as odds ratios for smoking
ps_or <- data.table(term = names(fit$ps_fit$beta), coef = fit$ps_fit$beta,
                    se = sqrt(diag(fit$ps_fit$vcov)))
ps_or[, `:=`(OR = exp(coef), OR_lo = exp(coef - 1.96 * se), OR_hi = exp(coef + 1.96 * se))]
fwrite(ps_or, "output/ps_model_odds_ratios.csv")

## ---------------------------------------------------------------------------
## 2. Positivity: how extreme are the estimated propensity scores?
## ---------------------------------------------------------------------------
cat("\n== Propensity-score distribution (positivity check) ==\n")
qs <- quantile(fit$e, c(0, .001, .01, .05, .25, .5, .75, .95, .99, .999, 1))
print(round(qs, 5))
cat("Smokers:     median PS", round(median(fit$e[t == 1]), 4), " range", signif(range(fit$e[t == 1]), 3), "\n")
cat("Non-smokers: median PS", round(median(fit$e[t == 0]), 4), " range", signif(range(fit$e[t == 0]), 3), "\n")
cat("Share of smokers with PS < 0.01:", round(100 * mean(fit$e[t == 1] < 0.01), 2), "%\n")
cat("Share of non-smokers with PS < 0.01:", round(100 * mean(fit$e[t == 0] < 0.01), 2), "%\n")

## ---------------------------------------------------------------------------
## 3. All estimators, primary truncation = clamp PS at its 1st and 99th percentiles
## ---------------------------------------------------------------------------
cat("\n== Effect estimates (PS clamped at 1st/99th percentile) ==\n")
res <- estimate_all(d, fit, 0.01, 0.99)
tab <- res$table

## Conditional odds ratio from an ordinary adjusted logistic regression
## (the "textbook" model: preterm ~ smoker + confounders).  Reported separately
## because it is a CONDITIONAL odds ratio, not the marginal (population) one.
adj <- fit_logit(update(rhs_main, ~ smoker + .), d, "preterm")
b <- adj$beta["smoker"]; s <- sqrt(adj$vcov["smoker", "smoker"])
cond_or <- data.table(estimator = "Adjusted logistic regression (conditional OR)",
                      OR = exp(b), OR_lo = exp(b - 1.96 * s), OR_hi = exp(b + 1.96 * s))

tab[, estimator := factor(estimator, levels = estimator)]
options(width = 200, datatable.print.nrows = 50)
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\nConditional OR (adjusted logistic regression):\n"); print(round(cond_or[, -1], 3))
cat("\nMantel-Haenszel conditional OR across PS quintiles:\n"); print(res$mh[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\nAIPTW effect among smokers (ATT, risk difference in percentage points):\n"); print(res$att[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\nPS-quintile strata (2x2 counts: a=smoker&preterm, b=smoker&term, c=nonsmoker&preterm, d=nonsmoker&term):\n"); print(res$strata)

fwrite(tab, "output/results_main.csv")
fwrite(rbind(cond_or, res$mh, fill = TRUE), "output/results_conditional_or.csv")
fwrite(res$att, "output/results_att.csv")
fwrite(res$strata, "output/ps_strata_counts.csv")

## ---------------------------------------------------------------------------
## 4. Weights: size, stability and effective sample size
## ---------------------------------------------------------------------------
p   <- mean(t)
e   <- res$e_used
w   <- ifelse(t == 1, p / e, (1 - p) / (1 - e))          # stabilized IPTW weights
ess <- function(w) sum(w)^2 / sum(w^2)                   # effective sample size
wsum <- data.table(
  group = c("Smokers", "Non-smokers"),
  n = c(sum(t == 1), sum(t == 0)),
  min = c(min(w[t == 1]), min(w[t == 0])), median = c(median(w[t == 1]), median(w[t == 0])),
  mean = c(mean(w[t == 1]), mean(w[t == 0])), max = c(max(w[t == 1]), max(w[t == 0])),
  effective_n = c(ess(w[t == 1]), ess(w[t == 0])))
cat("\n== Stabilized IPTW weights (after clamping) ==\n"); print(wsum[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(wsum, "output/weight_summary.csv")

## ---------------------------------------------------------------------------
## 5. Covariate balance before and after IPTW (standardized mean differences)
## ---------------------------------------------------------------------------
cat("\n== Covariate balance ==\n")
bal <- list()
for (v in c("age_grp", "race", "educ", "precare", "parity", "bmi_grp")) {
  for (lv in levels(d[[v]])) {
    x <- as.numeric(d[[v]] == lv)
    bal[[length(bal) + 1]] <- data.table(covariate = paste0(v, ": ", lv),
      smd_before = smd_binary(x, t), smd_after = smd_binary(x, t, w))
  }
}
for (v in c("wic", "ppterm", "pdiab", "phype")) {
  x <- as.numeric(d[[v]] == "Yes")
  bal[[length(bal) + 1]] <- data.table(covariate = paste0(v, ": Yes"),
    smd_before = smd_binary(x, t), smd_after = smd_binary(x, t, w))
}
bal <- rbindlist(bal)
bal[, `:=`(abs_before = abs(smd_before), abs_after = abs(smd_after))]
fwrite(bal, "output/balance_table.csv")
cat("Largest |SMD| before weighting:", round(max(bal$abs_before), 3), "(", bal$covariate[which.max(bal$abs_before)], ")\n")
cat("Largest |SMD| after  weighting:", round(max(bal$abs_after), 3),  "(", bal$covariate[which.max(bal$abs_after)], ")\n")
cat("Covariate levels with |SMD| > 0.10 before:", sum(bal$abs_before > 0.1), "of", nrow(bal), "; after:", sum(bal$abs_after > 0.1), "\n")

## ---------------------------------------------------------------------------
## 6. E-values (how strong would an unmeasured confounder have to be?)
## ---------------------------------------------------------------------------
aip <- tab[grepl("^AIPTW", estimator)]
ev <- data.table(estimate = c("AIPTW risk ratio", "AIPTW risk ratio, lower 95% limit"),
                 value = c(aip$RR, aip$RR_lo),
                 E_value = c(evalue_rr(aip$RR), evalue_rr(aip$RR_lo)))
cat("\n== E-values ==\n"); print(ev[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(ev, "output/evalues.csv")

## ---------------------------------------------------------------------------
## 7. Figures
## ---------------------------------------------------------------------------
theme_proj <- theme_minimal(base_size = 16) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),   # figure titles are given in the paper/slide captions
        plot.title.position = "plot",
        legend.position = "top", panel.grid.minor = element_blank())

## 7a. Propensity-score overlap (log scale because the exposure is rare)
set.seed(1234)
ix <- c(which(t == 1), sample(which(t == 0), 150000))        # all smokers + a random 150k non-smokers (for speed)
pd <- data.table(ps = fit$e[ix], group = factor(ifelse(t[ix] == 1, "Smoked in pregnancy", "Did not smoke"),
                                                levels = c("Did not smoke", "Smoked in pregnancy")))
q <- quantile(fit$e, c(.01, .99))
p1 <- ggplot(pd, aes(ps, fill = group, colour = group)) +
  geom_density(alpha = 0.45, adjust = 1.3, linewidth = 0.6) +
  scale_x_log10(breaks = c(.001, .003, .01, .03, .1, .3), labels = c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3")) +
  geom_vline(xintercept = q, linetype = "dashed", colour = GRAY) +
  scale_fill_manual(values = c(GRAY, ORANGE), name = NULL) + scale_colour_manual(values = c(GRAY, ORANGE), name = NULL) +
  labs(x = "Estimated propensity score (log scale)", y = "Density",
       title = "Overlap is limited for non-smokers",
       subtitle = "Dashed lines: 1st and 99th percentiles used to clamp the score") + theme_proj
ggsave("figures/fig_ps_overlap.png", p1, width = 8, height = 4.8, dpi = 200)

## 7b. Love plot: balance before and after weighting
bl <- melt(bal[, .(covariate, Before = abs_before, `After IPTW` = abs_after)], id.vars = "covariate",
           variable.name = "Sample", value.name = "absSMD")
bl[, covariate := factor(covariate, levels = bal[order(abs_before), covariate])]
p2 <- ggplot(bl, aes(absSMD, covariate, colour = Sample)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = GRAY) +
  geom_point(size = 2.3) +
  scale_colour_manual(values = c(Before = ORANGE, `After IPTW` = NAVY), name = NULL) +
  labs(x = "Absolute standardized mean difference", y = NULL,
       title = "Balance before and after weighting", subtitle = "Dashed line = 0.10 imbalance threshold") +
  theme_proj + theme(axis.text.y = element_text(size = 10.5))
ggsave("figures/fig_love_plot.png", p2, width = 8, height = 7, dpi = 200)

## 7c. Forest plot of the estimators (marginal risk difference and odds ratio)
fp <- copy(tab)
fp[, estimator := factor(estimator, levels = rev(estimator))]
p3 <- ggplot(fp, aes(RD_pp, estimator)) +
  geom_vline(xintercept = 0, colour = GRAY) +
  geom_errorbarh(aes(xmin = RD_lo, xmax = RD_hi), height = 0.2, colour = NAVY) +
  geom_point(aes(colour = grepl("AIPTW", estimator)), size = 3.2) +
  scale_colour_manual(values = c(`FALSE` = NAVY, `TRUE` = ORANGE), guide = "none") +
  labs(x = "Risk difference in preterm birth (percentage points)", y = NULL,
       title = "Adjustment lowers the crude risk difference by one-third to two-fifths") + theme_proj +
  theme(plot.margin = margin(6, 24, 6, 6))   # room for the right end of the axis title
ggsave("figures/fig_forest_rd.png", p3, width = 9, height = 4.2, dpi = 200)
p4 <- ggplot(fp, aes(OR, estimator)) +
  geom_vline(xintercept = 1, colour = GRAY) +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.2, colour = NAVY) +
  geom_point(aes(colour = grepl("AIPTW", estimator)), size = 3.2) +
  scale_colour_manual(values = c(`FALSE` = NAVY, `TRUE` = ORANGE), guide = "none") +
  labs(x = "Marginal odds ratio for preterm birth", y = NULL,
       title = "Odds ratio by estimator") + theme_proj
ggsave("figures/fig_forest_or.png", p4, width = 9, height = 4.2, dpi = 200)

saveRDS(list(e = fit$e, mu1 = fit$mu1, mu0 = fit$mu0), "output/nuisance_predictions_main.rds")
cat("\nDone. Results saved in output/ and figures/\n")
