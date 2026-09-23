suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
options(width=120)
cc <- readRDS("cc.rds")

## Full analytic cohort is 3.25M rows; fitting glm() on all of it exhausts this
## container's memory (confirmed: OOM-killed after ~41s). An 11-covariate logistic
## model is precisely estimable from a much smaller sample, so we fit the PS model
## on a random 15% subsample (stratified by exposure to preserve prevalence) and
## then SCORE the full cohort with the fitted coefficients (prediction is cheap;
## no refitting needed). This is a standard, fully documented efficiency step.
n_full <- nrow(cc)
frac <- 0.15
idx1 <- which(cc$smoker == 1); idx0 <- which(cc$smoker == 0)
set.seed(79001)
samp_idx <- c(sample(idx1, size = round(frac*length(idx1))),
              sample(idx0, size = round(frac*length(idx0))))
cc_fit <- cc[samp_idx, ]
cat("PS-model fitting subsample:", nrow(cc_fit), "of", n_full,
    sprintf("(%.1f%% , seed=79001, stratified by exposure)\n", 100*nrow(cc_fit)/n_full))
cat("  treated in subsample:", sum(cc_fit$smoker==1), " control in subsample:", sum(cc_fit$smoker==0), "\n")

t0 <- Sys.time()
ps_formula <- smoker ~ mager + mrace6 + meduc_f + precare_f + wic01 + parity0 +
                        ppterm01 + pdiab01 + phype01 + bmi + pay_f
ps_fit <- glm(ps_formula, data = cc_fit, family = binomial())
cat("PS model fit time:", round(as.numeric(Sys.time()-t0, units="secs"),1), "sec\n")
rm(cc_fit); gc(FALSE)

cc$ps <- predict(ps_fit, newdata = cc, type = "response")
cc$logit_ps <- predict(ps_fit, newdata = cc, type = "link")

cat("\n=== Propensity score model coefficients ===\n")
print(round(summary(ps_fit)$coefficients, 4))

cat("\nPS range treated:", round(range(cc$ps[cc$smoker==1]),4), "\n")
cat("PS range control:", round(range(cc$ps[cc$smoker==0]),4), "\n")

## overlap / positivity plot
p_overlap <- ggplot(cc, aes(x = ps, fill = factor(smoker, labels=c("Non-smoker","Smoker")))) +
  geom_density(alpha = 0.5, adjust = 1.2) +
  labs(x = "Estimated propensity score  P(smoker=1 | X)", y = "Density", fill = "Group",
       title = "Overlap of propensity scores by smoking status (unmatched, n = 3,251,950)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave("fig_overlap_unmatched.png", p_overlap, width = 7.2, height = 4.4, dpi = 220)

## simple empirical positivity check: proportion of each covariate stratum combo with both groups present
cat("\n% of common-support region (ps between 0.01 and 0.99):",
    round(100*mean(cc$ps > 0.01 & cc$ps < 0.99), 2), "%\n")

saveRDS(cc, "cc_ps.rds")
saveRDS(ps_fit, "ps_fit.rds")
cat("\nStep 2 (propensity score model + overlap) complete.\n")
