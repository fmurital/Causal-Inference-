suppressPackageStartupMessages({library(dplyr); library(sandwich); library(lmtest)})
options(width=120)
cc_all <- readRDS("cc_ps.rds")
cc_w   <- readRDS("cc_weighted.rds")
flushcat <- function(...) { cat(...); flush(stdout()) }

## All multivariable (many-covariate) model fits are capped at a stratified random
## subsample (same rationale/logic as the propensity-score model: an 11-covariate
## logistic model does not need millions of rows to be estimated precisely, and
## unrestricted glm() on the full 3.25M/1.9M-row frames repeatedly exhausted this
## container's memory). Simple 1-covariate / tabulation models use the full cohort.
subsample <- function(d, max_n = 400000, seed = 79003) {
  if (nrow(d) <= max_n) return(d)
  set.seed(seed)
  i1 <- which(d$smoker==1); i0 <- which(d$smoker==0)
  frac <- max_n/nrow(d)
  idx <- c(sample(i1, round(frac*length(i1))), sample(i0, round(frac*length(i0))))
  d[idx, ]
}

or_ci <- function(fit, term, label, robust=TRUE) {
  if (robust) {
    co <- coeftest(fit, vcov. = sandwich(fit))
  } else co <- coeftest(fit)
  b <- co[term,"Estimate"]; se <- co[term,"Std. Error"]
  data.frame(model=label, log_or=b, se=se, OR=exp(b), CI_low=exp(b-1.96*se), CI_high=exp(b+1.96*se))
}

results <- list()

flushcat("Fitting model 1 (crude, full cohort n=", nrow(cc_all), ")...\n")
f1 <- glm(preterm ~ smoker, data = cc_all, family = binomial())
results[["1"]] <- or_ci(f1, "smoker", "1. Crude (unadjusted, full cohort)")
write.csv(bind_rows(results), "outcome_results.csv", row.names=FALSE)
flushcat("  done. OR=", round(results[["1"]]$OR,3), "\n")

flushcat("Fitting model 2 (adjusted, subsample)...\n")
d2 <- subsample(cc_all)
flushcat("  subsample n=", nrow(d2), " treated=", sum(d2$smoker==1), "\n")
f2 <- glm(preterm ~ smoker + mager + mrace6 + meduc_f + precare_f + wic01 + parity0 +
                     ppterm01 + pdiab01 + phype01 + bmi + pay_f,
          data = d2, family = binomial())
results[["2"]] <- or_ci(f2, "smoker", "2. Multivariable-adjusted (DAG set, subsample)")
write.csv(bind_rows(results), "outcome_results.csv", row.names=FALSE)
flushcat("  done. OR=", round(results[["2"]]$OR,3), "\n")
rm(d2, f2); gc(FALSE)

flushcat("Fitting model 3 (IPTW only, subsample of trimmed cohort)...\n")
d3 <- subsample(cc_w)
f3 <- glm(preterm ~ smoker, data = d3, family = quasibinomial(), weights = iptw)
results[["3"]] <- or_ci(f3, "smoker", "3. IPTW (trimmed, stabilized weights, subsample)")
write.csv(bind_rows(results), "outcome_results.csv", row.names=FALSE)
flushcat("  done. OR=", round(results[["3"]]$OR,3), "\n")

flushcat("Fitting model 4 (IPTW + adjusted, subsample)...\n")
f4 <- glm(preterm ~ smoker + mager + mrace6 + meduc_f + precare_f + wic01 + parity0 +
                     ppterm01 + pdiab01 + phype01 + bmi + pay_f,
          data = d3, family = quasibinomial(), weights = iptw)
results[["4"]] <- or_ci(f4, "smoker", "4. IPTW + covariate-adjusted (doubly robust, subsample)")
write.csv(bind_rows(results), "outcome_results.csv", row.names=FALSE)
flushcat("  done. OR=", round(results[["4"]]$OR,3), "\n")
rm(d3, f3, f4); gc(FALSE)

flushcat("Computing model 5 (PS-quintile Mantel-Haenszel, full trimmed cohort)...\n")
mh_tabs <- lapply(levels(cc_w$ps_strata), function(s) {
  d <- cc_w[cc_w$ps_strata==s,]
  a <- sum(d$smoker==1 & d$preterm==1); b <- sum(d$smoker==1 & d$preterm==0)
  cc_ <- sum(d$smoker==0 & d$preterm==1); dd <- sum(d$smoker==0 & d$preterm==0)
  matrix(c(a,cc_,b,dd), nrow=2)
})
mh <- mantelhaen.test(array(unlist(mh_tabs), dim=c(2,2,5)))
results[["5"]] <- data.frame(model="5. PS-quintile stratified (Mantel-Haenszel, full cohort)",
    log_or=log(unname(mh$estimate)), se=NA,
    OR=unname(mh$estimate), CI_low=mh$conf.int[1], CI_high=mh$conf.int[2])
write.csv(bind_rows(results), "outcome_results.csv", row.names=FALSE)
flushcat("  done. OR=", round(results[["5"]]$OR,3), "\n")

res <- bind_rows(results)
cat("\n=== Odds ratio for preterm birth: smokers vs. non-smokers ===\n")
print(res[,c("model","OR","CI_low","CI_high")], digits=4)

rd_crude <- 100*(mean(cc_all$preterm[cc_all$smoker==1]) - mean(cc_all$preterm[cc_all$smoker==0]))
rd_w1 <- weighted.mean(cc_w$preterm[cc_w$smoker==1], cc_w$iptw[cc_w$smoker==1])
rd_w0 <- weighted.mean(cc_w$preterm[cc_w$smoker==0], cc_w$iptw[cc_w$smoker==0])
flushcat(sprintf("\nCrude risk: smoker %.2f%% vs non-smoker %.2f%% (diff %.2f pts)\n",
            100*mean(cc_all$preterm[cc_all$smoker==1]), 100*mean(cc_all$preterm[cc_all$smoker==0]), rd_crude))
flushcat(sprintf("IPTW-weighted risk: smoker %.2f%% vs non-smoker %.2f%% (diff %.2f pts)\n",
            100*rd_w1, 100*rd_w0, 100*(rd_w1-rd_w0)))

sink("outcome_summary.txt")
cat("n_full_cohort=", nrow(cc_all), "\n")
cat("n_trimmed_cohort=", nrow(cc_w), "\n")
cat("rd_crude_pp=", round(rd_crude,3), "\n")
cat("rd_iptw_pp=", round(100*(rd_w1-rd_w0),3), "\n")
cat("risk_smoker_crude_pct=", round(100*mean(cc_all$preterm[cc_all$smoker==1]),3), "\n")
cat("risk_nonsmoker_crude_pct=", round(100*mean(cc_all$preterm[cc_all$smoker==0]),3), "\n")
cat("risk_smoker_iptw_pct=", round(100*rd_w1,3), "\n")
cat("risk_nonsmoker_iptw_pct=", round(100*rd_w0,3), "\n")
sink()
flushcat("\nStep 4 (outcome models) complete.\n")
