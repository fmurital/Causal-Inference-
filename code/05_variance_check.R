###############################################################################
## 05_variance_check.R
## STEP N: does the fast "influence-function" standard error agree with a
## brute-force bootstrap?
##
## Why this matters: with 3.3 million births we cannot afford to bootstrap the
## whole analysis, so the main results use closed-form (influence-function)
## standard errors.  Here we check those formulas on a smaller random sample
## (n = 200,000) where a full bootstrap IS affordable.  If both methods give
## similar standard errors on the small sample, the formulas can be trusted on
## the large one.
##
## Input  : output/analytic_cohort.rds
## Output : output/variance_check.csv
## HOW TO RUN:  Rscript 05_variance_check.R     (about 30-40 minutes: 200 bootstrap fits)
###############################################################################

set.seed(1234)
source("pipeline_functions.R")
d_all <- readRDS("output/analytic_cohort.rds")
N_SUB <- 200000L          # size of the random sample
B     <- 200L             # number of bootstrap replications

sub <- d_all[sample(nrow(d_all), N_SUB)]        # simple random sample, seed 1234
rm(d_all)
cat("Subsample:", N_SUB, "births,", sum(sub$smoker), "smokers\n")

## One full pass of the pipeline on data frame `dd`; returns three risk differences
## (in percentage points): AIPTW, IPTW, and outcome regression (g-computation).
one_pass <- function(dd) {
  rhs <- specs$main$rhs
  e   <- predict_logit(fit_logit(rhs, dd, "smoker", chunk = 1e6L), dd, chunk = 1e6L)
  f1  <- fit_logit(rhs, dd, "preterm", rows = which(dd$smoker == 1), chunk = 1e6L)
  f0  <- fit_logit(rhs, dd, "preterm", rows = which(dd$smoker == 0), chunk = 1e6L)
  mu1 <- predict_logit(f1, dd, chunk = 1e6L); mu0 <- predict_logit(f0, dd, chunk = 1e6L)
  e   <- pmin(pmax(clamp_ps(e, 0.01, 0.99), 1e-6), 1 - 1e-6)
  y <- dd$preterm; t <- dd$smoker
  c(aiptw = est_aiptw(y, t, e, mu1, mu0)$RD_pp,
    iptw  = est_iptw(y, t, e)$RD_pp,
    gcomp = 100 * (mean(mu1) - mean(mu0)))
}

## (1) Closed-form influence-function standard errors on the subsample
rhs <- specs$main$rhs
smk <- which(sub$smoker == 1); non <- which(sub$smoker == 0)
e   <- predict_logit(fit_logit(rhs, sub, "smoker", chunk = 1e6L), sub, chunk = 1e6L)
f1  <- fit_logit(rhs, sub, "preterm", rows = smk, chunk = 1e6L)
f0  <- fit_logit(rhs, sub, "preterm", rows = non, chunk = 1e6L)
mu1 <- predict_logit(f1, sub, chunk = 1e6L); mu0 <- predict_logit(f0, sub, chunk = 1e6L)
c1  <- gcomp_correction(f1, sub, mu1, "preterm", smk, chunk = 1e6L)
c0  <- gcomp_correction(f0, sub, mu0, "preterm", non, chunk = 1e6L)
ec  <- pmin(pmax(clamp_ps(e, 0.01, 0.99), 1e-6), 1 - 1e-6)
y <- sub$preterm; t <- sub$smoker
if_tab <- rbind(est_aiptw(y, t, ec, mu1, mu0), est_iptw(y, t, ec),
                est_gcomp(y, t, mu1, mu0, c1, c0))
if_se <- (if_tab$RD_hi - if_tab$RD_lo) / (2 * qnorm(0.975))   # back out the SE in percentage points

## (2) Bootstrap: resample births with replacement and repeat the whole pipeline
cat("Running", B, "bootstrap replications...\n")
boot <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("aiptw", "iptw", "gcomp")))
tic <- Sys.time()
for (b in seq_len(B)) {
  boot[b, ] <- one_pass(sub[sample(N_SUB, replace = TRUE)])
  if (b %% 20 == 0) cat(sprintf("  replication %d of %d (%.1f min elapsed)\n", b, B,
                                as.numeric(Sys.time() - tic, units = "mins")))
}

out <- data.table(estimator = c("AIPTW", "IPTW", "Outcome regression"),
                  RD_pp_point = if_tab$RD_pp,
                  SE_influence_function = if_se,
                  SE_bootstrap = apply(boot, 2, sd),
                  ratio_IF_to_bootstrap = if_se / apply(boot, 2, sd))
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(out, "output/variance_check.csv")
cat("\nDone.\n")
