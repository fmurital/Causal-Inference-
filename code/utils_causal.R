###############################################################################
## utils_causal.R
## Helper functions shared by every script in this project.
##
## Project : Causal effect of maternal smoking during pregnancy on preterm birth
## Author  : Faruk Muritala (Kennesaw State University)
## Data    : 2023 U.S. NVSS natality public-use file (NCHS)
##
## HOW TO USE:  you never run this file directly.  Each analysis script calls
##   source("utils_causal.R")
## at the top.  Everything here is plain base R plus the data.table package.
##
## WHY THESE FUNCTIONS EXIST
##   The analytic cohort has about 3.3 million births.  R's built-in glm() copies
##   the full design matrix several times and can run out of memory on a laptop.
##   fit_logit() below produces the SAME maximum-likelihood logistic regression
##   as glm(), but it reads the data in chunks of 400,000 rows, so memory use
##   stays small.  (Section "self-test" in 03_main_analysis.R proves it matches
##   glm() to many decimal places.)
##
## SEED: every script that uses randomness calls set.seed(1234).
###############################################################################

suppressPackageStartupMessages(library(data.table))

## ---------------------------------------------------------------------------
## 1. Small math helpers
## ---------------------------------------------------------------------------
logit <- function(p) log(p / (1 - p))          # log-odds
expit <- function(x) 1 / (1 + exp(-x))         # inverse of logit (probability)

## ---------------------------------------------------------------------------
## 2. Chunked logistic regression (exact maximum likelihood, IRLS algorithm)
## ---------------------------------------------------------------------------
## IRLS = iteratively re-weighted least squares, the algorithm glm() itself uses.
## At every iteration we compute, from a "working response" z and weights w,
##      beta_new = (X'WX)^(-1) X'Wz
## X'WX and X'Wz are SUMS over rows, so we can add them up chunk by chunk.
##
## Arguments
##   rhs   : one-sided formula for the design matrix, e.g.  ~ age_grp + educ
##   data  : data.table that contains every variable named in rhs and yvar
##   yvar  : name (character) of the 0/1 response column
##   rows  : optional integer vector; fit only on these rows (e.g. smokers only)
##   chunk : number of rows processed at a time
##   ridge : a tiny number added to the diagonal of X'WX so that a design column
##           that is all zeros inside a subgroup (an empty cell) cannot make the
##           matrix singular.  1e-8 is far too small to change any estimate.
## Returns a list: beta (coefficients), vcov (model-based covariance), rhs,
##   n (rows used), iterations, deviance.
fit_logit <- function(rhs, data, yvar, rows = NULL, chunk = 400000L,
                      maxit = 30L, tol = 1e-9, ridge = 1e-8) {
  if (is.null(rows)) rows <- seq_len(nrow(data))
  chunks   <- split(rows, ceiling(seq_along(rows) / chunk))
  beta     <- NULL
  dev_old  <- Inf
  for (it in seq_len(maxit)) {
    XtWX <- 0; XtWz <- 0; dev <- 0
    for (ix in chunks) {
      d <- data[ix]                              # copy of this chunk only
      X <- model.matrix(rhs, d)                  # design matrix for the chunk
      y <- d[[yvar]]
      if (is.null(beta)) {                       # first pass: crude starting values
        mu  <- (y + 0.5) / 2
        eta <- logit(mu)
      } else {
        eta <- drop(X %*% beta)
        mu  <- expit(eta)
      }
      mu  <- pmin(pmax(mu, 1e-10), 1 - 1e-10)    # keep away from exactly 0 or 1
      w   <- mu * (1 - mu)                       # IRLS weights
      z   <- eta + (y - mu) / w                  # IRLS working response
      XtWX <- XtWX + crossprod(X * sqrt(w))      # accumulate X'WX
      XtWz <- XtWz + crossprod(X, w * z)         # accumulate X'Wz
      dev  <- dev - 2 * sum(y * log(mu) + (1 - y) * log(1 - mu))
    }
    beta <- solve(XtWX + diag(ridge, ncol(XtWX)), XtWz)   # Newton update
    if (abs(dev - dev_old) / (abs(dev) + 0.1) < tol) break   # converged?
    dev_old <- dev
  }
  list(beta = drop(beta), vcov = solve(XtWX + diag(ridge, ncol(XtWX))),
       rhs = rhs, n = length(rows), iterations = it, deviance = dev)
}

## Predicted probability P(Y=1 | X) from a fit_logit() object, for ALL rows of
## `data` (not just the rows the model was fitted on).  Also chunked.
predict_logit <- function(fit, data, chunk = 400000L) {
  n   <- nrow(data)
  out <- numeric(n)
  for (ix in split(seq_len(n), ceiling(seq_len(n) / chunk))) {
    X <- model.matrix(fit$rhs, data[ix])
    out[ix] <- expit(drop(X %*% fit$beta))
  }
  out
}

## ---------------------------------------------------------------------------
## 3. Turning influence-function values into estimates, SEs and 95% CIs
## ---------------------------------------------------------------------------
## Every estimator below can be written as an average of per-person values psi_i:
##       estimate = mean(psi_i)
## and its variance as  var(psi_i - mean(psi)) / n   (the "influence function"
## variance, a.k.a. the sandwich variance).  We apply that to the risk under
## smoking (psi1) and under no smoking (psi0), then use the delta method to get
## the risk difference, risk ratio and odds ratio with confidence intervals.
## `psi1` and `psi0` are vectors with one value per birth.
summarize_if <- function(psi1, psi0, label) {
  n  <- length(psi1)
  m1 <- mean(psi1); m0 <- mean(psi0)
  e1 <- psi1 - m1;  e0 <- psi0 - m0              # centred influence values
  se <- function(phi) sqrt(var(phi) / n)
  rd  <- m1 - m0;              se_rd  <- se(e1 - e0)
  lrr <- log(m1 / m0);         se_lrr <- se(e1 / m1 - e0 / m0)
  lor <- logit(m1) - logit(m0); se_lor <- se(e1 / (m1 * (1 - m1)) - e0 / (m0 * (1 - m0)))
  z <- qnorm(0.975)
  data.table(
    estimator = label,
    risk_smoke = 100 * m1, risk_nosmoke = 100 * m0,       # percent
    RD_pp = 100 * rd, RD_lo = 100 * (rd - z * se_rd), RD_hi = 100 * (rd + z * se_rd),
    RR = exp(lrr), RR_lo = exp(lrr - z * se_lrr), RR_hi = exp(lrr + z * se_lrr),
    OR = exp(lor), OR_lo = exp(lor - z * se_lor), OR_hi = exp(lor + z * se_lor)
  )
}

## ---------------------------------------------------------------------------
## 4. The four effect estimators (all target the ATE = E[Y(1)] - E[Y(0)])
## ---------------------------------------------------------------------------
## t = 1 if mother smoked, y = 1 if preterm, e = estimated propensity score,
## mu1 / mu0 = predicted risk of preterm birth if smoking / not smoking.

## (a) AIPTW  = augmented inverse-probability-of-treatment weighting (doubly robust)
##     psi1_i = mu1_i + t_i (y_i - mu1_i) / e_i
##     psi0_i = mu0_i + (1 - t_i)(y_i - mu0_i) / (1 - e_i)
##   Read psi1 as: "outcome-model prediction, plus a weighted correction for how
##   wrong that prediction is among the people who actually smoked."
est_aiptw <- function(y, t, e, mu1, mu0, label = "AIPTW (doubly robust)") {
  psi1 <- mu1 + t * (y - mu1) / e
  psi0 <- mu0 + (1 - t) * (y - mu0) / (1 - e)
  summarize_if(psi1, psi0, label)
}

## (b) IPTW (Hajek / stabilised form): weighted average outcome in each arm,
##     each person weighted by 1 / P(the treatment they actually got | X).
est_iptw <- function(y, t, e, label = "IPTW (stabilized)") {
  w1 <- t / e;  w0 <- (1 - t) / (1 - e)
  m1 <- sum(w1 * y) / sum(w1)
  m0 <- sum(w0 * y) / sum(w0)
  psi1 <- m1 + w1 * (y - m1) / mean(w1)          # influence values of a Hajek mean
  psi0 <- m0 + w0 * (y - m0) / mean(w0)
  summarize_if(psi1, psi0, label)
}

## (c) Regression standardisation (g-computation): average the outcome model's
##     predictions over EVERYONE, once pretending all smoke, once none smoke.
##     The influence values include the extra term that accounts for the
##     uncertainty in the estimated outcome-model coefficients.
est_gcomp <- function(y, t, mu1, mu0, corr1, corr0, label = "Outcome regression (g-computation)") {
  summarize_if(mu1 + corr1, mu0 + corr0, label)
}

## Parameter-uncertainty correction used by est_gcomp().
## fit  : arm-specific outcome model from fit_logit()
## mu   : its predictions for ALL rows
## arm  : integer row indices of the arm the model was fitted on
## Derivation: m = mean_j mu_j(beta).  A first-order expansion gives
##   m_hat - m = h' (beta_hat - beta),  h = mean_j mu_j (1 - mu_j) X_j,
##   beta_hat - beta ~ V * sum_(i in arm) X_i (y_i - mu_i)
## so each arm member i contributes  n * h' V X_i (y_i - mu_i)  to the influence.
gcomp_correction <- function(fit, data, mu, yvar, arm, chunk = 400000L) {
  n <- nrow(data)
  h <- 0
  for (ix in split(seq_len(n), ceiling(seq_len(n) / chunk))) {
    X <- model.matrix(fit$rhs, data[ix])
    h <- h + crossprod(X, mu[ix] * (1 - mu[ix]))
  }
  h <- h / n
  g <- n * (fit$vcov %*% h)                      # p x 1 vector
  corr <- numeric(n)
  for (ix in split(arm, ceiling(seq_along(arm) / chunk))) {
    X <- model.matrix(fit$rhs, data[ix])
    corr[ix] <- drop(X %*% g) * (data[[yvar]][ix] - mu[ix])
  }
  corr
}

## (d) Propensity-score STRATIFICATION: cut the cohort into PS quintiles, take
##     the smoker/non-smoker difference inside each stratum, then average the
##     strata weighted by their size.
est_strata <- function(y, t, e, k = 5, label = "PS-quintile stratification") {
  brk <- unique(quantile(e, probs = seq(0, 1, length.out = k + 1)))
  s   <- cut(e, breaks = brk, include.lowest = TRUE, labels = FALSE)
  n   <- length(y)
  m1 <- 0; m0 <- 0; v1 <- 0; v0 <- 0
  tab <- vector("list", max(s))
  for (j in seq_len(max(s))) {
    ii <- s == j
    n1 <- sum(t[ii] == 1); n0 <- sum(t[ii] == 0); w <- sum(ii) / n
    p1 <- mean(y[ii & t == 1]); p0 <- mean(y[ii & t == 0])
    m1 <- m1 + w * p1;  m0 <- m0 + w * p0
    v1 <- v1 + w^2 * p1 * (1 - p1) / n1
    v0 <- v0 + w^2 * p0 * (1 - p0) / n0
    ## 2x2 cell counts as DOUBLE (not integer) to avoid the 32-bit overflow that
    ## makes base R's mantelhaen.test() return NA for very large strata.
    tab[[j]] <- c(a = as.double(sum(ii & t == 1 & y == 1)), b = as.double(sum(ii & t == 1 & y == 0)),
                  c = as.double(sum(ii & t == 0 & y == 1)), d = as.double(sum(ii & t == 0 & y == 0)))
  }
  z <- qnorm(0.975)
  rd <- m1 - m0; se_rd <- sqrt(v1 + v0)
  lor <- logit(m1) - logit(m0)
  se_lor <- sqrt(v1 / (m1 * (1 - m1))^2 + v0 / (m0 * (1 - m0))^2)
  lrr <- log(m1 / m0)
  se_lrr <- sqrt(v1 / m1^2 + v0 / m0^2)
  ## Mantel-Haenszel pooled (conditional) odds ratio with the Robins-Breslow-
  ## Greenland variance of log(OR_MH).
  num <- den <- PR <- PS <- QR <- QS <- 0
  for (tb in tab) {
    nn <- sum(tb); R <- tb["a"] * tb["d"] / nn; S <- tb["b"] * tb["c"] / nn
    P <- (tb["a"] + tb["d"]) / nn; Q <- (tb["b"] + tb["c"]) / nn
    num <- num + R; den <- den + S
    PR <- PR + P * R; PS <- PS + P * S + Q * R; QS <- QS + Q * S
  }
  or_mh <- unname(num / den)
  se_mh <- unname(sqrt(PR / (2 * num^2) + PS / (2 * num * den) + QS / (2 * den^2)))
  list(
    marginal = data.table(estimator = label,
      risk_smoke = 100 * m1, risk_nosmoke = 100 * m0,
      RD_pp = 100 * rd, RD_lo = 100 * (rd - z * se_rd), RD_hi = 100 * (rd + z * se_rd),
      RR = exp(lrr), RR_lo = exp(lrr - z * se_lrr), RR_hi = exp(lrr + z * se_lrr),
      OR = exp(lor), OR_lo = exp(lor - z * se_lor), OR_hi = exp(lor + z * se_lor)),
    mantel_haenszel = data.table(estimator = "PS-quintile Mantel-Haenszel (conditional OR)",
      OR = or_mh, OR_lo = exp(log(or_mh) - z * se_mh), OR_hi = exp(log(or_mh) + z * se_mh)),
    strata_table = data.table(stratum = seq_along(tab), do.call(rbind, tab))
  )
}

## (e) AIPTW for the ATT (effect among mothers who actually smoked).  Secondary
##     estimand.  psi_i = [ t(y - mu0) - (1-t) e/(1-e) (y - mu0) ] / p,
##     p = P(T = 1);  ATT = mean(psi) and its standard error follows from the
##     centred values.
est_aiptw_att <- function(y, t, e, mu0, label = "AIPTW (ATT, effect among smokers)") {
  p   <- mean(t)
  num <- t * (y - mu0) - (1 - t) * (e / (1 - e)) * (y - mu0)
  att <- mean(num) / p
  phi <- (num - t * att) / p
  se  <- sqrt(var(phi) / length(y))
  z <- qnorm(0.975)
  data.table(estimator = label, ATT_RD_pp = 100 * att,
             lo = 100 * (att - z * se), hi = 100 * (att + z * se))
}


## (f) Generic weighted-mean estimator: any pair of weights (w1 for smokers,
##     w0 for non-smokers).  Used for overlap weights (ATO) below.
est_wmean <- function(y, t, w1, w0, label) {
  m1 <- sum(w1 * y) / sum(w1)
  m0 <- sum(w0 * y) / sum(w0)
  summarize_if(m1 + w1 * (y - m1) / mean(w1), m0 + w0 * (y - m0) / mean(w0), label)
}

## ---------------------------------------------------------------------------
## 5. Propensity-score truncation and the E-value
## ---------------------------------------------------------------------------
## Clamp (do NOT drop) propensity scores to the [lo, hi] quantiles of e.
## Clamping keeps every mother in the analysis, so the estimand stays the ATE
## for the whole cohort, while stopping a handful of extreme scores from
## producing huge weights.  lo = 0, hi = 1 means "no truncation".
clamp_ps <- function(e, lo = 0.01, hi = 0.99) {
  q <- quantile(e, c(lo, hi))
  pmin(pmax(e, q[1]), q[2])
}

## E-value (VanderWeele & Ding 2017): the minimum strength of association, on
## the risk-ratio scale, that an UNMEASURED confounder would need with BOTH
## smoking and preterm birth to fully explain away the observed risk ratio.
evalue_rr <- function(rr) {
  rr <- ifelse(rr < 1, 1 / rr, rr)
  rr + sqrt(rr * (rr - 1))
}

## ---------------------------------------------------------------------------
## 6. Weighted standardized mean difference (covariate balance)
## ---------------------------------------------------------------------------
## For a 0/1 covariate x: SMD = (p1 - p0) / sqrt((p1(1-p1) + p0(1-p0)) / 2),
## where p1, p0 are the (weighted) proportions among smokers / non-smokers.
smd_binary <- function(x, t, w = rep(1, length(x))) {
  p1 <- sum(w[t == 1] * x[t == 1]) / sum(w[t == 1])
  p0 <- sum(w[t == 0] * x[t == 0]) / sum(w[t == 0])
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
}
