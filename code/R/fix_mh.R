suppressPackageStartupMessages(library(dplyr))
cc_w <- readRDS("cc_weighted.rds")
strata <- levels(cc_w$ps_strata)
num <- 0; den <- 0
## Robins-Breslow-Greenland variance components
PR <- 0; PS <- 0; QR <- 0; QS <- 0
for (s in strata) {
  d <- cc_w[cc_w$ps_strata==s,]
  a <- as.double(sum(d$smoker==1 & d$preterm==1))
  b <- as.double(sum(d$smoker==1 & d$preterm==0))
  cc_ <- as.double(sum(d$smoker==0 & d$preterm==1))
  dd <- as.double(sum(d$smoker==0 & d$preterm==0))
  n <- a+b+cc_+dd
  num <- num + a*dd/n
  den <- den + b*cc_/n
  R <- a*dd/n; S <- b*cc_/n
  P <- (a+dd)/n; Q <- (b+cc_)/n
  PR <- PR + P*R; PS <- PS + P*S + Q*R; QS <- QS + Q*S
}
mh_or <- num/den
var_log_mh <- PR/(2*num^2) + PS/(2*num*den) + QS/(2*den^2)
se_log <- sqrt(var_log_mh)
log_or <- log(mh_or)
ci <- exp(c(log_or-1.96*se_log, log_or+1.96*se_log))
cat(sprintf("MH pooled OR = %.4f  (95%% CI %.4f - %.4f)\n", mh_or, ci[1], ci[2]))

res <- read.csv("outcome_results.csv")
res[res$model=="5. PS-quintile stratified (Mantel-Haenszel, full cohort)",
    c("log_or","OR","CI_low","CI_high")] <- c(log_or, mh_or, ci[1], ci[2])
write.csv(res, "outcome_results.csv", row.names=FALSE)
print(res[,c("model","OR","CI_low","CI_high")], digits=4)
