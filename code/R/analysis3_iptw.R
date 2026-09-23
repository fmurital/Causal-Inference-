suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr)})
options(width=120)
cc <- readRDS("cc_ps.rds")

cat("Common-support (ps in .01-.99) among TREATED:",
    round(100*mean(cc$ps[cc$smoker==1] > .01 & cc$ps[cc$smoker==1] < .99),2), "%\n")
cat("Common-support (ps in .01-.99) among CONTROL:",
    round(100*mean(cc$ps[cc$smoker==0] > .01 & cc$ps[cc$smoker==0] < .99),2), "%\n")

## trim extreme PS (standard positivity safeguard) before weighting
cc <- cc %>% filter(ps > 0.01, ps < 0.99)
cat("After trimming ps to [0.01,0.99]:", nrow(cc), "rows (",
    sum(cc$smoker==1), "treated /", sum(cc$smoker==0), "control )\n")

## ---- stabilized IPTW weights ----
p_treat <- mean(cc$smoker)
cc$iptw <- ifelse(cc$smoker==1, p_treat/cc$ps, (1-p_treat)/(1-cc$ps))
cat("\nStabilized IPTW weight summary:\n"); print(summary(cc$iptw))
cat("Weights > 10:", sum(cc$iptw>10), " (", round(100*mean(cc$iptw>10),3), "% )\n")
## cap at 1st/99th pctile as robustness (Cole & Hernan practice)
qs <- quantile(cc$iptw, c(.01,.99))
cc$iptw_trim <- pmin(pmax(cc$iptw, qs[1]), qs[2])

## ---- PS quintile strata (secondary approach) ----
cc$ps_strata <- cut(cc$ps, breaks = quantile(cc$ps, probs = seq(0,1,0.2)),
                     include.lowest = TRUE, labels = paste0("Q",1:5))
cat("\nStrata sizes:\n"); print(table(cc$ps_strata, cc$smoker))

saveRDS(cc, "cc_weighted.rds")

## ---- weighted balance table (SMD using IPTW weights) ----
wtd_mean <- function(x, w) sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE)
wtd_var  <- function(x, w) { m<-wtd_mean(x,w); sum(w*(x-m)^2, na.rm=TRUE)/sum(w, na.rm=TRUE) }
smd_w <- function(x, g, w) {
  m1<-wtd_mean(x[g==1], w[g==1]); m0<-wtd_mean(x[g==0], w[g==0])
  s1<-wtd_var(x[g==1], w[g==1]);  s0<-wtd_var(x[g==0], w[g==0])
  (m1-m0)/sqrt((s1+s0)/2)
}
smd_u <- function(x, g) {
  m1<-mean(x[g==1],na.rm=TRUE); m0<-mean(x[g==0],na.rm=TRUE)
  s1<-var(x[g==1],na.rm=TRUE);  s0<-var(x[g==0],na.rm=TRUE)
  (m1-m0)/sqrt((s1+s0)/2)
}
num_covs <- c("mager","bmi","wic01","parity0","ppterm01","pdiab01","phype01")
fac_covs <- c("mrace6","meduc_f","precare_f","pay_f")
rows <- list()
for (v in num_covs) {
  rows[[length(rows)+1]] <- data.frame(covariate=v,
      smd_before=smd_u(cc[[v]], cc$smoker),
      smd_after =smd_w(cc[[v]], cc$smoker, cc$iptw))
}
for (v in fac_covs) for (lv in levels(cc[[v]])) {
  xb <- as.numeric(cc[[v]]==lv)
  rows[[length(rows)+1]] <- data.frame(covariate=paste0(v,":",lv),
      smd_before=smd_u(xb, cc$smoker),
      smd_after =smd_w(xb, cc$smoker, cc$iptw))
}
bal <- bind_rows(rows); bal$abs_before<-abs(bal$smd_before); bal$abs_after<-abs(bal$smd_after)
cat("\n=== IPTW-weighted balance table ===\n"); print(bal, digits=3)
write.csv(bal, "balance_table_iptw.csv", row.names=FALSE)
cat("\nMax |SMD| before:", round(max(bal$abs_before),4), " after weighting:", round(max(bal$abs_after),4), "\n")
cat("N covariates |SMD|>0.1 before:", sum(bal$abs_before>0.1), " after:", sum(bal$abs_after>0.1), "\n")

bal_long <- bal %>% select(covariate, Before=abs_before, After=abs_after) %>%
  pivot_longer(c(Before,After), names_to="Sample", values_to="AbsSMD")
bal_long$Sample <- factor(bal_long$Sample, levels=c("Before","After"))
ord <- bal %>% arrange(abs_before) %>% pull(covariate)
bal_long$covariate <- factor(bal_long$covariate, levels=ord)
p_love <- ggplot(bal_long, aes(x=AbsSMD,y=covariate,color=Sample)) +
  geom_point(size=2.2) + geom_vline(xintercept=0.1, linetype="dashed", color="grey40") +
  scale_color_manual(values=c(Before="#d95f02", After="#1b9e77")) +
  labs(x="Absolute standardized mean difference", y=NULL,
       title="Covariate balance before vs. after IPTW weighting") +
  theme_minimal(base_size=11) + theme(legend.position="top")
ggsave("fig_love_plot.png", p_love, width=7.5, height=6, dpi=220)

p_overlap_w <- ggplot(cc, aes(x=ps, fill=factor(smoker, labels=c("Non-smoker","Smoker")))) +
  geom_density(alpha=0.5, adjust=1.2) +
  labs(x="Estimated propensity score", y="Density", fill="Group",
       title=sprintf("Propensity score overlap after trimming (n = %s)", format(nrow(cc), big.mark=","))) +
  theme_minimal(base_size=13) + theme(legend.position="top")
ggsave("fig_overlap_trimmed.png", p_overlap_w, width=7.2, height=4.4, dpi=220)

cat("\nStep 3 (IPTW weighting + stratification + balance) complete.\n")
