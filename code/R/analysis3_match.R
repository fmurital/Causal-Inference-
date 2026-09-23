suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(MatchIt)})
options(width=120)
cc <- readRDS("cc_ps.rds")

## Overlap by group before matching (relevant to positivity)
cat("Common-support (ps in .01-.99) among TREATED:",
    round(100*mean(cc$ps[cc$smoker==1] > .01 & cc$ps[cc$smoker==1] < .99),2), "%\n")
cat("Common-support (ps in .01-.99) among CONTROL:",
    round(100*mean(cc$ps[cc$smoker==0] > .01 & cc$ps[cc$smoker==0] < .99),2), "%\n")

caliper_sd <- 0.2 * sd(cc$logit_ps)
cat("Caliper (0.2 * SD of logit PS):", round(caliper_sd, 4), "\n")

## Matching all 96,698 treated against the full 3.15M-row control pool with
## MatchIt's formula interface did not finish in 10 minutes (model-frame
## construction + greedy search over 3.25M rows). An 8:1 control pool gives the
## matching algorithm ample choice per treated unit while making the search
## tractable; controls are drawn at random (same seed) before matching, so the
## discarded controls are missing completely at random with respect to the
## matching step.
set.seed(79002)
ctrl_idx <- which(cc$smoker == 0)
pool_n <- min(length(ctrl_idx), 8 * sum(cc$smoker == 1))
ctrl_pool <- sample(ctrl_idx, pool_n)
cc_m <- cc[c(which(cc$smoker == 1), ctrl_pool), ]
cat("Matching input:", nrow(cc_m), "rows (", sum(cc_m$smoker==1), "treated /",
    sum(cc_m$smoker==0), "candidate controls, 8:1 random pool, seed=79002 )\n")

t0 <- Sys.time()
m.out <- matchit(smoker ~ mager + mrace6 + meduc_f + precare_f + wic01 + parity0 +
                           ppterm01 + pdiab01 + phype01 + bmi + pay_f,
                  data = cc_m, method = "nearest", distance = cc_m$ps,
                  caliper = 0.2, std.caliper = TRUE, ratio = 1, replace = FALSE)
cat("Matching time:", round(as.numeric(Sys.time()-t0, units="secs"),1), "sec\n")
print(summary(m.out, un = FALSE)$nn)

md <- match.data(m.out)
cat("\nMatched sample size:", nrow(md), " (", sum(md$smoker==1), "treated /", sum(md$smoker==0), "control )\n")
saveRDS(md, "matched_data.rds")
saveRDS(m.out, "matchit_obj.rds")

## ---- balance table: standardized mean differences before/after matching ----
smd <- function(x, g) {
  m1 <- mean(x[g==1], na.rm=TRUE); m0 <- mean(x[g==0], na.rm=TRUE)
  s1 <- var(x[g==1], na.rm=TRUE);  s0 <- var(x[g==0], na.rm=TRUE)
  (m1 - m0) / sqrt((s1+s0)/2)
}
num_covs <- c("mager","bmi","wic01","parity0","ppterm01","pdiab01","phype01")
fac_covs <- c("mrace6","meduc_f","precare_f","pay_f")

bal_rows <- list()
for (v in num_covs) {
  bal_rows[[length(bal_rows)+1]] <- data.frame(
    covariate = v,
    smd_before = smd(cc[[v]], cc$smoker),
    smd_after  = smd(md[[v]], md$smoker))
}
for (v in fac_covs) {
  for (lv in levels(cc[[v]])) {
    xb <- as.numeric(cc[[v]] == lv); xa <- as.numeric(md[[v]] == lv)
    bal_rows[[length(bal_rows)+1]] <- data.frame(
      covariate = paste0(v, ":", lv),
      smd_before = smd(xb, cc$smoker),
      smd_after  = smd(xa, md$smoker))
  }
}
bal <- bind_rows(bal_rows)
bal$abs_before <- abs(bal$smd_before)
bal$abs_after  <- abs(bal$smd_after)
cat("\n=== Balance table (Standardized Mean Difference) ===\n")
print(bal, digits=3)
write.csv(bal, "balance_table.csv", row.names = FALSE)
cat("\nMax |SMD| before:", round(max(bal$abs_before),4), " after:", round(max(bal$abs_after),4), "\n")
cat("N covariates with |SMD|>0.1 before:", sum(bal$abs_before>0.1), " after:", sum(bal$abs_after>0.1), "\n")

## love plot
bal_long <- bal %>% select(covariate, Before=abs_before, After=abs_after) %>%
  tidyr::pivot_longer(c(Before, After), names_to="Sample", values_to="AbsSMD")
bal_long$Sample <- factor(bal_long$Sample, levels=c("Before","After"))
ord <- bal %>% arrange(abs_before) %>% pull(covariate)
bal_long$covariate <- factor(bal_long$covariate, levels=ord)
p_love <- ggplot(bal_long, aes(x=AbsSMD, y=covariate, color=Sample)) +
  geom_point(size=2.2) +
  geom_vline(xintercept=0.1, linetype="dashed", color="grey40") +
  scale_color_manual(values=c(Before="#d95f02", After="#1b9e77")) +
  labs(x="Absolute standardized mean difference", y=NULL,
       title="Covariate balance before vs. after 1:1 nearest-neighbor PS matching") +
  theme_minimal(base_size = 11) + theme(legend.position="top")
ggsave("fig_love_plot.png", p_love, width=7.5, height=6, dpi=220)

## overlap after matching
p_overlap_m <- ggplot(md, aes(x=ps, fill=factor(smoker, labels=c("Non-smoker","Smoker")))) +
  geom_density(alpha=0.5, adjust=1.2) +
  labs(x="Estimated propensity score", y="Density", fill="Group",
       title=sprintf("Propensity score overlap after matching (n = %s matched pairs)",
                      format(sum(md$smoker==1), big.mark=","))) +
  theme_minimal(base_size=13) + theme(legend.position="top")
ggsave("fig_overlap_matched.png", p_overlap_m, width=7.2, height=4.4, dpi=220)

cat("\nStep 3 (matching + balance) complete.\n")
