###############################################################################
## 07_paper_numbers.R
## Extra descriptive numbers quoted in the paper and slides (overlap and weights).
##
## Input  : output/analytic_cohort.rds, output/nuisance_predictions_main.rds
## Output : output/paper_numbers.csv
## HOW TO RUN:  Rscript 07_paper_numbers.R      (about 1 minute)
###############################################################################

set.seed(1234)
suppressPackageStartupMessages(library(data.table))
source("utils_causal.R")
d   <- readRDS("output/analytic_cohort.rds")
nu  <- readRDS("output/nuisance_predictions_main.rds")
e   <- nu$e; t <- d$smoker; y <- d$preterm

q01 <- quantile(e, 0.01); q99 <- quantile(e, 0.99)
ec  <- pmin(pmax(e, q01), q99)                         # the clamped propensity score used in the main analysis
w   <- ifelse(t == 1, mean(t) / ec, (1 - mean(t)) / (1 - ec))   # stabilized ATE weights

out <- data.table(
  quantity = c("PS 1st percentile (clamp lower bound)", "PS 99th percentile (clamp upper bound)",
               "PS median among smokers", "PS median among non-smokers",
               "% of non-smokers with PS < 0.01", "% of non-smokers with PS < 0.05",
               "% of smokers with PS < 0.01", "% of smokers with PS > 0.10",
               "% of smokers whose PS was clamped to the lower bound",
               "% of non-smokers whose PS was clamped to the lower bound",
               "max stabilized weight, smokers", "effective n, smokers (of 97,187)",
               "share of total smoker weight carried by top 1% of smokers"),
  value = c(q01, q99, median(e[t == 1]), median(e[t == 0]),
            100 * mean(e[t == 0] < 0.01), 100 * mean(e[t == 0] < 0.05),
            100 * mean(e[t == 1] < 0.01), 100 * mean(e[t == 1] > 0.10),
            100 * mean(e[t == 1] < q01), 100 * mean(e[t == 0] < q01),
            max(w[t == 1]), sum(w[t == 1])^2 / sum(w[t == 1]^2),
            { ws <- sort(w[t == 1], decreasing = TRUE); 100 * sum(ws[1:ceiling(0.01 * length(ws))]) / sum(ws) }))
fwrite(out, "output/paper_numbers.csv")
print(out[, .(quantity, value = round(value, 4))])
cat("\nDone.\n")
