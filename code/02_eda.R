###############################################################################
## 02_eda.R
## STEP N: exploratory data analysis (what is in the data?).
##
## Input  : output/analytic_cohort.rds (from 01_prepare_data.R), data/nvss2023_extract.csv
## Output : output/table1_by_smoking.csv, output/eda_numbers.csv, figures/fig_eda_*.png
## HOW TO RUN:  Rscript 02_eda.R      (about 2 minutes)
## Packages  :  data.table, ggplot2
###############################################################################

set.seed(1234)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
dir.create("figures", showWarnings = FALSE)
NAVY <- "#1F3864"; ORANGE <- "#C1611A"; GRAY <- "#6B7F99"
theme_proj <- theme_minimal(base_size = 16) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),   # figure titles are given in the paper/slide captions
        panel.grid.minor = element_blank(), legend.position = "top")

d <- readRDS("output/analytic_cohort.rds")
d[, group := factor(ifelse(smoker == 1, "Smoked in pregnancy", "Did not smoke"),
                    levels = c("Did not smoke", "Smoked in pregnancy"))]
n <- nrow(d)

## ---------------------------------------------------------------------------
## 0. What does the NCHS exposure flag CIG_REC mean?  Check it against the
##    cigarettes-before-pregnancy item (CIG_0) in the raw extract.
## ---------------------------------------------------------------------------
cig <- fread("data/nvss2023_extract.csv", select = c("cig_0", "cig_rec"),
             colClasses = list(character = "cig_rec"))
cig <- cig[cig_0 < 99]                                   # 99 = unknown
cat("CIG_REC versus cigarettes smoked daily BEFORE pregnancy (CIG_0), rows with known CIG_0:\n")
print(cig[, .(n = .N), by = .(cig_rec, smoked_before = cig_0 > 0)][order(cig_rec, smoked_before)])
cat("Share of CIG_REC = Y mothers with CIG_0 = 0 (smoking began or resumed after conception):",
    round(100 * cig[cig_rec == "Y", mean(cig_0 == 0)], 1), "%\n")
cat("Share of mothers with CIG_0 > 0 who have CIG_REC = N (quit before or during pregnancy):",
    round(100 * cig[cig_0 > 0, mean(cig_rec == "N")], 1), "%\n\n")
rm(cig)

## ---------------------------------------------------------------------------
## 1. Table 1: characteristics by smoking status
## ---------------------------------------------------------------------------
pct <- function(v) {                                      # one categorical variable -> rows of n (%)
  tab <- d[, .N, by = .(level = get(v), smoker)]
  tab <- dcast(tab, level ~ smoker, value.var = "N", fill = 0)
  setnames(tab, c("0", "1"), c("n_nonsmoker", "n_smoker"))
  tab[, `:=`(variable = v, pct_nonsmoker = 100 * n_nonsmoker / sum(n_nonsmoker),
             pct_smoker = 100 * n_smoker / sum(n_smoker))]
  pt <- d[, .(pt = 100 * mean(preterm)), by = .(level = get(v), smoker)]
  tab[, preterm_pct_nonsmoker := pt[smoker == 0][match(tab$level, level), pt]]
  tab[, preterm_pct_smoker    := pt[smoker == 1][match(tab$level, level), pt]]
  tab[, level := as.character(level)]
  tab[]
}
t1 <- rbindlist(lapply(c("age_grp", "race", "educ", "precare", "wic", "parity", "ppterm", "pdiab", "phype", "bmi_grp"), pct))
setcolorder(t1, c("variable", "level"))
fwrite(t1, "output/table1_by_smoking.csv")
options(width = 200)
cat("Table 1 (first rows):\n"); print(t1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)][1:14])

## Headline numbers used in the paper
num <- data.table(
  quantity = c("N analytic cohort", "N smokers", "N non-smokers", "smoking prevalence %",
               "preterm % overall", "preterm % smokers", "preterm % non-smokers",
               "mean gestational age smokers (wk)", "mean gestational age non-smokers (wk)",
               "mean maternal age smokers", "mean maternal age non-smokers",
               "mean BMI smokers", "mean BMI non-smokers"),
  value = c(n, sum(d$smoker), sum(1 - d$smoker), 100 * mean(d$smoker), 100 * mean(d$preterm),
            100 * mean(d$preterm[d$smoker == 1]), 100 * mean(d$preterm[d$smoker == 0]),
            mean(d$ga_weeks[d$smoker == 1]), mean(d$ga_weeks[d$smoker == 0]),
            mean(d$mager[d$smoker == 1]), mean(d$mager[d$smoker == 0]),
            mean(d$bmi[d$smoker == 1]), mean(d$bmi[d$smoker == 0])))
fwrite(num, "output/eda_numbers.csv"); print(num[, .(quantity, value = round(value, 2))])

## ---------------------------------------------------------------------------
## 2. Figures
## ---------------------------------------------------------------------------
## 2a. Gestational-age distribution by smoking status
ga <- d[ga_weeks >= 24, .N, by = .(group, ga_weeks)]
ga[, share := 100 * N / sum(N), by = group]
p1 <- ggplot(ga, aes(ga_weeks, share, fill = group)) +
  geom_col(position = "dodge", width = 0.85) +
  geom_vline(xintercept = 36.5, linetype = "dashed", colour = "black") +
  annotate("text", x = 36.3, y = max(ga$share) * 0.95, label = "Preterm\n(< 37 weeks)", hjust = 1, size = 3.8) +
  scale_fill_manual(values = c(GRAY, ORANGE), name = NULL) +
  labs(x = "Gestational age at delivery (completed weeks)", y = "% of births within group",
       title = "Smokers deliver earlier: the whole distribution shifts left") + theme_proj
ggsave("figures/fig_eda_gestational_age.png", p1, width = 9, height = 4.8, dpi = 200)

## 2b. Who smokes? Prevalence of smoking by education, race, age
prev <- function(v, lab) d[, .(prev = 100 * mean(smoker), n = .N), by = .(level = get(v))][, var := lab][]
pv <- rbind(prev("educ", "Education"), prev("age_grp", "Maternal age"), prev("race", "Race/ethnicity"))
pv[, level := factor(level, levels = unlist(lapply(c("educ", "age_grp", "race"), function(v) levels(d[[v]]))))]
p2 <- ggplot(pv, aes(level, prev)) +
  geom_col(fill = NAVY, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", prev)), vjust = -0.4, size = 4.6) +
  facet_wrap(~var, scales = "free_x", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "% who smoked in pregnancy", title = "Smoking in pregnancy is strongly patterned by education, race and age") +
  theme_proj + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("figures/fig_eda_who_smokes.png", p2, width = 13, height = 5.2, dpi = 200)

## 2c. Crude preterm risk by smoking status within levels of key confounders
cr <- function(v, lab) d[, .(risk = 100 * mean(preterm), n = .N), by = .(level = get(v), group)][, var := lab][]
cc <- rbind(cr("educ", "Education"), cr("age_grp", "Maternal age"), cr("race", "Race/ethnicity"))
cc[, level := factor(level, levels = unlist(lapply(c("educ", "age_grp", "race"), function(v) levels(d[[v]]))))]
p3 <- ggplot(cc, aes(level, risk, fill = group)) +
  geom_col(position = "dodge", width = 0.75) +
  facet_wrap(~var, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c(GRAY, ORANGE), name = NULL) +
  labs(x = NULL, y = "% preterm", title = "Preterm risk is higher in smokers at every level, and also varies with the confounders") +
  theme_proj + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave("figures/fig_eda_preterm_by_group.png", p3, width = 13, height = 5.2, dpi = 200)

## 2d. Missing data
ms <- fread("output/missingness.csv")
ms[, variable := factor(variable, levels = ms[order(pct_missing), variable])]
p4 <- ggplot(ms, aes(pct_missing, variable)) + geom_col(fill = NAVY, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", pct_missing)), hjust = -0.15, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "% of births with the item missing or unknown", y = NULL,
       title = "Missing data are limited") + theme_proj
ggsave("figures/fig_eda_missing.png", p4, width = 8, height = 4.4, dpi = 200)
cat("\nDone.\n")
