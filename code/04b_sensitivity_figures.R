###############################################################################
## 04b_sensitivity_figures.R
## Draws the robustness and subgroup figures from the CSV files written by
## 04_sensitivity.R (no models are fitted here, so it runs in a few seconds).
##
## Input  : output/sens_truncation.csv, output/sens_specifications.csv,
##          output/subgroup_aiptw.csv
## Output : figures/fig_robustness.png, figures/fig_subgroups.png
## HOW TO RUN:  Rscript 04b_sensitivity_figures.R
###############################################################################

set.seed(1234)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
NAVY <- "#1F3864"; ORANGE <- "#C1611A"; GRAY <- "#6B7F99"
dir.create("figures", showWarnings = FALSE)
trunc_tab <- fread("output/sens_truncation.csv")
spec_tab  <- fread("output/sens_specifications.csv")
sg        <- fread("output/subgroup_aiptw.csv")
sg[, level := as.character(level)]
trunc_tab[, scenario := gsub("^DROP ", "Delete ", scenario)]   # older files used the word DROP

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------
theme_proj <- theme_minimal(base_size = 16) +
  theme(plot.title = element_blank(), panel.grid.minor = element_blank())   # titles are in the captions

rb <- rbind(
  trunc_tab[estimator %in% c("AIPTW (doubly robust)", "Overlap-weighted (ATO)"),
            .(group = "Propensity\nscore\nhandling", label = scenario, RD_pp, RD_lo, RD_hi)],
  spec_tab[estimator == "AIPTW (doubly robust)",
           .(group = "Covariate\nset and\nmodel form", label = specification, RD_pp, RD_lo, RD_hi)])
rb[, label := factor(label, levels = rev(unique(label)))]
p <- ggplot(rb, aes(RD_pp, label)) +
  geom_vline(xintercept = 0, colour = GRAY) +
  geom_errorbarh(aes(xmin = RD_lo, xmax = RD_hi), height = 0.25, colour = NAVY) +
  geom_point(size = 2.8, colour = ORANGE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Risk difference (percentage points)", y = NULL) + theme_proj +
  theme(strip.text.y = element_text(face = "bold", size = 12, angle = 0))
ggsave("figures/fig_robustness.png", p, width = 11, height = 6.4, dpi = 200)

sgp <- copy(sg); sgp[, lab := paste0(level)]
sgp[, lab := factor(lab, levels = rev(unique(lab)))]
sgp[, variable := factor(variable, levels = c("age_grp", "race", "educ", "bmi_grp", "parity", "ppterm"),
                         labels = c("Age", "Race/\nethnicity", "Education", "BMI", "Parity", "Prior\npreterm"))]
p2 <- ggplot(sgp, aes(RD_pp, lab)) +
  geom_vline(xintercept = 0, colour = GRAY) +
  geom_errorbarh(aes(xmin = RD_lo, xmax = RD_hi), height = 0.25, colour = NAVY) +
  geom_point(size = 2.4, colour = ORANGE) +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = "AIPTW risk difference (percentage points)", y = NULL) + theme_proj +
  theme(strip.text.y = element_text(size = 12, angle = 0))
ggsave("figures/fig_subgroups.png", p2, width = 9, height = 10, dpi = 200)
cat("Saved figures/fig_robustness.png and figures/fig_subgroups.png\n")
