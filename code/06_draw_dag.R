###############################################################################
## 06_draw_dag.R
## Draws the causal diagram (DAG) used to choose the adjustment set.
##
## Output : figures/fig_dag.png
## HOW TO RUN:  Rscript 06_draw_dag.R      (a few seconds; needs only ggplot2)
##
## READING THE DIAGRAM
##   * Smoking in pregnancy (T) -> Preterm birth (Y) is the causal effect we want.
##   * Each of the ten confounders has one arrow into T and one arrow into Y.
##     That gives ten "back-door" paths  T <- X -> Y  (gray arrows).
##   * Blue dashed arrows are assumed relationships AMONG the confounders. They do
##     not change which variables must be adjusted for (see Section 3 of the paper):
##     every confounder is a direct parent of both T and Y, so conditioning on all
##     ten blocks every back-door path, whatever the arrows among them.
##   * Insurance type is deliberately NOT in the diagram (see the paper).
###############################################################################

set.seed(1234)
suppressPackageStartupMessages(library(ggplot2))
dir.create("figures", showWarnings = FALSE)

NAVY <- "#1F3864"; ORANGE <- "#C1611A"; GRAY <- "#8A97A8"; STEEL <- "#2F6FA8"

## Left-to-right order chosen so that related confounders sit next to each other
conf <- data.frame(
  id    = c("race", "educ", "wic", "pnc", "age", "parity", "ppt", "bmi", "diab", "hyp"),
  label = c("Race /\nethnicity", "Maternal\neducation", "WIC\nparticipation", "Prenatal-care\ntiming",
            "Maternal\nage", "Parity\n(prior births)", "Prior\npreterm birth", "BMI\n(pre-preg.)",
            "Diabetes\n(pre-preg.)", "Hypertension\n(pre-preg.)"),
  stringsAsFactors = FALSE)
conf$x <- seq(1, by = 1.58, length.out = nrow(conf))
conf$y <- 5

T_node <- data.frame(label = "Smoking during\npregnancy (T)", x = mean(conf$x) - 3.4, y = 1.2)
Y_node <- data.frame(label = "Preterm birth\n< 37 weeks (Y)", x = mean(conf$x) + 3.4, y = 1.2)

## Arrows: every confounder -> T and -> Y. Arrow tips are spread along the top edge of the
## T and Y boxes so that the ten arrows can be told apart.
spread <- seq(-1.0, 1.0, length.out = nrow(conf))
arr_T <- data.frame(x = conf$x, y = 4.62, xend = T_node$x + spread, yend = T_node$y + 0.40)
arr_Y <- data.frame(x = conf$x, y = 4.62, xend = Y_node$x + spread, yend = Y_node$y + 0.40)

## Assumed arrows among confounders (from, to)
inter <- data.frame(from = c("race", "educ", "age", "age", "parity", "bmi", "bmi"),
                    to   = c("educ", "wic", "parity", "bmi", "ppt", "diab", "hyp"))
inter$x    <- conf$x[match(inter$from, conf$id)] + 0.15
inter$xend <- conf$x[match(inter$to,   conf$id)] - 0.15
inter$y    <- 5.38;  inter$yend <- 5.38

p <- ggplot() +
  geom_segment(data = arr_T, aes(x, y, xend = xend, yend = yend), colour = GRAY, linewidth = 0.45,
               arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  geom_segment(data = arr_Y, aes(x, y, xend = xend, yend = yend), colour = GRAY, linewidth = 0.45,
               arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  geom_curve(data = inter[inter$from != "age" | inter$to != "bmi", ],
             aes(x, y, xend = xend, yend = yend), curvature = -0.9, colour = STEEL, linewidth = 0.6,
             linetype = "dashed", arrow = arrow(length = unit(0.11, "inches"), type = "closed")) +
  geom_curve(data = inter[inter$from == "age" & inter$to == "bmi", ],
             aes(x, y, xend = xend, yend = yend), curvature = -0.45, colour = STEEL, linewidth = 0.6,
             linetype = "dashed", arrow = arrow(length = unit(0.11, "inches"), type = "closed")) +
  geom_segment(aes(x = T_node$x + 1.5, y = T_node$y, xend = Y_node$x - 1.5, yend = Y_node$y),
               colour = ORANGE, linewidth = 1.8, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  geom_label(data = conf, aes(x, y, label = label), fill = "#EEF2F7", colour = "#1A2233",
             label.size = 0.5, label.r = unit(0.25, "lines"), size = 3.7, lineheight = 0.95,
             label.padding = unit(0.45, "lines")) +
  geom_label(data = T_node, aes(x, y, label = label), fill = ORANGE, colour = "white", fontface = "bold",
             size = 5.0, label.size = 0, label.padding = unit(0.6, "lines")) +
  geom_label(data = Y_node, aes(x, y, label = label), fill = NAVY, colour = "white", fontface = "bold",
             size = 5.0, label.size = 0, label.padding = unit(0.6, "lines")) +
  annotate("text", x = mean(conf$x), y = 1.55, label = "Causal effect of interest", colour = ORANGE,
           fontface = "bold", size = 5) +
  annotate("text", x = mean(conf$x), y = 0.5, size = 4.1, colour = "#444444",
           label = "Gray arrows: each confounder causes both smoking and preterm birth (ten back-door paths)\npre-preg. = pre-pregnancy; adjusting for all ten blocks every back-door path") +
  annotate("text", x = mean(conf$x), y = 6.95, size = 4.1, colour = STEEL,
           label = "Blue dashed arrows: assumed relationships among confounders (they do not change the adjustment set)") +
  coord_cartesian(xlim = c(min(conf$x) - 0.9, max(conf$x) + 0.9), ylim = c(-0.1, 7.2)) +
  theme_void() + theme(plot.margin = margin(8, 8, 8, 8), plot.background = element_rect(fill = "white", colour = NA))

ggsave("figures/fig_dag.png", p, width = 14.5, height = 6.6, dpi = 200)
cat("Saved figures/fig_dag.png\n")
