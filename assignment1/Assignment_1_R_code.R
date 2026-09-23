# STAT 7900 - Assignment 1: R code for Problems 4 and 5
# Run in RStudio; needs dplyr and ggplot2. Expected: Problem 4 crude 3.242, adjusted 1.985; negative-effect run 0.242 and -1.015.

# Problem 4
set.seed(7901)
n <- 5000
C <- rnorm(n)
T <- rbinom(n, 1, plogis(C))
Y <- 2*T + 1.5*C + rnorm(n)
dat <- data.frame(Y, T, C)

crude <- lm(Y ~ T, data = dat)
adj   <- lm(Y ~ T + C, data = dat)
cat("== 4(a)(b): true effect = 2 ==\n")
print(round(summary(crude)$coefficients, 4))
print(round(summary(adj)$coefficients, 4))
cat("crude T:", round(coef(crude)["T"], 3), " adjusted T:", round(coef(adj)["T"], 3), "\n")
cat("cor(C,T) =", round(cor(dat$C, dat$T), 3), " mean C | T=1:", round(mean(C[T==1]),3), " mean C | T=0:", round(mean(C[T==0]),3), "\n")

# 4(e): rerun with negative effect
set.seed(7901)
n <- 5000
C <- rnorm(n)
T <- rbinom(n, 1, plogis(C))
Y <- -1*T + 1.5*C + rnorm(n)
dat2 <- data.frame(Y, T, C)
crude2 <- lm(Y ~ T, data = dat2)
adj2   <- lm(Y ~ T + C, data = dat2)
cat("\n== 4(e): true effect = -1 ==\n")
print(round(summary(crude2)$coefficients, 4))
print(round(summary(adj2)$coefficients, 4))
cat("crude T:", round(coef(crude2)["T"], 3), " adjusted T:", round(coef(adj2)["T"], 3), "\n")


suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
data(UCBAdmissions)
ucb <- as.data.frame(UCBAdmissions)
print(head(ucb))

a <- ucb %>% group_by(Gender) %>%
  summarise(applicants = sum(Freq),
            admitted = sum(Freq[Admit == 'Admitted']),
            admission_rate = admitted / applicants)
cat("\n== 5(a) ==\n"); print(as.data.frame(a), digits = 4)

dept_rates <- ucb %>% group_by(Dept, Gender) %>%
  summarise(applicants = sum(Freq),
            admitted = sum(Freq[Admit == 'Admitted']),
            admission_rate = admitted / applicants,
            .groups = 'drop')
cat("\n== 5(b) ==\n"); print(as.data.frame(dept_rates), digits = 4)

# extra diagnostics for 5(d)
w <- dept_rates %>% group_by(Dept) %>%
  summarise(total_apps = sum(applicants),
            dept_rate = sum(admitted)/sum(applicants),
            male_apps = applicants[Gender=="Male"], female_apps = applicants[Gender=="Female"])
w <- w %>% mutate(male_share = male_apps/total_apps, female_share = female_apps/total_apps)
cat("\n== dept-level: selectivity & applicant mix ==\n"); print(as.data.frame(w), digits = 3)
g <- dept_rates %>% group_by(Gender) %>% mutate(share = applicants/sum(applicants))
cat("\n== share of each gender's applicants by dept ==\n"); print(as.data.frame(g %>% select(Dept, Gender, share)), digits=3)
cat("\nAdmitted/rejected total: ", sum(ucb$Freq), "\n")

# Standardization: apply common (total) dept mix
tot <- w %>% mutate(wt = total_apps/sum(total_apps)) %>% select(Dept, wt)
std <- dept_rates %>% left_join(tot, by="Dept") %>% group_by(Gender) %>% summarise(std_rate = sum(admission_rate*wt))
cat("\n== Standardized (common dept mix) ==\n"); print(as.data.frame(std), digits=4)

p <- ggplot(dept_rates, aes(x = Dept, y = admission_rate, fill = Gender)) +
  geom_col(position = 'dodge') +
  labs(x = 'Department', y = 'Admission Rate',
       title = 'UC Berkeley Admission Rates by Department and Gender')
print(p)
cat("\nplot saved\n")
