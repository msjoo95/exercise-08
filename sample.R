library(tidyverse)
library(skimr)

# Step 1
f<-"https://raw.githubusercontent.com/difiore/ada-datasets/main/Street_et_al_2017.csv"
d<-read_csv(f, col_names=TRUE)

skim(d)

# Step 2
plot_1 <- d %>% 
  ggplot(aes(x=Group_size, y=ECV)) +
  geom_point(shape = 21, colour = "red", fill = "black", size=1, na.rm = TRUE) +
  ggtitle("ECV vs Group Size") +
  theme_minimal()
plot_2 <- d %>% 
  ggplot(aes(x=Longevity, y=ECV)) +
  geom_point(shape = 21, colour = "red", fill = "black", size=1, na.rm = TRUE) +
  ggtitle("ECV vs Longevity") +
  theme_minimal()
plot_3 <- d %>% 
  ggplot(aes(x=Weaning, y=ECV)) +
  geom_point(shape = 21, colour = "red", fill = "black", size=1, na.rm = TRUE) +
  ggtitle("ECV vs Weaning") +
  theme_minimal()
plot_4 <- d %>% 
  ggplot(aes(x=Repro_lifespan, y=ECV)) +
  geom_point(shape = 21, colour = "red", fill = "black", size=1, na.rm = TRUE) +
  ggtitle("ECV vs Repro_lifespan") +
  theme_minimal()
library(cowplot)
plot_grid(plot_1, plot_2, plot_3, plot_4)

# Step 3
d_clean <- d %>% 
  filter(!is.na(ECV) & !is.na(Group_size))
(beta_1 <- cov(d_clean$Group_size, d_clean$ECV)/var(d_clean$Group_size))
(beta_0 <- mean(d_clean$ECV) - beta_1*mean(d_clean$Group_size))

# Step 4
(m <- lm(ECV~Group_size, data=d_clean))
summary(m)

# Step 5 - similar slope among different groups
d_Catarrhini <- d %>% 
  filter(Taxonomic_group=='Catarrhini') %>% 
  filter(!is.na(ECV) & !is.na(Group_size))
d_Platyrrhini <- d %>% 
  filter(Taxonomic_group=='Platyrrhini') %>% 
  filter(!is.na(ECV) & !is.na(Group_size))
d_Strepsirhini <- d %>% 
  filter(Taxonomic_group=='Strepsirhini') %>% 
  filter(!is.na(ECV) & !is.na(Group_size))

lm(ECV~Group_size, data=d_Catarrhini)
lm(ECV~Group_size, data=d_Platyrrhini)
lm(ECV~Group_size, data=d_Strepsirhini)

# Step 6
SSY <- sum((m$model$ECV - mean(m$model$ECV))^2) # y - y_mean
SSR <- sum((m$fitted.values - mean(m$model$ECV))^2) # y_predicted - y_mean
SSE <- sum((m$model$ECV - m$fitted.values)^2) # or SSE <- sum(m$residuals^2) -> y - y_predicted
SSY == SSE + SSR # just for checking
(df_regression <- 1)  # p = 1 (the number of predictor variables)
(df_error <- nrow(d_clean) - df_regression - 1)  # n - p - 1
(df_y <- nrow(d_clean) - df_regression)  # n - p

MSR <- SSR/df_regression  # mean variance explained by the regression equation
MSE <- SSE/df_error  # mean remaining variance
MSY <- SSY/df_y  # mean overall variance

SSX <- sum((m$model$Group_size - mean(m$model$Group_size))^2)  # how much x variation there is
(SEbeta1 <- sqrt(MSE/SSX)) # the standard error for the slope coefficient!!!

beta_1 + qt(p = c(0.025, 0.975), df=df_error)*SEbeta1 # 95% of CI!!!

t.calc <- beta_1 / SEbeta1 # t-statistic
(p.calc = 2 * (1 - pt(abs(t.calc), df = nrow(m$model) - 2))) # p-value!!!!

summary(m) # same results
confint(m) # same results

# Step 7
# permutation approach to calculate p value
m <- lm(ECV ~ Group_size, data = d_clean)
beta1 <- coef(m)[2]

n_permutations <- 1000
perm_beta1 <- numeric(n_permutations)

for (i in 1:n_permutations) {
  permuted_ECV <- sample(d_clean$ECV)
  perm_model <- lm(permuted_ECV ~ d_clean$Group_size)
  perm_beta1[i] <- coef(perm_model)[2]
}

SE_permutation <- sd(perm_beta1)
t_stat_perm <- beta1 / SE_permutation
p_value <- 2 * (1 - pt(abs(t_stat_perm), df = nrow(d_clean) - 2))
p_value

# Step 8
m <- lm(ECV ~ Group_size, data = d_clean)
beta1 <- coef(m)[2]

n_bootstrap <- 1000
bootstrap_beta1 <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  boot_sample <- d_clean[sample(nrow(d_clean), replace = TRUE), ]
  boot_model <- lm(ECV ~ Group_size, data = boot_sample)
  bootstrap_beta1[i] <- coef(boot_model)[2]
}
#quantile method
CI_quantile <- quantile(bootstrap_beta1, c(0.025, 0.975))
CI_quantile
#theory-based method
SE_bootstrap <- sd(bootstrap_beta1)
t_critical <- qt(0.975, df = nrow(d_clean) - 2)
CI_theory <- beta1 + c(-1, 1) * t_critical * SE_bootstrap
CI_theory
#The slope (β₁) is different from zero, suggesting that the relationship is likely to be significant.
