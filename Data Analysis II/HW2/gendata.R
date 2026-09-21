library(tidyverse)
library(lme4)

set.seed(123)

n_subjects <- 5
n_per_subject <- 25

# 1. Create Subject-level offsets
# We purposely make subjects with higher x-means have higher y-intercepts
subject_data <- data.frame(
  Subject = as.factor(1:n_subjects),
  x_center = seq(-4, 4, length.out = n_subjects),  # Centers of x move right
  intercept = seq(2, 18, length.out = n_subjects), # Intercepts move UP
  slope = rnorm(n_subjects, -1, 0.2)             # All slopes are NEGATIVE
)

# 2. Generate individual observations
df <- expand.grid(Subject = as.factor(1:n_subjects), 
                  obs = 1:n_per_subject) %>%
  left_join(subject_data, by = "Subject") %>%
  mutate(
    # x varies around the subject's center
    x = rnorm(n(), mean = x_center, sd = 2),
    # y = individual intercept + (negative slope * x) + noise
    y = intercept + (slope * x) + rnorm(n(), sd = 0.7)
  )

# 3. Compare the Models
# OLS (ignores Subject) -> Expect Positive Slope
model_ols <- lm(y ~ x, data = df)

# LMM (accounts for Subject) -> Expect Negative Slope
model_lmer <- lmer(y ~ x + (1 + x | Subject), data = df)

# Check results
print(paste("OLS Slope:", round(coef(model_ols)[2], 3)))
print(fixef(model_lmer))

ggplot(df, aes(x = x, y = y, color = Subject)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) + # Subject-specific (Negative)
  geom_smooth(aes(group = 1), method = "lm", color = "black", 
              linetype = "dashed", size = 1) + # Global OLS (Positive)
  theme_minimal() +
  labs(title = "Simpson's Paradox: Global Positive vs. Local Negative")


df <- tibble(x=df$x, y=df$y, subject=df$Subject)
write.csv(df, "HW2data2.csv")
