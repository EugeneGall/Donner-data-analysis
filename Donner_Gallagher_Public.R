# Donner_Gallagher
# Written by Eugene.Gallagher@umb.edu 7/10/23, last revised 7/19/23
# Assisted by GPT-4
# Analysis of data from Grayson (1990, Table 1; Grayson 1994)
# References
# Grayson, D. K. 1990. Donner Party Deaths: a Demographic Assessment. J.
#    Anthropological Research 46: 223-242.
# Grayson, D. K. 1994. Differential Mortality and the Donner Party Disaster. 
#    Evolutionary Anthropology 2: 151-159.
# Ramsey, F. L. and D. W. Schafer. 2013. The Statistical Sleuth: a Course in
#    Methods of Data Analysis, 3rd Edition. Brooks/Cole Cengage Learning,
#    Boston MA, 760 pp.
# Approach: added the under 15 data to the Donner data from Statistical Sleuth
#    3rd edition, Changed age of Patrick Breen to 51 (Grayson, 1994, p 155)
# Grayson (1990) argued Family Group Size, Age, and Sex control survival. This
# analysis will test all three variables.

# Use data imputation to fill in the 2 missing ages for the 2 Wolfingers
# Have R determine family size by the numbers of individuals with the same
# last name (not used here). But, also analyze Grayson's Family Group Size.
# Code aided by GPT-4 with many dozens iterative prompts.

# Install and load packages
library(plotly)
library(rms)
library(tidyverse) # contains dplyr and ggplot2

# Read the data
file_path <- "../data/Donner.csv"
# Read the CSV file
Donner <- read.csv(file_path)
str(Donner)

# Calculate family size based on Last_Name, but Family_Group_Size from Grayson
# (1990) Table 1 will be used in this code's analyses.
Donner$Family_Size <- as.integer(ave(Donner$Last_Name, Donner$Last_Name, FUN =length))

# This produced slightly different Family Sizes than Grayson (1990 Table 1) in
# that he merged some families using data that he didn't describe and that
# Gallagher does not have. For example, he lists the Donner party as having 16
# members, but the Donner surname has just 14 members.

# Convert Status to binary
Donner$Status <- ifelse(Donner$Status == "Survived", 1, 0)

# Impute missing Age values for the two Wolfingers using median imputation
Donner$Age[is.na(Donner$Age)] <- median(Donner$Age, na.rm = TRUE)
# Both are assigned the median age of 18.
# I had GPT-4 write a dplyr pipe to replace NA's by the medians for each sex,
# but the median age for Females was 13, so I opted to using median imputation
# producing an imputed median age of 18 for both Woffingers.

# Prepare data for rms, Harrell's 'Regression Modeling Strategies'
ddist <- datadist(Donner)
options(datadist = "ddist")

# Fit the model with Age, Sex, andd Grayson's (1990) Family_Group_Size
mod  <- Glm(Status ~ rcs(Age,3) * Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod1 <- Glm(Status ~ rcs(Family_Group_Size,3), data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod2 <- Glm(Status ~ Age + Sex + rcs(Family_Group_Size,3), data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod3 <- Glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,3), data = Donner, family = binomial(), x = TRUE, y = TRUE)

# Summary and ANOVA for mod
summary(mod)
anova(mod)
# Sex amd Age with restricted cubic spline 3-knot curve and interaction all with
# p < 0.05

# Summary and ANOVA for mod1, Family Group Size alone
summary(mod1)
anova(mod1)
# Very strong nonlinear effect of group size

# Summary and ANOVA for mod2, an additive model, Age, Sex and 
# rcs(Family_Group Size)
summary(mod2)
anova(mod2)
# Strong effect of Family_Group_Size and Sex, but not Age (p=0.24)

# Fit a restricted cubic spline regression (3 knots) for Age and 
# Family Group Size with an interaction
summary(mod3)
anova(mod3)

# Generate 4 ggplot2 graphics
# Define new levels for the predictors
new_age <- seq(min(Donner$Age), max(Donner$Age), length.out = 100)
new_family_group_size <- seq(min(Donner$Family_Group_Size), 
                             max(Donner$Family_Group_Size), length.out = 100)
new_sex <- unique(Donner$Sex)

# Make predictions on the link scale (logit scale)
link_pred  <- Predict(mod, Age = new_age, Sex = new_sex)
link_pred1 <- Predict(mod1, Family_Group_Size = new_family_group_size)
link_pred2 <- Predict(mod2, Age = new_age, Sex = new_sex,
                      Family_Group_Size = new_family_group_size)
link_pred3 <- Predict(mod3, Age = new_age, Sex = new_sex, 
                      Family_Group_Size = new_family_group_size)

# Transform predictions back to the original scale (probability scale)
pred <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred$yhat),
  lower = plogis(link_pred$lower),
  upper = plogis(link_pred$upper)
)

pred1 <- data.frame(
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = new_family_group_size,
  fit = plogis(link_pred1$yhat),
  lower = plogis(link_pred1$lower),
  upper = plogis(link_pred1$upper)
)

# Transform predictions back to the original scale (probability scale)
pred2 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = rep(new_family_group_size, each = length(new_family_group_size)),
  fit = plogis(link_pred2$yhat),
  lower = plogis(link_pred2$lower),
  upper = plogis(link_pred2$upper)
)

pred3 <-  data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = rep(new_family_group_size, each = length(new_family_group_size)),
  fit = plogis(link_pred3$yhat),
  lower = plogis(link_pred3$lower),
  upper = plogis(link_pred3$upper)
)
# Adjust y-values for jittered points
Donner$AdjustedStatus <- ifelse(Donner$Status == 0, -0.03, 1.03)

# Plot the results for the Family_Group_Size Model
ggplot(pred, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus,
                                 color = Sex), width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "rcs(Age, 3) * Sex with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

ggplot(pred1, aes(x = Family_Group_Size, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Family_Group_Size, y = AdjustedStatus,
              color = Sex), width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Family Group Size", y = "Estimated Probability of Survival",
       title = "rcs(Family Group Size, 3) with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the additive model, mod2
plot_3d_mod2 <- plot_ly(data = subset(pred2, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred2, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
         title = "Age + Sex + rcs(Family_Group_Size,3)")
plot_3d_mod2


# Plot the results for the rcs(Age,3) * Sex model
plot_3d_mod3 <- plot_ly(data = subset(pred3, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred3, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
         title = "rcs(Age,3) * Sex + rcs(Family_Group_Size,3)")
plot_3d_mod3

### These 3 graphs plot confidence regions as wedges: too busy
# Plot the results for the additive model, mod2
plot_3d_mod2 <- plot_ly(data = pred2, x = ~Age, y = ~Family_Group_Size, z = ~fit, color = ~Sex, 
                        type = "mesh3d") %>% 
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
         title = "Age + Sex + rcs(Family_Group_Size,3)")
plot_3d_mod2

# Plot the results for the rcs(Age,3) * Sex model
plot_3d_mod3 <- plot_ly(data = pred3, x = ~Age, y = ~Family_Group_Size, z = ~fit, color = ~Sex, 
                        type = "mesh3d") %>% 
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
         title = "rcs(Age,3) * Sex + rcs(Family_Group_Size,3)")
plot_3d_mod3

####### Redo the Statistical Sleuth Analysis using cases with Age>=15 ##########

# Create a new data frame with cases where Age is greater than or equal to 15
Donner_15up <- Donner[Donner$Age >= 15, ]

mod5 <- Glm(Status ~ Age + Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod6 <- Glm(Status ~ Age * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod7 <- Glm(Status ~ rcs(Age,3) * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod8 <- Glm(Status ~ rcs(Age,3) + Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod9 <- Glm(Status ~ rcs(Age,3),       data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod10 <- Glm(Status ~ rcs(Family_Group_Size,3), data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)

# Summary and ANOVA for mod5
summary(mod5)
anova(mod5)
# Age (p=0.060), Sex (p = 0.027)

# Summary and ANOVA for mod6
summary(mod6)
anova(mod6)
# Age (p=0.06), Sex (p = 0.07), Age:Sex (p=0.07)

# Summary and ANOVA for mod7
summary(mod7)
anova(mod7)
# All p's > 0.15

# Summary and ANOVA for mod8
summary(mod8)
anova(mod8)
# Only sex important (p = 0.02)

# Summary and ANOVA for mod9
summary(mod9)
anova(mod9)
# No rcs(Age, 3) effect

# Any effect of rcs (Family Group Size)
summary(mod10)
anova(mod10)
# p = 0.005

# Check whether Age is important in a Wilks drop in deviance test:
# requires glm, not Harrell's Glm
mod5g <- glm(Status ~ Age + Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod6g <- glm(Status ~ Age * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
anova(mod5g, mod6g, test = "Chi")
# Good evidence for an interaction effect p =0.03

# Should we spend the extra df on the restricted cubic spline for Age: NO!
mod7g <- glm(Status ~ rcs(Age,3) * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
anova(mod6g, mod7g, test = "Chi")
# p = 0.575, no extra explanatory value of rcs(Age,3)

# Plot the data
# Define new levels for the predictors
new_age <- seq(min(Donner_15up$Age), max(Donner_15up$Age), length.out = 100)
new_family_group_size <- seq(min(Donner_15up$Family_Group_Size), 
                             max(Donner_15up$Family_Group_Size), length.out = 100)
new_sex <- unique(Donner_15up$Sex)

# Make predictions on the link scale (logit scale)
link_pred5 <- Predict(mod5, Age = new_age, Sex = new_sex)
link_pred6 <- Predict(mod6, Age = new_age, Sex = new_sex)
link_pred7 <- Predict(mod7, Age = new_age, Sex = new_sex)
link_pred8 <- Predict(mod8, Age = new_age, Sex = new_sex)
link_pred9 <- Predict(mod9, Age = new_age)
link_pred10 <- Predict(mod10, Family_Group_Size = new_family_group_size)

# Transform predictions back to the original scale (probability scale)
pred5<- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred5$yhat),
  lower = plogis(link_pred5$lower),
  upper = plogis(link_pred5$upper)
)

# Transform predictions back to the original scale (probability scale)
pred6 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred6$yhat),
  lower = plogis(link_pred6$lower),
  upper = plogis(link_pred6$upper)
)

# Transform predictions back to the original scale (probability scale)
pred7 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred7$yhat),
  lower = plogis(link_pred7$lower),
  upper = plogis(link_pred7$upper)
)

# Transform predictions back to the original scale (probability scale)
pred8 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred8$yhat),
  lower = plogis(link_pred8$lower),
  upper = plogis(link_pred8$upper)
)

# Transform predictions back to the original scale (probability scale)
pred9 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred9$yhat),
  lower = plogis(link_pred9$lower),
  upper = plogis(link_pred9$upper)
)

# Transform predictions back to the original scale (probability scale)
pred10 <- data.frame(
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = new_family_group_size,
  fit = plogis(link_pred1$yhat),
  lower = plogis(link_pred1$lower),
  upper = plogis(link_pred1$upper)
)

# Adjust y-values for jittered points
Donner_15up$AdjustedStatus <- ifelse(Donner_15up$Status == 0, -0.03, 1.03)

# Plot the results for the additive model: Not as informative as intxn model
ggplot(pred5, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Age, y = AdjustedStatus, color = Sex),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "Age + Sex model with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the interaction model
ggplot(pred6, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Age, y = AdjustedStatus, color = Sex),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "Age * Sex model with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the rcs(Age,3) * sex model
ggplot(pred7, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Age, y = AdjustedStatus, color = Sex),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "rcs(Age,3) * Sex, with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the rcs(Age,3) + sex model
ggplot(pred10, aes(x = Family_Group_Size, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Family_Group_Size, y = AdjustedStatus,
                                 color = Sex), width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Family Group Size", y = "Estimated Probability of Survival",
       title = "rcs(Family Group Size, 3) with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Overall conclusion:
# With Age>= 15, there is a poor fit with the rcs(Age,3), but as with the
# Sleuth3 analysis, the Age * Sex interaction is important as determined by the
# Wilks drop in deviance test (p=0.018), indicating the need for an interaction
# term, just not a cubic spline for age.


