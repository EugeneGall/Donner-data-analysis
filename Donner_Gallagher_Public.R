# Donner_Gallagher_Public
# Written by Eugene.Gallagher@umb.edu 7/10/23, last revised 8/13/23
# Analysis of Donner data from Grayson (1990, Table 1 &  Grayson 1994)
# Aided by OpenAI GPT-4
# References
# Burnham, K. P. and Anderson, D.R. (2004), “Multimodal inference:
#   understanding AIC and BIC in model selection,” Sociological Methods and
# Grayson, D. K. 1990. Donner Party Deaths: a Demographic Assessment. J.
#    Anthropological Research 46: 223-242.
# Grayson, D. K. 1994. Differential Mortality and the Donner Party Disaster. 
#    Evolutionary Anthropology 2: 151-159.
# Grayson, D. K. (1997) “The timing of Donner Party deaths” Appendix 3 in
#    Hardesty, D. L. (1997) The Archaeology of the Donner Party, Reno:
#    University of Nevada Press.
# Ramsey, F. L. and D. W. Schafer. 2013. The Statistical Sleuth: a Course in
#    Methods of Data Analysis, 3rd Edition. Brooks/Cole Cengage Learning,
#    Boston MA, 760 pp.
# Rarek, E. 2008. Desperate Passage: The Donner Party’s Perilous Journey West.
#    Oxford University Press, Oxford. 304 pp.
# Stewart, G. E. 1960. Ordeal by Hunger: the Ordeal of the Donner Party.
#    Houghton Mifflin, Boston.392 pp.
# 
# Approach: 1) add the under 15 Age data to the Donner data from Statistical
# Sleuth 3rd edition, 2) change age of Patrick Breen to 51 (Grayson, 1994, 
# p 155). 3) Grayson (1990, 1997) argued Family Group Size, Age, and Gender
# control survival and Rarek emphasized the poor survivorship of teamsters and
# servants. This R code will analyze the effects of all four variables.
# 4) Reviewed above books to find death dates for surviving travelers and
# performed survival analyses confirming Grayson (1997) on death timing.

# Used data imputation to fill in the missing age for Mr. Wolfinger
# Have R determine family size by the numbers of individuals with the same
# last name (not used here). But, also analyzed Grayson's (1990) Family Group
# Size, which incorporates information from Stewart's (1960) roster on who was
# traveling with each Family Travel Group.

# Optional statement needed to see full output in a word processor:
sink("my_output.txt")   # Optional Redirect output to a file, make sure that
                        # that the file readme.txt is not open in another app.

# Install and load packages
library(boot) # for cv.glm function
library(caret) # For GAM cross-validation
library(forestplot)  # for the Cox PH plots
library(mgcv) # Wood's package for GAMs
library(plotly) # for 3-d graphics
library(rms)    # Harrell's regression modeling strategies pacakge
library(survival)  # for Cox Proportional Hazards Modeling
library(survminer) # for Cox Proportional Hazards Modeling
library(tidyverse) # contains dplyr and ggplot2

# Read the data from Gallagher's github site (public has read, not edit access)

Donner <- read.csv("https://raw.githubusercontent.com/EugeneGall/donner-data-analysis/main/Donner.csv")

# Calculate family size based on Last_Name, but Family_Group_Size from Grayson
# (1990) Table 1 will be used in this code's analyses.
Donner$Family_Size <- as.integer(ave(Donner$Last_Name, Donner$Last_Name, FUN =length))

str(Donner)

# This produced slightly different Family Sizes than Grayson (1990 Table 1) in
# that Grayson (1990) included family links using data from Stewart's (1960, 
# p 363-364) Donner Party Roster.

# Convert Status to binary
Donner$Status <- ifelse(Donner$Status == "Survived", 1, 0)

# Calculate the Survival_Time which is Death Date or last rescue (Lewis 
# Keseberg) - First_Snow, October 28 1846.
# Convert the character strings to Date objects
Donner$First_Snow_Date <- as.Date(Donner$First_Snow, format="%Y-%m-%d")
Donner$Last_Date_Date <- as.Date(Donner$Last_Date, format="%Y-%m-%d")

# Calculate the difference in days for individual survivorship
Donner$Survival_Time <- as.numeric(Donner$Last_Date_Date - Donner$First_Snow_Date)

# View the first few rows of the data frame to confirm the results
head(Donner)

# Impute missing Age values for Mr. Wolfinger, using median imputation
# Rarek (2008, xi) provides Dora Wolfinger's age as 20. This section will
# Assign Mr. Wolfinger an age of 23, which seems reasonable.
# Compute the median age for each sex
medians <- Donner %>%
  group_by(Sex) %>%
  summarize(median_age = median(Age, na.rm = TRUE))
# Replace NAs in Age based on the median age of each sex
Donner <- Donner %>%
  left_join(medians, by = "Sex") %>%
  mutate(Age = ifelse(is.na(Age), median_age, Age)) %>%
  select(-median_age)
Donner$Age[is.na(Donner$Age)] <- median(Donner$Age, na.rm = TRUE)

### Optional: Delete 8 cases for individuals who died before the first Snowstorm
# on 1846-10-28. This will reduce the number of cases to 87 - 8 = 79

Donner <- Donner %>%
  filter(Delete != 1)
# View the structure to confirm the filtering
str(Donner)

# Prepare data for rms, Harrell's 'Regression Modeling Strategies' These
# two statements are required for the summary of the Harrell's Glm
ddist <- datadist(Donner)
options(datadist = "ddist")

# Fit the model with Age, Sex, and Grayson's (1990) Family_Group_Size
# In deciding on knot size, I used the rule that the minimum knot size should
# be chosen (rcs has a min of 3) unless there is another larger number of knots
# with an AIC which is more than than 4 AIC units lower (Burnham & Anderson
# 2004, p 271)

# Two methods will be coded for finding the appropriate number of knots:
# 1) Brute force by manually plugging in different knot numbers in Glm
# 2) An optimization routine based on AIC

##### First approach: brute force fitting to find minimal AICs.
##### Don't use for pared data, too time-consuming.
#     3 is the minimum knot size permitted with Harrell's rcs function
mod  <- Glm(Status ~ rcs(Age,3) * Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod1 <- Glm(Status ~ rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod2 <- Glm(Status ~ Age + Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
mod3 <- Glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)

# Summary and ANOVA for mod
summary(mod)
anova(mod)
AIC(mod)
# Sex and Age with restricted cubic spline 3-knot curve and interaction all with
# p < 0.05
#   mod     AIC
# 3 knots  92.89261 * Chosen because within 4 of lowest AIC
# 4 knots  91.25976
# 5 knots  92.63563 ** Not chosen because only 0.26957 less than 3 knot AIC
# 6 knots Apparently Singular Matrix, no estimation possible

# Summary and ANOVA for mod1, Family Group Size alone
summary(mod1)
anova(mod1)
AIC(mod1)
# AIC = 106.51
# Very strong nonlinear effect of group size
#   mod2    AIC
# 3 knots 104.1358
# 4 knots 101.3932
# 5 knots 100.6443 * Chosen because within 4 AIC of minimum of 97.10704
# 6 knots 99.11325
# 7 knots 97.10704
# 8 knots 97.10704
# 9 knots 97.10704
# Summary and ANOVA for mod2, an additive model, Age, Sex and 
# rcs(Family_Group Size)
summary(mod2)
anova(mod2)
AIC(mod2)
#   mod2    AIC
# 3 knots 101.5496
# 4 knots 100.6285
# 5 knots  97.38427 * Chosen because within 4 AIC of minimum
# 6 knots  94.06806
# 7 knots  94.03242
# 8 knots  94.03242
# 9 knots  94.03242

# Strong effect of Sex, but not Age (p=0.24) or Family Group Size (p= 0.29)

# Fit a restricted cubic spline regression for Age and Family Group Size with an
# rcs(Age,3)* Sex interaction effect
summary(mod3)
anova(mod3)
AIC(mod3)
#   mod3    AIC
# 3 knots 94.07029
# 4 knots 78.58391 *
# 5 knots 79.2632
# 6 knots Singular Matrix no estimation possible
##### End of brute force approach

##### Second approach: GPT4 Streamlined the above code block with functions:

# Function to fit models with one rcs term
fit_best_rcs <- function(formula_template, data, knot_range) {
  best_aic <- Inf
  best_model <- NULL
  best_knots <- knot_range[1]
  results <- list()
  
  for (k in knot_range) {
    formula <- as.formula(gsub("KNOTS", k, formula_template))
    model <- glm(formula, data = data, family = binomial())
    current_aic <- AIC(model)
    results[[as.character(k)]] <- current_aic
    
    if (current_aic < best_aic && (best_aic - current_aic > 4 || k < best_knots)) {
      best_aic <- current_aic
      best_model <- model
      best_knots <- k
    }
  }
  
  return(list(model = best_model, knots = best_knots, results = results))
}

# Function to fit model with two rcs terms
fit_best_double_rcs <- function(formula, data, knot_range, family = binomial()) {
  best_aic <- Inf
  best_k_age <- max(knot_range)
  best_k_fgs <- max(knot_range)
  
  for (k_age in knot_range) {
    for (k_fgs in knot_range) {
      tryCatch({
        current_formula <- as.formula(
          gsub("knots_age", as.character(k_age), 
               gsub("knots_fgs", as.character(k_fgs), formula)))
        current_model <- glm(current_formula, data = data, family = family)
        current_aic <- AIC(current_model)
        
        if (current_aic < best_aic) {
          best_aic <- current_aic
          best_k_age <- k_age
          best_k_fgs <- k_fgs
        }
      }, error = function(e) {}) # Handle potential errors
    }
  }
  
  best_formula <- as.formula(
    gsub("knots_age", as.character(best_k_age), 
         gsub("knots_fgs", as.character(best_k_fgs), formula)))
  best_model <- glm(best_formula, data = data, family = family)
  
  return(list(model = best_model, knots_age = best_k_age, knots_fgs = best_k_fgs))
}

# Extract and print details
print_details <- function(model_result, model_name) {
  cat("\n---", model_name, "---\n")
  if ("knots_age" %in% names(model_result)) {
    cat("Knots: ", model_result$knots_age, " & ", model_result$knots_fgs, "\n")
  } else {
    cat("Knots: ", model_result$knots, "\n")
  }
  cat("AIC: ", AIC(model_result$model), "\n")
  cat("\nModel Summary:\n")
  print(summary(model_result$model))
  cat("\nModel Anova:\n")
  print(anova(model_result$model))
}

# Fitting models
# Note that 3 is the minimum knot size for rcs, dozens of warnings with knots>5
knot_range <- 3:5
mod_results <- fit_best_rcs("Status ~ rcs(Age, KNOTS) * Sex", Donner, knot_range)
mod1_results <- fit_best_rcs("Status ~ rcs(Family_Group_Size, KNOTS)", Donner, knot_range)
mod2_results <- fit_best_rcs("Status ~ rcs(Family_Group_Size, KNOTS) * Sex", Donner, knot_range)
mod3_results <- fit_best_double_rcs("Status ~ rcs(Age, knots_age) * Sex + rcs(Family_Group_Size, knots_fgs)", Donner, knot_range)

# Printing the details
print_details(mod_results, "mod")
print_details(mod1_results, "mod1")
print_details(mod2_results, "mod2")
print_details(mod3_results, "mod3")

###### End of GPT-4 AIC optimization code ######

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

# Plot the results for the rcs(Age,3) * Sex model
set.seed(8) 
ggplot(pred, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus,
            shape = Sex, color = Sex), width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
#       title = "rcs(Age, 3) * Sex with 95% confidence intervals",
  title = "79 travelers, rcs(Age, 3) * Sex with 95% confidence intervals",
        color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the rcs(Family_Group_Size,5) model
ggplot(pred1, aes(x = Family_Group_Size, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Family_Group_Size, y = AdjustedStatus,
              shape = Sex, color = Sex), width = 0.3, height = 0.03,
              size = 1.5) +
  labs(x = "Family Group Size", y = "Estimated Probability of Survival",
#      title = "rcs(Family Group Size, 5) with 95% confidence intervals",
title = "79 travelers, rcs(Family Group Size, 5) with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the additive model, mod2 (graphic not used in manuscript)
plot_3d_mod2 <- plot_ly(data = subset(pred2, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred2, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
#         title = "Age + Sex + rcs(Family_Group_Size,5)")
       title = "79 travelers, Age + Sex + rcs(Family_Group_Size,5)")
plot_3d_mod2

# Plot the results for the rcs(Age,3) * Sex * rcs(Family_Group-Size,5) 3-d model
plot_3d_mod3 <- plot_ly(data = subset(pred3, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred3, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
#          title = "rcs(Age,3) * Sex + rcs(Family_Group_Size,5)")
      title = "79 travelers, rcs(Age,3) * Sex + rcs(Family_Group_Size,5)")
plot_3d_mod3

##### GAM analysis of the 79 Traveler #####
######## k-fold cross validation to find optimal k for the GAM

# Set a seed for reproducibility
set.seed(123)

# Create a function to perform cross-validation for a given k
cv_gam <- function(k_value){
  
  # Number of folds
  n_folds <- 5
  
  # Create a vector to store CV errors
  cv_errors <- numeric(n_folds)
  
  # Create fold indices
  fold_indices <- sample(rep(1:n_folds, length.out = nrow(Donner)))
  
  for(i in 1:n_folds){
    
    # Split the data into training and test sets
    training_data <- Donner[fold_indices != i, ]
    test_data <- Donner[fold_indices == i, ]
    
    # Fit the GAM model on the training data
    gam_model <- gam(Status ~ s(Age, by = as.numeric(Sex == "Male"), k=k_value) +
                       s(Age, by = as.numeric(Sex == "Female"), k=k_value),
                     data = training_data, family = binomial())
    
    # Predict on the test data
    predictions <- predict(gam_model, newdata = test_data, type = "response")
    
    # Compute the binomial deviance (log loss) for the current fold and store it
    cv_errors[i] <- -2 * sum(test_data$Status * log(predictions) + (1 - test_data$Status) * log(1 - predictions))
    
  }
  
  # Return the mean CV error
  return(mean(cv_errors))
}

# Test a range of k values
# This statement produced many warnings, reduced max to 8 and increased min to 3
# k_values <- 2:10
k_values <- 3:8
cv_results <- sapply(k_values, cv_gam)
cv_results
# Determine the optimal k based on minimum CV error
optimal_k <- k_values[which.min(cv_results)]

print(optimal_k)

## optimal_k is 2, so redo the model

mod_gam_k2 <- gam(Status ~ s(Age, by = as.numeric(Sex == "Male"), k=2) +
                    s(Age, by = as.numeric(Sex == "Female"), k=2), 
                  data = Donner, family = binomial())
mod_gam_k2

# Check the summary
summary(mod_gam_k2)

# Create a prediction dataset
new_age <- seq(min(Donner$Age), max(Donner$Age), length.out = 100)
new_sex <- unique(Donner$Sex)
pred_data <- expand.grid(Age = new_age, Sex = new_sex)

# Predict using the GAM model
predictions <- predict(mod_gam_k2, newdata = pred_data, type = "link", se.fit = TRUE)
pred_data$fit <- plogis(predictions$fit)
pred_data$lower <- plogis(predictions$fit - 1.96 * predictions$se.fit)
pred_data$upper <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# plot the GAM analysis (Figure 3 in manuscript)
# Adjust y-values for jittered points
Donner$AdjustedStatus <- ifelse(Donner$Status == 0, -0.03, 1.03)
set.seed(8)
ggplot(pred_data, aes(x = Age, y = fit, color = Sex, shape = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
#       title = "GAM (k=2) with Age smooths by Sex with 95% confidence intervals") +
 title = "79 travelers, GAM (k=2) with Age smooths by Sex with 95% confidence intervals") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "pink2"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  guides(color = guide_legend(override.aes = list(shape = c(17, 16))))

### Family Group Size, code from GPT-4 #########################################

# Define the k-fold cross-validation manually for GAM:
set.seed(123) # for reproducibility
folds <- createFolds(Donner$Status, k = 10)

cv_errors <- data.frame(k = integer(), RMSE = numeric(), Rsquared = numeric())

for (k_value in 3:8) {
  rmse_values <- numeric()
  rsq_values <- numeric()
  
  # Attempt to fit models and catch errors
  for (fold in folds) {
    train_data <- Donner[-fold, ]
    test_data <- Donner[fold, ]
    success <- TRUE
    
    model <- tryCatch({
      gam(Status ~ s(Family_Group_Size, k = k_value), data = train_data, family = binomial())
    }, error = function(e) {
      message(paste("Error with k =", k_value, ":", e$message))
      success <- FALSE
      return(NULL)
    }, warning = function(w) {
      message(paste("Warning with k =", k_value, ":", w$message))
    })
    
    if(success && !is.null(model)){
      preds <- predict(model, newdata = test_data, type = "response")
      
      rmse <- sqrt(mean((test_data$Status - preds)^2))
      rsq <- 1 - (sum((test_data$Status - preds)^2) / sum((test_data$Status - mean(test_data$Status))^2))
      
      rmse_values <- c(rmse_values, rmse)
      rsq_values <- c(rsq_values, rsq)
    }
  }
  
  if (length(rmse_values) > 0 && length(rsq_values) > 0) {
    cv_errors <- rbind(cv_errors, data.frame(k = k_value, RMSE = mean(rmse_values), Rsquared = mean(rsq_values)))
  }
}

# Check the results
cv_errors

# The value of k with the smallest RMSE would be the optimal choice
best_k <- cv_errors[which.min(cv_errors$RMSE),]$k
best_k

# Fit the GAM with optimal k

mod_gam <- gam(Status ~ s(Family_Group_Size, k = best_k), data = Donner, family = binomial())

summary(mod_gam)

# Generate ggplot graphics for GAM
# Create a prediction dataset
new_family_group_size <- seq(min(Donner$Family_Group_Size), 
                             max(Donner$Family_Group_Size), length.out = 100)
# Create prediction data
# Predict on the link (logit) scale
link_predictions <- predict(mod_gam, newdata = data.frame(Family_Group_Size = new_family_group_size), type = "link", se.fit = TRUE)

# Compute confidence intervals on the link scale
link_lower <- link_predictions$fit - 1.96 * link_predictions$se.fit
link_upper <- link_predictions$fit + 1.96 * link_predictions$se.fit

# Back-transform to the probability scale
pred_data <- data.frame(
  Family_Group_Size = new_family_group_size,
  fit = plogis(link_predictions$fit),
  lower = plogis(link_lower),
  upper = plogis(link_upper)
)

# Plot the GAM (k = best_k = 5) model.
ggplot(pred_data, aes(x = Family_Group_Size, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Family_Group_Size, y = AdjustedStatus,
                                 shape = Sex, color = Sex), width = 0.3, height = 0.03,
              size = 1.5) +
  labs(x = "Family Group Size", y = "Estimated Probability of Survival",
#      title = paste0("GAM (k=", best_k, ") with Family Group Size with 95% confidence intervals"),
 title = paste0("79 travelers, GAM (k=", best_k, ") with Family Group Size with 95% confidence intervals"),
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

##### Cox Proportional Hazards Model

##    Code & analysis from Harrell (2015) Chapter 20 (p 475-519)
#    Couldn't get Harrell code to work

Donner$Death <- abs(Donner$Status-1)
# Code from GPT-4
# 1. Load necessary packages
# install.packages("survival")
# install.packages("survminer")

# Assuming the data frame is already loaded as `Donner`
# 2. Fit the Cox proportional hazards model
cox_model <- coxph(Surv(Survival_Time, Death) ~ Sex, data = Donner)

# Print the summary of the model
print(summary(cox_model))

# Check the proportional hazards assumption with 
ph_test <- cox.zph(cox_model)

# Print the results
print(ph_test)

# Plot the Schoenfeld residuals
plot(ph_test)

# 3. Visualize the results 
# Kaplan-Meier survival curve
km_fit <- survfit(Surv(Survival_Time, Death) ~ Sex, data = Donner)
ggsurvplot(km_fit, data = Donner, risk.table = TRUE,
           pval = TRUE, pval.coord = c(0.8, 0.25),
           legend = "right", legend.title = "Sex",
           legend.labs = c("Female", "Male"), conf.int=TRUE,
           risk.table.y.text = TRUE)

# Hazard ratio plot
# Convert model coefficients to hazard ratios
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))

# Create a data frame for plotting
df <- data.frame(
  Variable = names(hr),
  HR = hr,
  lower = ci[, 1],
  upper = ci[, 2]
)

# Organize data for forestplot
labeltext <- list(
  c("", df$Variable),
  c("HR", as.character(round(df$HR, 2))),
  c("Lower 95% CI", as.character(round(df$lower, 2))),
  c("Upper 95% CI", as.character(round(df$upper, 2)))
)

# Plot using forestplot
forestplot(
  labeltext = labeltext,
  graph.pos = 3,
  mean = c(NA, df$HR),
  lower = c(NA, df$lower),
  upper = c(NA, df$upper),
  xlog = TRUE, # because hazard ratios are typically plotted on a log scale
  clip = c(0.5, 2), # you can adjust these values if necessary
  xticks = c(0.5, 1, 2),
  title = "Hazard Ratios (95% CI)"
)

# 4. Teamster/Servant analyses

cox_model_2 <- coxph(Surv(Survival_Time, Death) ~ Teamster_Hired_Hands, data = Donner)

# Print the summary of the model
print(summary(cox_model_2))

# Check the proportional hazards assumption
ph_test_2 <- cox.zph(cox_model_2)

# Print the results
print(ph_test_2)

# Plot the Schoenfeld residuals
plot(ph_test_2)

# 5. Visualize the results 
# Kaplan-Meier survival curve
km_fit <- survfit(Surv(Survival_Time, Death) ~ Teamster_Hired_Hands, data = Donner)
ggsurvplot(km_fit, data = Donner, risk.table = TRUE,
           pval = TRUE, pval.coord = c(0.8, 0.25),
           legend = "right", legend.title = "Sex",
           legend.labs = c("Family Member", "Teamster"), conf.int=TRUE,
           risk.table.y.text = TRUE)

# Hazard ratio plot
# Convert model coefficients to hazard ratios
hr <- exp(coef(cox_model_2))
ci <- exp(confint(cox_model_2))

# Create a data frame for plotting
df <- data.frame(
  Variable = names(hr),
  HR = hr,
  lower = ci[, 1],
  upper = ci[, 2]
)

# Organize data for forestplot
labeltext <- list(
  c("", df$Variable),
  c("HR", as.character(round(df$HR, 2))),
  c("Lower 95% CI", as.character(round(df$lower, 2))),
  c("Upper 95% CI", as.character(round(df$upper, 2)))
)

# Plot using forestplot
forestplot(
  labeltext = labeltext,
  graph.pos = 3,
  mean = c(NA, df$HR),
  lower = c(NA, df$lower),
  upper = c(NA, df$upper),
  xlog = TRUE, # because hazard ratios are typically plotted on a log scale
  clip = c(0.5, 2), # you can adjust these values if necessary
  xticks = c(0.5, 1, 2.5),
  title = "Hazard Ratios (95% CI)"
)

# The proportional hazards assumption is clearly violated, so GPT-4 suggested
# an alternative:
Donner$SurvTime_Teamster <- with(Donner, Survival_Time * Teamster_Hired_Hands)
cox_time_dep <- coxph(Surv(Survival_Time, Death) ~ Teamster_Hired_Hands + SurvTime_Teamster, data = Donner)
summary(cox_time_dep)


# Print the summary of the model
print(summary(cox_time_dep))

# Check the proportional hazards assumption

# 5. Visualize the results 
# Kaplan-Meier survival curve
km_fit <- survfit(Surv(Survival_Time, Death) ~ Teamster_Hired_Hands, data = Donner)
ggsurvplot(km_fit, data = Donner, risk.table = TRUE,
           pval = TRUE, pval.coord = c(0.8, 0.25),
           legend = "right", legend.title = "Sex",
           legend.labs = c("Family Member", "Teamster"), conf.int=TRUE,
           risk.table.y.text = TRUE)

# Hazard ratio plot
# Convert model coefficients to hazard ratios
hr <- exp(coef(cox_time_dep))
ci <- exp(confint(cox_time_dep))

# Create a data frame for plotting
df <- data.frame(
  Variable = names(hr),
  HR = hr,
  lower = ci[, 1],
  upper = ci[, 2]
)

# Organize data for forestplot
labeltext <- list(
  c("", df$Variable),
  c("HR", as.character(round(df$HR, 2))),
  c("Lower 95% CI", as.character(round(df$lower, 2))),
  c("Upper 95% CI", as.character(round(df$upper, 2)))
)

# Plot using forestplot
forestplot(
  labeltext = labeltext,
  graph.pos = 3,
  mean = c(NA, df$HR),
  lower = c(NA, df$lower),
  upper = c(NA, df$upper),
  xlog = TRUE, # because hazard ratios are typically plotted on a log scale
  clip = c(0.5, 2), # you can adjust these values if necessary
  xticks = c(0.5, 1, 2.5),
  title = "Hazard Ratios (95% CI)"
)

### Forlorn Hope Fisher's exact test:
# Create the data matrix
# The matrix format is:
#        Survived Died
# Women      5     0
# Men        2     8

donner_matrix <- matrix(c(5, 2, 0, 8), nrow=2, byrow=TRUE)
colnames(donner_matrix) <- c("Survived", "Died")
rownames(donner_matrix) <- c("Women", "Men")
print(donner_matrix)

# Fisher's exact test
test_result <- fisher.test(donner_matrix)

# Print the results
print(test_result)

# Odds ratio
print(paste("Odds Ratio:", test_result$estimate))

# Confidence interval for the odds ratio
print(paste("95% Confidence Interval:", round(test_result$conf.int[1], 3), "to", round(test_result$conf.int[2], 3)))

sink()   # Optional Turn off redirection
 