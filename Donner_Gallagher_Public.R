# Donner_Gallagher_Public
# Written by Eugene.Gallagher@umb.edu 7/10/23, last revised 7/31/23 with help
# from ChatGPT-4
# Analysis of data from Grayson (1990, Table 1; Grayson 1994)
# References
# Burnham, K. P. and Anderson, D.R. (2004), “Multimodal inference:
#   understanding AIC and BIC in model selection,” Sociological Methods and
# Grayson, D. K. 1990. Donner Party Deaths: a Demographic Assessment. J.
#    Anthropological Research 46: 223-242.
# Grayson, D. K. 1994. Differential Mortality and the Donner Party Disaster. 
#    Evolutionary Anthropology 2: 151-159.
# Ramsey, F. L. and D. W. Schafer. 2013. The Statistical Sleuth: a Course in
#    Methods of Data Analysis, 3rd Edition. Brooks/Cole Cengage Learning,
#    Boston MA, 760 pp.
# Stewart, G. E. 1960. Ordeal by Hunger: the Ordeal of the Donner Party.
#    Houghton Mifflin, Boston.
# Approach: added the under 15 data to the Donner data from Statistical Sleuth
#    3rd edition, Changed age of Patrick Breen to 51 (Grayson, 1994, p 155)
# Grayson (1990) argued Family Group Size, Age, and Sex control survival. This
#    R code will analyze the effects of all three variables.

# Use data imputation to fill in the 2 missing ages for the 2 Wolfingers
# Have R determine family size by the numbers of individuals with the same
# last name (not used here). But, also analyze Grayson's (1990) Family Group
# Size, which incorporates information from Stewart's (1960) roster on who was
# traveling with each Family Travel Group.

# Optional Needed to see more output, say in a word processor
# sink("my_output.txt")   # Optional Redirect output to a file, make sure that
                          # that the file readme.txt is not open in another app.

# Install and load packages
library(boot) # for cv.glm function
library(caret) # For GAM cross-validation
library(mgcv)
library(plotly) # for 3-d graphics
library(rms)
library(tidyverse) # contains dplyr and ggplot2


# Read the data from Gallagher's github site (public has read, not edit access)

Donner <- read.csv("https://raw.githubusercontent.com/EugeneGall/donner-data-analysis/main/Donner.csv")


# Calculate family size based on Last_Name, but Family_Group_Size from Grayson
# (1990) Table 1 will be used in this code's analyses.
Donner$Family_Size <- as.integer(ave(Donner$Last_Name, Donner$Last_Name, FUN =length))

str(Donner)

# This produced slightly different Family Sizes than Grayson (1990 Table 1) in
# that he merged some families using data from Stewart's (1960) Roster.

# Convert Status to binary
Donner$Status <- ifelse(Donner$Status == "Survived", 1, 0)

# Impute missing Age values for the two Wolfingers, using median imputation
Donner$Age[is.na(Donner$Age)] <- median(Donner$Age, na.rm = TRUE)
# Both are assigned the median age of 18.
# I had GPT-4 write a dplyr pipe to replace NA's by the medians for each sex,
# but the median age for Females was 13, so I opted to using median imputation
# producing an imputed median age of 18 for both Wolfingers.

# Prepare data for rms, Harrell's 'Regression Modeling Strategies' These
# two statements are required for the summary of the Harrell's Glm
ddist <- datadist(Donner)
options(datadist = "ddist")

# Fit the model with Age, Sex, and Grayson's (1990) Family_Group_Size
# In deciding on knot size, I used the rule that the minimum knot size should
# be chosen (rcs has a min of 3) unless there is another minimum less than 4
# AIC units lower (Burnham & Anderson 2004, p 271)

# Two methods will be shown for finding the appropriate number of knots:

##### First approach: brute force fitting to find minimal AICs. 
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
# 3 knots 104.8571 * Chosen because within 4 of lowest AIC
# 4 knots 104.8716
# 5 knots 103.1862
# 6 knots Singular Matrix, no estimation possible

# Summary and ANOVA for mod1, Family Group Size alone
summary(mod1)
anova(mod1)
AIC(mod1)
# AIC = 106.51
# Very strong nonlinear effect of group size
#   mod2    AIC
# 3 knots 112.2558
# 4 knots 111.5511
# 5 knots 106.5101 *
# 6 knots 109.9358
# 7 knots 105.819
# 8 knots 105.819
# 9 knots 105.819

# Summary and ANOVA for mod2, an additive model, Age, Sex and 
# rcs(Family_Group Size)
summary(mod2)
anova(mod2)
AIC(mod2)
#   mod2    AIC
# 3 knots 110.3669
# 4 knots 111.1772
# 5 knots 102.851 *
# 6 knots 106.1738
# 7 knots 102.8313
# 8 knots 102.8313
# 9 knots 102.8313

# Strong effect of Family_Group_Size and Sex, but not Age (p=0.24)

# Fit a restricted cubic spline regression for Age and Family Group Size with an
# interaction effect
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
       title = "rcs(Age, 3) * Sex with 95% confidence intervals",
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
       title = "rcs(Family Group Size, 5) with 95% confidence intervals",
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
         title = "Age + Sex + rcs(Family_Group_Size,5)")
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
         title = "rcs(Age,3) * Sex + rcs(Family_Group_Size,5)")
plot_3d_mod3

##### GAM analysis of the full data #####

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

# plot the GAM analysis (Figure 4 in manuscript)
# Adjust y-values for jittered points
Donner$AdjustedStatus <- ifelse(Donner$Status == 0, -0.03, 1.03)
set.seed(8)
ggplot(pred_data, aes(x = Age, y = fit, color = Sex, shape = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "GAM (k=2) with Age smooths by Sex with 95% confidence intervals") +
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

for (k_value in 3:9) {
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
       title = paste0("GAM (k=", best_k, ") with Family Group Size with 95% confidence intervals"),
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

####### Redo the Statistical Sleuth Analysis using cases with Age>=15 ##########

# Note Statistical Sleuth (all 3 editions) used only Age 15 and older data and
# dropped the 2 Wolfingers because Grayson (1990) didn't provide ther age. 
# This analysis will analyze those age-pared data.
# Only mod5 & mod6 will be plotted with ggplot

# Create a new data frame with cases where Age is greater than or equal to 15
Donner_15up <- Donner[Donner$Age >= 15, ]

mod5 <- Glm(Status ~ Age + Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod6 <- Glm(Status ~ Age * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod7 <- Glm(Status ~ rcs(Age,3) * Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod8 <- Glm(Status ~ rcs(Age,3) + Sex, data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod9 <- Glm(Status ~ rcs(Age,3),       data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)
mod10 <- Glm(Status ~ rcs(Family_Group_Size,5), data = Donner_15up, family = binomial(), x = TRUE, y = TRUE)

# Summary and ANOVA for mod5
summary(mod5)
mod5$coefficients
anova(mod5)
# Age (p=0.060), Sex (p = 0.027)

# Summary and ANOVA for mod6
summary(mod6)
mod6$coefficients
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
# p = 0.005, so strong Family_Group_Size Effect

# Check whether Age is important in a Wilks drop in deviance chi square test:
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

# Adjust y-values for jittered points
Donner_15up$AdjustedStatus <- ifelse(Donner_15up$Status == 0, -0.03, 1.03)

# Plot the results for the additive model: Not as informative as intxn model
set.seed(8)
ggplot(pred5, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Age, y = AdjustedStatus, 
                                      shape = Sex, color = Sex),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "Age + Sex model with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the interaction model
set.seed(8)
ggplot(pred6, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner_15up, aes(x = Age, y = AdjustedStatus, 
                                      shape = Sex, color = Sex),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
       title = "Age * Sex model with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# New code for a GAM 3-d plot. Sex has to be numeric.########################

# Convert 'Sex' from character to numeric (0 for Male, 1 for Female)
Donner$Sex_numeric <- ifelse(Donner$Sex == "Male", 0, 1)

# Setting up cross-validation
set.seed(123) # for reproducibility
k_values <- c(2, 3, 4, 5, 6) # potential k values
cv_errors <- data.frame()

# k-fold CV

for (k_value in k_values) {
  mod_gam <- gam(Status ~ s(Age, k = k_value, by = Sex_numeric) + s(Family_Group_Size, k = k_value), 
                 data = Donner, family = binomial())
  
  cv_result <- cv.glm(Donner, mod_gam, K = 10) # 10-fold CV
  cv_errors <- rbind(cv_errors, data.frame(k = k_value, CVError = cv_result$delta[1]))
}

# Optimal k value
best_k <- cv_errors[which.min(cv_errors$CVError),]$k
best_k

# Fitting GAM with best k
mod_gam <- gam(Status ~ s(Age, k = best_k, by = Sex_numeric) + s(Family_Group_Size, k = best_k), 
               data = Donner, family = binomial())
summary(mod_gam)
anova(mod_gam)

# Generate predictions
pred_grid <- expand.grid(Age = new_age, 
                         Sex = new_sex, 
                         Family_Group_Size = new_family_group_size)
pred_grid$Sex_numeric <- ifelse(pred_grid$Sex == "Male", 0, 1)
pred_gam <- predict(mod_gam, newdata = pred_grid, type = "response")

# Combine predictions with grid for plotting
pred_data <- cbind(pred_grid, fit = pred_gam)

# 3D plot
plot_3d_gam <- plot_ly(data = subset(pred_data, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                       type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred_data, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
         title = paste0("GAM with k = ", best_k))
plot_3d_gam

# Overall conclusion on the 15-up analysis
# With Age>= 15, there is a poor fit with the rcs(Age,3), but as with the
# Sleuth3 analysis, the Age * Sex interaction is important, as determined by the
# Wilks drop in deviance test (p=0.018), indicating the need for an interaction
# term,
# There is no justification for using a restricted cubic spline fit on the
# age-pared Donner data.

# sink()   # Optional Turn off redirection