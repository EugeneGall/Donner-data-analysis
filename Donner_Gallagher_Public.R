# Donner_Gallagher_Public
# Written by Eugene Gallagher, Associate Professor, School for the Environment
# UMass Boston, Eugene.Gallagher@umb.edu 7/10/23, last revised 8/23/23
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
# Grayson, D. K. (2018), Sex and Death on the Western Emigrant Trail: The
#    Biology of Three American Tragedies. Salt Lake City: University of Utah
#    Press.
# Ramsey, F. L. and D. W. Schafer. 2013. The Statistical Sleuth: a Course in
#    Methods of Data Analysis, 3rd Edition. Brooks/Cole Cengage Learning,
#    Boston MA, 760 pp.
# Rarick, E. 2008. Desperate Passage: The Donner Party’s Perilous Journey West.
#    Oxford University Press, Oxford. 304 pp.
# Stewart, G. E. 1960. Ordeal by Hunger: the Ordeal of the Donner Party.
#    Houghton Mifflin, Boston.392 pp.
# 
# Approach: 1) add the under 15 Age data to the Donner data from Statistical
# Sleuth 3rd edition, 2) change age of Patrick Breen to 51 (Grayson, 1994, 
# p 155). 3) Grayson (1990, 1997, 2018) argued Family Group Size, Age, and Sex
# control survival and Rarick emphasized the poor survivorship of teamsters and
# servants, combined here as Employees. This R code will analyze the effects of
# all four variables.
# 4) Reviewed above books to find death dates for surviving travelers and
# performed survival analyses confirming Grayson (1997) on death timing.
# Updated Demographic data in Donner.csv to conform to Grayson (2018)
# 5) Renumbered the models and ran both glm and Glm on each model, the latter
#    to get the data on rcs effect sizes and to plot the data Harrell-style.
# 6) It appears that the 3-d plots are badly overfit due to the AIC optimization
#    I'll reduce the number of knots and then have GPT-4 develop a k-fold
#    cross validation to objectively choose appropriate knot numbers.

# Used data imputation to fill in the missing age for Mr. Wolfinger
# Have R determine family size by the numbers of individuals with the same
# last name (not used in manuscript) but used Grayson's (1990) Family Group
# Size, which incorporates information from Stewart's (1960) roster on who was
# traveling with each Family Travel Group.

# Optional statement needed to see full output in a word processor:
sink("my_output.txt")   # Optional Redirect output to a file, make sure that
                        # that the file readme.txt is not open in another app.

# Install and load packages
library(boot) # for cv.glm function
library(caret) # For GAM cross-validation
library(forestplot)  # for the Cox PH plots
library(gmodels)  # For Fisher tests of the Forlorn Hope data
library(grid) # for the Risk Table beneath the KM plots
library(gridExtra) # for the Risk Table beneath the KM plots
library(gtable) # for the Risk Table beneath the KM plots
library(mgcv) # Wood's package for GAMs
library(patchwork) # for the Risk Table below the KM plots.
library(plotly) # for 3-d graphics
library(rms)    # Harrell's regression modeling strategies pacakge
library(survival)  # for Cox Proportional Hazards Modeling
library(survminer) # for Cox Proportional Hazards Modeling
library(tidyverse) # contains dplyr and ggplot2

# Read the data from Gallagher's github site (public has read, not edit access)

Donner <- read.csv("https://raw.githubusercontent.com/EugeneGall/donner-data-analysis/main/Donner.csv")

# Calculate family size based on Last_Name, but Family_Group_Size from Grayson
# (1990, 2018) Table will be used in this code's analyses.
Donner$Family_Size <- as.integer(ave(Donner$Last_Name, Donner$Last_Name, FUN =length))

str(Donner)

# This produced slightly different Family Sizes than Grayson (1990 Table 1, 
# 2018 Table 2.1. Used Grayson (1990, 2018) Family Group Sizes with included
# family links using data from Stewart's (1960, p 363-364) Donner Party Roster.

# Convert Status to binary
Donner$Status <- ifelse(Donner$Status == "Survived", 1, 0)

# Calculate the Survival_Time which is Death Date or last rescue (Lewis 
# Keseberg 4/29/1847) - First_Snow, October 28 1846.
# Convert the character strings to Date objects
Donner$First_Snow_Date <- as.Date(Donner$First_Snow, format="%Y-%m-%d")
Donner$Last_Date_Date <- as.Date(Donner$Last_Date, format="%Y-%m-%d")

# Calculate the difference in days for individual survivorship
Donner$Survival_Time <- as.numeric(Donner$Last_Date_Date - Donner$First_Snow_Date)

# View the first few rows of the data frame to confirm the results
head(Donner)

# Impute missing Age values for Mr. Wolfinger, using median imputation
# Rarick (2008, xi) provides Dora Wolfinger's age as 20. This section will
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

### Delete 8 cases for individuals who died before the first Snowstorm
# on 1846-10-28. This will reduce the number of cases to 87 - 8 = 79

Donner <- Donner %>%
  filter(Delete != 1)
# View the structure to confirm the filtering
str(Donner)

# Prepare data for rms, Harrell's 'Regression Modeling Strategies' These
# two statements are required for the summary of the Harrell's Glm
ddist <- datadist(Donner)
options(datadist = "ddist")

# Fit the model with Age, Sex, and Grayson's (1990, 2018) Family_Group_Size
# In deciding on knot size, I used the rule that the minimum knot size should
# be chosen (rcs has a min of 3) unless there is another larger number of knots
# with an AIC which is more than than 4 AIC units lower (Burnham & Anderson
# 2004, p 271)

# Two methods will be coded for finding the appropriate number of knots:
# 1) Brute force by manually plugging in different knot numbers in Glm
# 2) An optimization routine based on AIC

##### First approach: brute force fitting to find minimal AICs.
# 3 is the minimum knot size permitted with Harrell's rcs function

# Fit Model 1:    Status ~ rcs(Age, 3) * Sex
Mod1  <- glm(Status ~ rcs(Age,3) * Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)

# Need Harrell's Glm for effect sizes and graphic with rms::Predict
mod1  <- Glm(Status ~ rcs(Age,3) * Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)
# Where are the knots located?
k <- attr (rcs(Donner$Age,3) ,'parms')
k
# knots located at Ages 1, 14, and 46
# Null model
mod_null <- glm(Status ~ 1, data = Donner, family = binomial())
# Additive model
mod2  <- Glm(Status ~ rcs(Age,3) + Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)
Mod2  <- glm(Status ~ rcs(Age,3) + Sex, data = Donner, family = binomial(), x = TRUE, y = TRUE)
# Test the models
chi_sq_n1 <- anova(mod_null, Mod1, test="Chisq")
# Print the results
print(chi_sq_n1)
# Summary and ANOVA for Mod1
# Test whether the interaction model adds value.
chi_sq_21 <- anova(Mod2, Mod1, test="Chisq")
# Print the results
print(chi_sq_21)
# Summaries and effect sizes of the Age(rcs4) * Sex model
summary(Mod1)
anova(Mod1)
AIC(Mod1)

#Summary of Brute force fitting of AICs for interaction model
#   mod1     AIC
# 3 knots  95.40752 * Chosen because within 4 of lowest AIC
# 4 knots  92.3582 ** Not chosen because only 3.04932 less than 3 knot solution
# 5 knots  92.88625 
# 6 knots Apparently Singular Matrix, no estimation possible

# Fit Family_Group_Size fit with base glm and Harrell's rms::Glm
Mod3 <- glm(Status ~ rcs(Family_Group_Size,5), data = Donner, 
            family = binomial(), x = TRUE, y = TRUE)
mod3 <- Glm(Status ~ rcs(Family_Group_Size,5), data = Donner, 
            family = binomial(), x = TRUE, y = TRUE)
# Null model
# Test the models
chi_sq_n2 <- anova(mod_null, Mod3, test="Chisq")
# Print the results
print(chi_sq_n2)
# Summary and ANOVA for mod
summary(Mod3)
anova(Mod3)
summary(mod3)  # For effect size
AIC(mod3)
AIC(Mod3)

mod4 <- Glm(Status ~ Age + Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
Mod4 <- glm(Status ~ Age + Sex + rcs(Family_Group_Size,4), data = Donner, family = binomial(), x = TRUE, y = TRUE)
# for mod5, knot sizes of 3 & 5 determined by k-fold cross validation, code below
mod5 <- Glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
Mod5 <- glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)

##### Second approach: GPT4 Streamlined the above code block with functions:

# Function to fit models with one rcs term, updated on 8/17/23 by GPT-4
fit_best_rcs <- function(formula_template, data, knot_range) {
  
  aic_values <- numeric(length(knot_range))
  models <- vector("list", length(knot_range))
  
  # Calculate AIC for each knot size
  for (k in knot_range) {
    formula <- as.formula(gsub("KNOTS", k, formula_template))
    model <- glm(formula, data = data, family = binomial())
    aic_values[k - knot_range[1] + 1] <- AIC(model)
    models[[k - knot_range[1] + 1]] <- model
  }
  
  min_aic <- min(aic_values)
  best_knots <- knot_range[which.min(aic_values)]
  
  # Iterate through the aic_values and apply the rules
  for (i in seq_along(aic_values)) {
    if (aic_values[i] - min_aic <= 4 && knot_range[i] < best_knots) {
      best_knots <- knot_range[i]
    }
  }
  
  return(list(model = models[[which(knot_range == best_knots)]], knots = best_knots, results = setNames(aic_values, knot_range)))
}

# Function to fit model with two rcs terms, updated by GPT-4 on 8/17/23
fit_best_double_rcs <- function(formula, data, knot_range, family = binomial()) {
  best_aic <- Inf
  best_k_age <- max(knot_range)
  best_k_fgs <- max(knot_range)
  aic_matrix <- matrix(NA, nrow = length(knot_range), ncol = length(knot_range), dimnames = list(knot_range, knot_range))
  
  for (k_age in knot_range) {
    for (k_fgs in knot_range) {
      tryCatch({
        current_formula <- as.formula(
          gsub("knots_age", as.character(k_age), 
               gsub("knots_fgs", as.character(k_fgs), formula)))
        current_model <- glm(current_formula, data = data, family = family)
        current_aic <- AIC(current_model)
        
        aic_matrix[as.character(k_age), as.character(k_fgs)] <- current_aic
        
        # Update if current AIC is smaller than the best so far
        if (current_aic < best_aic) {
          best_aic <- current_aic
          best_k_age <- k_age
          best_k_fgs <- k_fgs
        }
      }, error = function(e) {}) # Handle potential errors
    }
  }
  
  min_aic <- best_aic  # best AIC among all combinations
  best_sum_knots <- best_k_age + best_k_fgs  # initialize with the sum of best knots
  
  # Find combinations with the smallest sum of knots but AIC within 4 of min_aic
  for (k_age in knot_range) {
    for (k_fgs in knot_range) {
      current_sum_knots <- k_age + k_fgs
      if ((aic_matrix[as.character(k_age), as.character(k_fgs)] - min_aic <= 4) && current_sum_knots < best_sum_knots) {
        best_k_age <- k_age
        best_k_fgs <- k_fgs
        best_sum_knots <- current_sum_knots
      }
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
# Note that 3 is the minimum knot size for rcs, dozens of warnings with knots>6
knot_range <- 3:7
mod1_results <- fit_best_rcs("Status ~ rcs(Age, KNOTS) * Sex", Donner, knot_range)
mod3_results <- fit_best_rcs("Status ~ rcs(Family_Group_Size, KNOTS)", Donner, knot_range)
mod4_results <- fit_best_rcs("Status ~ rcs(Family_Group_Size, KNOTS) * Sex", Donner, knot_range)
mod5_results <- fit_best_double_rcs("Status ~ rcs(Age, knots_age) * Sex + rcs(Family_Group_Size, knots_fgs)", Donner, knot_range)

# Printing the details
print_details(mod1_results, "mod1")
print_details(mod3_results, "mod3")
print_details(mod4_results, "mod4")
print_details(mod5_results, "mod5")

###### End of GPT-4 AIC rcs optimization code ######

# Generate 4 ggplot2 graphics
# Define new levels for the predictors
new_age <- seq(min(Donner$Age), max(Donner$Age), length.out = 100)
new_family_group_size <- seq(min(Donner$Family_Group_Size), 
                             max(Donner$Family_Group_Size), length.out = 100)
new_sex <- unique(Donner$Sex)

# Make predictions on the link scale (logit scale)
# The Predict function require fit by Glm not glm
link_pred1 <- Predict(mod1, Age = new_age, Sex = new_sex)
link_pred3 <- Predict(mod3, Family_Group_Size = new_family_group_size)
link_pred4 <- Predict(mod4, Age = new_age, Sex = new_sex,
                      Family_Group_Size = new_family_group_size)
link_pred5 <- Predict(mod5, Age = new_age, Sex = new_sex, 
                      Family_Group_Size = new_family_group_size)

# Transform predictions back to the original scale (probability scale)
pred1 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  fit = plogis(link_pred1$yhat),
  lower = plogis(link_pred1$lower),
  upper = plogis(link_pred1$upper)
)

pred3 <- data.frame(
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = new_family_group_size,
  fit = plogis(link_pred3$yhat),
  lower = plogis(link_pred3$lower),
  upper = plogis(link_pred3$upper)
)

# Transform predictions back to the original scale (probability scale)
pred4 <- data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = rep(new_family_group_size, each = length(new_family_group_size)),
  fit = plogis(link_pred4$yhat),
  lower = plogis(link_pred4$lower),
  upper = plogis(link_pred4$upper)
)

pred5 <-  data.frame(
  Age = rep(new_age, times = length(new_sex)),
  Sex = rep(new_sex, each = length(new_age)),
  Family_Group_Size = rep(new_family_group_size, each = length(new_family_group_size)),
  fit = plogis(link_pred5$yhat),
  lower = plogis(link_pred5$lower),
  upper = plogis(link_pred5$upper)
)
# Adjust y-values for jittered points
Donner$AdjustedStatus <- ifelse(Donner$Status == 0, -0.03, 1.03)

# Plot the results for Figure 1 rcs(Age,3) * Sex model
set.seed(13) 
ggplot(pred1, aes(x = Age, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus,
            shape = Sex, color = Sex), width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
  title = "Fig. 1. rcs(Age, 3) * Sex with 95% confidence intervals",
        color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the rcs(Family_Group_Size,4) model
ggplot(pred3, aes(x = Family_Group_Size, y = fit, color = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Family_Group_Size, y = AdjustedStatus,
              shape = Sex, color = Sex), width = 0.3, height = 0.03,
              size = 1.5) +
  labs(x = "Family Group Size", y = "Estimated Probability of Survival",
  title = "Fig. 3. rcs(Family Group Size, 5) with 95% confidence intervals",
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

# Plot the results for the additive model, mod2 (graphic not used in manuscript)
plot_3d_mod4 <- plot_ly(data = subset(pred4, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred4, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
#         title = "Age + Sex + rcs(Family_Group_Size,6)")
       title = "Age + Sex + rcs(Family_Group_Size,6)")
plot_3d_mod4

# Plot the results for the rcs(Age,3) * Sex * rcs(Family_Group-Size,5) 3-d model
# Graphic not used in the manuscript (too curvy)
plot_3d_mod5 <- plot_ly(data = subset(pred5, Sex == "Male"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
                        type = "mesh3d", opacity = 0.6, name = "Male", showscale = FALSE) %>%
  add_trace(data = subset(pred5, Sex == "Female"), x = ~Age, y = ~Family_Group_Size, z = ~fit, 
            type = "mesh3d", opacity = 0.6, name = "Female", showscale = FALSE) %>%
  layout(scene = list(zaxis = list(range = c(0, 1)),
                      xaxis = list(title = "Age"),
                      yaxis = list(title = "Family Group Size"),
                      zaxis = list(title = "Estimated Probability of Survival")),
                      title = "rcs(Age,3) * Sex + rcs(Family_Group_Size,5)")
plot_3d_mod5

### 8/20/23 addendum. Fitting rcs knots with k-fold cross validation ##########
# The AIC approach produces a badly overfit mod5 showing
# a 3-d plot that was far too flexible. I sent a 3-p prompt on the problem to
# GPT-4 and it provided this solution

cv_rcs <- function(k_age, k_fgs, data, n_folds = 10){
  
  # Create a vector to store CV errors
  cv_errors <- numeric(n_folds)
  
  # Create fold indices
  fold_indices <- sample(rep(1:n_folds, length.out = nrow(data)))
  
  for(i in 1:n_folds){
    
    # Split the data into training and test sets
    training_data <- data[fold_indices != i, ]
    test_data <- data[fold_indices == i, ]
    
    # Fit the GLM model on the training data with rcs terms
    glm_model <- glm(Status ~ rcs(Age, k_age) * Sex + rcs(Family_Group_Size, k_fgs), 
                     data = training_data, family = binomial())
    
    # Predict on the test data
    predictions <- predict(glm_model, newdata = test_data, type = "response")
    
    # Compute the binomial deviance (log loss) for the current fold and store it
    cv_errors[i] <- -2 * sum(test_data$Status * log(predictions) + (1 - test_data$Status) * log(1 - predictions))
    
  }
  
  # Return the mean CV error
  return(mean(cv_errors))
}

# Perform k-fold CV for a range of knot sizes
k_values <- 3:7
cv_results_matrix <- matrix(NA, nrow = length(k_values), ncol = length(k_values),
                            dimnames = list(k_values, k_values))

for(k_age in k_values){
  for(k_fgs in k_values){
    cv_results_matrix[as.character(k_age), as.character(k_fgs)] <- cv_rcs(k_age, k_fgs, Donner)
  }
}

# Determine the optimal k_age and k_fgs based on minimum CV error
optimal_indices <- which(cv_results_matrix == min(cv_results_matrix), arr.ind = TRUE)
optimal_k_age <- as.numeric(rownames(cv_results_matrix)[optimal_indices[1,1]])
optimal_k_fgs <- as.numeric(colnames(cv_results_matrix)[optimal_indices[1,2]])

print(paste("Optimal k for Age: ", optimal_k_age))
print(paste("Optimal k for Family_Group_Size: ", optimal_k_fgs))

# Fit the optimal model
best_glm_model <- glm(Status ~ rcs(Age, optimal_k_age) * Sex + 
                      rcs(Family_Group_Size, optimal_k_fgs), 
                      data = Donner, family = binomial())
summary(best_glm_model)
## This optimization produces Age: 3 knots, Family Group Size 5 knots.

# Summary and ANOVA for mod5 with knot sizes determined by k-fold cross validation
# Just repeating thise calculations for clarity.
mod5 <- Glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
Mod5 <- glm(Status ~ rcs(Age,3) * Sex + rcs(Family_Group_Size,5), data = Donner, family = binomial(), x = TRUE, y = TRUE)
summary(Mod5)
anova(Mod5)
summary(mod5)  # For effect size
AIC(mod5)

### end of 8/20/23 GPT-4 code. #################################################

##### GAM analysis of the 79 Travelers #####
######## k-fold cross validation to find optimal k for the GAMs

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

# plot the GAM analysis for Sex (Figure 2 in manuscript)
# Adjust y-values for jittered points
Donner$AdjustedStatus <- ifelse(Donner$Status == 0, -0.03, 1.03)

# Plot Fig. 2, GAM analysis of Sex
set.seed(3)
ggplot(pred_data, aes(x = Age, y = fit, color = Sex, shape = Sex)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = Donner, aes(x = Age, y = AdjustedStatus),
              width = 0.3, height = 0.03, size = 1.5) +
  labs(x = "Age (Years)", y = "Estimated Probability of Survival",
   title = "Fig. 2. GAM (k=2) with Age smooths by Sex with 95% confidence intervals") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "pink2"), name = "Sex") +
  scale_shape_manual(values = c("Male" = 17, "Female" = 16), name = "Sex") + 
  guides(color = guide_legend(override.aes = list(shape = c(17, 16))))

### GAMs for Family Group Size, code from GPT-4 ################################

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
# title = paste0("GAM (k=", best_k, ") with Family Group Size with 95% confidence intervals"),
 title = paste0("Fig. 4. GAM (k=", best_k, ") with Family Group Size with 95% confidence intervals"),
       color = "Sex") +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.06, 1.06), breaks = seq(0, 1, 0.2))

##### Cox Proportional Hazards Model for Sex

Donner$Death <- abs(Donner$Status-1)
# Code initially from GPT-4

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

# Fit Kaplan-Meier survival curves for Sex
km_fit <- survfit(Surv(Survival_Time, Death) ~ Sex, data = Donner)

# Define the days of arrival of the relief parties
relief_days <- c(83, 124, 136, 171)
colors <- c("red", "blue", "magenta", "purple")
linetypes <- c("dashed", "dotted", "dotdash", "longdash")
labels <- c("1st", "2nd", "3rd", "4th")

# Plot Fig. 6 the KM curve & Risk Table for Sex

# Generate KM plot without the risk table
p1 <- ggsurvplot(km_fit, data = Donner, risk.table = FALSE,
                 pval = TRUE, pval.coord = c(0.8, 0.25),
                 legend = "right", legend.title = "Sex",
                 title = "Fig. 6. Kaplan-Meier Survivorship by Sex",
                 legend.labs = c("Female", "Male"), conf.int = TRUE,
                 break.x.by = 20)

# Add the relief party vertical lines and annotations
for (i in 1:4) {
  p1$plot <- p1$plot + 
    annotate("segment", x = relief_days[i], xend = relief_days[i], 
             y = 0, yend = 1, colour = colors[i], linetype = linetypes[i]) +
    annotate("text", x = relief_days[i], y = 0.2, 
             label = labels[i], angle = 90, vjust = 0, color = colors[i])
}

# Generate risk table separately
p2 <- ggsurvplot(km_fit, data = Donner, risk.table = TRUE,
                 pval = FALSE, legend = "none",
                 break.x.by = 20)

# Use patchwork to combine the two plots

# Use patchwork to combine the two plots with specified heights
combined_plot <- p1$plot / p2$table + plot_layout(heights = c(5, 1))

# Print the combined plot
combined_plot

#### Hazard ratio plot
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
  title = "Hazard Ratios by Sex (95% CI)"
)

# 4. Employee (Teamsters & 2 Servants) Survival Analysis
# 4.1 Cox PH model

cox_model_2 <- coxph(Surv(Survival_Time, Death) ~ Employee, data = Donner)

# Print the summary of the model
print(summary(cox_model_2))

# Check the proportional hazards assumption
ph_test_2 <- cox.zph(cox_model_2)

# Print the results
print(ph_test_2)
# Some evidence that the constant proportional hazard assumption is violated.
# p = 0.041 for Employee and Global chisq 4.17 with 1 df

# Plot the Schoenfeld residuals
plot(ph_test_2)

# 4.2 Kaplan-Meier survival curve for Employment

km_fit <- survfit(Surv(Survival_Time, Death) ~ Employee, data = Donner)

# Define the days of arrival of the relief parties
relief_days <- c(83, 124, 136, 171)
colors <- c("red", "blue", "magenta", "purple")
linetypes <- c("dashed", "dotted", "dotdash", "longdash")
labels <- c("1st", "2nd", "3rd", "4th")

# Generate Fig. 7 KM plot for employment with the risk table
p1 <- ggsurvplot(km_fit, data = Donner, risk.table = FALSE,
                 pval = TRUE, pval.coord = c(0.8, 0.25),
                 legend = "right", legend.title = "Sex",
                 title = "Fig. 7. Kaplan-Meier Survivorship by Employment",
                 legend.labs = c("Family Member", "Employee"), conf.int = TRUE,
                 break.x.by = 20)

# Add the relief party vertical lines and annotations
for (i in 1:4) {
  p1$plot <- p1$plot + 
    annotate("segment", x = relief_days[i], xend = relief_days[i], 
             y = 0, yend = 1, colour = colors[i], linetype = linetypes[i]) +
    annotate("text", x = relief_days[i], y = 0.05, 
             label = labels[i], angle = 90, vjust = 0, color = colors[i])
}

# Generate risk table separately
p2 <- ggsurvplot(km_fit, data = Donner, risk.table = TRUE,
                 pval = FALSE, legend = "none",
                 break.x.by = 20)

# Use patchwork to combine the two plots

# Use patchwork to combine the two plots with specified heights
combined_plot <- p1$plot / p2$table + plot_layout(heights = c(5, 1))

# Print the combined plot
combined_plot

# GPT-4 suggested that I try stratified KM plots to show the change in
# Proportional Hazard ratio with time. Here is the codefrom GPT-4

# Split the Survival_Time into intervals

Donner$TimeGroup <- cut(Donner$Survival_Time, breaks = c(0, 50, 100, max(Donner$Survival_Time)), 
                        labels = c("Early", "Middle", "Late"), include.lowest = TRUE)

# KM fit using Employee and stratified by TimeGroup
km_fit_time <- survfit(Surv(Survival_Time, Death) ~ Employee + strata(TimeGroup), data = Donner)

# Plot Stratified Kaplan_Meier curves.
ggsurvplot(km_fit_time, data = Donner, risk.table = FALSE, 
           legend = "right", legend.title = "Sex",
           title = "Fig. 8. Kaplan-Meir Employee Survivorship Stratified by Time Intervals", conf.int = TRUE,
           break.x.by = 20)

# 4.3 Calculate and plot the hazard ratio for Employment

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
Donner$SurvTime_Employee <- with(Donner, Survival_Time * Employee)
cox_time_dep <- coxph(Surv(Survival_Time, Death) ~ Employee + SurvTime_Employee, data = Donner)
summary(cox_time_dep)

# Print the summary of the model
print(summary(cox_time_dep))

# Check the proportional hazards assumption
ph_test_3 <- cox.zph(cox_time_dep)

# Print the results
print(ph_test_3)
# Evidence for the violation of the constant hazard ratio assumption for Employee
# But not for Global

# Plot the Schoenfeld residuals
plot(ph_test_3)
# Fairly level distribution for Betas for SurvTime_Employees

### 6. Forlorn Hope Fisher's exact test:
# Create the data matrix
# The matrix format is:
#        Survived Died
# Women      5     0
# Men        2     8

donner_matrix <- matrix(c(5, 0, 2, 8), nrow=2, byrow=TRUE)
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

# Produces an SPSS-like table.  Note that with the Essex data, warnings are
# issued because the minimum sample sizes for 2 x 2 have expectations less
# than 5. The proportions and Fisher's exact test are correct, hence CrossTable
CrossTable(donner_matrix,digits=3,fisher = TRUE, chisq = TRUE, expected = TRUE,
           format = "SPSS")

sink()   # Optional Turn off redirection
