library(tidyverse)
library(brms)         
library(broom.mixed)
library(tidybayes)
library(ggplot2)
library(ggthemes)
library(scales)
library(compositions)
library(rstan)
library(corrplot)
library(bayesplot)

set.seed(1234)
data <- read.csv("finalbis.csv")

data1 <- data[] # Remove rows containing N/A - none apply, remove none

# store the mean and sd for tariff before scaling
avg_entry_tariff_mean <- mean(as.numeric(data1$avg_entry_tariff), na.rm = TRUE)
avg_entry_tariff_sd <- sd(as.numeric(data1$avg_entry_tariff), na.rm = TRUE)

# Create dataframes for ethnic and POLAR proportions
eth_comp <- as.data.frame(cbind(data1$White.ethnic.group, 
  data1$Black.ethnic.group, data1$Asian.ethnic.group, 
  data1$Mixed.ethnic.group, data1$Other.ethnic.group))

pol_comp <- as.data.frame(cbind(data1$POLAR4.Q1, data1$POLAR4.Q2, 
     data1$POLAR4.Q3, data1$POLAR4.Q4, data1$POLAR4.Q5))

# Apply the ILR transformations to those variables representing 
# compositional data 
ilr_eth <- ilr(eth_comp)
ilr_pol <- ilr(pol_comp)

# store results as data frames
ilr_eth_df <- as.data.frame(ilr_eth)
ilr_pol_df <- as.data.frame(ilr_pol)

# return column names for identification
colnames(ilr_eth_df) <- paste0(c("White", "Black", "Asian", "Mixed"))
colnames(ilr_pol_df) <- paste0("POLAR4.Q", 1:ncol(ilr_pol_df))

# Scale the response so it is between 0 and 1
data1$career_after_15_month  <- as.numeric(data1$career_after_15_month) / 100

# Standardise the covariates
data1$avg_entry_tariff <- scale(as.numeric(data1$avg_entry_tariff))
data1$satisfied_feedback <- scale(data1$satisfied_feedback)
data1$satisfied_teaching <- scale(data1$satisfied_teaching)
data1$students_staff_ratio <- scale(data1$students_staff_ratio)
data1$added_value <- scale(data1$added_value)
data1$Women <- scale(data1$Women)
data1$spent_per_student <- scale(data1$spent_per_student)

# Create combined dataframe with ILR transformed variables for more complex fits
polar.data <- cbind(data1$career_after_15_month,data1$avg_entry_tariff,
                    data1$satisfied_feedback, data1$satisfied_teaching, 
                    data1$spent_per_student, data1$added_value, 
                    data1$students_staff_ratio,data1$Women,
                    ilr_eth_df,ilr_pol_df)

# Remove data1$ from column headings
colnames(polar.data) <- gsub("data1\\$", "", colnames(polar.data))


#correlation plot
predictors <- data1 %>%
  # Select the desired columns and convert them to numeric
  select(satisfied_teaching, satisfied_feedback, students_staff_ratio, 
         spent_per_student, avg_entry_tariff, added_value, 
         career_after_15_month, continuation, POLAR4.Q1) %>%
  mutate(across(everything(), as.numeric))

# Calculate the correlation matrix:
cor_matrix <- cor(predictors, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black",  # text label colour
         tl.srt = 35, 
         tl.cex = 0.9,      # var names text size
         addCoef.col = "black", # add coefficient numbers in black
         number.cex = 0.8)  # correlation coefficients text size



# Data visualisation plot
#create a scottish variable

data1$Scottish <- ifelse(is.na(data1$SIMD.2016.Q3), "Non-Scottish", "Scottish")
#just for corrplot
data3 <- data1[-c(21,40),]
ggplot(data3, aes(x = avg_entry_tariff, y = POLAR4.Q1, color = Scottish)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.2) +
  labs(x = "Average Entry Tariff",
       y = "Proportion from Disadvantaged Background (POLAR.Q1)",
       title = "Average Entry Tariff vs. Proportion of Disadvantaged Students",
       color = "University Type") +
  theme_minimal()



################################################
# Gaussian vs Beta Model Comparison
################################################
model_formula <- bf(
  career_after_15_month ~ avg_entry_tariff*POLAR4.Q1
)

# Gaussian bayesian model
model_gaus_bayes6 <- brm(
  formula = model_formula,
  data = data1,
  brmsfamily("gaussian"),
  chains = 4, iter = 8000, warmup = 2000,
  cores = 4, seed = 1234, 
  backend = "cmdstanr",
  file = "model_gaus_bayes6",  # Save this so it doesn't have to always rerun
  # set normal prior on regression coefficients (mean of 0, location of 3)
  prior =  c(prior(normal(0, 10), "b"), 
             # set normal prior on intercept (mean of 0, location of 3)             
             prior(normal(0, 10), "Intercept")) 
)

pp_check(model_gaus_bayes6)
#plot includes draws outside the range

#Choose BETA
# Simplest bayesian beta regression model
model_beta_bayes <- brm(
  bf(career_after_15_month ~ avg_entry_tariff*POLAR4.Q1,
     phi ~ avg_entry_tariff*POLAR4.Q1),
  data = data1,
  family = Beta(),
  prior = c(  # Prior for intercept (logit scale)
    set_prior("normal(1.5, 1)", class = "Intercept"),
    set_prior("normal(0, 1)", class = "b"),
    
    # Prior for the dispersion parameter (phi)
    set_prior("normal(0, 1)", class = "Intercept", dpar = "phi")),
  chains = 4, iter = 8000, warmup = 2000,
  cores = 4, seed = 1234, 
  backend = "cmdstanr",
  file = "model_beta_bayes"  # Save this so it doesn't have to always rerun
)

# Plot the posterior coefficients to assess significance, 
posterior_beta <- model_beta_bayes %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", "Mean"),
         intercept = str_detect(.variable, "Intercept"))

ggplot(posterior_beta, aes(x = .value, y = fct_rev(.variable), 
                           fill = component)) +
  geom_vline(xintercept = 0) +
  stat_halfeye(aes(slab_alpha = intercept), 
               .width = c(0.8, 0.95), point_interval = "median_hdi") +
  scale_fill_viridis_d(option = "viridis", end = 0.6) +
  scale_slab_alpha_discrete(range = c(1, 0.4)) +
  guides(fill = "none", slab_alpha = "none") +
  labs(x = "Coefficient", y = "Variable",
       caption = "80% and 95% credible intervals shown in black") +
  facet_wrap(vars(component), ncol = 1, scales = "free_y") +
  theme_clean()

# Use a dataset where POLAR4.Q1 is average and avg_ucas_tariff is a range
beta_bayes_pred_1 <- model_beta_bayes %>% 
  epred_draws(newdata = expand_grid(POLAR4.Q1 = mean(data1$POLAR4.Q1),
               avg_entry_tariff = seq(-1.5, 2.8, by = 0.05)))

# Plot how career outcomes change as entry tariff of the university increases
ggplot(beta_bayes_pred_1, aes(x = avg_entry_tariff * 
              avg_entry_tariff_sd + avg_entry_tariff_mean, y = .epred)) +
  stat_lineribbon() + 
  scale_y_continuous(labels = label_percent()) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Average UCAS entry tariff", 
       y = "Predicted Score for Career outcome after 15 months",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

# Use a dataset where tariff is average and POLAR4.Q1 is a range
beta_bayes_pred_1 <- model_beta_bayes %>% 
  epred_draws(newdata = expand_grid(avg_entry_tariff = 0,
                                    POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)))

# Plot how career outcomes change for universities with higher proportions of 
# quintile 1 students

ggplot(beta_bayes_pred_1, aes(x = POLAR4.Q1, y = .epred)) +
  stat_lineribbon() + 
  scale_y_continuous(labels = label_percent()) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "POLAR4", y = "Predicted Score for Career outcome after 15 months",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

# Posterior predictive check
pp_check(model_beta_bayes, ndraws = 100)
ppd <- posterior_predict(model_beta_bayes)
observed <- data1$career_after_15_month
pred_lower <- apply(ppd, 2, quantile, probs= 0.025)
pred_upper <- apply(ppd, 2, quantile, probs= 0.975)
coverage <- mean(observed >= pred_lower &observed <= pred_upper)
cat("Coverage Probability for 95% Credible Interval", coverage, "\n")


# LOOCV for this model
loo1 <- loo(model_beta_bayes, moment_match = TRUE)

# Create a grid for prediction: 
# For illustration, we use a range of avg_entry_tariff values 
#and fix POLAR4.Q1 at a particular value.
newdata <- expand_grid(
  avg_entry_tariff = seq(-1.5, 2.8, by = 0.05),
  POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)  # or a specific value of interest
)

# Extract posterior draws for the relevant coefficients from the fitted model.
post_draws <- spread_draws(model_beta_bayes, 
                           b_Intercept, 
                           b_avg_entry_tariff, 
                           b_POLAR4.Q1, 
                           `b_avg_entry_tariff:POLAR4.Q1`)



# Combine newdata with the posterior draws
predictions <- newdata %>%
  crossing(post_draws) %>%
  mutate(
    # Base linear predictor (using full tariff effect in interaction)
    linpred_base = b_Intercept +
      b_avg_entry_tariff * avg_entry_tariff +
      b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * avg_entry_tariff * POLAR4.Q1,
    
    # Hypothetical: increase avg_entry_tariff by 10% for students.
    # This multiplies avg_entry_tariff by 1.1 
    # in both the main effect and the interaction.
    # We must transform and then transform back this standardised variable
    linpred_cf = b_Intercept +
      b_avg_entry_tariff * (((((avg_entry_tariff * avg_entry_tariff_sd) +
  avg_entry_tariff_mean) * 1.1) - avg_entry_tariff_mean) / 
    avg_entry_tariff_sd) + b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * ((((((avg_entry_tariff * 
     avg_entry_tariff_sd) + avg_entry_tariff_mean) * 1.1) - 
       avg_entry_tariff_mean) / avg_entry_tariff_sd) * POLAR4.Q1)
  ) %>%

  # apply the inverse link (logistic function)
  mutate(
    epred_base = plogis(linpred_base),
    epred_cf   = plogis(linpred_cf)
  )

# Now, summarize the predictions for each scenario.
pred_summary <- predictions %>%
  group_by(avg_entry_tariff) %>%
  summarise(
    base_median = median(epred_base),
    base_lower  = quantile(epred_base, 0.025),
    base_upper  = quantile(epred_base, 0.975),
    cf_median   = median(epred_cf),
    cf_lower    = quantile(epred_cf, 0.025),
    cf_upper    = quantile(epred_cf, 0.975)
  )

# Plotting the two scenarios using ggplot2.

ggplot(pred_summary, aes(x = avg_entry_tariff * 
                           avg_entry_tariff_sd + avg_entry_tariff_mean)) +
  # Base scenario ribbon and line
  geom_ribbon(aes(ymin = base_lower, ymax = base_upper, fill = "Base"), 
              alpha = 0.3) +
  geom_line(aes(y = base_median, color = "Base"), linewidth = 1) +
  # Hypothetical scenario ribbon and line
  geom_ribbon(aes(ymin = cf_lower, ymax = cf_upper, 
                  fill = "10% Raise in tariff"), alpha = 0.3) +
  geom_line(aes(y = cf_median, color = "10% Raise in tariff"), size = 1) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Average UCAS Entry Tariff",
       y = "Pred. Career Outcome after 15 Mths",
       title = "Base vs. Hypothetical Predictions",
       fill = "Scenario",
       color = "Scenario") +
  theme_minimal() +
  theme(
    # Increase the axis title sizes
    axis.title = element_text(size = 18, face = "bold"),
    # Increase the axis tick label sizes
    axis.text = element_text(size = 14),
    # Increase the legend title and text sizes (scenario labels)
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    # Optionally increase the plot title size and center it
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    # Keep legend at the bottom
    legend.position = "bottom"
  )
# Create a grid for prediction: 
# For illustration, we use a range of avg_entry_tariff values and 
# fix POLAR4.Q1 at a particular value.
newdata <- expand_grid(
  avg_entry_tariff = seq(-1.5, 2.8, by = 0.05),
  POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)  # or a specific value of interest
)

# Extract posterior draws for the relevant coefficients from the fitted model.
post_draws <- spread_draws(model_beta_bayes, 
                           b_Intercept, 
                           b_avg_entry_tariff, 
                           b_POLAR4.Q1, 
                           `b_avg_entry_tariff:POLAR4.Q1`)



# Combine newdata with the posterior draws
predictions <- newdata %>%
  crossing(post_draws) %>%
  mutate(
    # Base linear predictor (using full tariff effect in interaction)
    linpred_base = b_Intercept +
      b_avg_entry_tariff * avg_entry_tariff +
      b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * avg_entry_tariff * POLAR4.Q1,
    
    # Hypothetical: increase avg_entry_tariff by 10% for Q1 students
    # This multiplies avg_entry_tariff by 1.1 in the interaction.
    linpred_cf = b_Intercept +
      b_avg_entry_tariff * avg_entry_tariff +
      b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * ((((((avg_entry_tariff * 
        avg_entry_tariff_sd) + avg_entry_tariff_mean) * 1.1) - 
          avg_entry_tariff_mean) / avg_entry_tariff_sd) * POLAR4.Q1)
  ) %>%
  # For a beta regression, we used logit link , 
  #so apply the inverse link (logistic function)
  mutate(
    epred_base = plogis(linpred_base),
    epred_cf   = plogis(linpred_cf)
  )

# Now, summarize the predictions for each scenario.
pred_summary <- predictions %>%
  group_by(POLAR4.Q1) %>%
  summarise(
    base_median = median(epred_base),
    base_lower  = quantile(epred_base, 0.025),
    base_upper  = quantile(epred_base, 0.975),
    cf_median   = median(epred_cf),
    cf_lower    = quantile(epred_cf, 0.025),
    cf_upper    = quantile(epred_cf, 0.975)
  )

# Plotting the two scenarios using ggplot2.
ggplot(pred_summary, aes(x = POLAR4.Q1)) +
  # Base scenario ribbon and line
  geom_ribbon(aes(ymin = base_lower, ymax = base_upper, 
                  fill = "Base"), alpha = 0.3) +
  geom_line(aes(y = base_median, color = "Base"), linewidth = 1) +
  # Hypothetical scenario ribbon and line
  geom_ribbon(aes(ymin = cf_lower, ymax = cf_upper, 
                  fill = "10% Raise in Tariff in Quintile 1"), alpha = 0.3) +
  geom_line(aes(y = cf_median, color = "10% Raise in Tariff in Quintile 1"), 
            size = 1) + scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Proportion belonging to the first quintile",
       y = "Predicted Career Outcome in 15 Mths",
       title = "Base vs. Hypothetical Predictions",
       fill = "Scenario",
       color = "Scenario") +
  theme_minimal() +
  theme(
    # Increase the axis title sizes
    axis.title = element_text(size = 18, face = "bold"),
    # Increase the axis tick label sizes
    axis.text = element_text(size = 14),
    # Increase the legend title and text sizes (scenario labels)
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    # Optionally increase the plot title size and center it
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    # Keep legend at the bottom
    legend.position = "bottom"
  )



##############################################################################
# Fit a more complex that fits all ILR transform POLAR4 quintiles along with 
# their interactions with UCAS Tariff
##############################################################################

model_beta_bayes_polar <- brm(
  bf(career_after_15_month ~ avg_entry_tariff+POLAR4.Q1+POLAR4.Q2+
       POLAR4.Q3+POLAR4.Q4+avg_entry_tariff:POLAR4.Q1+
       avg_entry_tariff:POLAR4.Q2+avg_entry_tariff:POLAR4.Q3+
       avg_entry_tariff:POLAR4.Q4,
     phi ~ avg_entry_tariff+POLAR4.Q1+POLAR4.Q2+POLAR4.Q3+POLAR4.Q4+
       avg_entry_tariff:POLAR4.Q1+avg_entry_tariff+POLAR4.Q2+
       avg_entry_tariff:POLAR4.Q3+avg_entry_tariff:POLAR4.Q4),
  data = polar.data,
  family = Beta(),
  prior = c(# Prior for intercept (logit scale)
    set_prior("normal(1.5, 1)", class = "Intercept"),  
    set_prior("normal(0, 1)", class = "b"),
    
    # Prior for the dispersion parameter (phi)
    set_prior("normal(0, 1)", class = "Intercept", dpar = "phi")),
  chains = 4, iter = 8000, warmup = 2000,
  cores = 4, seed = 1234, 
  backend = "cmdstanr",
  file = "model_beta_bayes_polar"  # Save this so it doesn't have to always rerun
)

# LOOCV for the more complex model
loo2 <- loo(model_beta_bayes_polar, moment_match=TRUE)

# Compare loo with simple model
loo::loo_compare(loo1, loo2)
################################################
# New model that includes ILR transformed racial categories
################################################
model_beta_bayes_race <- brm(
  bf(career_after_15_month ~ avg_entry_tariff+POLAR4.Q1+
       POLAR4.Q2+POLAR4.Q3+POLAR4.Q4+avg_entry_tariff:POLAR4.Q1+
       avg_entry_tariff:POLAR4.Q2+avg_entry_tariff:POLAR4.Q3+
       avg_entry_tariff:POLAR4.Q4 +  White+ Black+ Asian+ Mixed,
     phi ~ avg_entry_tariff+POLAR4.Q1+POLAR4.Q2+POLAR4.Q3+
       POLAR4.Q4+avg_entry_tariff:POLAR4.Q1+avg_entry_tariff:POLAR4.Q2+
       avg_entry_tariff:POLAR4.Q3+avg_entry_tariff:POLAR4.Q4+ 
       White+ Black+ Asian+ Mixed),
  data = polar.data,
  family = Beta(),
  prior = c(# Prior for intercept (logit scale)
    set_prior("normal(1.5, 1)", class = "Intercept"),  
    set_prior("normal(0, 1)", class = "b"),
    
    # Prior for the dispersion parameter (phi)
    set_prior("normal(0, 1)", class = "Intercept", dpar = "phi")),
  chains = 4, iter = 8000, warmup = 2000,
  cores = 4, seed = 1234, 
  backend = "cmdstanr",
  file = "model_beta_bayes_race"  # Save this so it doesn't have to always rerun
)

loo_race <- loo(model_beta_bayes_race, moment_match = TRUE)

loo::loo_compare(loo1, loo_race)


################################################
# Model assessing Continuation as the result
################################################
data1$continuation <- as.numeric(data1$continuation )/100
data2<-data1[-c(40),]
# Simplest bayesian beta regression model
model_beta_bayes_continuation <- brm(
  bf(continuation ~ avg_entry_tariff*POLAR4.Q1,
     phi ~ avg_entry_tariff*POLAR4.Q1),
  data = data2,
  family = Beta(),
  chains = 4, iter = 8000, warmup = 2000,
  cores = 4, seed = 1234, 
  backend = "cmdstanr",
  # Save this so it doesn't have to always rerun
  file = "model_beta_bayes_continuation"  
)

# Plot the posterior coefficients to assess significance, 
posterior_beta <- model_beta_bayes_continuation %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", "Mean"),
         intercept = str_detect(.variable, "Intercept"))

ggplot(posterior_beta, aes(x = .value, y = fct_rev(.variable), 
  fill = component)) + geom_vline(xintercept = 0) +
  stat_halfeye(aes(slab_alpha = intercept), 
               .width = c(0.8, 0.95), point_interval = "median_hdi") +
  scale_fill_viridis_d(option = "viridis", end = 0.6) +
  scale_slab_alpha_discrete(range = c(1, 0.4)) +
  guides(fill = "none", slab_alpha = "none") +
  labs(x = "Coefficient", y = "Variable",
       caption = "80% and 95% credible intervals shown in black") +
  facet_wrap(vars(component), ncol = 1, scales = "free_y") +
  theme_clean()

# Use a dataset where POLAR4.Q1 is average and avg_ucas_tariff is a range
beta_bayes_pred_cont_1 <- model_beta_bayes_continuation %>% 
  epred_draws(newdata = expand_grid(POLAR4.Q1 = mean(data2$POLAR4.Q1),
                                avg_entry_tariff = seq(-1.5, 2.8, by = 0.05)))

# Plot how contiuation outcomes change as entry tariff of the university increases
ggplot(beta_bayes_pred_cont_1, aes(x = avg_entry_tariff * 
                    avg_entry_tariff_sd + avg_entry_tariff_mean, y = .epred)) +
  stat_lineribbon() + 
  scale_y_continuous(labels = label_percent()) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Average UCAS entry tariff", y = "Predicted Continuation Outcome",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

# Use a dataset where tariff is average and POLAR4.Q1 is a range
beta_bayes_pred_cont_1 <- model_beta_bayes_continuation %>% 
  epred_draws(newdata = expand_grid(avg_entry_tariff = 0,
                                    POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)))

# Plot how Continuatio changes for universities with higher proportions of 
#quintile 1 students

ggplot(beta_bayes_pred_cont_1, aes(x = POLAR4.Q1, y = .epred)) +
  stat_lineribbon() + 
  scale_y_continuous(labels = label_percent()) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "POLAR4", y = "Predicted Continuation outcome",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

# Posterior predictive check
pp_check(model_beta_bayes_continuation, ndraws = 100)
ppd <- posterior_predict(model_beta_bayes_continuation)
observed <- data2$continuation
pred_lower <- apply(ppd, 2, quantile, probs= 0.025)
pred_upper <- apply(ppd, 2, quantile, probs= 0.975)
coverage <- mean(observed >= pred_lower &observed <= pred_upper)
cat("Coverage Probability for 95% Credible Interval", coverage, "\n")


# LOOCV for this model
loo4 <- loo(model_beta_bayes_continuation)
#mostly good values

# Use a dataset where tariff is average and POLAR4.Q1 is a range
beta_bayes_pred_cont_2 <- model_beta_bayes_continuation %>% 
  epred_draws(newdata = expand_grid(avg_entry_tariff = 0,
                                    POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)))

# Plot how Coninuation change for universities with higher proportions of 
#quintile 1 students

ggplot(beta_bayes_pred_cont_2, aes(x = POLAR4.Q1, y = .epred)) +
  stat_lineribbon() + 
  scale_y_continuous(labels = label_percent()) +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "POLAR4", y = "Predicted Coninuation outcome",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")


# Create a grid for prediction: 
# For illustration, we use a range of avg_entry_tariff values and fix 
#POLAR4.Q1 at a particular value.
newdata <- expand_grid(
  avg_entry_tariff = seq(-1.5, 2.8, by = 0.05),
  POLAR4.Q1 = seq(0.033, 0.342, by = 0.003)  # or a specific value of interest
)

# Extract posterior draws for the relevant coefficients from the fitted model.
post_draws <- spread_draws(model_beta_bayes_continuation, 
                           b_Intercept, 
                           b_avg_entry_tariff, 
                           b_POLAR4.Q1, 
                           `b_avg_entry_tariff:POLAR4.Q1`)



# Combine newdata with the posterior draws
predictions <- newdata %>%
  crossing(post_draws) %>%
  mutate(
    # Base linear predictor (using full tariff effect in interaction)
    linpred_base = b_Intercept +
      b_avg_entry_tariff * avg_entry_tariff +
      b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * avg_entry_tariff * POLAR4.Q1,
    
    # Hypothetical: increase avg_entry_tariff by 10% for Q1 students
    # This multiplies avg_entry_tariff by 1.1 in the interaction.
    linpred_cf = b_Intercept +
      b_avg_entry_tariff * avg_entry_tariff +
      b_POLAR4.Q1 * POLAR4.Q1 +
      `b_avg_entry_tariff:POLAR4.Q1` * ((((((avg_entry_tariff * 
  avg_entry_tariff_sd) + avg_entry_tariff_mean) * 1.1) - 
    avg_entry_tariff_mean) / avg_entry_tariff_sd) * POLAR4.Q1)
  ) %>%
  # For a beta regression, the link used was logit, 
  #so apply the inverse link (logistic function)
  mutate(
    epred_base = plogis(linpred_base),
    epred_cf   = plogis(linpred_cf)
  )

# Now, summarize the predictions for each scenario.
pred_summary <- predictions %>%
  group_by(POLAR4.Q1) %>%
  summarise(
    base_median = median(epred_base),
    base_lower  = quantile(epred_base, 0.025),
    base_upper  = quantile(epred_base, 0.975),
    cf_median   = median(epred_cf),
    cf_lower    = quantile(epred_cf, 0.025),
    cf_upper    = quantile(epred_cf, 0.975)
  )

# Plotting the two scenarios using ggplot2.
ggplot(pred_summary, aes(x = POLAR4.Q1)) +
  # Base scenario ribbon and line
  geom_ribbon(aes(ymin = base_lower, ymax = base_upper, fill = "Base"), 
              alpha = 0.3) +
  geom_line(aes(y = base_median, color = "Base"), linewidth = 1) +
  # Hypothetical scenario ribbon and line
  geom_ribbon(aes(ymin = cf_lower, ymax = cf_upper, 
                  fill = "10% Raise in Tariff in Quintile 1"), alpha = 0.3) +
  geom_line(aes(y = cf_median, color = "10% Raise in Tariff in Quintile 1"), 
            size = 1) +scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Proportion belonging to the first quintile",
       y = "Predicted Continuation",
       title = "Base vs. Hypothetical Predictions",
       fill = "Scenario",
       color = "Scenario") +
  theme_minimal() +
  theme(legend.position = "bottom")
