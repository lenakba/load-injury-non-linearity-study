
# This script is intended as a guide for modeling the relationship between
# a measure of training load and injury in sports research, using Restricted Cubic Splines (RCS).
# The examples will go through standard logistic regression, 
# and later a mixed effects logistic regression model.However, the steps to model
# training load using splines in a Poisson model and other regression models are the same.
# The best-practice information on splines is based on Frank Harrell's "Regression Modeling Strategies"
# https://link.springer.com/book/10.1007/978-3-319-19425-7

library(ggplot2) # for creating figures
library(rms) # Harrell's rms package includes the functions we need for splines
library(lme4) # package for mixed model functions
library(merTools) # a sister-package to lme4 for extra functions on mixed models

# load your data. At the minimum, your data should have 
# - one column for load, measured in any metric (such as sRPE, GPS measures or ACWR)
# - one column for injury, must either be a logical variable (TRUE/FALSE) or coded (0 for no injury, 1 for injury)
#   to work in our logistic regression model. 
# - one column for athlete ID, coupling the load values and injuries to the right person

# here, we used the d_example_guide.rds data, simulated from football data, for the example to be reproducible, but go ahead and 
# replace the object "d" with your own data. 
d = readRDS("d_example_guide.rds")

#-------------------------------------------first steps, standard logistic regression

# for a regular logistic regression model, we can write
fit_logistic = glm(injury ~ load, data = d)

# running fit_logistic or summary(fit_logistic) will provide us the results and information from the fit
# For restricted cubic splines, 3-5 knots are sufficient in the vast majority of cases. 
# Choosing the number of knots can also be determined using Akaike's information Criterion. In this guide, we choose 3.
# for a restricted cubic splines with 3 knots, we simply write
fit_splines = glm(injury ~ rcs(load, 3), data = d)

# Here, the placement is determined by the default settings for the rcs function in the rms package.
# according to the documentation, knots are placed based on data quartiles and tertiles etc.
# if we run
hist(d$load)
# we see that the example data is fairly skewed. Partitioning by quartiles may not be the best fit.
# running a histogram to check the data should be done before running the first model (unlike in our example).
# In the code below, instead of the number of knots, we feed the argument with a vector of where along our load 
# variable we wish our knots to be placed
fit_splines_loc = glm(injury ~ rcs(load, c(500, 1500, 2500)), data = d)

# we can check which was better using the model with the lowest AIC
AIC(fit_splines)
AIC(fit_splines_loc)
#looks like placing the knots had a better fit.

# if we run
summary(fit_splines_loc)
# we'll see our coefficients (Estimate), standard error and p-values
# Even should we transform our estimates to an Odds Ratio, they make little sense.
# The p-values can be used and understood as usual.
# The best way to interpret splines is by visualizing predictions. 

# in a case where multiple variables have been included in the dataset, 
# a new dataset needs to be created where the other variables have been 
# set to a constant parameter, i.e. using the mean age
# we make a new object name for this dataset here, so we don't overwrite our original data object
pred_data = d
pred_data$age = 17
# we then add our data to the argument, newdata. type = "response" means 
# we will receive probability of injury instead of logodds of injury.
pred_load = predict(fit_splines_loc, type = "response", newdata = pred_data)

# to create a figure from our predictions, 
# let's add the predictions to our dataset with the load values used for predictions
pred_data$yhat = pred_load

# we loaded the tidyverse package, which includes the ggplot2 package.
# ggplot of the simplest, default form
ggplot(pred_data, aes(x = load, y = yhat)) +
  geom_line()

# Mapping the load values to the predicted probabilities showed that the splines
# modeled a U-shaped relationship in the example data.

#--------------------------------------------mixed model with splines

# a mixed effects logistic regression model is a bit more complicated
# most notably, because glmm isn't part of base R, so not everything is easily compatible

# for a mixed model with binomial distribution and random intercept per athlete
fit_mixed = glmer(injury ~ load + (1 | p_id), family = "binomial", data = d)
# for a mixed model with binomial distribution, random intercept and random slope per athlete
fit_mixed_slope = glmer(injury ~ load + (load | p_id), family = "binomial", data = d)
# finally, our splines model with knots placed where we think they should be
# note that glmms can take some time to run
# in this example we use an intercept-only model, as the random slope model failed to converge
fit_mixed_splines = glmer(injury ~ rcs(load, c(500, 1500, 2500)) + (1 | p_id), family = "binomial", data = d)

# now for predictions
# since we now have the random effect, we must set the exampledata to a fixed athlete as our example
# we make a new object name for this dataset here, so we don't overwrite our original data object
pred_data_mixed = d
pred_data_mixed$p_id = 1

# this time, we use the predictInterval() function that also predict the confidence intervals. 
# this is why we loaded the merTools package. predictInterval() can't be used on anything but a glmer()-created objected
# note that estimating the confidence intervals might take some time
pred_load_mixed_ci = predictInterval(fit_mixed_splines, ignore.fixed.terms = 1, type = "probability", newdata = pred_data_mixed)

# we add the load data we used for predictions to our predicted values
pred_load_mixed_ci$load = pred_data_mixed$load

# we plot, now with confidence intervals
ggplot(pred_load_mixed_ci, aes(x = load, y = fit, min = lwr, max = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) # alpha is for transperancy

### if you don't want multiple predictions per load value, you can predict based on a distinct set of values
pred_data_distinct = as.data.frame(unique(d$load))
pred_data_distinct$p_id = 1
# the names need to be exactly the same as the dataset used to fit the model
names(pred_data_distinct)[1] = "load"
pred_load_mixed_distinct = predictInterval(fit_mixed_splines, ignore.fixed.terms = 1, type = "probability", newdata = pred_data_distinct)
pred_load_mixed_distinct$load = pred_data_distinct$load

# we plot, now with confidence intervals
ggplot(pred_load_mixed_distinct, aes(x = load, y = fit, min = lwr, max = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.3)

