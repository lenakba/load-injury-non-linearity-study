# This script is intended as a guide for modeling the relationship between
# a measure of training load and injury in sports research, using fractional polynomials (FP).
# The examples will go through standard logistic regression, 
# and later a mixed effects logistic regression model.However, the steps to model
# training load using FP in a Poisson model and other regression models are the same.

library(tidyverse) # for creating figures, and for a self-made function we'll need later on
library(rlang) # for creating functions with tidyverse syntax
library(mfp) # the mfp package has functions for automatic determination of FP
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
# fitting load with fractional polynomial terms, we can run
# (note that it may take a little time)
fit_fp = mfp(injury ~ fp(load), data = d, family = "binomial")

# the model searches for the best fit of a number of possible polynomial transformations
# it uses a backwards selection process, so it might take some time

# running summary() of the fit will provide information of how many polynomials terms were added (of 1 or 2), 
# and to what power they were chosen. Estimate is the beta-value and represents the logodds. 
# The p-value indicates significance of each term.
summary(fit_fp)
# by saving the object like this
fit_summary = summary(fit_fp)
# it's possible to access different data from the fit
# for example, the aic, or the coefficients, to use them later
fit_summary$aic
fit_summary$coefficients

# we can interpret our results numerically, but it may still be useful to
# create a figure of the predictions to see the relationship shape

# in a case where multiple variables have been included in the dataset, 
# a new dataset needs to be created where the other variables have been 
# set to a constant parameter, i.e. using the mean age
# we make a new object name for this dataset here, so we don't overwrite our original data object
pred_data = d
pred_data$age = 17
# we then add our data to the argument, newdata. type = "response" means 
# we will receive probability of injury instead of logodds of injury.
pred_load = predict(fit_fp, type = "response", newdata = pred_data)

# to create a figure from our predictions, 
# let's add the predictions to our dataset with the load values used for predictions
pred_data$yhat = pred_load

# we loaded the tidyverse package, which includes the ggplot2 package.
# ggplot of the simplest, default form
ggplot(pred_data, aes(x = load, y = yhat)) +
  geom_line()

# Mapping the load values to the predicted probabilities showed that the polynomials
# modeled a U-shaped relationship in the example data.

#--------------------------------------------mixed model with polynomials

# a mixed effects logistic regression model is a bit more complicated to model
# glmm isn't part of base R, so not everything is easily compatible
# most notably, the lme4 package (or nlme package if you prefer) is not
# compatible with the mfp package that determine the polynomials for us in a easy-to-use function

# for standard mixed models without fractional polynomials:
# a mixed model with binomial distribution and random intercept per athlete
fit_mixed = glmer(injury ~ load + (1 | p_id), family = "binomial", data = d)
# a mixed model with binomial distribution, random intercept and random slope per athlete
fit_mixed_slope = glmer(injury ~ load + (load | p_id), family = "binomial", data = d)

# running the mfp function won't work. We need another solution
# one solution would be manually model every possible incarnation of an FP2 model and manually determine the best fit
# but here is an automatized solution

# Below is a function, glmer_fp() that searches for the best FP fit in a standard regression model 
# with the mfp-function we used earlier. The FP-formula that was found to be the best one 
# is then extracted from the standard model and run in a random effects model using glmer().
# the default is a logistic regression (binomial distribution of outcome measure)

# The function was created using tidyverse, rlang (which is why we loaded these at the top), 
# the lme4 and mfp package.

# the function has the following arguments:
# d                 The dataset with the load and injury variable
# injury            The variable used to denote injury. Can be heealth problems or any other definition of injury.
# load              The load variable in the dataset. 
# rdm_effect        The random effect term. Must be surrounded by quotes "". Examples are "(1|your_id_variable)"
#                   for a random intercept and "(load|your_id_variable)" for a random slope + random intercept. 
# family            Determine the model family. The defualt is binomial, meaning the models will run logistic regression.
#                   can be sett to "poisson" or other alternatives, see ?glmer()
# for a multivariable model, extra covariates will have to be added
# using + after fp(!!load), i.e. !!injury ~ fp(!!load) + sex + age, but the names must match those in the data
# and the function will no longer be general
glmer_fp = function(d, injury, load, rdm_effect, family = "binomial"){
  injury = enexpr(injury)
  load = enexpr(load) 
  
  # run mfp to find best FP terms
  prox_fit = eval_bare(expr(mfp(!!injury ~ fp(!!load), data = d, family = family)))
  fp_form = prox_fit$formula 
  
  # automatically use that formula in a random effects model 
  formula_start = paste0("!!",fp_form[2], " ", fp_form[1])
  formula = paste0(formula_start, " ", fp_form[3], " + ",rdm_effect,"") 
  glmm_fit = eval_bare(expr(glmer(as.formula(formula), family = family, data = d)))
  glmm_fit
}

# an example of using the function with a random intercept and random slope
fit_mixed = glmer_fp(d, injury, load, "(load|p_id)")

# an example of using the function with a random intercept only
fit_mixed_intercept = glmer_fp(d, injury, load, "(1|p_id)")

fit_mixed = glmer_fp(d, injury, load, "(1|p_id) + (load|p_id)")

# if both were able to converge, we can determine best fit with AIC
AIC(fit_mixed)
AIC(fit_mixed_intercept)

# in the example data, the random intercept model had better fit

# now for predictions
# since we now have the random effect, we must set the exampledata to a fixed athlete as our example
# we make a new object name for this dataset here, so we don't overwrite our original data object
pred_data_mixed = d
pred_data_mixed$p_id = 1

# this time, we use the predictInterval() function that also predict the confidence intervals. 
# this is why we loaded the merTools package. predictInterval() can't be used on anything but a glmer()-created objected
# note that estimating the confidence intervals might take some time
pred_load_mixed_ci = predictInterval(fit_mixed_intercept, ignore.fixed.terms = 1, type = "probability", newdata = pred_data_mixed)

# we add the load data we used for predictions to our predicted values
pred_load_mixed_ci$load = pred_data_mixed$load

# we plot, now with confidence intervals
ggplot(pred_load_mixed_ci, aes(x = load, y = fit, min = lwr, max = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) + # alpha is for transparency
  theme_light() # a simple way to clean up

### if you don't want multiple predictions per load value, you can predict based on a distinct set of values
pred_data_distinct = as.data.frame(unique(d$load))
pred_data_distinct$p_id = 1
# the names need to be exactly the same as the dataset used to fit the model
names(pred_data_distinct)[1] = "load"
pred_load_mixed_distinct = predictInterval(fit_mixed_intercept, ignore.fixed.terms = 1, type = "probability", newdata = pred_data_distinct)
pred_load_mixed_distinct$load = pred_data_distinct$load

# we plot, now with confidence intervals
ggplot(pred_load_mixed_distinct, aes(x = load, y = fit, min = lwr, max = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  theme_light()

