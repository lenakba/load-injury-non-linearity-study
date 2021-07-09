# This script is intended as a guide for modeling the relationship between
# a measure of training load and injury in sports research, using fractional polynomials (FP).
# The examples will go through standard logistic regression, 
# and later a mixed effects logistic regression model.However, the steps to model
# training load using FP in a Poisson model and other regression models are the same.

library(tidyverse) # for creating figures, and for a self-made function we'll need later on
library(rlang) # for creating functions with tidyverse syntax
library(mfp) # the mfp package has functions for automatic determination of FP
library(lme4) # package for mixed model functions
library(sjPlot) # for plotting predicted values with splines
library(ggeffects) # for model predictions with splines
library(clubSandwich) # for cluster-robust confidence intervals

# load your data. At the minimum, your data should have 
# - one column for load, measured in any metric (such as sRPE, GPS measures or ACWR)
# - one column for injury, must either be a logical variable (TRUE/FALSE) or coded (0 for no injury, 1 for injury)
#   to work in our logistic regression model. 
# - one column for athlete ID, coupling the load values and injuries to the right person

# here, we used the d_example_guide.rds data, simulated from football data, for the example to be reproducible, but go ahead and 
# replace the object "d" with your own data. 
d = readRDS("d_example_guide.rds")

# For model predictions to be calculated by
# the ggeffects package, categorical variables
# must be factors, not the character class
d = d %>% mutate(p_id = as.factor(p_id))
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
# the sjPlot package has a handy function 
# which plots the model predictions with confidence intervals
# unfortunately, the syntax is not compatible with that in the mfp-package
# we therefore have to write out our model manually
# by checking the formula in the object "fit_fp"
fit_fp_manual = glm(injury ~ load + I(((load + 1)/100)^3), data = d, family = "binomial")
plot_model(fit_fp_manual, type = "pred", terms = "load [all]")

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

# if both were able to converge, we can determine best fit with AIC
AIC(fit_mixed)
AIC(fit_mixed_intercept)

# in the example data, the random intercept model had better fit

# We can visualize in the same manner as the model without random effects
# Again, we have to write out our fit manually to do so,
# as ggeffects is not compatible with the mfp syntax
fit_mixed_manual = glmer(injury ~ load + I(((load + 1)/100)^3) + (1|p_id), data = d,  family = "binomial")

# Since this is a mixed model, we can calculate cluster-robust
# confidence intervals which take into account the uncertainty stemming from random effects variance
# we can do this by calculating the predictions with ggeffects
# and specifying the function for obtaining the variance-covariance matrix from the model
# We use vcovCR() from the clubSandwich package and specify the clusters in a list in the vcov.args argument
preds = ggpredict(
  fit_mixed_manual, 
  "load [all]", 
  vcov.fun = "vcovCR", 
  vcov.type = "CR0", 
  vcov.args = list(p_id = d$p_id),
  type = "re.zi"
)
plot(preds)
