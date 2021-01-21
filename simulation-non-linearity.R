
# all functions used to perform the simulations in the study is in this script
# if data are loaded correctly, running the rest of the script should reproduce results in the study

# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors
options(scipen = 17, 
        stringsAsFactors = FALSE)

# Packages
library(tidyverse) # for everything
library(rlang) # for tidy evaluation
library(rms) # for restriced cubic splines
library(slider) # calculations on sliding windows
library(mfp) # fractional polynomials
library(lme4) # mixed models with random effects
library(SimCorMultRes) # for simulating longitudinal data
library(merTools) # tools for working with merMod objects from lme4 package
library(parameters) # easy functions for obtaining model parameters
library(DescTools) # for brier score and other model fits


# data is available in the GitHub repository. 
# includes only the U19 load and ACWR values, nothing else.
# Here, we assume working directory is set to "source file location", otherwise, the path must we added
# into the loading of the data
d_load = readRDS("d_load_selected.rds")
d_acwr = readRDS("d_acwr_selected.rds")

# Step 1) Create 2 datasets with 2 levels of n load values: 
# 1) original load values
# 3) 3 teams (3*25 football players) followed 1 season. Soccer season is 10 months (300 days) = 75*300 = 22500

# simulated cluster size and sample size
clsize = 300
n_athletes = 25

# cluster size and sample size in original data
clsize_orig = nrow(d_load)/62
n_athletes_orig = nrow(d_load)/clsize_orig
clsize_orig_acwr = 332
n_athletes_orig_acwr = 19

#----------------------- function for sampling Load-measures
sample_nonmissing = function(tl, n){
  sample_tl = sample(tl, n, replace = TRUE)
  sample_tl
}

# function for sampling load measueres based on number of athletes and cluster size
sample_longitudinal = function(load, n_athletes, clsize){
  d = enframe(sample_nonmissing(load, clsize*n_athletes), name = NULL, value = "load") %>% 
    mutate(p_id = as.character(rep(1:n_athletes, each = clsize)))
  d 
}

#---------------------------logreg and logit function used more than once
log_reg = function(tl_coef){
  res = 1 / (1 + exp(-tl_coef))
  res
}
logit = function(x) log(x/(1-x))

###---------------------------------------Creating risk profiles--------------------------------------

# Step 2) For each dataset created, create datasets with fake injuries by pre-determined probability distributions: 
# - Injuries from a linear relationship with load
# - Injuries from a quadratic relationship with load
# - No association (flat), injuries are added from a uniform distribution
# - For ACWR-load only, injuries from Gabbett’s J-shape-function

#--------------- sRPE functions
lin_function = function(load){
  y = log_reg(-0.5 + 0.001*load)
  y}
u_function = function(load){
  y = log_reg((0.000002*(load-1500)^2) - 1)
  y}

#---------------- ACWR function
# create risk shapes
# coefficients based on the risk profiles in Carey et al. 2018
coefs_acwr = function(d, acwr){
  
  # create categories based on cutoffs
  # using cutoffs from Carey et al. 2018
  # based on figure in Gabbett 2016
  acwr = enexpr(acwr)
  d = d %>% mutate(acwr_j =  case_when(!!acwr < 1 ~ 1,
                                       !!acwr >= 1 & !!acwr < 1.75 ~ 2, 
                                       !!acwr >= 1.75 ~ 3)
  )
  
  # add coefficients as column to dataset
  d = d %>% 
    mutate(coefs_j_acwr = log_reg(case_when(acwr_j == 1 ~ (-3.4 + 2*(1-!!acwr)^2),
                                            acwr_j == 2 ~ -3.4 + ((1-!!acwr)^2),
                                            acwr_j == 3 ~ (1.5*!!acwr)-5.4))
    )
  d
}

#--------------------simulating fake injuries based on coefficients
# we assume correlation between measures in a longitudinal study before adding injuries

# function for creating a covariance matrix with autoregressive correlation
# https://i.stack.imgur.com/I3uwR.jpg 
# values closer in time are more correlated than those further away in time
ar1_cor = function(n, rho = 0.8) {
  exponent = abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                   (1:n - 1))
  rho^exponent
}

# function that simulates injuries based on simulated longitudinal correlations
sim_long_injury = function(d, clsize, formula){# Simulation of correlated binary responses
  formula = enquo(formula)
  # Define the marginal risk [Log-odds scale]
  d_formula = d %>% mutate(formula_logit = logit(!!formula))
  FUN = d_formula %>% pull(formula_logit)
  
  #autoregressive covariance matrix
  matrix = ar1_cor(clsize)
  
  d_sim_long = SimCorMultRes::rbin(
    # Number of repeated obs. per athlete
    clsize = clsize,
    # Formula for the marginal risk model
    xformula = ~FUN, 
    # Intercept for marginal risk
    intercepts = 0,
    # Coefficents for marginal risk
    betas = 1,
    # Correlation matrix for response
    cor.matrix = matrix,
    # Link function
    link = "logit"
  )
  d_sim = d_sim_long$simdata %>% tibble() %>% dplyr::select(y)
  d_sim
}
# test
#sim_long_injury(l_loads[[1]], clsize_orig, coefs_lin) %>% rename(injury_u = y)

#--------------------------------helper function for sampling, creating coefficients and simulating injuries all in one

# function which bootstraps data, adds injuries based on all the different risk functions
simulate_injury = function(d, rep = 1, n_athletes = 25, clsize = 300, load = "load", sim = TRUE){
  
  if(sim){
    # resample load values
    d_sample = sample_longitudinal(d$load, n_athletes, clsize)
  } else {
    d_sample = d
  }
  
  if(load == "acwr"){
    # create risk shape and simulate injuries
    d_sample = d_sample %>% coefs_acwr(., load) 
    d_sample = d_sample %>% mutate(injury_j = sim_long_injury(d_sample, clsize, coefs_j_acwr) %>% pull(),
                                   load_cat = acwr_j)
  } else if(load == "load"){
    d_sample = d_sample %>% mutate(coefs_u = u_function(.$load),
                                   coefs_lin = lin_function(.$load),
                                   coefs_flat = 0.5)
    d_u = sim_long_injury(d_sample, clsize, coefs_u) 
    d_lin = sim_long_injury(d_sample, clsize, coefs_lin) 
    d_flat = sim_long_injury(d_sample, clsize, coefs_flat) 
    d_sample = d_sample %>% mutate(injury_u = d_u$y, injury_lin = d_lin$y, injury_flat = d_flat$y,
                                   load_cat = case_when(load < 500 ~ 1,
                                                        load >= 500 & load < 1500 ~ 2,
                                                        load >= 1500 & load < 2500 ~ 3,
                                                        load >= 2500 ~ 4))
  }
  
  # create load variable with some noise in it
  # z = max(x) - min(x) (assuming the usual case). 
  # The amount a to be added is either provided as positive argument amount or otherwise computed from z, as follows:
  # If amount == 0, we set a <- factor * z/50 (same as S). Meaning:
  # 1*((max(acwr_7_21) - min(acwr_7_21))/50)
  d_sample = d_sample %>% mutate(load_noised = jitter(load, 1, 0), load_noised = ifelse(load_noised < 0, 0, load_noised))
  d_sample
  
  # adding n label
  d_sample = d_sample %>% mutate(label = paste0("n = ",n_athletes*clsize,""))
  
  # create categories for categorical analysis and exponential for quadratic regression
  d_sample = d_sample %>% mutate(load_quarts = Hmisc::cut2(load_noised, g=4))
  d_sample
}

# for test run:
# l_acwr = list(
#   simulate_injury(d_acwr, load = "acwr",  sim = FALSE, clsize = clsize_orig_acwr, n_athletes = n_athletes_orig_acwr),
#   simulate_injury(d_acwr, load = "acwr",  sim = TRUE, clsize = 3*clsize)
# )


#-----------------------------------------Model fitting
# The following functions are used inside the final function for the simulation. See simulate_and_calc() below.

# There are no packages developed to automatically calculate fractional polynomials in a mixed model
# Here's a function that first determines the best fp values (which powers of x, load, have the best fit)
# in a non-mixed model. These fps are then used in a mixed model with a random intercept.
glmer_fp = function(d, injury){
  injury = enexpr(injury)
  
  prox_fit = eval_bare(expr(mfp(!!injury ~ fp(load), data = d, family = "binomial")))
  fp_form = prox_fit$formula 
  
  formula_start = paste0("!!",fp_form[2], " ", fp_form[1])
  formula = paste0(formula_start, " ", fp_form[3], " + (1 | p_id)")
  glmm_fit = eval_bare(expr(glmer(as.formula(formula), family = "binomial", data = d)))
  glmm_fit
}

#----------- Predictions
# function for adding confidence interval to predicted values, and the values used to predict
predict_and_combine = function(list_fit, list_example, method, load = load){
  load = enquo(load)
  
  l_pred = list_fit %>% map2(.x = ., .y = list_example, ~predict(.x, type = "response", newdata = .y))
  l_conf = list_fit %>% map2(.x = ., .y = list_example, ~merTools::predictInterval(.x, newdata = .y, ignore.fixed.terms = 1, type = "probability"))
  l_pred = l_pred %>% map(. %>% enframe(name = NULL, value = "yhat"))
  l_pred = l_pred %>% map2(.x = ., .y = l_conf, ~mutate(.x, upper = .y$upr, lower = .y$lwr))
  l_pred = l_pred %>% map2(.x = ., .y = list_example, ~mutate(.x, load = .y %>% dplyr::select(!!load) %>% pull(), label = .y$label))
  l_pred = l_pred %>% map(. %>% mutate(method = method))
  # unwrap our list to a table
  d_pred = do.call(bind_rows, l_pred) %>% as_tibble()
  d_pred
}

# to calculate coverage and other scores for the categorical data, we will first create fits using the categories,
# but then in the predicted data, categories are made into continuous. This way, they can be plotted in the same figure.
continuous_cats = function(d_cat, load = "acwr", subjective = FALSE){
  d_cat = tibble(d_cat)
  
  if(load == "acwr"){
    regex_last = ",[0-9][.][0-9][0-9][0-9]"
    regex_first = "[0-9]\\.[0-9][0-9][0-9],"
    regex_symbol = ","
  } else if(load == "load"){
    regex_last = ",.{1,}"
    regex_first = "[0-9]{1,},"
    regex_symbol = "[^[:alnum:]=\\.]"
  }
  
  last_nums = str_extract(unique(d_cat$load), regex(regex_last)) 
  last_nums_numeric = str_replace_all(last_nums, regex(regex_symbol), "") %>% str_replace_all(., regex("[:blank:]"), "") %>% as.numeric() %>% round()
  d_last_nums_numeric = sort(last_nums_numeric) %>% enframe(name = "load_quarts_code")
  
  first_nums = str_extract(unique(d_cat$load), regex(regex_first)) 
  first_nums_numeric = str_replace_all(first_nums, regex(regex_symbol), "") %>% str_replace_all(., regex("[:blank:]"), "") %>% as.numeric() %>% round()
  d_first_nums_numeric = sort(first_nums_numeric) %>% enframe(name = "load_quarts_code")
  
  d_cat = d_cat %>% mutate(load_quarts_code = as.numeric(factor(load)))
  size = nrow(d_cat)/2
  d_split1 = d_cat %>% slice(1:size) 
  d_split2 = d_cat %>% slice(size+1:n())
  
  if(!subjective){
    d_split1 = d_split1 %>% left_join(d_last_nums_numeric, by = "load_quarts_code")
    d_split2 = d_split2 %>% left_join(d_first_nums_numeric, by = "load_quarts_code")
  } else if(subjective == TRUE & load == "load"){
    d_split1 = d_split1 %>% mutate(value = case_when(load == 1 ~ 0,
                                                     load == 2 ~ 500,
                                                     load == 3 ~ 1499,
                                                     load == 4 ~ 2499))
    d_split2 = d_split2 %>% mutate(value = case_when(load == 1 ~ 499,
                                                     load == 2 ~ 1500,
                                                     load == 3 ~ 2500,
                                                     load == 4 ~ 3280))
  } else if(subjective == TRUE & load == "acwr"){
    d_split1 = d_split1 %>% mutate(value = case_when(load == 1 ~ 0,
                                                     load == 2 ~ 1,
                                                     load == 3 ~ 1.75))
    
    d_split2 = d_split2 %>% mutate(value = case_when(load == 1 ~ 1,
                                                     load == 2 ~ 1.75,
                                                     load == 3 ~ 3))
    
  }
  d_full = bind_rows(d_split1, d_split2) %>% dplyr::select(-load) %>% rename(load = value)
  d_full
}

#function for obtaining root mean squared error RMSE
rmse = function(fit){
  fit_sum = summary(fit)
  rmse = sqrt(mean(fit_sum$residuals^2))
  rmse
}

# function for obtaining parameters from any model fit
# specify method for a column with the model-type
get_params = function(fit, method){
  d_params = parameters::parameters(fit) %>% tibble()
  d_params = d_params %>% mutate(method = method, 
                                 rmse = rmse(fit))
  d_params
}

#function for obtaining c-statistics and brier score
validation_stats = function(l_fit, injury, l_data){
  injury = enexpr(injury)
  brier = map2(.x = l_fit, .y = l_data, ~DescTools::BrierScore(pred=predict(.x, type="response"), resp = .y %>% dplyr::select(!!injury) %>% pull())) %>% 
    map(., . %>% enframe(name = NULL, value = "brier"))
  c_statistics = map2(.x = l_fit, .y = l_data, ~DescTools::Cstat(x=predict(.x, type="response"), resp = .y %>% dplyr::select(!!injury) %>% pull())) %>% 
    map(., . %>% enframe(name = NULL, value = "c_stat"))
  validation_stats = map2(.x = brier, .y = c_statistics, ~.x %>% mutate(c_stat = .y$c_stat))
  validation_stats
}

# function for calculating coverage
coverage = function(d_preds, example_list, coefs, rep, cat = FALSE, load_type = "load", subjective = FALSE){
  coefs = enquo(coefs)
  rep = enquo(rep)
  
  if(cat){
    l_preds = d_preds %>% group_by(label) %>% nest()
    if(subjective){
      l_preds$data = l_preds$data %>% map(., ~continuous_cats(., load = load_type, subjective = TRUE))     
    } else {
      l_preds$data = l_preds$data %>% map(., ~continuous_cats(., load = load_type))      
    }
    d_preds = unnest(l_preds, cols = data) %>% ungroup()
  }
  
  
  # make dataset with coefficients to compare
  d_collated_load = example_list %>% bind_rows()
  l_example_labelled_load = list(d_collated_load)
  methods = d_preds %>% distinct(method) %>% dplyr::select(method) %>% pull()
  for(i in 1:length(methods)){
    l_example_labelled_load = map_at(l_example_labelled_load, i, ~mutate(., method = methods[i]))
  }
  
  d_coefs_load = l_example_labelled_load[[1]]
  d_coefs_load = d_coefs_load %>% dplyr::select(label, load, !!coefs, method, load_cat)  
  
  if(cat){
    d_preds_coefs = d_preds %>% full_join(d_coefs_load, by = c("method", "label", "load_quarts_code" = "load_cat"))
  } else {
    d_preds_coefs = d_preds %>% full_join(d_coefs_load, by = c("method", "label", "load"))
  }
  remove(d_preds, d_coefs_load, l_example_labelled_load)
  
  d_coverage = d_preds_coefs %>% 
    mutate(coverage = ifelse((upper >= !!coefs & lower <= !!coefs) & !is.na(yhat), 1, 0), 
           missing_pred = ifelse(!is.na(!!coefs) & is.na(yhat), 1, 0), 
           missing_and_coverage = ifelse((upper < !!coefs | lower > !!coefs) | is.na(yhat), 1, 0),
           rep = !!rep)
  
  d_prop_coverage = d_coverage %>% group_by(label, method, !!rep) %>% 
    summarise(n_coverage = sum(coverage==1, na.rm = TRUE), 
              denominator = n(), 
              prop = n_coverage/denominator,
              n_missing_pred = sum(missing_pred==1, na.rm = TRUE),
              prop_missing = n_missing_pred/denominator,
              n_missing_cov = sum(missing_and_coverage==1, na.rm = TRUE),
              prop_missing_cov = n_missing_cov/denominator) %>% 
    ungroup() %>% mutate(label = as.integer(str_replace_all(label, "n = ", "")))
  d_prop_coverage
}

#-----------------------------------------------------Final helper function

# simulate calc uses all of the functions above
# it simulates injuries with the simulate_injury function, given an injury type as an argument
# - fits all the 7 models
# - runs predictions on the original data, thereby coverage can be calculated
# - calculates all other performance metrics.
# Required arguments:
# load_type                       Either "load", the default for sRPE, or "acwr" to indicate what the load metric used
# injury                          Either injury_u, injury_j, injury_lin, injury_flat
# coefs                           Either coefs_u, coefs_j_acwr, coefs_lin, coefs_flat
# rep                             The permutation. rep = 1 by default. In the for-loop, adding rep = i will then
#                                 add the needed permutation number to each dataset.
simulate_and_calc = function(load_type = "load", injury, coefs, rep = 1){
  
  injury = enexpr(injury)
  coefs = enquo(coefs)
  
  # simulate injuries and create a list of datasets of unequal sizes
  if(load_type == "acwr"){
    l_load = list(
      simulate_injury(d_acwr, load = "acwr",  sim = FALSE, clsize = clsize_orig_acwr, n_athletes = n_athletes_orig_acwr),
      simulate_injury(d_acwr, load = "acwr",  sim = TRUE, clsize = 3*clsize)
    )
  } else if(load_type == "load"){
    l_load = list(
      simulate_injury(d_load, load = load_type, sim = FALSE, clsize = clsize_orig, n_athletes = n_athletes_orig),
      simulate_injury(d_load, load = load_type,  sim = TRUE, clsize = 3*clsize)
    ) 
  }
  
  # buffer for convergence issues
  optctrl = list(maxfun=100000)
  
  # fit models
  fit_splines = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ rcs(load_noised, 3) + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  
  if(load_type == "load"){
    fit_splines_loc = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ rcs(load_noised, c(500, 1500, 2500)) + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  } else if(load_type == "acwr"){
    fit_splines_loc = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ rcs(load_noised, c(1, 1.75, 2)) + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  }
  fit_fp = l_load %>% map(~eval_bare(expr(glmer_fp(., !!injury))))
  fit_lin =  l_load %>% map(~eval_bare(expr(glmer(!!injury ~ load_noised + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  fit_cat = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ load_quarts + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  fit_cat_loc = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ as.factor(load_cat) + (1 | p_id), family = "binomial", data = ., control=glmerControl(optCtrl=optctrl)))))
  
  # # quadratic regression on regular load values caused these problems: https://rpubs.com/jimsavage/scale_issues
  if(load_type == "load"){
    fit_quad = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ I(load_noised/500) + I((load_noised/500)^2) + (1 | p_id), family = binomial("logit"), data = ., control=glmerControl(optCtrl=optctrl)))))
  } else if(load_type == "acwr"){
    fit_quad = l_load %>% map(~eval_bare(expr(glmer(!!injury ~ I(load_noised) + I((load_noised)^2) + (1 | p_id), family = binomial("logit"), data = ., control=glmerControl(optCtrl=optctrl)))))
  }
  
  # example datasets used for predictions
  list_example = l_load[[2]] %>% mutate(p_id = "1")
  p_id_list = l_load %>% map(. %>% distinct(p_id) %>% slice(1) %>% pull())
  list_example = l_load %>% map(. %>% dplyr::select(p_id, label, starts_with("load"), starts_with("coefs"))) %>% map2(.x = ., .y = p_id_list, ~.x %>% mutate(p_id = .y))
  
  # calculate coverage
  d_preds_splines = predict_and_combine(fit_splines, list_example, method = "Restricted Cubic Splines (Data-driven)")
  d_coverage_splines = coverage(d_preds_splines, list_example, coefs = !!coefs, rep = rep)
  remove(d_preds_splines)
  d_preds_splines_loc = predict_and_combine(fit_splines_loc, list_example, method = "Restricted Cubic Splines (Subjectively)")
  d_coverage_splines_loc = coverage(d_preds_splines_loc, list_example, coefs = !!coefs, rep = rep)  
  remove(d_preds_splines_loc)
  d_preds_fp = predict_and_combine(fit_fp, list_example, method = "Fractional Polynomials")
  d_coverage_fp = coverage(d_preds_fp, list_example, coefs = !!coefs, rep = rep)  
  remove(d_preds_fp)
  d_preds_lin = predict_and_combine(fit_lin, list_example, method = "Linear Regression")
  d_coverage_lin = coverage(d_preds_lin, list_example, coefs = !!coefs, rep = rep) 
  remove(d_preds_lin)
  d_preds_cat = predict_and_combine(fit_cat, list_example, method = "Categorized (Quartiles)", load = load_quarts)
  d_coverage_cat = coverage(d_preds_cat, list_example, coefs = !!coefs, rep = rep, cat = TRUE, load_type = load_type)
  remove(d_preds_cat)
  d_preds_cat_loc = predict_and_combine(fit_cat_loc, list_example, method = "Categorized (Subjectively)", load = load_cat)
  d_coverage_cat_loc = coverage(d_preds_cat_loc, list_example, coefs = !!coefs, rep = rep, cat = TRUE, load_type = load_type, subjective = TRUE)
  remove(d_preds_cat_loc)
  d_preds_quad = predict_and_combine(fit_quad, list_example, method = "Quadratic Regression")
  d_coverage_quad = coverage(d_preds_quad, list_example, coefs = !!coefs, rep = rep) 
  remove(d_preds_quad)
  d_prop_coverage = bind_rows(d_coverage_splines, d_coverage_splines_loc, d_coverage_fp, d_coverage_lin, d_coverage_cat, d_coverage_cat_loc, d_coverage_quad)
  
  # obtain c-statistics and brier score
  v_splines = validation_stats(fit_splines, injury = !!injury, l_data = l_load)
  v_splines_loc = validation_stats(fit_splines_loc, injury = !!injury, l_data = l_load)
  v_fp = validation_stats(fit_fp, injury = !!injury, l_data = l_load)
  v_lin = validation_stats(fit_lin, injury = !!injury, l_data = l_load)
  v_cat = validation_stats(fit_cat, injury = !!injury, l_data = l_load)
  v_cat_loc = validation_stats(fit_cat_loc, injury = !!injury, l_data = l_load)
  v_quad = validation_stats(fit_quad, injury = !!injury, l_data = l_load)
  l_v_stats = list(v_splines, v_splines_loc, v_fp, v_lin, v_cat, v_cat_loc, v_quad)
  l_v_stats = flatten(l_v_stats)
  remove(v_splines, v_splines_loc, v_fp, v_lin, v_cat, v_cat_loc, v_quad)
  
  # obtain parameters
  l_splines = fit_splines %>% map(., ~get_params(., "Restricted Cubic Splines (Data-driven)"))
  l_splines_loc = fit_splines_loc %>% map(., ~get_params(., "Restricted Cubic Splines (Subjectively)"))
  l_fp = fit_fp %>% map(., ~get_params(., "Fractional Polynomials"))
  l_lin = fit_lin %>% map(., ~get_params(., "Linear Regression"))
  l_cat = fit_cat %>% map(., ~get_params(., "Categorized (Quartiles)"))
  l_cat_loc = fit_cat_loc %>% map(., ~get_params(., "Categorized (Subjectively)"))
  l_quad = fit_quad %>% map(., ~get_params(., "Quadratic Regression"))
  remove(fit_splines, fit_splines_loc, fit_fp, fit_lin, fit_cat, fit_cat_loc, fit_quad)
  
  # combine to one big listC
  l_params = list(l_splines, l_splines_loc, l_fp, l_lin, l_cat, l_cat_loc, l_quad)
  remove(l_splines, l_splines_loc, l_fp, l_cat, l_cat_loc, l_quad)
  n_methods = length(l_params)
  l_params = flatten(l_params)
  l_params = l_params %>% map2(.x = ., .y = l_v_stats, ~cbind(., .y) %>% as_tibble()) # add validation stats
  
  # simulation data size
  l_n = l_load %>% map(~nrow(.))
  l_n = rep(l_n, n_methods)
  remove(l_load)
  
  # add sample size to datasets and collapse to one big dataframe
  d_params = l_params %>% map2(.x = ., .y = l_n, ~mutate(.x, n = .y)) %>% do.call(bind_rows, .) %>% as_tibble()
  d_params = d_params %>% mutate(sig = ifelse(p < 0.05, 1, 0), rep = rep) # 1 for sig, 0 for insig
  d_params = d_params %>% left_join(d_prop_coverage, by = c("method", "n" = "label", "rep"))
  d_params
}

# for a test, run
# simulate_and_calc(load_type = "acwr", injury = injury_j, coefs = coefs_j_acwr)

#-----------------for loop running simulations
n_sim = 1900 # set to the number of permutations. 
set.seed = 1234

# In the for-loop, each permutation of results is saved to 1 datafile
# this is to ensure that, should the simulation meet an error or something else should happen
# the simulations performed up to that point are still saved.
# Here, we have created a folder named "u" that the data are saved in. 
# We assume that the directory is "source file location" and that the u, j, lin and flat folders are created in that location. 
options(warn=-1)
rel = "u"
for(i in 1:n_sim){
  # capture seeds
  cat(capture.output(.Random.seed), file=paste0("seeds\\",i,"_random_seed.txt"), sep="\n")
  # run a permutation of simulate and calc
  d_u = simulate_and_calc(injury = injury_u, coefs = coefs_u, rep = i)
  # save 1 permutation worth of results to 1 file. We've made a folder before hand with the name "u"
  saveRDS(d_u, file=paste0("",rel,"\\",i,"_d_",rel,".rds"))
}
options(warn=0)

# Repeat for other risk shapes
# J
options(warn=-1)
rel = "j"
for(i in 4:n_sim){
  # capture seeds
  cat(capture.output(.Random.seed), file=paste0("seeds\\",i,"_random_seed.txt"), sep="\n")
  d_j = simulate_and_calc(load_type = "acwr", injury = injury_j, coefs = coefs_j_acwr, rep = i)
  saveRDS(d_j, file=paste0("",rel,"\\",i,"_d_",rel,".rds"))
}
options(warn=0)

options(warn=-1)
rel = "lin"
for(i in 4:n_sim){
  # capture seeds
  cat(capture.output(.Random.seed), file=paste0("seeds\\",i,"_random_seed.txt"), sep="\n")
  d_lin = simulate_and_calc(load_type = "acwr", injury = injury_j, coefs = coefs_j_acwr, rep = i)
  saveRDS(d_lin, file=paste0("",rel,"\\",i,"_d_",rel,".rds"))
}
options(warn=0)

options(warn=-1)
rel = "flat"
for(i in 4:n_sim){
  # capture seeds
  cat(capture.output(.Random.seed), file=paste0("seeds\\",i,"_random_seed.txt"), sep="\n")
  d_flat = simulate_and_calc(load_type = "acwr", injury = injury_j, coefs = coefs_j_acwr, rep = i)
  saveRDS(d_flat, file=paste0("",rel,"\\",i,"_d_",rel,".rds"))
}
options(warn=0)



#--------------------------reading RDS files 
# to check a single dataset:
readRDS("my\\file_location\\u\\1_d_u.rds") 

# to read all the simulated files at once
# we first make an empty dataset with the right dimensions
n_d_col = 23 # number of columns
n_d_rows_u = 43 # number of rows

make_dframe = function(n_d_col, n_d_rows){
  d_estimates = data.frame(matrix(ncol = n_d_col, nrow = (n_sim*n_d_rows))) # 60 = number of rows generated with 1 rep of simulate_and_fit()
  col_names = c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "z", "df_error", "p", "method", "rmse", "brier", "c_stat", "n", "sig", "rep", "n_coverage", "denominator", "prop", "n_missing_pred", "prop_missing", "n_missing_cov", "prop_missing_cov", "mse")
  colnames(d_estimates) = col_names
  d_estimates
}

d_estimates_u = make_dframe(n_d_col, n_d_rows_u)

# object with the file location where all the simulated dataset files are
# then save an object with the list of filenames
folder_u = "my\\file\\location\\u"
u_files = list.files(path = folder_u)

# for-loop that fills in the empty dataset with all of the results
for(i in 1:length(u_files)){
  pos_end = i*n_d_rows_u
  pos_start =  pos_end - (n_d_rows_u-1)
  
  d_estimates_u[pos_start:pos_end,] = readRDS(paste0(folder_u, "\\",i,"_d_u.rds"))
}

d_estimates_u = d_estimates_u %>% as_tibble()

# to save the results for use (calculating mean RMSE and so on)
saveRDS(d_estimates_u, file=paste0("D:\\r skript\\simulation_nonlinearity\\results_combined\\d_u.rds"))

