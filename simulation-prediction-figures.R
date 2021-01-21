# This R script performs a single simulation and creates figures visualizing the predictions
# against the (true) simulated coefficients
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
library(devEMF) # for exporting figures as .emf files
library(SimCorMultRes) # for simulating longitudinal data
library(merTools) # tools for working with merMod objects from lme4 package

# data is available in the GitHub repository. 
# includes only the U19 load and ACWR values, nothing else.
# if data are loaded correctly, running the rest of the script should reproduce results in the study
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

# sample data
set.seed(1234) # same seed as in Carey et al. 2018
load_3t = sample_longitudinal(d_load$load, 3*n_athletes, clsize)
l_loads = list(d_load, load_3t)

acwr_3t = sample_longitudinal(d_acwr$load, 3*n_athletes, clsize)
l_acwr = list(d_acwr, acwr_3t)

#---------------------------logreg and logit function used more than once
log_reg = function(tl_coef){
  res = 1 / (1 + exp(-tl_coef))
  res
}
logit = function(x) log(x/(1-x))
###---------------------------------------Creating risk profiles--------------------------------------

# Step 2) For each dataset created in dataset 4, create datasets with fake injuries by pre-determined probability distributions: 
# - Injuries from a linear relationship with load
# - Injuries from a quadratic relationship with load
# - No association (flat), injuries are added from a uniform distribution
# - One or more probability functions discovered in the non-simulated analysis
# - For ACWR-load only, injuries from Gabbett’s J-shape-function

#--------------- Load functions
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
# for test run
# sim_long_injury(l_loads[[1]], clsize_orig, coefs_lin) %>% rename(injury_u = y)

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

# create lists with datasets we use for our single simulation
# acwr list
l_acwr = list(
  simulate_injury(d_acwr, load = "acwr",  sim = FALSE, clsize = clsize_orig_acwr, n_athletes = n_athletes_orig_acwr),
  simulate_injury(d_acwr, load = "acwr",  sim = TRUE, clsize = 3*clsize)
)

# load list
l_load = list(
  simulate_injury(d_load, load = "load", sim = FALSE, clsize = clsize_orig, n_athletes = n_athletes_orig),
simulate_injury(d_load, load = "load",  sim = TRUE, clsize = 3*clsize)
)

#------------------------------------Visualization----------------------------------------------------------


#------------------------------------Objects to base predictions on
# example datasets used for predictions
d_acwr_example = l_acwr[[2]] %>% mutate(p_id = "1665")
d_load_example = l_load[[2]] %>% mutate(p_id = "1665")

p_id_list = l_acwr %>% map(. %>% distinct(p_id) %>% slice(1) %>% pull())
l_acwr_example = l_acwr %>% map(. %>% dplyr::select(p_id, label, starts_with("load"), starts_with("coefs"))) %>% map2(.x = ., .y = p_id_list, ~.x %>% mutate(p_id = .y))
l_load_example = l_load %>% map(. %>% dplyr::select(p_id, label, starts_with("load"), starts_with("coefs"))) %>% map2(.x = ., .y = p_id_list, ~.x %>% mutate(p_id = .y))

# values used in prediction
acwr_7_21_values = d_acwr_example$load
load_values = d_load_example$load

text_size = 14
# check distributions of load values
plot_hist_rpe = ggplot(enframe(load_values), aes(x = load_values)) + 
  geom_histogram() + theme_light() + xlab("sRPE (AU)") + ylab("Count")  +
    theme(axis.text = element_text(size=text_size),
          axis.title =  element_text(size=text_size))
plot_hist_acwr = ggplot(enframe(acwr_7_21_values), aes(x = acwr_7_21_values)) + 
  geom_histogram() + theme_light() + xlab("ACWR 7:21")  + ylab("Count") + 
  theme(axis.text = element_text(size=text_size),
        axis.title =  element_text(size=text_size))

#---------------------------------- Dataset with true coefficients for use in figures
d_collated_load = l_load_example %>% bind_rows()
l_example_labelled_load = list(d_collated_load, d_collated_load, d_collated_load, d_collated_load, d_collated_load, d_collated_load, d_collated_load)
d_collated_acwr = l_acwr_example %>% bind_rows()
l_example_labelled_acwr = list(d_collated_acwr, d_collated_acwr, d_collated_acwr, d_collated_acwr, d_collated_acwr, d_collated_acwr, d_collated_acwr)
methods = c("Restricted Cubic Splines\n(Data-driven)", "Restricted Cubic Splines\n(Subjectively)", "Fractional Polynomials", "Linear Model", "Quadratic Model", "Categorized (Quartiles)", "Categorized (Subjectively)")
for(i in 1:length(methods)){
  l_example_labelled_load = map_at(l_example_labelled_load, i, ~mutate(., method = methods[i]))
  l_example_labelled_acwr = map_at(l_example_labelled_acwr, i, ~mutate(., method = methods[i]))
}

d_coefs_load = l_example_labelled_load %>% bind_rows()
d_coefs_acwr = l_example_labelled_acwr %>% bind_rows()
#------------------------------------------------------model fitting---------------------------
# •	Logistic regression (linearity assumed)
# •	Logistic regression with load categorized (by quartiles)
# •	Logistic regression with load categorized (subjectively)
# •	Logistic regression with load as a quadratic term (quadratic regression)
# •	Logistic regression with load using splines (default placement of knots)
# •	Logistic regression with load using splines (subjective placement of knots)
# •	Logistic regression with load using fractional polynomials

# Function for finding the optimal number of knots in a splines model given AIC
find_optimal_knot = function(d, tl, injury){
  tl = enexpr(tl)
  injury = enexpr(injury)
  aic = rep(NA, 3)
  knots = 3:5
  for (i in seq_along(knots)) {
    knot = knots[i]
    fit_splines = eval_bare(expr(glmer(!!injury ~ rcs(!!tl, knot) + (1 | p_id), data = d, family = "binomial"))) 
    aic[i] = AIC(fit_splines)
  }
  d_knots = tibble(knots, aic)
  knot_optimal = d_knots %>% filter(aic == min(aic)) %>% pull("knots")
  knot_optimal
}

# function for fitting an lrm with splines given a load variable, dataset and injury variable, 
# where the optimal knot is chosen with AIC
fit_load_splines = function(d, load, injury){
  load = enexpr(load)
  injury = enexpr(injury)
  
  knot = find_optimal_knot(d, !!load, !!injury)
  knot = enexpr(knot)
  fit_splines = eval_bare(expr(glmer(!!injury ~ rcs(!!load, !!knot) + (1 | p_id), data = d, family = "binomial"))) 
  fit_splines
}

#------------- fit models
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

# to test, one can run:
#preds_splines_j_acwr = predict_and_combine(fit_splines_j_acwr, l_acwr_example, method = "Restricted Cubic Splines")
#preds_fp_j_acwr = predict_and_combine(fit_fp_j_acwr, l_acwr_example, method = "Fractional Polynomials")

# to plot with the other methods, we need to pretend our categories are continuous
# this function extracts the categories and makes them into a continuous variable
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

# helper function that fits all the models and performs all the predictions
fit_and_predict = function(list_loads, list_example, injury, load_type = "load"){
  
  injury = enexpr(injury)
  
  # fit all different models
  fit_splines = list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ rcs(load_noised, 3) + (1 | p_id), family = "binomial", data = .))))
  fit_splines_loc = list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ rcs(load_noised, c(500,1500,2500)) + (1 | p_id), family = "binomial", data = .))))
  fit_fp = list_loads %>% map(~eval_bare(expr(glmer_fp(., !!injury))))
  fit_lin =  list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ load_noised + (1 | p_id), family = "binomial", data = .))))
  fit_cat = list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ load_quarts + (1 | p_id), family = "binomial", data = .))))
  fit_cat_loc = list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ load_cat + (1 | p_id), family = "binomial", data = .))))
  # quadratic regression on regular load values caused these problems: https://rpubs.com/jimsavage/scale_issues
  fit_quad = list_loads %>% map(~eval_bare(expr(glmer(!!injury ~ I(load_noised/500) + I((load_noised/500)^2) + (1 | p_id), family = binomial("logit"), data = .))))

  # make all the predictions
  d_preds_splines = predict_and_combine(fit_splines, list_example, method = "Restricted Cubic Splines\n(Data-driven)")
  d_preds_splines_loc = predict_and_combine(fit_splines_loc, list_example, method = "Restricted Cubic Splines\n(Subjectively)")
  d_preds_fp = predict_and_combine(fit_fp, list_example, method = "Fractional Polynomials")
  d_preds_lin = predict_and_combine(fit_lin, list_example, method = "Linear Model")
  d_preds_cat = predict_and_combine(fit_cat, list_example, method = "Categorized (Quartiles)", load = load_quarts) 
  d_preds_cat_loc = predict_and_combine(fit_cat_loc, list_example, method = "Categorized (Subjectively)", load = load_cat) 
  d_preds_quad = predict_and_combine(fit_quad, list_example, method = "Quadratic Model")
  
  # forcing categorized values into continuous
  l_preds_cat = d_preds_cat %>% group_by(label) %>% nest()
  l_preds_cat$data = l_preds_cat$data %>% map(., ~continuous_cats(., load = load_type))
  d_preds_cat = unnest(l_preds_cat, cols = data)
  l_preds_cat_loc = d_preds_cat_loc %>% group_by(label) %>% nest()
  l_preds_cat_loc$data = l_preds_cat_loc$data %>% map(., ~continuous_cats(., load = load_type, subjective = TRUE))
  d_preds_cat_loc = unnest(l_preds_cat_loc, cols = data)
  
  # combine to single dataset
  d_preds = bind_rows(d_preds_splines, d_preds_splines_loc, d_preds_fp, d_preds_lin, d_preds_quad, d_preds_cat, d_preds_cat_loc)
  d_preds
}

# run predictions for all models
d_preds_j = fit_and_predict(l_acwr, l_acwr_example, injury = injury_j, load_type = "acwr")
d_preds_u = fit_and_predict(l_load, l_load_example, injury = injury_u, load_type = "load")
d_preds_flat = fit_and_predict(l_load, l_load_example, injury = injury_flat, load_type = "load")
d_preds_lin = fit_and_predict(l_load, l_load_example, injury = injury_lin, load_type = "load")

# function for creating the same type of plot
plot_sim_shapes = function(d, d_coefs = d_coefs_load %>% filter(label == "n = 8494"), coefs, x_label = "sRPE (AU)"){
  coefs = enexpr(coefs)
  
  text_size = 14
  
  ggplot(data = d, aes(x = load)) +
    facet_wrap(~method, ncol = 4) + 
    geom_line(data = d_coefs, aes(y = !!coefs, color = "True Risk"), size = 0.5) +
    geom_ribbon(aes(min = lower, max = upper), alpha = 0.3, fill = nih_distinct[3]) + 
    geom_line(aes(y = yhat, color = "Estimation"), size = 0.5) + 
    xlab(x_label) + ylab('Probability\nof Injury') + 
    scale_color_manual(values = c("Green", "Black")) + 
    theme(legend.title = element_blank(), legend.position = "bottom",
          axis.text = element_text(size=text_size),
          strip.text.x = element_text(size = text_size),
          axis.title =  element_text(size=text_size),
          legend.text = element_text(size = text_size))
  
}

# figure used in article below
# it won't be 100% reproduced, because colors and themes in the figures 
# used in the article are from local self-made functions
# but the results should be the same
method_arrangement = methods %>% enframe() %>% mutate(arrangement = c(6, 7,5 ,1 ,4 ,2 ,3)) %>% arrange(arrangement)
d_preds_u = d_preds_u %>% mutate(method = factor(method, levels = method_arrangement$value))
d_coefs_load = d_coefs_load %>% mutate(method = factor(method, levels = method_arrangement$value))

png("u-shape-presentation.png", width = 12, height = 6, units="in", res=600)
plot_sim_shapes(d_preds_u %>% filter(label == "n = 8494"), coefs = coefs_u)
dev.off()

# For visualizing all shapes
# U shape
plot_sim_shapes(d_preds_u %>% filter(label == "n = 8494"), coefs = coefs_u)
# J shape
plot_sim_shapes(d_preds_j %>% filter(label == "n = 8494"), 'Gabbett 2016 J shape', d_coefs_acwr, coefs_j_acwr, x_label= "ACWR 7:21")
# linear shape
plot_sim_shapes(d_preds_lin %>% filter(label == "n = 8494"), 'Linear relationship', coefs = coefs_lin)
# flat shape
plot_sim_shapes(d_preds_flat %>% filter(label == "n = 8494"), 'Flat shape (no relationship) sRPE', coefs = coefs_flat)

