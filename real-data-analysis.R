
# This script includes everything performed in the real data
# analysis of the non-linearity study
# script is irreproducible, as data could not be
# made available in necessary format due to anonymisation
# The script is uploaded for the sake of transparency.

# Packages
library(tidyverse) # for everything
library(rlang) # for tidy evaluation
library(ostrc) # self-made micro cycle functions
library(lmisc) # for NIH colors
library(rms) # for restriced cubic splines
library(lme4) # for random effects modeling
library(slider) # for calculating mean and ewma over sliding windows (of equal or unequal size)
library(mice) # for pooling an imputed model
library(merTools) # tools for working with merMod objects from lme4 package


# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors:
options(scipen = 30, 
        stringsAsFactors = FALSE)

# loading imputed datasets
# code is only for illustration
# l_u19, l_handball and l_strom are lists of imputated datasets that have been imputed beforehand
folder_base = paste0("our\\data\\folder\\")
l_u19 = readRDS(paste0(folder_base, "u19_imputed.RDS"))
l_handball = readRDS(paste0(folder_base, "handball_imputed.RDS"))
l_strom = readRDS(paste0(folder_base, "strom_imputed.RDS"))


# make sure sex is treated as categorical
# we also create it as a dummy variable 0 = female, 1 = male
l_u19 = l_u19 %>% map(. %>% mutate(sex = ifelse(sex == 2, 0, 1), sex_fac = factor(sex)))
l_handball = l_handball %>% map(. %>% mutate(sex_fac = factor(sex)))

# list with just covariates we append later
l_covariates_u19 = l_u19 %>% map(. %>% distinct(p_id, .keep_all = TRUE) %>% dplyr::select(p_id, sex_fac, age))
l_covariates_handball = l_handball %>% map(. %>% distinct(p_id, .keep_all = TRUE) %>% dplyr::select(p_id, sex_fac, age))
l_covariates_strom = l_strom %>% map(. %>% distinct(p_id, .keep_all = TRUE) %>% dplyr::select(p_id, position, age))



#---------------------------------------Training load measures
# Now that we've loaded our imputed data, we calculate training load measures
#------------- 1 load unmodified. Check! This ones done. 

#------------- 2 ACWR over days. 

# function calculates ewma on a sliding window of 21 days
chronic_ewma = function(x){
  l =  slide(x, ~ewma(., 21, 1), .before = 20, step = 1, .complete =TRUE) %>% map(last)
  l = compact(l)
  l = unlist(l)
  l
}

# function calculates mean on a sliding window of 7 days
acute_mean = function(x){
  l =  slide(x, ~mean(.), .before = 6, step = 1, .complete =TRUE)
  l = compact(l)
  l = unlist(l)
  l
}

# function sum of injuries in a 4 day period
inj_sum = function(x){
  l =  slide(x, ~sum(.), .before = 3, step = 1, .complete =TRUE) 
  l = compact(l)
  l = unlist(l)
  l 
}

# function to run a user-specified function on each dataset in a list
function_on_list = function(nested_list, FUN = acute_mean, load = TRUE){
  
  for(i in 1:length(nested_list)){
    l_nest = pluck(nested_list, i)
    
    if(load){
      l_nest$data = l_nest$data %>% map(., ~FUN(.$load))
    } else {
      l_nest$data = l_nest$data %>% map(., ~FUN(.$injury))    
    }
    
    l_unnest = unnest(l_nest, cols = c(data)) %>% mutate(index = 1:n()) 
    nested_list[[i]] = l_unnest 
  }
  nested_list 
}

# nesting datasets for acute and chronic calculations, respectiveley
# the first acute window is from 22-28 days
l_u19_ewma_chronic = l_u19 %>% map(. %>% group_by(p_id) %>% nest) 
# The first 21 days need to be discarded for the first acute window
l_u19_acute_times = l_u19 %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% filter(index > 21))
l_u19_mean_acute = l_u19_acute_times %>% map(. %>% group_by(p_id) %>% nest) 
# The first 29 daysof injury discarded for each individual. 
# We want to look at the injury in the next 4 days, not on the same day + 
# ACWR requires 28 measures before the first injury that can be analyzed
l_u19_inj = l_u19 %>% map(. %>% group_by(p_id) %>% slice(-1:-28) %>% nest) 
# we also want injuries on the same day
l_u19_inj_0 = l_u19 %>% map(. %>% group_by(p_id) %>% slice(-1:-27) %>% dplyr::select(p_id, injury0 = injury, load) %>% mutate(index = 1:n()) %>% ungroup) 

# running each function on each player in each imputed dataset
l_u19_mean_acute = function_on_list(l_u19_mean_acute, FUN = acute_mean) %>% map(. %>% rename(acute_mean = data)) 
l_u19_ewma_chronic = function_on_list(l_u19_ewma_chronic, FUN = chronic_ewma) %>% map(. %>% rename(chronic_ewma = data)) 
l_u19_inj = function_on_list(l_u19_inj, FUN = inj_sum, load = FALSE) %>% map(. %>% rename(injury_1_4 = data) %>% mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0))) 

# combine our lists and calculate ACWR
l_u19_acute_chronic = map2(l_u19_mean_acute, l_u19_ewma_chronic, ~left_join(.x, .y, by = c("index", "p_id")))
l_u19_acute_chronic_inj = map2(l_u19_acute_chronic, l_u19_inj, ~left_join(.x, .y, by = c("index", "p_id")))
l_acwr_u19 = l_u19_acute_chronic_inj %>% map(. %>% mutate(acwr_7_21 = ifelse(chronic_ewma == 0, 3, acute_mean/chronic_ewma)))
l_acwr_u19 = map2(l_acwr_u19, l_u19_inj_0, ~left_join(.x, .y, by = c("index", "p_id")))

# Hulin et al. 2014: 
# Weekly workloads that were below 1 SD for the individual's chronic workloads were removed from the analysis. 
# This was performed so that the analysis would not consider small absolute increases of acute workload at 
# low chronic workloads (ie, if a fast bowler had a chronic external workload of six deliveries, a 300% 
# increase would be an acute workload of 18 deliveries).
# code to remove
#high_acwr_pos_u19 = l_acwr_u19 %>% map(~which(.$acwr_7_21 >= 3))
#l_acwr_u19 = map2(l_acwr_u19, high_acwr_pos_u19, ~slice(.x, -.y))
# code to replace
l_acwr_u19 = l_acwr_u19 %>% map(. %>% mutate(acwr_7_21 = ifelse(acwr_7_21 > 3, 3, acwr_7_21)))

# add covariates
l_acwr_u19 = map2(l_acwr_u19, l_covariates_u19, ~left_join(.x, .y, by = c("p_id")))
l_acwr_u19 = l_acwr_u19 %>% map(. %>% ungroup)

#---Handball
# nesting datasets for acute and chronic calculations, respectiveley
l_handball_ewma_chronic = l_handball %>% map(. %>% group_by(p_id) %>% nest) 
l_handball_acute_times = l_handball %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% filter(index > 21))
l_handball_mean_acute = l_handball_acute_times %>% map(. %>% group_by(p_id) %>% nest) 
l_handball_inj= l_handball %>% map(. %>% group_by(p_id) %>% slice(-1:-28) %>% nest) 
l_handball_inj_0 = l_handball %>% map(. %>% group_by(p_id) %>% slice(-1:-27) %>% dplyr::select(p_id, injury0 = injury, load) %>% mutate(index = 1:n()) %>% ungroup) 

# running each function on each player in each imputed dataset
l_handball_mean_acute = function_on_list(l_handball_mean_acute, FUN = acute_mean) %>% map(. %>% rename(acute_mean = data)) 
l_handball_ewma_chronic = function_on_list(l_handball_ewma_chronic, FUN = chronic_ewma) %>% map(. %>% rename(chronic_ewma = data)) 
l_handball_inj = function_on_list(l_handball_inj, FUN = inj_sum, load = FALSE) %>% map(. %>% rename(injury_1_4 = data) %>% mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0))) 

# combine our lists and calculate ACWR
l_handball_acute_chronic = map2(l_handball_mean_acute, l_handball_ewma_chronic, ~left_join(.x, .y, by = c("index", "p_id")))
l_handball_acute_chronic_inj = map2(l_handball_acute_chronic, l_handball_inj, ~left_join(.x, .y, by = c("index", "p_id")))
l_acwr_handball = l_handball_acute_chronic_inj %>% map(. %>% mutate(acwr_7_21 = ifelse(chronic_ewma == 0, 3, acute_mean/chronic_ewma)))
l_acwr_handball = map2(l_acwr_handball, l_handball_inj_0, ~left_join(.x, .y, by = c("index", "p_id")))
l_acwr_handball = l_acwr_handball %>% map(. %>% ungroup)

# add covariates
l_acwr_handball = map2(l_acwr_handball, l_covariates_handball, ~left_join(.x, .y, by = c("p_id")))
l_acwr_handball = l_acwr_handball %>% map(. %>% ungroup)

# fix insane outliers
l_acwr_handball = l_acwr_handball %>% map(. %>% mutate(acwr_7_21 = ifelse(acwr_7_21 > 3, 3, acwr_7_21)))


#---stromsgodset
# nesting datasets for acute and chronic calculations, respectiveley
l_strom_ewma_chronic = l_strom %>% map(. %>% group_by(p_id) %>% nest) 
l_strom_acute_times = l_strom %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% filter(index > 21))
l_strom_mean_acute = l_strom_acute_times %>% map(. %>% group_by(p_id) %>% nest) 
l_strom_inj= l_strom %>% map(. %>% group_by(p_id) %>% slice(-1:-28) %>% nest) 
l_strom_inj_0 = l_strom %>% map(. %>% group_by(p_id) %>% slice(-1:-27) %>% dplyr::select(p_id, injury0 = injury, load) %>% mutate(index = 1:n()) %>% ungroup) 

# running each function on each player in each imputed dataset
l_strom_mean_acute = function_on_list(l_strom_mean_acute, FUN = acute_mean) %>% map(. %>% rename(acute_mean = data)) 
l_strom_ewma_chronic = function_on_list(l_strom_ewma_chronic, FUN = chronic_ewma) %>% map(. %>% rename(chronic_ewma = data)) 
l_strom_inj = function_on_list(l_strom_inj, FUN = inj_sum, load = FALSE) %>% map(. %>% rename(injury_1_4 = data) %>% mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0))) 

# combine our lists and calculate ACWR
l_strom_acute_chronic = map2(l_strom_mean_acute, l_strom_ewma_chronic, ~left_join(.x, .y, by = c("index", "p_id")))
l_strom_acute_chronic_inj = map2(l_strom_acute_chronic, l_strom_inj, ~left_join(.x, .y, by = c("index", "p_id")))
l_acwr_strom = l_strom_acute_chronic_inj %>% map(. %>% mutate(acwr_7_21 = ifelse(chronic_ewma == 0, 3, acute_mean/chronic_ewma)))
l_acwr_strom = map2(l_acwr_strom, l_strom_inj_0, ~left_join(.x, .y, by = c("index", "p_id")))
l_acwr_strom = l_acwr_strom %>% map(. %>% ungroup)

# add covariates
l_acwr_strom = map2(l_acwr_strom, l_covariates_strom, ~left_join(.x, .y, by = c("p_id")))
l_acwr_strom = l_acwr_strom %>% map(. %>% ungroup)

# fix insane outliers
l_acwr_strom = l_acwr_strom %>% map(. %>% mutate(acwr_7_21 = ifelse(acwr_7_21 > 3, 3, acwr_7_21)))


####---------- 3 ACWR over micro cycles.
l_u19 = l_u19 %>% map(. %>% mutate(match_cycle = as.numeric(match_cycle)))
l_handball = l_handball %>% map(. %>% mutate(match_cycle = as.numeric(match_cycle)))
l_strom = l_strom %>% map(. %>% mutate(match_cycle = as.numeric(match_cycle)))

# calculating ACWR over micro cycles per player using micro_cycle_loads from ostrc-package
# could be done much easier with a group_by(p_id), but the for-loop ensures we keep our p_id variable in the end
l_u19_acwr_mc = l_u19 %>% map(. %>% group_by(p_id) %>% nest)
for(i in 1:length(l_u19_acwr_mc)){
  l_u19_nest = pluck(l_u19_acwr_mc, i)
  l_u19_nest$data = l_u19_nest$data %>% map(., ~micro_cycle_loads(., load, match_cycle))
  l_u19_unnest = unnest(l_u19_nest, cols = c(data))
  l_u19_acwr_mc[[i]] = l_u19_unnest
}

l_handball_acwr_mc = l_handball %>% map(. %>% group_by(p_id) %>% nest)
for(i in 1:length(l_handball_acwr_mc)){
  l_handball_nest = pluck(l_handball_acwr_mc, i)
  l_handball_nest$data = l_handball_nest$data %>% map(., ~micro_cycle_loads(., load, match_cycle))
  l_handball_unnest = unnest(l_handball_nest, cols = c(data))
  l_handball_acwr_mc[[i]] = l_handball_unnest
}

l_strom_acwr_mc = l_strom %>% map(. %>% group_by(p_id) %>% nest)
for(i in 1:length(l_strom_acwr_mc)){
  l_strom_nest = pluck(l_strom_acwr_mc, i)
  l_strom_nest$data = l_strom_nest$data %>% map(., ~micro_cycle_loads(., load, match_cycle))
  l_strom_unnest = unnest(l_strom_nest, cols = c(data))
  l_strom_acwr_mc[[i]] = l_strom_unnest
}


# we need to add injuries to do any modeling
# this required a bit of code, made function to easier replicate over more than 1 list of datasets
obtain_injury_nextmc = function(list, n_cycles){
  
  l_inj = list %>% map(. %>% dplyr::select(p_id, date_training, match_cycle, injury))
  
  # run a distinct on match_cycles and make sure those with injury have the injury denoted
  l_mc_distinct = l_inj %>% map(. %>% group_by(p_id, match_cycle) %>% arrange(desc(injury), date_training) %>% 
                                  distinct(match_cycle, .keep_all = TRUE) %>% 
                                  arrange(p_id, date_training, match_cycle) %>% 
                                  ungroup()) 
  
  # in this case, ACWR is 1:3 match cycles, where 3 is uncoupled with the 4th. 
  # Therefore, the first acute load is the 4th match_cycle. 
  # adding injuries for the next match_cycle in the first
  l_inj_nextmc = l_mc_distinct %>% 
    map(. %>% filter(!match_cycle %in% 1:n_cycles) %>% group_by(p_id) %>%  
          arrange(p_id, date_training, match_cycle) %>% 
          mutate(injury_nextmc = lead(injury),
                 index = 1:n()) %>% 
          ungroup() %>% 
          dplyr::select(-date_training, -match_cycle, -injury))
  l_inj_nextmc
}
l_inj_nextmc_u19 = obtain_injury_nextmc(l_u19, 3)
l_inj_nextmc_handball = obtain_injury_nextmc(l_handball, 3)
l_inj_nextmc_strom = obtain_injury_nextmc(l_strom, 3)

# last step: using map2 to combine two lists.
l_acwr_mp_u19 = map2(l_u19_acwr_mc, l_inj_nextmc_u19, ~left_join(.x, .y, by = c("p_id", "index")))
l_acwr_mp_handball = map2(l_handball_acwr_mc, l_inj_nextmc_handball, ~left_join(.x, .y, by = c("p_id", "index")))
l_acwr_mp_strom = map2(l_strom_acwr_mc, l_inj_nextmc_strom, ~left_join(.x, .y, by = c("p_id", "index")))

# add covariates
l_acwr_mp_u19 = map2(l_acwr_mp_u19, l_covariates_u19, ~left_join(.x, .y, by = c("p_id"))) %>% map(. %>% ungroup)
l_acwr_mp_handball = map2(l_acwr_mp_handball, l_covariates_handball, ~left_join(.x, .y, by = c("p_id"))) %>%  map(. %>% ungroup)
l_acwr_mp_strom = map2(l_acwr_mp_strom, l_covariates_strom, ~left_join(.x, .y, by = c("p_id"))) %>%  map(. %>% ungroup)


# We add lagged injuries in the time period 1 to 4 days after training
l_u19_inj_1_4 = l_u19 %>% map(. %>% group_by(p_id) %>% slice(-1) %>% nest)
l_u19_inj_1_4_names = function_on_list(l_u19_inj_1_4, FUN = inj_sum, load = FALSE) %>%
  map(. %>% rename(injury_1_4 = data) %>%
        mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0)))
l_u19 = l_u19 %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% ungroup)
l_u19 = map2(l_u19, l_u19_inj_1_4_names, ~left_join(.x, .y, by = c("p_id", "index")))

l_handball_inj_1_4 = l_handball %>% map(. %>% group_by(p_id) %>% slice(-1) %>% nest)
l_handball_inj_1_4_names = function_on_list(l_handball_inj_1_4, FUN = inj_sum, load = FALSE) %>%
  map(. %>% rename(injury_1_4 = data) %>%
        mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0)))
l_handball = l_handball %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% ungroup)
l_handball = map2(l_handball, l_handball_inj_1_4_names, ~left_join(.x, .y, by = c("p_id", "index")))

l_strom_inj_1_4 = l_strom %>% map(. %>% group_by(p_id) %>% slice(-1) %>% nest)
l_strom_inj_1_4_names = function_on_list(l_strom_inj_1_4, FUN = inj_sum, load = FALSE) %>%
  map(. %>% rename(injury_1_4 = data) %>%
        mutate(injury_1_4 = ifelse(injury_1_4 >= 1, 1, 0)))
l_strom = l_strom %>% map(. %>% group_by(p_id) %>% mutate(index = 1:n()) %>% ungroup)
l_strom = map2(l_strom, l_strom_inj_1_4_names, ~left_join(.x, .y, by = c("p_id", "index")))

#-------------------------------------------------- Remove days the players did not train (no possible injury risk)


l_u19_training_only = l_u19 %>% map(. %>% filter(load != 0))
l_acwr_u19_training_only = l_acwr_u19  %>% map(. %>% filter(load != 0) %>% dplyr::select(-load))
l_acwr_u19 = l_acwr_u19  %>% map(. %>% dplyr::select(-load))

l_handball_training_only = l_handball %>% map(. %>% filter(load != 0))
l_acwr_handball_training_only = l_acwr_handball  %>% map(. %>% filter(load != 0) %>% dplyr::select(-load))
l_acwr_handball = l_acwr_handball  %>% map(. %>% dplyr::select(-load))

l_strom_training_only = l_strom %>% map(. %>% filter(load != 0))
l_acwr_strom_training_only = l_acwr_strom  %>% map(. %>% filter(load != 0) %>% dplyr::select(-load))
l_acwr_strom = l_acwr_strom  %>% map(. %>% dplyr::select(-load))

#-------------------------------------------------- checking for overly influential values
d_test = l_u19$`5`
f  <- lrm(injury ~ rcs(age, 3) + rcs(load, 3), data = d_test, x = TRUE, y = TRUE)
w <- which.influence(f, .2)
nam <- names(w)
for(i in 1:length(nam)) {
  print(paste("Influential observations for effect of ",nam[i]),quote=FALSE)
  print(d_test[w[[i]],])
}

d_test = l_handball$`5`
f  <- lrm(injury ~ rcs(age, 3) + rcs(load, 3), data = d_test, x = TRUE, y = TRUE)
w <- which.influence(f, .2)
nam <- names(w)
for(i in 1:length(nam)) {
  print(paste("Influential observations for effect of ",nam[i]),quote=FALSE)
  print(d_test[w[[i]],])
}

l_handball = l_handball %>% map(. %>% mutate(load = ifelse(load > 3500, 3500, load)))

d_test = l_strom$`5`
f  <- lrm(injury ~ rcs(age, 3) + rcs(load, 3), data = d_test, x = TRUE, y = TRUE)
w <- which.influence(f, .2)
nam <- names(w)
for(i in 1:length(nam)) {
  print(paste("Influential observations for effect of ",nam[i]),quote=FALSE)
  print(d_test[w[[i]],])
}

#-------------------------------------------------- Checking whether the ACWR is a fesible ratio
# l_acwr_handball and l_acwr_u19 has daily ACWR. l_acwr_mp_u19 and l_acwr_mp_handball have micro cycle ACWR.

#The numerator Y (left) and the ratio Y/X (right) against the denominator X for four models (A–D). 
# Only when a straight-line relationship between Y and X goes through the origin (model A) 

d_acwr_measures = bind_rows(l_acwr_u19$`5` %>% mutate(data = "Football U-19", type = "Daily") %>% rename(chronic_load = chronic_ewma, acwr = acwr_7_21), 
                            l_acwr_handball$`5` %>% mutate(data = "Handball", type = "Daily") %>% rename(chronic_load = chronic_ewma, acwr = acwr_7_21),
                            l_acwr_strom$`5` %>% mutate(data = "Football Elite", type = "Daily") %>% rename(chronic_load = chronic_ewma, acwr = acwr_7_21),
                            l_acwr_mp_u19$`5` %>% mutate(data = "Football U-19", type = "Microcycle"), 
                            l_acwr_mp_handball$`5` %>% mutate(data = "Handball", type = "Microcycle"),
                            l_acwr_mp_strom$`5` %>% mutate(data = "Football Elite", type = "Microcycle"))
# note that the figure won't be reproduced in the same colors
library(devEMF)
devEMF::emf("acwr_ratio_assumption.emf", width = 8, height = 6)
ggplot(d_acwr_measures, aes(x = chronic_load, y = acwr)) +
  facet_wrap(c("type", "data")) +
  geom_point() + 
  ylab("ACWR") +
  xlab("Chronic Load") +
  geom_hline(yintercept=1) + 
  theme_base()
dev.off()


#--------------------------------------------------------------------Modeling-------------------------------------------------------------------

# Remember our objects. l_handball and l_u19 have unmodified loads, and difference-adjusted loads.
# l_acwr_handball and l_acwr_u19 has daily ACWR. l_acwr_mp_u19 and l_acwr_mp_handball have micro cycle ACWR.

# Function for finding the logliks for each knot value between 3-6. 
# reports the logliks for each
find_optimal_knot = function(d, tl, injury){
  tl = enexpr(tl)
  injury = enexpr(injury)
  logliks = rep(NA, 3)
  knots = 3:7
  for (i in seq_along(knots)) {
    knot = knots[i]
    fit_splines = eval_bare(expr(lrm(!!injury ~ rcs(!!tl, knot) + age  + position , data = d, method = "lrm.fit"))) 
    logliks[i] = logLik(fit_splines)
  }
  d_knots = tibble(knots, logliks)
  knot_optimal = d_knots %>% filter(logliks == min(logliks)) %>% pull("knots")
  knot_optimal
}

#Adjusting for age and sex in all models, perform the following:
# •	GLMM with binomial distribution with load unmodified (linearity assumed) with random intercept and random slope
# •	GLMM with binomial distribution with load unmodified (linearity assumed) with random intercept only
# •	GLMM with binomial distribution with load modified by splines, with random intercept and random slope (k number of knots decided by AIC) 
# •	GLMM with binomial distribution with load modified by splines, with random intercept only

# function for obtaining model fits for 5 model levels
# needs a dataset and the training load variable in the dataset
find_opt_model = function(d, tl, injury){
  tl = enexpr(tl)
  injury = enexpr(injury)
  n_knot = find_optimal_knot(d, !!tl, !!injury)
  
  # base models without training load measures. 
  fit_mm_base = eval_bare(expr(glmer(!!injury ~ rcs(age, 3) + position + (1 | p_id), data = d, family = "binomial")))
  fit_mm_base_slope = eval_bare(expr(glmer(!!injury ~  rcs(age, 3) + position + (!!tl | p_id), data = d, family = "binomial")))
  
  # without splines
  fit_mm_intercept = eval_bare(expr(update(fit_mm_base, . ~ . + !!tl)))
  fit_mm_slope = eval_bare(expr(update(fit_mm_base_slope, . ~ . + !!tl)))
  
  # with splines
  fit_mm_splines_intercept = eval_bare(expr(update(fit_mm_base, . ~ . + rcs(!!tl, n_knot))))
  fit_mm_splines_slope = eval_bare(expr(update(fit_mm_base_slope, . ~ . + rcs(!!tl, n_knot))))
  
  anova = anova(fit_mm_splines_slope, fit_mm_splines_intercept, fit_mm_slope, fit_mm_intercept)
  anova = rlang::set_attrs(anova, knot = paste0("Modelled with ",n_knot," knots."))
  anova
}

# find best models

models_load_injury0_u19 = l_u19_training_only %>% map(~find_opt_model(., load, injury))
models_acwr_injury0_u19 = l_acwr_u19_training_only %>% map(~find_opt_model(., acwr_7_21, injury0))
models_load_u19 = l_u19 %>% map(~find_opt_model(., load, injury_1_4))
models_acwr_u19 = l_acwr_u19 %>% map(~find_opt_model(., acwr_7_21, injury_1_4))
models_mpacwr_u19 = l_acwr_mp_u19 %>% map(~find_opt_model(., acwr, injury_nextmc))
attributes(models_load_injury0_u19$`5`)
attributes(models_acwr_injury0_u19$`5`)
attributes(models_load_u19$`5`)
attributes(models_acwr_u19$`5`)
attributes(models_mpacwr_u19$`5`)

models_load_injury0_handball = l_handball_training_only %>% map(~find_opt_model(., load, injury))
models_acwr_injury0_handball = l_acwr_handball_training_only %>% map(~find_opt_model(., acwr_7_21, injury0))
models_load_handball = l_handball %>% map(~find_opt_model(., load, injury_1_4))
models_acwr_handball = l_acwr_handball %>% map(~find_opt_model(., acwr_7_21, injury_1_4))
models_mpacwr_handball = l_acwr_mp_handball %>% map(~find_opt_model(., acwr, injury_nextmc))
attributes(models_load_injury0_handball$`5`)
attributes(models_acwr_injury0_handball$`5`)
attributes(models_load_handball$`5`)
attributes(models_acwr_handball$`5`)
attributes(models_mpacwr_handball$`5`)

models_load_injury0_strom = l_strom_training_only %>% map(~find_opt_model(., load, injury))
models_acwr_injury0_strom = l_acwr_strom_training_only %>% map(~find_opt_model(., acwr_7_21, injury0))
models_load_strom = l_strom %>% map(~find_opt_model(., load, injury_1_4))
models_acwr_strom = l_acwr_strom %>% map(~find_opt_model(., acwr_7_21, injury_1_4))
models_mpacwr_strom = l_acwr_mp_strom %>% map(~find_opt_model(., acwr, injury_nextmc))
attributes(models_load_injury0_strom$`5`)
attributes(models_acwr_injury0_strom$`5`)
attributes(models_load_strom$`5`)
attributes(models_acwr_strom$`5`)
attributes(models_mpacwr_strom$`5`)

#---------------- perform best model on each load-value

# u19, injury on the same day as training load
fit_load_injury0_u19 = l_u19_training_only %>% map(glmer, formula = injury ~ rcs(load, 3) + sex_fac + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_injury0_u19 = l_acwr_u19_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~  rcs(load, 3) + sex_fac + age + (1 | p_id), family = "binomial") 

# u19 injury in the next 1-4 days after training, or in the next micro cycle
fit_load_u19 = l_u19 %>% map(glmer, formula = injury_1_4 ~ rcs(load, 3) + sex_fac + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_u19 = l_acwr_u19 %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~ rcs(load, 3) + sex_fac + age + (load | p_id), family = "binomial") 
fit_acwr_mp_u19 = l_acwr_mp_u19 %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~ rcs(load, 4) + sex_fac + age + (1 | p_id), family = "binomial") 

# handball, injury0
fit_load_injury0_handball = l_handball_training_only %>% map(glmer, formula = injury ~ rcs(load, 3) + sex_fac  + age + (1 | p_id), family = "binomial")
fit_acwr_7_21_injury0_handball = l_acwr_handball_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~ rcs(load, 3) + sex_fac + age + (1 | p_id), family = "binomial")

# handball, injury 1-4
fit_load_handball = l_handball %>% map(glmer, formula = injury_1_4 ~ rcs(load,3) + sex_fac  + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_handball = l_acwr_handball %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~  rcs(load, 3) + sex_fac  + age + (load | p_id), family = "binomial") 
fit_acwr_mp_handball = l_acwr_mp_handball %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~  rcs(load, 3) + sex_fac + age + (1 | p_id), family = "binomial") 

# stroml, injury0
fit_load_injury0_strom = l_strom_training_only %>% map(glmer, formula = injury ~ rcs(load, 3) + age + (1 | p_id), family = "binomial")
fit_acwr_7_21_injury0_strom = l_acwr_strom_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~ rcs(load, 3) + age + (1 | p_id), family = "binomial")

# strom, injury 1-4
fit_load_strom = l_strom %>% map(glmer, formula = injury_1_4 ~ rcs(load,3) + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_strom = l_acwr_strom %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~  age  + rcs(load,3) + (1 | p_id), family = "binomial") 
fit_acwr_mp_strom = l_acwr_mp_strom %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~  age + rcs(load,3) + (1 | p_id), family = "binomial") 



# perform models assuming linearity on each load-value

# u19, injury on the same day as training load
fit_load_injury0_u19_linear = l_u19_training_only %>% map(glmer, formula = injury ~ load + sex_fac + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_injury0_u19_linear = l_acwr_u19_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~  load + sex_fac + age + (1 | p_id), family = "binomial") 

# u19 injury in the next 1-4 days after training, or in the next micro cycle
fit_load_u19_linear = l_u19 %>% map(glmer, formula = injury_1_4 ~ load + sex_fac + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_u19_linear = l_acwr_u19 %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~ load + sex_fac + age + (load | p_id), family = "binomial") 
fit_acwr_mp_u19_linear = l_acwr_mp_u19 %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~ load + sex_fac + age + (1 | p_id), family = "binomial") 

# handball, injury0
fit_load_injury0_handball_linear = l_handball_training_only %>% map(glmer, formula = injury ~ load + sex_fac + age + (1 | p_id), family = "binomial")
fit_acwr_7_21_injury0_handball_linear = l_acwr_handball_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~ load + sex_fac  + age + (1 | p_id), family = "binomial")

# handball, injury 1-4
fit_load_handball_linear = l_handball %>% map(glmer, formula = injury_1_4 ~ load + sex_fac + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_handball_linear = l_acwr_handball %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~  load + sex_fac + age + (load | p_id), family = "binomial") 
fit_acwr_mp_handball_linear = l_acwr_mp_handball %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~  load + sex_fac + age + (1 | p_id), family = "binomial") 

# strom, injury0
fit_load_injury0_strom_linear = l_strom_training_only %>% map(glmer, formula = injury ~ load + age + (1 | p_id), family = "binomial")
fit_acwr_7_21_injury0_strom_linear = l_acwr_strom_training_only %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury0 ~ load + age + (1 | p_id), family = "binomial")

# strom, injury 1-4
fit_load_strom_linear = l_strom %>% map(glmer, formula = injury_1_4 ~ load + age + (1 | p_id), family = "binomial") 
fit_acwr_7_21_strom_linear = l_acwr_strom %>% map(. %>% rename(load = acwr_7_21)) %>% map(glmer, formula = injury_1_4 ~  load  + age + (1 | p_id), family = "binomial") 
fit_acwr_mp_strom_linear = l_acwr_mp_strom %>% map(. %>% rename(load = acwr)) %>% map(glmer, formula = injury_nextmc ~  load + age + (1 | p_id), family = "binomial") 


# We want predicted values

# we can either use the "predict than combine"-method https://journals.sagepub.com/doi/pdf/10.1177/0049124115610345
# or cheese coefficients to a regular merMod object https://stackoverflow.com/questions/52713733/how-to-use-predict-function-with-my-pooled-results-from-mice
# This is the predict-than-combine method
# function for predicting values on a list of fitted models om imputed datasets
# the output is the mean predicted value of each prediction
predict_on_imputed = function(fit, d_example){
  
  l_pred = fit %>% map(predict, type = "response", newdata = d_example)
  l_conf = fit %>% map(predictInterval, newdata = d_example, ignore.fixed.terms = 1, type = "probability")
  
  # unwrap our list to a table
  d_pred = do.call(bind_cols, l_pred) %>% as_tibble()
  d_conf = do.call(bind_cols, l_conf) %>% as_tibble()
  
  # calulate mean predicted value
  d_pred_means = d_conf %>% transmute(load_values = d_example$load,
                                      lower = rowMeans(dplyr::select(., starts_with("lwr"))),
                                      upper = rowMeans(dplyr::select(., starts_with("upr"))),
                                      yhat = rowMeans(d_pred))
  d_pred_means
  
}

# predicted values adjusted for age and sex. Age is set to mean age. sex = 1 male.
d_base =  tibble(
  age = rep(17.22, 25),
  sex_fac = rep("0", 25),
  p_id = rep(101, 25)
)

d_base_strom =  tibble(
  age = rep(17.22, 25),
  sex_fac = rep("0", 25),
  p_id = rep(100, 25),
  position = rep("Striker", 25)
)

fake_loads = function(x){
  seq = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out=25)
  seq
}

# make sample data
d_sample_loads_u19 = d_base %>% mutate(load = seq(min(l_u19$`5`$load, na.rm = TRUE), 3000, length.out=25))
d_sample_loads_handball = d_base %>% mutate(load = fake_loads(l_handball$`5`$load))
d_sample_loads_strom = d_base_strom %>% mutate(load = fake_loads(l_strom$`5`$load) )
d_sample_acwr = d_base %>% mutate(load = seq(min(l_u19$`5`$load, na.rm = TRUE), 2.5, length.out=25))
d_sample_acwr_strom = d_base_strom %>% mutate(load = seq(min(l_u19$`5`$load, na.rm = TRUE), 2.5, length.out=25))

# u19, injury 0
d_preds_load_injury0_u19 = predict_on_imputed(fit_load_injury0_u19, d_sample_loads_u19) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_u19 = predict_on_imputed(fit_acwr_7_21_injury0_u19, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")

# u19, injury 1-4 or next mp
d_preds_load_u19 = predict_on_imputed(fit_load_u19, d_sample_loads_u19) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_u19 = predict_on_imputed(fit_acwr_7_21_u19, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_u19 = predict_on_imputed(fit_acwr_mp_u19, d_sample_acwr) %>% mutate(load_type = "ACWR 1:3 micro cycles")

# handball, injury 0
d_preds_load_injury0_handball = predict_on_imputed(fit_load_injury0_handball, d_sample_loads_handball) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_handball = predict_on_imputed(fit_acwr_7_21_injury0_handball, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")

# handball, injury 1-4
d_preds_load_handball = predict_on_imputed(fit_load_handball, d_sample_loads_handball) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_handball = predict_on_imputed(fit_acwr_7_21_handball, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_handball = predict_on_imputed(fit_acwr_mp_handball, d_sample_acwr) %>% mutate(load_type = "ACWR 1:3 micro cycles")

# strom, injury 0
d_preds_load_injury0_strom = predict_on_imputed(fit_load_injury0_strom, d_sample_loads_strom) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_strom = predict_on_imputed(fit_acwr_7_21_injury0_strom, d_sample_acwr_strom)  %>% mutate(load_type = "ACWR 7:21")

# strom, injury 1-4
d_preds_load_strom = predict_on_imputed(fit_load_strom, d_sample_loads_strom) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_strom = predict_on_imputed(fit_acwr_7_21_strom, d_sample_acwr_strom)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_strom = predict_on_imputed(fit_acwr_mp_strom, d_sample_acwr_strom) %>% mutate(load_type = "ACWR 1:3 micro cycles")



# combine u19 and handball data
d_preds_u19 = bind_rows(d_preds_load_u19, d_preds_acwr_7_21_u19, d_preds_acwr_mp_u19) %>% mutate(pop = "Football U19")
d_preds_handball = bind_rows(d_preds_load_handball, d_preds_acwr_7_21_handball, d_preds_acwr_mp_handball) %>% mutate(pop = "Handball")
d_preds_strom = bind_rows(d_preds_load_strom, d_preds_acwr_7_21_strom, d_preds_acwr_mp_strom) %>% mutate(pop = "Football\nStrømsgodset")
d_preds = bind_rows(d_preds_u19, d_preds_handball, d_preds_strom)

d_preds_injury0_u19 = bind_rows(d_preds_load_injury0_u19, d_preds_acwr_7_21_injury0_u19) %>% mutate(pop = "Football U19")
d_preds_injury0_handball = bind_rows(d_preds_load_injury0_handball, d_preds_acwr_7_21_injury0_handball) %>% mutate(pop = "Handball")
d_preds_injury0_strom = bind_rows(d_preds_load_injury0_strom, d_preds_acwr_7_21_injury0_strom) %>% mutate(pop = "Football\nStrømsgodset")
d_preds_injury0 = bind_rows(d_preds_injury0_u19, d_preds_injury0_handball)



#--------------------------------------------linear models

# u19, injury 0
d_preds_load_injury0_u19_linear = predict_on_imputed(fit_load_injury0_u19_linear, d_sample_loads_u19) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_u19_linear = predict_on_imputed(fit_acwr_7_21_injury0_u19_linear, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")

# u19, injury 1-4 or next mp
d_preds_load_u19_linear = predict_on_imputed(fit_load_u19_linear, d_sample_loads_u19) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_u19_linear = predict_on_imputed(fit_acwr_7_21_u19_linear, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_u19_linear = predict_on_imputed(fit_acwr_mp_u19_linear, d_sample_acwr) %>% mutate(load_type = "ACWR 1:3 micro cycles")

# handball, injury 0
d_preds_load_injury0_handball_linear = predict_on_imputed(fit_load_injury0_handball_linear, d_sample_loads_handball) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_handball_linear = predict_on_imputed(fit_acwr_7_21_injury0_handball_linear, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")

# handball, injury 1-4
d_preds_load_handball_linear = predict_on_imputed(fit_load_handball_linear, d_sample_loads_handball) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_handball_linear = predict_on_imputed(fit_acwr_7_21_handball_linear, d_sample_acwr)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_handball_linear = predict_on_imputed(fit_acwr_mp_handball_linear, d_sample_acwr) %>% mutate(load_type = "ACWR 1:3 micro cycles")

# strom, injury 0
d_preds_load_injury0_strom_linear = predict_on_imputed(fit_load_injury0_strom_linear, d_sample_loads_strom) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_injury0_strom_linear = predict_on_imputed(fit_acwr_7_21_injury0_strom_linear, d_sample_acwr_strom)  %>% mutate(load_type = "ACWR 7:21")

# strom, injury 1-4
d_preds_load_strom_linear = predict_on_imputed(fit_load_strom_linear, d_sample_loads_strom) %>% mutate(load_type = "sRPE")
d_preds_acwr_7_21_strom_linear = predict_on_imputed(fit_acwr_7_21_strom_linear, d_sample_acwr_strom)  %>% mutate(load_type = "ACWR 7:21")
d_preds_acwr_mp_strom_linear = predict_on_imputed(fit_acwr_mp_strom_linear, d_sample_acwr_strom) %>% mutate(load_type = "ACWR 1:3 micro cycles")

# combine football and handball data
d_preds_u19_linear = bind_rows(d_preds_load_u19_linear, d_preds_acwr_7_21_u19_linear, d_preds_acwr_mp_u19_linear) %>% mutate(pop = "Football U-19")
d_preds_handball_linear = bind_rows(d_preds_load_handball_linear, d_preds_acwr_7_21_handball_linear, d_preds_acwr_mp_handball_linear) %>% mutate(pop = "Handball")
d_preds_strom_linear = bind_rows(d_preds_load_strom_linear, d_preds_acwr_7_21_strom_linear, d_preds_acwr_mp_strom_linear) %>% mutate(pop = "Football Elite")
d_preds_linear = bind_rows(d_preds_u19_linear, d_preds_handball_linear, d_preds_strom_linear)

d_preds_injury0_u19_linear = bind_rows(d_preds_load_injury0_u19_linear, d_preds_acwr_7_21_injury0_u19_linear) %>% mutate(pop = "Football U-19")
d_preds_injury0_handball_linear = bind_rows(d_preds_load_injury0_handball_linear, d_preds_acwr_7_21_injury0_handball_linear) %>% mutate(pop = "Handball")
d_preds_injury0_strom_linear = bind_rows(d_preds_load_injury0_strom_linear, d_preds_acwr_7_21_injury0_strom_linear) %>% mutate(pop = "Football Elite")
d_preds_injury0_linear = bind_rows(d_preds_injury0_u19_linear, d_preds_injury0_handball_linear, d_preds_injury0_strom_linear)




#------------------------Figures with all results (also flat results) for appendix
# renaming labels to be consistent with terms in the article
d_preds_full = d_preds %>% filter(pop == "Football U-19" & load_values <= 2000 | pop == "Handball" | pop == "Football Elite")
d_preds_full = d_preds_full %>% mutate(load_type = case_when(load_type == "ACWR 1:3 Match periods" ~ "Micro-cycle ACWR 1:3",
                                                             load_type == "ACWR 7:21" ~ "Daily ACWR 7:21",
                                                             TRUE ~ load_type))

d_preds_linear_full = d_preds_linear %>% filter(pop == "Football U-19" & load_values <= 2000 | pop == "Handball" | pop == "Football Elite")
d_preds_linear_full = d_preds_linear_full %>% mutate(load_type = case_when(load_type == "ACWR 1:3 Match periods" ~ "Micro-cycle ACWR 1:3",
                                                                           load_type == "ACWR 7:21" ~ "Daily ACWR 7:21",
                                                                           TRUE ~ load_type))

d_preds_injury0_full = d_preds_injury0 %>% filter(pop == "Football U-19" & load_values <= 2000 | pop == "Handball" | pop == "Football Elite")
d_preds_injury0_full= d_preds_injury0_full %>% mutate(load_type = case_when(load_type == "ACWR 1:3 Match periods" ~ "Micro-cycle ACWR 1:3",
                                                                            load_type == "ACWR 7:21" ~ "Daily ACWR 7:21",
                                                                            TRUE ~ load_type))
d_preds_injury0_linear_full = d_preds_injury0_linear %>% filter(pop == "Football U-19" & load_values <= 2000 | pop == "Handball" | pop == "Football Elite")
d_preds_injury0_linear_full= d_preds_injury0_linear_full %>% mutate(load_type = case_when(load_type == "ACWR 1:3 Match periods" ~ "Micro-cycle ACWR 1:3",
                                                                                          load_type == "ACWR 7:21" ~ "Daily ACWR 7:21",
                                                                                          TRUE ~ load_type))

# note that the figures won't be reproduced with the same colors 
# and theme, but the results are the same
# (the theme used for the article were from self-made R theme functions)
text_size = 12
ostrc_yellow = "#FF9900"
d_preds_full = d_preds %>% filter(pop == "Football U19" & load_values <= 2000 | pop == "Handball" | pop == "Football\nStrømsgodset")
d_preds_linear_full = d_preds_linear %>% filter(pop == "Football U19" & load_values <= 2000 | pop == "Handball" | pop == "Football\nStrømsgodset")
ggplot(d_preds_full, aes(x = load_values, y = yhat)) +
  facet_wrap(vars(load_type, pop), scales = "free", ncol = 3) +
  geom_line(size = 1) + 
  geom_line(data = d_preds_linear_full, aes(x = load_values, y = yhat)) + 
  geom_ribbon(aes(min = lower, max = upper), alpha = 0.25, fill = ostrc_yellow) +
  ylab("Probability of injury\Next 4 days or\nNext Microcycle") +
  xlab("Load value")  +
  scale_y_continuous(limits = c(0, 1.0)) +
  theme(axis.text = element_text(size=text_size),
        strip.text.x = element_text(size = text_size),
        axis.title =  element_text(size=text_size),
        legend.position = "bottom",
        legend.title = element_blank())

d_preds_injury0_full = d_preds_injury0 %>% filter(pop == "Football U19" & load_values <= 2000 | pop == "Handball" | pop == "Football\nStrømsgodset")
d_preds_injury0_linear_full = d_preds_injury0_linear %>% filter(pop == "Football U19" & load_values <= 2000 | pop == "Handball" | pop == "Football\nStrømsgodset")
ggplot(d_preds_injury0_full, aes(x = load_values, y = yhat)) +
  facet_wrap(vars(load_type, pop), scales = "free", ncol = 3) +
  geom_line(size = 1) + 
  geom_line(data = d_preds_injury0_linear_full, aes(x = load_values, y = yhat)) + 
  geom_ribbon(aes(min = lower, max = upper), alpha = 0.25, fill = ostrc_yellow) +
  ylab("Probability of injury\non same day") +
  xlab("Load value")  +
  scale_y_continuous(limits = c(0, 1.0))  +
  theme(axis.text = element_text(size=text_size),
        strip.text.x = element_text(size = text_size),
        axis.title =  element_text(size=text_size),
        legend.position = "bottom",
        legend.title = element_blank())

#------------------------- Figure made for article with handball only
d_preds_load_injury0_handball = d_preds_load_injury0_handball %>% mutate(injury = "(A) Current Day (472 Injuries)")
d_preds_load_handball = d_preds_load_handball %>% mutate(injury = "(B) Next 4 Days (1 136 Injuries)")

d_preds_load_injury0_handball_linear = d_preds_load_injury0_handball_linear %>% mutate(injury = "(A) Current Day (472 Injuries)")
d_preds_load_handball_linear = d_preds_load_handball_linear %>% mutate(injury = "(B) Next 4 Days (1 136 Injuries)")

d_preds_sig = bind_rows(d_preds_load_injury0_handball, d_preds_load_handball) 
d_preds_sig_lin = bind_rows(d_preds_load_injury0_handball_linear, d_preds_load_handball_linear)

# note that the figures won't be reproduced with the same colors 
# and theme, but the results are the same
# (the theme used for the article were from self-made R theme functions)
text_size = 12
# there are so few values to go on beyond 2000 that the confidence intervals are incredibly broad
ggplot(d_preds_sig, aes(x = load_values, y = yhat)) +
  facet_wrap(~injury, ncol = 2) +
  geom_line(size = 1) + 
  geom_line(data = d_preds_sig_lin, aes(x = load_values, y = yhat)) + 
  geom_ribbon(aes(min = lower, max = upper), alpha = 0.25, fill = ostrc_yellow) +
  ylab("Probability\nof Injury") +
  xlab("sRPE (AU)")  +
  scale_y_continuous(limits = c(0, 1.0)) +
  theme(axis.text = element_text(size=text_size),
        strip.text.x = element_text(size = text_size),
        axis.title =  element_text(size=text_size),
        legend.position = "bottom",
        legend.title = element_blank())
