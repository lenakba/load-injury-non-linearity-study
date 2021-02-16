# this scripts show  how the multiple imputation
# was performed in the non-linearity study
# The analyses cannot be reproduced as the data
# are only available as anonymised versions
# where columns of background information etc. are removed
# this information was used to impute the data used in the article
# therefore, this script is mostly available for transparency.
#  the anonymised datasets are not in the GitHub repository, only available
# as supplementary to the non-linearity study

library(tidyverse) # for datareading, wrangling and handling
library(mice) # multiple imputation

# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors
options(scipen = 17, 
        stringsAsFactors = FALSE)

# loading unimputed datasets
# code is only for illustration
folder_base = paste0("our\\data\\folder\\")
d_u19 = read_delim(paste0(folder_base, "data_football_u19_anonymised.csv"), delim = ";")
d_handball = read_delim(paste0(folder_base, "data_handball_anonymised.csv"), delim = ";")
d_strom = read_delim(paste0(folder_base, "data_football_premierdiv_anonymised.csv"), delim = ";")

#------------------------imputation

# the imputation is performed in the default manner of mice, which is predicted mean matching
# note that the data will note replicate the imputed datasets available with the study
# the data had to be anonymized before uploading, this included removing background variables
# that can identify an athlete. All background variables were used in imputing the missing load values
# for the data used in the study.
# most notable, for the premier division data, since they were few, and all are members of the same team
# age could be used to identify an athlete, and had to be removed. Multiple imputation won't work unless
# other variables can be used to predict the load. Only training date can be used for this data.

set.seed(123)
l_imputed_u19 = d_u19 %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") 

l_imputed_handball = d_handball %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") 

l_imputed_strom = d_strom %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") 

# imputation can also be validated using mice::densityplot 
l_imputed_u19 = d_match_cycles_u19 %>%
  mice(seed = 123, print = FALSE)
mice::densityplot(l_imputed_u19, ~load)
bwplot(l_imputed_u19, ~load)


# for each dataset, remove other imputed variables and 
# insert just the variables we wanted to impute: load and age
# this ensures that the injury variable would remain the same, too
tl_cols = c("load", "minutes", "intensity", "age")

d_woload_u19 = d_u19 %>% select(-tl_cols)
l_imputed_u19 = l_imputed_u19 %>% map(. %>% select(tl_cols) %>% bind_cols(d_woload_u19))

d_woload_handball = d_handball %>% select(-tl_cols)
l_imputed_handball = l_imputed_handball %>% map(. %>% select(tl_cols) %>% bind_cols(d_woload_handball))

d_woload_strom = d_strom %>% select(-load, -minutes, -intensity)
l_imputed_strom = l_imputed_strom %>% map(. %>% select(load, minutes, intensity) %>% bind_cols(d_woload_strom))

#------------------------------------------------------ validating imputation
# check that the imputed values are similarly distributed to the imputed values
imp_u19 = mice(d_u19  %>% rename(sRPE = load, Age = age), seed = 123)
densityplot_u19 = densityplot(x=imp_u19, data = ~sRPE + Age)

imp_handball = mice(d_handball  %>% rename(sRPE = load, Age = age), seed = 123)
densityplot_handball = densityplot(x=imp_handball, data = ~sRPE + Age)

imp_strom = mice(d_strom  %>% rename(sRPE = load, Age = age), seed = 123)
densityplot_strom = densityplot(x=imp_strom, data = ~sRPE)

# you'll need the package ggarrange installed for this code
# to recreate figures used in the study, one can run
ggpubr::ggarrange(densityplot_u19, densityplot_strom, densityplot_handball, labels= c("Football U19", "Football Premier", "Handball"))


