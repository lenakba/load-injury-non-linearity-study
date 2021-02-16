# This script is for reading simulation_results.csv and 
# simulation_results_u_categorization_fixed.csv
# and calculate the mean and percentages of performance measures
# should be able to reproduce results found in Table 1 and Table 2 in the article
# And Figure 3

# so we don't have to deal with scientific notations
# and strings aren't automaticcaly read as factors
options(scipen = 17,stringsAsFactors = FALSE)

# Packages
library(tidyverse) # for everything
library(ggpubr) # for multiple figures in 1 panel
folder = paste0("O:\\Prosjekter\\Bache-Mathiesen-Biostatistikk\\Lena Kristin Bache-Mathiesen\\r-code\\git\\non-linearity-study\\")

#----------------------------------------------------Table 1 Performance Measure Results

# reading simulations
# error in coverage for categorized subjectively for the U shape
d_rels = read_delim(paste0(folder, "simulation_results.csv"), delim = ";")

# calc mean or percent parameters for each relationship, method, and N
group_vars = c("relationship", "n", "method")
num_vars = names(d_rels %>% select(-sig, -rep, -all_of(group_vars), -Parameter))

d_res = d_rels %>% filter(relationship != "Flat") %>% 
  distinct(!!!syms(group_vars), rep, .keep_all = TRUE) %>% 
  group_by(!!!syms(group_vars)) %>% summarise_at(vars(all_of(num_vars)), ~mean(., na.rm = TRUE)) %>% 
  select(all_of(group_vars), rmse, brier, c_stat, coverage = prop, coverage_missing = prop_missing_cov) %>% ungroup()

d_sig = d_rels  %>% filter(Parameter != "(Intercept)") %>%
  group_by(!!!syms(group_vars)) %>% arrange(!!!syms(group_vars), desc(sig)) %>% distinct(!!!syms(group_vars), rep, .keep_all = TRUE) %>% 
  summarise(n_sig = sum(sig==1), denom = n(), prop = n_sig/denom) %>% ungroup()

d_res = d_res %>% left_join(d_sig, by = group_vars)

d_res = 
  d_res %>% mutate(method = case_when(method == "Restricted Cubic Splines (Data-driven)" ~ "Restricted Cubic Splines\n(Data-Driven)",
                                      method == "Restricted Cubic Splines (Subjectively)" ~ "Restricted Cubic Splines\n(Subjectively)",
                                      method == "Categorized" ~ "Categorized (Subjectively)",
                                      method == "Categorized by Quartiles" ~ "Categorized (Quartiles)",
                                      method == "Linear Regression" ~ "Linear Model",
                                      method == "Quadratic Regression" ~ "Quadratic Model",
                                      TRUE ~ method),
                   coverage_missing = 1-coverage_missing,
                   coverage = coverage*100, 
                   coverage_missing = coverage_missing*100)

# due to an error, coverage had to be calculated for categorized (subjectively) by itself
d_u_cat = read_delim(paste0(folder, "simulation_results_u_categorization_fixed.csv"), delim = ";")
d_u_cat = d_u_cat %>% rename(rep = Rep, method = Method, n = Label)

d_res_cat = d_u_cat %>% 
  distinct(!!!syms(group_vars), rep, .keep_all = TRUE) %>% 
  group_by(!!!syms(group_vars)) %>% summarise_at(vars("prop_coverage", "prop_coverage_missing"), ~mean(., na.rm = TRUE)) %>% 
  select(all_of(group_vars), coverage = prop_coverage, coverage_missing = prop_coverage_missing) %>% ungroup()

# calc percentages, and we now have the coverages missing in the other data
d_res_cat = d_res_cat %>% mutate(coverage_missing = 1-coverage_missing,
                       coverage = coverage*100, 
                       coverage_missing = coverage_missing*100) 

#----------------------------------------------------RMSE FIGURE, Figure 3 in the article
# For creating the RMSE figure used in the article
# function for making a dot plot with a consistent style
# figure won't have the same background and colors as the one used in article
# but will show the same results
fig_dots = function(d, x, y, title, x_lab = NULL, percent = FALSE){
  text_size = 10
  
  x = enquo(x)
  y = enquo(y)
  
  p = ggplot(d, aes(x = !!x, y = !!y)) +
    geom_point(size = 3) + 
    ggtitle(title) +
    xlab(x_lab) +
    ylab(NULL)
  
  if(percent){
    p = p + scale_x_continuous(labels=axis_percent)
  }
  p = p + theme(axis.text = element_text(size=text_size),
          strip.text.x = element_text(size = text_size),
          axis.title =  element_text(size=text_size))
  p
}


d_res_u = d_res %>% filter(relationship == "U", n == "n = 8494")
d_res_u = d_res_u %>% arrange(rmse) %>% mutate(method_fct = fct_inorder(method))

d_res_j = d_res %>% filter(relationship == "J ACWR", n == "n = 6308")
d_res_j = d_res_j %>% arrange(rmse) %>% mutate(method_fct = fct_inorder(method))

plot_rmse_u = fig_dots(d_res_u, x = rmse, y = method_fct, title = "(A) U-shaped Relationship", "Root-Mean-Square Error")
plot_rmse_j = fig_dots(d_res_j, x = rmse, y = method_fct, title = "(B) J-shaped Relationship", "Root-Mean-Square Error") 

# unquote to save file
#cairo_pdf("Figure 3 Mono Image.pdf", width = 9.9, height = 3)
ggarrange(plot_rmse_u, plot_rmse_j, ncol = 2)
#dev.off()

#----------------------------------------------------- False discovery rate, Table 2 in the article
d_flat = d_rels %>% filter(relationship == "Flat")

# calc percent significant for the flat shape
# note that if the method has multiple terms or categories, it counts as significant if at least 1 has a p-value < 0.05
d_false_disc = d_flat  %>% filter(Parameter != "(Intercept)") %>%
  group_by(!!!syms(group_vars)) %>% arrange(!!!syms(group_vars), desc(sig)) %>% distinct(!!!syms(group_vars), rep, .keep_all = TRUE) %>% 
  summarise(n_sig = sum(sig==1), denom = n(), prop = n_sig/denom, perc = round(100*prop, 1)) %>% ungroup()
