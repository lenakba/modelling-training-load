
# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors
options(scipen = 17, 
        stringsAsFactors = FALSE)

# loading packages
library(tidyverse) # for datawrangling

###################################Reading fits and results from simulations############################################
base_folder = "O:\\Prosjekter\\Bache-Mathiesen-003-modelling-training-load\\Data\\simulations\\"
# amount of load
folder_j_constant = paste0(base_folder, "amount_j_constant\\")
folder_j_decay = paste0(base_folder, "amount_j_decay\\")
folder_j_exponential_decay = paste0(base_folder, "amount_j_exponential_decay\\")
folder_lin_direction_flip = paste0(base_folder, "amount_lin_direction_flip\\")

# change in load
folder_lin_constant = paste0(base_folder, "change_lin_constant\\")
folder_lin_decay = paste0(base_folder, "change_lin_decay\\")
folder_lin_exponential_decay = paste0(base_folder, "change_lin_exponential_decay\\")

# vector of file types
file_types = c("fits", "res")

# output from the for-loops and functions below are saved as "simulation_results_perfparams.csv"
# read the csv file to save time, or run all the for-loops and functions again
# perf_estimates_all = read_delim("simulation_results_perfparams.csv", delim = ";")

#------------amount, j constant
# reading the simulated results from fits
files_fits = list.files(path = folder_j_constant)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_j_constant = list()
d_res_j_constant = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_j_constant, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_j_constant, "res_",i,"_.rds"))
  l_fits_j_constant = append(l_fits_j_constant, temp_data_fits)
  d_res_j_constant = rbind(d_res_j_constant, temp_data_res)
}

###################################Source functions needed to calc performance measures###########################################

# fetching functions for performance parameters
# we assume working directory is the same location as this script
#source("performance-measure-functions.R", encoding = "UTF-8")

source("functions-relationships.R", encoding = "UTF-8")

###################################Calculate performance############################################

tl_predvalues = d_res_j_constant %>% filter(rep == 1, method == "ra") %>% pull(t_load)
true_effect = calc_coefs(tl_predvalues, lag_seq, fjdecay)
truecumcoefs = rowSums(true_effect)

cumeff_preds_ra = cp_preds_ra$allfit
cumeff_preds_ewma = cp_preds_ewma$allfit
cumeff_preds_dlnm = cp_preds_dlnm$allfit

coverage = function(estimate, target, se){
  qn = qnorm(0.95)  
  coef_high = estimate + qn*se
  coef_low = estimate - qn*se
  coverage = target <= coef_high & target >= coef_low
  numerator = sum(coverage == TRUE)
  denominator = length(coverage)
  coverage_prop = numerator/denominator
  coverage_prop
}

coverage_ra = coverage(cumeff_preds_ra, truecumcoefs, cp_preds_ra$allse)
coverage_ewma = coverage(cumeff_preds_ewma, truecumcoefs, cp_preds_ewma$allse)
coverage_dlnm = coverage(cumeff_preds_dlnm, truecumcoefs, cp_preds_dlnm$allse)

average_width = function(estimate, target, se){
  qn = qnorm(0.95)  
  coef_high = estimate + qn*se
  coef_low = estimate - qn*se
  aw = mean(coef_high-coef_low)
  aw
}

aw_ra = average_width(cumeff_preds_ra, truecumcoefs, cp_preds_ra$allse)
aw_ewma = average_width(cumeff_preds_ewma, truecumcoefs, cp_preds_ewma$allse)
aw_dlnm = average_width(cumeff_preds_dlnm, truecumcoefs, cp_preds_dlnm$allse)

# Root-mean-squared-error
rmse = function(estimate, target){
  sqrt(mean((estimate - target)^2)) 
}

rmse_ra = rmse(cumeff_preds_ra, truecumcoefs)
rmse_ewma = rmse(cumeff_preds_ewma, truecumcoefs)
rmse_dlnm = rmse(cumeff_preds_dlnm, truecumcoefs)
