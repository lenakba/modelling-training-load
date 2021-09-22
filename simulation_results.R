
# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors
options(scipen = 17, 
        stringsAsFactors = FALSE)

# loading packages
library(tidyverse) # for datawrangling

################################### Reading fits and results from simulations############################################
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

################################### Amount
#------------amount, j constant
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
d_res_j_constant = d_res_j_constant %>% mutate(relationship = "Constant")

#------------amount, j decay
files_fits = list.files(path = folder_j_decay)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_j_decay = list()
d_res_j_decay = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_j_decay, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_j_decay, "res_",i,"_.rds"))
  l_fits_j_decay = append(l_fits_j_decay, temp_data_fits)
  d_res_j_decay = rbind(d_res_j_decay, temp_data_res)
}
d_res_j_decay = d_res_j_decay %>% mutate(relationship = "Decay")

#------------amount, j exponential decay
files_fits = list.files(path = folder_j_exponential_decay)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_j_exponential_decay = list()
d_res_j_exponential_decay = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_j_exponential_decay, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_j_exponential_decay, "res_",i,"_.rds"))
  l_fits_j_exponential_decay = append(l_fits_j_exponential_decay, temp_data_fits)
  d_res_j_exponential_decay = rbind(d_res_j_exponential_decay, temp_data_res)
}
d_res_j_exponential_decay = d_res_j_exponential_decay %>% mutate(relationship = "Exponential Decay")

# make list of datasets for amount results
l_res_amount = list(d_res_j_constant, d_res_j_decay, d_res_j_exponential_decay)

################################### Change

#------------change, lin constant
files_fits = list.files(path = folder_lin_constant)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_lin_constant = list()
d_res_lin_constant = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_lin_constant, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_lin_constant, "res_",i,"_.rds"))
  l_fits_lin_constant = append(l_fits_lin_constant, temp_data_fits)
  d_res_lin_constant = rbind(d_res_lin_constant, temp_data_res)
}
d_res_lin_constant = d_res_lin_constant %>% mutate(relationship = "Constant")

#------------amount, lin decay
files_fits = list.files(path = folder_lin_decay)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_lin_decay = list()
d_res_lin_decay = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_lin_decay, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_lin_decay, "res_",i,"_.rds"))
  l_fits_lin_decay = append(l_fits_lin_decay, temp_data_fits)
  d_res_lin_decay = rbind(d_res_lin_decay, temp_data_res)
}
d_res_lin_decay = d_res_lin_decay %>% mutate(relationship = "Decay")

#------------amount, lin exponential decay
files_fits = list.files(path = folder_lin_exponential_decay)
n_sim = length(files_fits)/length(file_types) # divide by the number of file types
l_fits_lin_exponential_decay = list()
d_res_lin_exponential_decay = data.frame()
for(i in 1:n_sim){
  temp_data_fits = readRDS(paste0(folder_lin_exponential_decay, "fits_",i,"_.rds"))
  temp_data_res = readRDS(paste0(folder_lin_exponential_decay, "res_",i,"_.rds"))
  l_fits_lin_exponential_decay = append(l_fits_lin_exponential_decay, temp_data_fits)
  d_res_lin_exponential_decay = rbind(d_res_lin_exponential_decay, temp_data_res)
}
d_res_lin_exponential_decay = d_res_lin_exponential_decay %>% mutate(relationship = "Exponential Decay")

# make list of datasets for change results
l_res_change = list(d_res_lin_constant, d_res_lin_decay, d_res_lin_exponential_decay)


###################################Source functions needed to calc performance measures###########################################

# fetching functions for performance parameters
# and for estimating the true relationship
# we assume working directory is the same location as this script
source("functions-performance.R", encoding = "UTF-8")
source("functions-relationships.R", encoding = "UTF-8")

# list of relationship functions
fw_funs = list(fjconst, fjdecay, fjexponential_decay)
fw_funs_change = list(flinconst, flindecay, flinexponential_decay)

###################################Calculate performance############################################
# lag set at 4 weeks (28) as is often used in tl studies
lag_min = 0
lag_max = 28
lag_seq = lag_min:lag_max # number of days before current day assumed to affect risk of injury
tl_predvalues = d_res_j_constant %>% filter(rep == 1, method == "ra") %>% pull(t_load)
tl_predvalues_change = d_res_lin_constant %>% filter(rep == 1, method == "dlnm") %>% pull(t_load_change)

# function for adding the true cumulative effect given the fw-function
add_true_coefs = function(d_res, tl_predvalues, lag_seq, fw){
  true_effect = calc_coefs(tl_predvalues, lag_seq, fw)
  truecumcoefs = rowSums(true_effect_j_constant)
  d_res = d_res %>% mutate(true_cumul_coefs = truecumcoefs)
  d_res
}

#---------------Amount

# group by method and simulation number (rep)
# then use map2 to match the right dataset to the right fw function
# and add the (simulated) true cumulative effect
# collapse to dataset
d_res_amount = l_res_amount %>% 
  map(. %>% group_by(method, rep)) %>% 
  map2(.x = .,
       .y = fw_funs,
       ~add_true_coefs(.x, tl_predvalues, lag_seq, .y)) %>% bind_rows()

nsims = max(d_res_amount$rep)
d_res_amount = d_res_amount %>% group_by(relationship, method, rep) %>% 
                  mutate(ci_high = calc_ci(cumul, se, "high"), 
                         ci_low = calc_ci(cumul, se, "low"),
                         coverage = coverage(ci_high, ci_low, true_cumul_coefs),
                         aw = average_width(ci_high, ci_low),
                         rmse = rmse(cumul, true_cumul_coefs),
                         mcse_rmse = mcse_rmse(cumul, true_cumul_coefs, nsims),
                         mcse_coverage = mcse_coverage(ci_low, ci_high, true_cumul_coefs, n(), nsims)) %>% 
                  ungroup()

perf_cols = c("aic", "coverage", "aw", "rmse")
d_perf_params_amount = d_res_amount %>% group_by(relationship, method) %>% summarise_at(vars(perf_cols, starts_with("mcse")), mean)

#---------------Change
d_res_change = l_res_change %>% bind_rows()
d_res_change = d_res_change %>% 
  group_by(relationship, method, rep) %>% 
  mutate(ci_high = calc_ci(cumul, se, "high"), 
         ci_low = calc_ci(cumul, se, "low"),
         aw = average_width(ci_high, ci_low)) %>% 
  ungroup()

perf_cols = c("aic", "aw")
d_perf_params_change = d_res_lin_constant %>% group_by(relationship, method) %>% summarise_at(vars(perf_cols), mean)

