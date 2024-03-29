
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

# output from the for-loops and functions below are saved as 
# "simulation_results_amount.csv" and "simulation_results_change.csv"
# read the csv file to save time, or run all the for-loops and functions again
d_res_amount = read_delim("simulation_results_amount.csv", delim = ";")
d_res_change = read_delim("simulation_results_change.csv", delim = ";", col_types = cols(
  t_load_acwr = col_double(),
  cumul = col_double(),
  se = col_double(),
  aic = col_double(),
  rmse_residuals = col_double(),
  method = col_character(),
  t_load_weekly_change = col_double(),
  t_load_change = col_double(),
  rep = col_double(),
  relationship = col_character(),
  true_cumul_coefs = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  aw = col_double()
))

################################### Amount
#------------amount, j constant
files_fits = list.files(path = folder_j_constant)
n_sim = 1900
d_res_j_constant = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_j_constant, "res_",i,"_.rds"))
  d_res_j_constant = rbind(d_res_j_constant, temp_data_res)
}
d_res_j_constant = d_res_j_constant %>% mutate(relationship = "Constant")

#------------amount, j decay
files_fits = list.files(path = folder_j_decay)
d_res_j_decay = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_j_decay, "res_",i,"_.rds"))
  d_res_j_decay = rbind(d_res_j_decay, temp_data_res)
}
d_res_j_decay = d_res_j_decay %>% mutate(relationship = "Decay")

#------------amount, j exponential decay
files_fits = list.files(path = folder_j_exponential_decay)
d_res_j_exponential_decay = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_j_exponential_decay, "res_",i,"_.rds"))
  d_res_j_exponential_decay = rbind(d_res_j_exponential_decay, temp_data_res)
}
d_res_j_exponential_decay = d_res_j_exponential_decay %>% mutate(relationship = "Exponential Decay")

#------------amount, lin direction_flip
files_fits = list.files(path = folder_lin_direction_flip)
d_res_lin_direction_flip = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_lin_direction_flip, "res_",i,"_.rds"))
  d_res_lin_direction_flip = rbind(d_res_lin_direction_flip, temp_data_res)
}
d_res_lin_direction_flip = d_res_lin_direction_flip %>% mutate(relationship = "Linear Direction Change")

# make list of datasets for amount results
l_res_amount = list(d_res_j_constant, d_res_j_decay, d_res_j_exponential_decay, d_res_lin_direction_flip)

################################### Change

#------------change, lin constant
files_fits = list.files(path = folder_lin_constant)
d_res_lin_constant = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_lin_constant, "res_",i,"_.rds"))
  d_res_lin_constant = rbind(d_res_lin_constant, temp_data_res)
}
d_res_lin_constant = d_res_lin_constant %>% mutate(relationship = "Constant")

#------------amount, lin decay
files_fits = list.files(path = folder_lin_decay)
d_res_lin_decay = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_lin_decay, "res_",i,"_.rds"))
  d_res_lin_decay = rbind(d_res_lin_decay, temp_data_res)
}
d_res_lin_decay = d_res_lin_decay %>% mutate(relationship = "Decay")

#------------amount, lin exponential decay
files_fits = list.files(path = folder_lin_exponential_decay)
d_res_lin_exponential_decay = data.frame()
for(i in 1:n_sim){
  temp_data_res = readRDS(paste0(folder_lin_exponential_decay, "res_",i,"_.rds"))
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
fw_funs = list(fjconst, fjdecay, fjexponential_decay, flindirection_flip)
fw_funs_change = list(flinconst, flindecay, flinexponential_decay)

###################################Calculate performance############################################
# lag set at 4 weeks (28) as is often used in tl studies
lag_min = 0
lag_max = 27
lag_seq = lag_min:lag_max # number of days before current day assumed to affect risk of injury
tl_predvalues = d_res_j_constant %>% filter(rep == 1, method == "ra") %>% pull(t_load)
tl_predvalues_change = d_res_lin_constant %>% filter(rep == 1, method == "dlnm %Δ") %>% pull(t_load_change)

# function for adding the true cumulative effect given the fw-function
add_true_coefs = function(d_res, tl_predvalues, lag_seq, fw){
  true_effect = calc_coefs(tl_predvalues, lag_seq, fw)
  truecumcoefs = rowSums(true_effect)
  d_res = d_res %>% mutate(true_cumul_coefs = truecumcoefs)
  d_res
}

# variables in both amount and change
perf_internal = c("aic", "rmse_residuals", "aw")

#---------------Amount

# group by method and simulation number (rep)
# then use map2 to match the right dataset to the right fw function
# and add the (simulated) true cumulative effect
# collapse to dataset
d_res_amount = l_res_amount %>% 
  map(. %>% group_by(method, rep)) %>% 
  map2(.x = .,
       .y = fw_funs,
       ~add_true_coefs(.x, tl_predvalues, lag_seq, .y)) %>% 
  bind_rows() %>% 
  ungroup()

nsims = max(d_res_amount$rep)
d_res_amount = d_res_amount %>% group_by(relationship, method, rep) %>% 
                  mutate(ci_high = calc_ci(cumul, se, "high"), 
                         ci_low = calc_ci(cumul, se, "low"),
                         coverage = coverage(ci_high, ci_low, true_cumul_coefs),
                         aw = average_width(ci_high, ci_low),
                         bias = raw_bias(cumul, true_cumul_coefs), 
                         rmse = rmse(cumul, true_cumul_coefs),
                         mcse_bias = mcse_bias(cumul, true_cumul_coefs, nsims),
                         mcse_rmse = mcse_rmse(cumul, true_cumul_coefs, nsims),
                         mcse_coverage = mcse_coverage(ci_low, ci_high, true_cumul_coefs, n(), nsims)) %>% 
                  ungroup()

# variables in amount only
perf_external = c("rmse", "coverage")
d_perf_params_amount = d_res_amount %>% group_by(relationship, method) %>% summarise_at(vars(all_of(perf_internal), all_of(perf_external), starts_with("mcse")), mean)

#---------------Change

# add true coefs (can only compare with DLNM)
d_res_change_dlnm = l_res_change %>% 
  map(. %>% filter(method == "dlnm %Δ") %>% group_by(rep)) %>% 
  map2(.x = .,
       .y = fw_funs_change,
       ~add_true_coefs(.x, tl_predvalues_change, lag_seq, .y)) %>% 
  bind_rows() %>% 
  ungroup()

d_res_change_other = l_res_change %>%  map(. %>% filter(method != "dlnm %Δ")) %>% bind_rows()
d_res_change = bind_rows(d_res_change_other, d_res_change_dlnm)

d_res_change = d_res_change %>% 
  group_by(relationship, method, rep) %>% 
  mutate(ci_high = calc_ci(cumul, se, "high"), 
         ci_low = calc_ci(cumul, se, "low"),
         aw = average_width(ci_high, ci_low)) %>% 
  ungroup()

d_perf_params_change = d_res_change %>% group_by(relationship, method) %>% summarise_at(vars(all_of(perf_internal)), mean)

#--------------- find the sample size after 500 sims
sample_size_needed = d_res_amount  %>% 
  group_by(relationship, method) %>% 
  summarise(variance_est = var(bias, na.rm = TRUE), n_sim = (variance_est^2)/0.25) %>% 
  ungroup() %>% summarise(n_sim_needed = max(n_sim)) %>% 
  pull(n_sim_needed)

############################################ Saving the results as datasets #################################################

# write_delim is preferable, but write_excel_csv is required for excel to understand
# that the file encoding is UTF-8
# write_excel_csv(d_res_amount, "simulation_results_amount.csv", delim = ";", na = "")
# write_excel_csv(d_res_change, "simulation_results_change.csv", delim = ";", na = "")
