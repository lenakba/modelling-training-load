library(tidyverse) # data preperation and modifitcation ++
library(dlnm) # distributed lag nonlinear models
library(PermAlgo) # simulated survival data dependent on conditions and covariates
library(survival) # for wrangling time-to-event data
library(splines) # for natural splines/cubic splines
library(slider) # for running functions on a sliding window of values, moving iteratively 1 step at a time
# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors:
options(scipen = 40, 
        stringsAsFactors = FALSE)

# assume the training load values are in the same folder
# as the r script, and that work repository is in the source file location
d_load = read_delim("norwegian_premier_league_football_td_vec.csv", delim = ";")

################################################################################
nsub = 100 # number of subjects
t_max = 300 # number of days (length of study)

# observed values
tl_observed = (d_load %>% filter(srpe <= 1200))$srpe
tl_observed_change = lead(tl_observed) - tl_observed
tl_observed_change = tl_observed_change[-length(tl_observed_change)]
tl_valid = min(tl_observed):max(tl_observed)
# lag set at 4 weeks (28) as is often used in tl studies
lag_min = 0
lag_max = 28
lag_seq = lag_min:lag_max # number of days before current day assumed to affect risk of injury

# vector of tl values used in visualizations of predictions
tl_predvalues = seq(min(tl_observed), max(tl_observed), 25)
tl_predvalues_change = seq(min(tl_observed_change), max(tl_observed_change), 25)

###################################Training load and lag structure functions###########################################
# FUNCTIONS TO COMPUTE THE LOGIT AND INVERSE LOGIT TRANSFORMATIONS
# logit = function(prob) log(prob/(1-prob))
# invlogit = function(linpred) exp(linpred)/(1+exp(linpred))

# functions for simulating the effect of the amount of training load on injury risk (J-shape)
# and the effect of amount of change of training load on injury risk (linear)
# the only exception is for the lag-function wdirection_flip, the linear function will be used instead of J-shape
fj = function(x)case_when(x < 600 ~ ((600-x)/200)^1.5/10,
                          x >= 600 ~ ((x-600)/200)^3/30)
flin = function(x) 0.0009*x

# functions to simulate the lag structure for the
# long-term time-varying effect of training load
# 4 scenarios: constant, decay, exponential decay, decay from positive to negative effect
wconst = function(lag)lag-lag+0.80
wdecay = function(lag)exp(-lag/100)
wexponential_decay = function(lag)exp(-lag/10)^2
wdirection_flip = function(lag)case_when(lag <= 6 ~ exp(-lag/10)^2,
                                         lag >  6 ~ -exp(lag/50)^2)

# create the f*w functions for amount of training load
fjconst = function(x, lag) fj(x) * wconst(lag)
fjdecay = function(x, lag) fj(x) * wdecay(lag)
fjexponential_decay = function(x, lag) fj(x) * wexponential_decay(lag)
flindirection_flip = function(x, lag) flin(x) * wdirection_flip(lag)

# f*w functions for change in load
flinconst = function(x, lag) flin(x) * wconst(lag)
flindecay = function(x, lag) flin(x) * wdecay(lag)
flinexponential_decay = function(x, lag) flin(x) * wexponential_decay(lag)

# calculate coefficients for each combination of f*w
calc_coefs = function(x, lag, FUN){
  coefs = outer(x, lag, FUN)
  dimnames(coefs) = list(x, paste(lag))
  coefs
}

#################################################Simulating exposure histories###################################
# simulate exposure histories
# for nsub number of subjects
# and t_max days under study
sim_tl_history = function(nsub, t_max, tl_values){
  l_tl_hist = vector(mode = "list", length = nsub)
  l_tl_hist = l_tl_hist %>% map(~sample(tl_values, t_max))
  # add IDs and calc the change in training load from previous day
  d_sim_tl_hist = l_tl_hist %>% 
    map(. %>% enframe(name = NULL)) %>% 
    bind_rows() %>% 
    rename(t_load = value) %>% 
    mutate(id = rep(1:nsub, each = t_max),
           day = rep(1:t_max, nsub),
           t_load_change = lead(t_load)-t_load) 
  d_sim_tl_hist %>% select(id, day, t_load, t_load_change)
}

#############################Calculate cumulative effect of training load given exposure histories#######################

# function to compute the cumulative effect of training load exposure, requires:
# tl_hist     the tl exposure history from t0 to max(t) for one participant
# fvar        the function for the effect of training load
# flag        the function for the time-lagged effect of training load
calc_cumeff = function(tl_hist, fvar, flag){
  l_effect_seqs = slide(tl_hist, ~pluck(.), .before = lag_max, step = 1, .complete = FALSE)
  l_effect_lags = l_effect_seqs %>% map(~(length(.)-1):0)
  l_flin = l_effect_seqs %>% map(~fvar(.))
  l_wdecay = l_effect_lags %>% map(~flag(.))
  combfun_effect = l_flin %>% map2(.x = ., .y = l_wdecay, ~ .x * .y)
  cumeffect = combfun_effect %>% map(., sum)
  cumeffect_mat = unlist(cumeffect)
  cumeffect_mat
}

# function for calculating the cumulative effect for each individual in a dataset
# with the choice of var-function (fvar) and lag-function (flag)
# and specification of whether it is the amount of training load "amount"
# or change in training load "change" the calculation is performed on.
calc_cumeffs_all = function(d_sim_hist, t_load_type = "amount", fvar, flag){
  l_sim_tl_hist = d_sim_hist %>% group_by(id) %>% nest() 
  if(t_load_type == "amount"){
    l_sim_tl_hist$data = l_sim_tl_hist$data %>% map(~calc_cumeff(.$t_load, fvar, flag))  
  } else if(t_load_type == "change"){
    l_sim_tl_hist$data = l_sim_tl_hist$data %>% map(~calc_cumeff(.$t_load_change, fvar, flag))
  }
  d_cumeffs = unnest(l_sim_tl_hist, cols = c(data)) %>% rename(cumeff = data) %>% ungroup()
  d_cumeffs
}

########################################Simulate injuries based on cumulative effect of training load#####################

# helper function which simulates survival data from training load values
# step 1 calculate cumulative effect of training load for each participant 
#        based on a function for trainig load fvar, and a lag-function flag.
#        Choice of t_load_type = "amount" or "change" depending on relationship.
# step 2 simulate injuries based on cumulative effects
sim_survdata = function(d_hist, nsub, t_max, t_load_type, fvar, flag){
  d_cumeff = calc_cumeffs_all(d_hist, t_load_type = t_load_type, fvar, flag)
  d_cumeff_mat = d_cumeff %>% mutate(day = rep(1:t_max, nsub)) %>% 
    pivot_wider(names_from = id, values_from = cumeff) %>% 
    select(-day) %>% as.matrix
  d_survival_sim = permalgorithm(nsub, t_max, Xmat = d_cumeff_mat, censorRandom = runif(nsub, 1, t_max*2), betas=1)
  d_survival_sim
} 

########################################Define the crossbasis for the DLNM method#####################

# function to arrange the simulated survival data in counting process form
counting_process_form = function(d_survival_sim){
  d_follow_up_times = d_survival_sim %>% distinct(Id, Fup)
  d_surv_lim = map2(.x = d_follow_up_times$Id,
                    .y = d_follow_up_times$Fup,
                    ~d_survival_sim %>% filter(Id == .x) %>% slice(.y)) %>% 
    bind_rows() %>% select(event = Event, exit = Stop, id = Id)
  
  # extracting timepoints in which an event happened
  ftime = d_surv_lim %>% filter(event == 1) %>% distinct(exit) %>% arrange(exit) %>% pull()
  
  # arrange the survival data so that, for each individual, we have an interval of enter and exit times
  # for each of the exit times above, with the information of whether or not they were injured at that time
  # meaning we will have the same time intervals per participant
  d_counting_process = survSplit(Surv(exit, event)~., d_surv_lim, cut = ftime, start="enter") %>% arrange(id)
  d_counting_process
}

# function for calculating the q matrix given the survival data in counting process form
# and the exposure history spread in wide format in a matrix
calc_q_matrix = function(d_counting_process, d_tl_hist_wide){
  # for each individual, for each of these exit times, we will extract the exposure history 
  # for the given lag-time which we are interested in
  # This is called the Q-matrix. The Q-matrix should be nrow(dataspl) X 0:lag_max dimensions.
  q = map2(.x = d_counting_process$id, 
           .y = d_counting_process$exit, 
           ~exphist(d_tl_hist_wide[.x,], .y, c(lag_min, lag_max))) %>% 
    do.call("rbind", .)
  q
  
}

# helper function to use both functions above in 1 step
# resulting in a q-matrix
from_sim_surv_to_q = function(d_survival_sim, d_tl_hist_wide){
  d_counting_process = counting_process_form(d_survival_sim)
  q = calc_q_matrix(d_counting_process, d_tl_hist_wide)
  q
}

####################################### Modify training load with different methods ####################################

# Perform the typical methods for handling training load amount: rolling average and EWMA

# Function for calculating rolling averages on a chooseable number of days
# Based on rollapplyr, not rollmean, as rollmean will only start calculating averages
# at n values, while rollapplyr allows the user to decide preliminary values.
ra = function(x, n_days, window = TRUE, ...){
  zoo::rollapplyr(x, n_days, mean, partial = window)
}

# function for calculating exponentially waited moving averages
# using similar syntax as the RA-function
ewma = function(x, n_days, window){
  TTR::EMA(x, n = window, ratio = 2/(1+n_days))
}

# Perform typical methods for handling change in training load, ACWR and week-to-week change
# calculate 7:28 coupled ACWR using RA on training load amount (this becomes, in theory, a measure of change)
# move 1 day at a time as advised in Carey et al. 2017
# function calculates RA on a sliding window of 21 days
chronic_ra = function(x){
  l = slide(x, ~ra(., 28), .before = 27, step = 1, .complete =TRUE) %>% map(last)
  l = compact(l)
  l = unlist(l)
  l
}

# function calculates mean on a sliding window of 7 days
acute_mean = function(x){
  l = slide(x, ~mean(.), .before = 6, step = 1, .complete =TRUE)
  l = compact(l)
  l = unlist(l)
  l
}

# function calculating sums on a sliding window of 7 days
weekly_sum = function(x){
  l = slide(x, ~sum(.), .before = 6, step = 1, .complete =TRUE)
  l = compact(l)
  l = unlist(l)
  l
}

# function to nest the exposure history data by each individual, 
# and run a user-specified function on each of their datasets in the list
function_on_list = function(d_sim_hist, FUN = NULL, day_start){
  nested_list = d_sim_hist %>% group_by(id) %>% nest()
  nested_list$data = nested_list$data %>% map(., ~FUN(.$t_load))
  l_unnest = unnest(nested_list, cols = c(data)) %>% ungroup() %>% mutate(day = rep(day_start:t_max, nsub)) 
  l_unnest
}

################################################## Calculate numeric performance measures ############################

aic_ra = AIC(fit_ra)
aic_ewma = AIC(fit_ewma)
aic_dlnm = AIC(fit_dlnm)  

tl_predvalues_ra = seq(min(d_survival_sim_cpform_mods$ra_t_load), max(d_survival_sim_cpform_mods$ra_t_load), 25) 
tl_predvalues_ewma = seq(min(d_survival_sim_cpform_mods$ewma_t_load), max(d_survival_sim_cpform_mods$ewma_t_load), 25)

cp_preds_ra = crosspred(ob_ra, fit_ra, at = tl_predvalues, cen = 600, cumul = TRUE)
cp_preds_ewma = crosspred(ob_ewma, fit_ewma, at = tl_predvalues, cen = 600, cumul = TRUE)
cp_preds_dlnm = crosspred(cb_dlnm, fit_dlnm, at = tl_predvalues, cen = 600, cumul = TRUE)

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

######################################################## Helper function which does all of the above

sim_fit_and_res = function(nsub, t_max, tl_values, t_load_type, tl_var, fvar, flag, predvalues, i, folder){
  tl_var = enquo(tl_var)

  # simulate exposure history
  d_tl_hist = sim_tl_history(nsub, t_max, tl_values)
  # simulate injuries based on a function on the tl variable (fvar) and the timelag (flag)
  d_survival_sim = sim_survdata(d_tl_hist, nsub, t_max, t_load_type = t_load_type, fvar = fvar, flag = flag) 
  
  # arrange the exposure history in wide format in a matrix
  # which is neeeded for calculating the q-matrix for the crossbasis
  d_sim_tl_hist_spread_day = 
    d_tl_hist %>% select(id, day, !!tl_var) %>% 
    pivot_wider(names_from = day, values_from = !!tl_var) %>% select(-id) %>% as.matrix
  
  d_survival_sim_cpform = counting_process_form(d_survival_sim)
  q_mat = from_sim_surv_to_q(d_survival_sim, d_sim_tl_hist_spread_day)

  if(t_load_type == "amount"){
  # calc rolling average and ewma on training load amount
    d_tl_hist = d_tl_hist %>% 
    mutate(ra_t_load = ra(d_tl_hist$t_load, 28),
           ewma_t_load = ewma(d_tl_hist$t_load, 28, 1))
  d_survival_sim_cpform_mods = d_survival_sim_cpform %>% left_join(d_tl_hist, by = c("id", "exit" = "day"))
  ob_ra = onebasis(d_survival_sim_cpform_mods$ra_t_load, "poly", degree = 2)
  ob_ewma = onebasis(d_survival_sim_cpform_mods$ewma_t_load, "poly", degree = 2)
  cb_dlnm = crossbasis(q_mat, lag=c(lag_min, lag_max), 
                       argvar = list(fun="poly", degree = 2),
                       arglag = list(fun="ns", knots = 3))

  fit_ra = coxph(Surv(enter, exit, event) ~ ob_ra, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
  fit_ewma = coxph(Surv(enter, exit, event) ~ ob_ewma, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
  fit_dlnm = coxph(Surv(enter, exit, event) ~ cb_dlnm, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
  list_fits = list(fit_ra, fit_ewma, fit_dlnm)

  cp_preds_ra = crosspred(ob_ra, fit_ra, at = predvalues, cen = 600, cumul = TRUE)
  cp_preds_ewma = crosspred(ob_ewma, fit_ewma, at = predvalues, cen = 600, cumul = TRUE)
  cp_preds_dlnm = crosspred(cb_dlnm, fit_dlnm, at = predvalues, cen = 600, cumul = TRUE)
  
  d_ra = bind_cols(t_load = predvalues, 
                   cumul = cp_preds_ra$allfit, 
                   se = cp_preds_ra$allse,
                   aic = AIC(fit_ra))
  
  d_ewma = bind_cols(t_load = predvalues, 
                   cumul = cp_preds_ewma$allfit, 
                   se = cp_preds_ewma$allse,
                   aic = AIC(fit_ewma))
  
  d_dlnm = bind_cols(t_load = predvalues, 
                   cumul = cp_preds_dlnm$allfit, 
                   se = cp_preds_dlnm$allse,
                   aic = AIC(fit_dlnm))
  list_res = list(d_ra, d_ewma, d_dlnm)
  
  } else if(t_load_type == "change"){
    
    # the first acute load can be calulated from day 22 to day 28
    # to have an equal number acute and chronic values
    d_sim_tl_hist_acwr_acute = d_tl_hist %>% filter(day >= 22)
    d_sim_hist_acute = function_on_list(d_sim_tl_hist_acwr_acute, FUN = acute_mean, 28) %>% rename(acute_load = data)
    d_sim_hist_chronic = function_on_list(d_tl_hist, FUN = chronic_ra, 28) %>% rename(chronic_load = data)
    
    # calculate the ACWR be acute/chronic
    d_sim_hist_acwr = d_sim_hist_acute %>% 
      left_join(d_sim_hist_chronic, by = c("id", "day")) %>% 
      mutate(acwr = acute_load/chronic_load)
    
    # week to week change requires 2 full weeks before calculation
    # the sliding window thereafter jumps one day at a time
    # but the calculation is "uncoupled"
    d_sim_hist_weekly = function_on_list(d_tl_hist, FUN = weekly_sum, 7) %>% 
      rename(week_sum = data) %>% 
      mutate(week_sum_lead = lead(week_sum, 7),
             weekly_change = week_sum_lead-week_sum)
    
    # the difference can't be measured until a week after the first week
    d_sim_hist_weekly = d_sim_hist_weekly %>% 
      select(-week_sum, -week_sum_lead) %>% 
      group_by(id) %>% 
      mutate(day = lead(day, 7)) %>% 
      filter(!is.na(day)) %>% ungroup()
    
    # couple survival data with change in load data
    d_survival_sim_cpform_mods = d_survival_sim_cpform %>% 
      left_join(d_sim_hist_acwr, by = c("id", "exit" = "day"))
    
    d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% 
      left_join(d_sim_hist_weekly, by = c("id", "exit" = "day"))
    
    ob_acwr = onebasis(d_survival_sim_cpform_mods$acwr, "lin")
    ob_weekly_change = onebasis(d_survival_sim_cpform_mods$weekly_change, "lin")
    cb_dlnm_change = crossbasis(q_mat, lag=c(lag_min, lag_max), 
                                argvar = list(fun="lin"),
                                arglag = list(fun="ns", knots = 3))
    
    fit_acwr = coxph(Surv(enter, exit, event) ~ ob_acwr, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
    fit_weekly_change = coxph(Surv(enter, exit, event) ~ ob_weekly_change, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
    fit_dlnm_change = coxph(Surv(enter, exit, event) ~ cb_dlnm_change, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
    list_fits = list(fit_acwr, fit_weekly_change, fit_dlnm_change)
    
    predvalues_acwr = seq(min(d_sim_hist_acwr$acwr, na.rm = TRUE), max(d_sim_hist_acwr$acwr, na.rm = TRUE), 25)
    predvalues_wchange = seq(min(d_sim_hist_weekly$weekly_change, na.rm = TRUE), max(d_sim_hist_weekly$weekly_change, na.rm = TRUE), 25)
    
    cp_preds_acwr = crosspred(ob_acwr, fit_acwr, at = predvalues_acwr, cen = 1, cumul = TRUE)
    cp_preds_weekly_change = crosspred(ob_weekly_change, fit_weekly_change, at = predvalues_wchange, cen = 0, cumul = TRUE)
    cp_preds_dlnm = crosspred(cb_dlnm_change, fit_dlnm_change, at = predvalues, cen = 0, cumul = TRUE)
    
    d_acwr = bind_cols(t_load_acwr = predvalues_acwr, 
                     cumul = cp_preds_acwr$allfit, 
                     se = cp_preds_acwr$allse,
                     aic = AIC(fit_acwr))
    
    d_ewma = bind_cols(t_load_weekly_change = predvalues_wchange, 
                       cumul = cp_preds_weekly_change$allfit, 
                       se = cp_preds_weekly_change$allse,
                       aic = AIC(fit_weekly_change))
    
    d_dlnm = bind_cols(t_load_change = predvalues, 
                       cumul = cp_preds_dlnm$allfit, 
                       se = cp_preds_dlnm$allse,
                       aic = AIC(fit_dlnm_change))
    list_res = list(d_acwr, d_ewma, d_dlnm)
  }
  saveRDS(list_fits, file = paste0(folder, "fits_",i,"_.rds"))
  saveRDS(list_res, file = paste0(folder, "res_",i,"_.rds"))
}

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


startsim = 1
nsim = 3
seqsim = startsim:nsim
set.seed(1234)
for(i in seqsim){
  # amount of training load
  sim_fit_and_res(nsub, t_max, tl_valid, "amount", tl_var = t_load, fvar = fj, flag = wconst, 
                  predvalues = tl_predvalues, i = i, folder = folder_j_constant)
  sim_fit_and_res(nsub, t_max, tl_valid, "amount", tl_var = t_load, fvar = fj, flag = wdecay, 
                  predvalues = tl_predvalues, i = i, folder = folder_j_decay)
  sim_fit_and_res(nsub, t_max, tl_valid, "amount", tl_var = t_load, fvar = fj, flag = wexponential_decay, 
                  predvalues = tl_predvalues, i = i, folder = folder_j_exponential_decay)
  
  # change in training load
  sim_fit_and_res(nsub, t_max, tl_valid, "change", tl_var = t_load_change, fvar = flin, flag = wconst, 
                  predvalues = tl_predvalues_change, i = i, folder = folder_lin_constant)
  sim_fit_and_res(nsub, t_max, tl_valid, "change", tl_var = t_load_change, fvar = flin, flag = wdecay, 
                  predvalues = tl_predvalues_change, i = i, folder = folder_lin_decay)
  sim_fit_and_res(nsub, t_max, tl_valid, "change", tl_var = t_load_change, fvar = flin, flag = wexponential_decay, 
                  predvalues = tl_predvalues_change, i = i, folder = folder_lin_exponential_decay)
}

#-----------------------Option: multiple cores
library(foreach)
library(doParallel)
set.seed(1234)
numCores = 4
registerDoParallel(numCores)

startsim = 1
nsim = 4
seqsim = startsim:nsim

options(warn=-1)
foreach (i = seqsim) %dopar% {
  library(tidyverse)
  library(dlnm)
  library(PermAlgo)
  library(survival)
  library(splines)
  library(slider) 
  sim_fit_and_res(nsub, t_max, tl_valid, "amount", tl_var = t_load, fvar = fj, flag = wdecay, 
                  predvalues = tl_predvalues, i = i, folder = folder_j_decay)
}
options(warn=0)



