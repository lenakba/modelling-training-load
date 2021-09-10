
library(tidyverse) # data preperation and modifitcation ++
library(dlnm) # distributed lag nonlinear models
library(PermAlgo) # simulated survival data dependent on conditions and covariates
library(survival) # for wrangling time-to-event data
library(splines) # for natural splines/cubic splines
library(lmisc) # NIH colors for figures
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


# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim = cbind(x_funs=rep(c("flin", "fj"), each=2),
                lag_funs=rep(c("wconst", "wdecay"), 2)) %>% 
  as_tibble()

# create the f*w functions for amount of training load
fjconst = function(x, lag) fj(x) * wconst(lag)
fjdecay = function(x, lag) fj(x) * wdecay(lag)
fjexponential_decay = function(x, lag) fj(x) * wexponential_decay(lag)
flindirection_flip = function(x, lag) flin(x) * wdirection_flip(lag)
fw_funs = list(fjconst, fjdecay, fjexponential_decay, flindirection_flip)

# f*w functions for change in load
flinconst = function(x, lag) flin(x) * wconst(lag)
flindecay = function(x, lag) flin(x) * wdecay(lag)
flinexponential_decay = function(x, lag) flin(x) * wexponential_decay(lag)
fw_funs_change = list(flinconst, flindecay, flinexponential_decay)

# calculate coefficients for each combination of f*w
calc_coefs = function(x, lag, FUN){
  coefs = outer(x, lag, FUN)
  dimnames(coefs) = list(x, paste(lag))
  coefs
}
l_coefs = fw_funs %>% map(., ~calc_coefs(tl_predvalues, lag_seq, .)) %>% map(., ~exp(.))
l_coefs_change = fw_funs_change %>% map(., ~calc_coefs(tl_predvalues_change, lag_seq, .)) %>% map(., ~exp(.))

# figures of f*w functions for amount of training load
persp(x = tl_predvalues, y = lag_seq, l_coefs[[1]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE (AU)")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE(AU)")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[3]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE (AU)")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[4]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

# figures for change in load
persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[1]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "ΔsRPE (AU)")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "ΔsRPE (AU)")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[3]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "ΔsRPE (AU)")


dag0 = l_coefs_change[[2]][,1]
plot(x = tl_predvalues_change, y = dag0)

x500 = l_coefs_change[[2]][51,1:29]
plot(x = lag_seq, y = x500)



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

set.seed(1234)
sim_tl_history(nsub, t_max, tl_valid)

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

d_cumeff = calc_cumeffs_all(d_sim_tl_hist, t_load_type = "amount", fj, wdecay)

########################################Simulate injuries based on cumulative effect of training load#####################
d_cumeff_mat = d_cumeff %>% mutate(day = rep(1:t_max, nsub)) %>% 
                            pivot_wider(names_from = id, values_from = cumeff) %>% select(-day) %>% as.matrix

set.seed(1234)
d_survival_sim = permalgorithm(nsub, t_max, Xmat = d_cumeff_mat, censorRandom = runif(nsub, 1, t_max*2), betas=1)

# helper function which simulates survival data from training load values
# step 1 simulate study with nsub participants with an training load 
#        exposure history lasting t_max days sampled from tl_values
# step 2 calculate cumulative effect of training load for each participant 
#        based on a function for trainig load fvar, and a lag-function flag.
#        Choice of t_load_type = "amount" or "change" depending on relationship.
# step 3 simulate injuries based on cumulative effects
sim_survdata = function(nsub, t_max, tl_values, t_load_type, fvar, flag){
  d_tl_hist = sim_tl_history(nsub, t_max, tl_valid)
  d_cumeff = calc_cumeffs_all(d_tl_hist, t_load_type = t_load_type, fvar, flag)
  d_cumeff_mat = d_cumeff %>% mutate(day = rep(1:t_max, nsub)) %>% 
                              pivot_wider(names_from = id, values_from = cumeff) %>% 
                              select(-day) %>% as.matrix
  d_survival_sim = permalgorithm(nsub, t_max, Xmat = d_cumeff_mat, censorRandom = runif(nsub, 1, t_max*2), betas=1)
  d_survival_sim
} 
set.seed(1234)
d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "amount", fvar = fj, flag = wdecay)
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

# arrange the exposure history in wide format in a matrix
# which is neeeded for calculating the q-matrix for the crossbasis
d_sim_tl_hist_spread_day = 
  d_sim_tl_hist %>% select(-t_load_change) %>% 
  pivot_wider(names_from = day, values_from = t_load) %>% select(-id) %>% as.matrix

d_sim_tl_hist_spread_day_change = 
  d_sim_tl_hist %>% select(-t_load) %>%
  pivot_wider(names_from = day, values_from = t_load_change) %>% select(-id) %>% as.matrix

# list of q-matrices
l_q_matrices = l_survival_sim %>% map(., ~from_sim_surv_to_q(., d_sim_tl_hist_spread_day))
l_q_matrices_change = l_survival_sim_change %>% map(., ~from_sim_surv_to_q(., d_sim_tl_hist_spread_day_change))

# need data on counting process form
l_counting_survival_sim = l_survival_sim %>% map(., ~counting_process_form(.))
l_counting_survival_sim_change = l_survival_sim_change %>% map(., ~counting_process_form(.))


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

# calc rolling average and ewma on training load amount
d_sim_tl_hist_mod = d_sim_tl_hist %>% 
  mutate(ra_t_load = ra(d_sim_tl_hist$t_load, 28),
         ewma_t_load = ewma(d_sim_tl_hist$t_load, 28, 1)) %>% select(-starts_with("t_load"))

l_survival_sim_basemethods = l_survival_sim %>% map(. %>% left_join(d_sim_tl_hist_mod, by = c("Id" = "id", "Stop" = "day")))

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

# the first acute load can be calulated from day 22 to day 28
# to have an equal number acute and chronic values
d_sim_tl_hist_acwr_acute = d_sim_tl_hist %>% filter(day >= 22)
d_sim_hist_acute = function_on_list(d_sim_tl_hist_acwr_acute, FUN = acute_mean, 28) %>% rename(acute_load = data)
d_sim_hist_chronic = function_on_list(d_sim_tl_hist, FUN = chronic_ra, 28) %>% rename(chronic_load = data)

# calculate the ACWR be acute/chronic
d_sim_hist_acwr = d_sim_hist_acute %>% 
  left_join(d_sim_hist_chronic, by = c("id", "day")) %>% 
  mutate(acwr = acute_load/chronic_load)

# week to week change requires 2 full weeks before calculation
# the sliding window thereafter jumps one day at a time
# but the calculation is "uncoupled"
d_sim_hist_weekly = function_on_list(d_sim_tl_hist, FUN = weekly_sum, 7) %>% 
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
l_survival_sim_change_acwr = l_survival_sim_change %>% 
  map(. %>% left_join(d_sim_hist_acwr %>% 
                        select(-ends_with("load")), 
                      by = c("Id" = "id", "Stop" = "day")))

l_survival_sim_change_basemethods = l_survival_sim_change_acwr %>% 
  map(. %>% left_join(d_sim_hist_weekly, 
                      by = c("Id" = "id", "Stop" = "day")))

################################################### Fit the models ##################################################
library(rlang)
map_cox = function(l_sim_surv, var){
  var = enexpr(var)
  eval_bare(expr(l_sim_surv %>% map(., ~coxph(Surv(Start, Stop, Event) ~ !!var, ., y = FALSE, ties = "efron"))))
}
l_fit_ra = map_cox(l_survival_sim_basemethods, ra_t_load)
l_fit_ewma = map_cox(l_survival_sim_basemethods, ewma_t_load)
l_fit_acwr = map_cox(l_survival_sim_change_basemethods, acwr)
l_fit_weekly_change = map_cox(l_survival_sim_change_basemethods, weekly_change)

# amount, 1 is j_constant, 2 j_decay, 3 j_exponential_decay, 4 lin_direction_flip
# change, 1 is lin_constant, 2 lin_decay, 3 lin_exponential_decay
arglist_final = list(list(fun = "lin"), list(fun = "lin"), list(fun="ns", knots = 3), list(fun="ns", knots = 6))
arglist_final_change = list(list(fun = "lin"), list(fun = "lin"), list(fun="ns", knots = 3))
l_crossbases_amount = map2(.x = l_q_matrices,
                           .y = arglist_final,
                           ~crossbasis(.x,
                                       lag=c(lag_min, lag_max),
                                       argvar = list(fun="ns", knots = 3, intercept = FALSE),
                                       arglag = list(fun="ns", knots = 3, intercept = FALSE)))

l_crossbases_change = map2(.x = l_q_matrices_change,
                           .y = arglist_final_change,
                           ~crossbasis(.x,
                                       lag=c(lag_min, lag_max),
                                       argvar = list(fun="lin", intercept = FALSE),
                                       arglag = .y))


l_fit_dlnm_amount = map2(.x = l_counting_survival_sim,
                         .y = l_crossbases_amount,
                         ~coxph(Surv(enter, exit, event) ~ .y, 
                                .x, y = FALSE, ties = "efron"))


l_fit_dlnm_change = map2(.x = l_counting_survival_sim_change,
                         .y = l_crossbases_change,
                         ~coxph(Surv(enter, exit, event) ~ .y, 
                                .x, y = FALSE, ties = "efron"))

# RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
# AIC, RMSE, coverage
aic_j_decay_ra = AIC(l_fit_ra[[2]])
aic_j_decay_ewma = AIC(l_fit_ewma[[2]])
aic_j_decay_dlnm = AIC(l_fit_dlnm_amount[[2]])

aic_lin_decay_acwr = AIC(l_fit_acwr[[2]])
aic_lin_decay_weekly_change = AIC(l_fit_weekly_change[[2]])
aic_lin_decay_dlnm = AIC(l_fit_dlnm_change[[2]])

# Root-mean-squared-error
rmse = function(estimate, target){
  sqrt(mean((estimate - target)^2)) 
}

#function for obtaining root mean squared error RMSE
rmse = function(fit){
  rmse = sqrt(mean(fit$residuals^2))
  rmse
}

rmse(l_fit_ra[[2]])
rmse(l_fit_ewma[[2]])
rmse(l_fit_dlnm_amount[[2]])

parameters::parameters(l_fit_ra[[2]])
parameters::parameters(l_fit_ewma[[2]])
parameters::parameters(l_fit_dlnm_amount[[2]])

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

# 3D GRAPHS OF PREDICTED VALUES FOR ASSESSING MODEL FIT
cb_j_decay = crossbasis(l_q_matrices[[2]],
                        lag=c(lag_min, lag_max),
                        argvar = list(fun="ns", knots = 3, intercept = FALSE),
                        arglag = list(fun="ns", knots = 3))

mod_j_decay = coxph(Surv(enter, exit, event) ~ l_crossbases_amount[[2]], l_counting_survival_sim[[2]], y = FALSE, ties = "efron")

# for a onebasis model (rolling average)
test_ra = l_survival_sim_basemethods[[2]]$ra_t_load
ob_ra = onebasis(test_ra,"lin")
mod_j_decay_ra = coxph(Surv(Start, Stop, Event) ~ ob_ra, l_survival_sim_basemethods[[2]], y = FALSE, ties = "efron")
pred_j_decay_ra = crosspred(ob_ra, mod_j_decay_ra, at = tl_predvalues, cen = 300, cumul = TRUE)

# for change in load
cb_lin_decay = crossbasis(l_q_matrices_change[[2]],
                          lag=c(lag_min, lag_max),
                          argvar = list(fun="ns", knots = 3, intercept = FALSE),
                          arglag = list(fun="ns", knots = 3))

mod_lin_decay = coxph(Surv(enter, exit, event) ~ cb_lin_decay, l_counting_survival_sim_change[[2]], y = FALSE, ties = "efron")

# OBTAIN THE PREDICTED RISK FOR A SEQUENCE OF TL LEVELS
cb_j_decay = l_crossbases_amount[[2]]

pred_j_decay = crosspred(cb_j_decay, mod_j_decay, at = tl_predvalues, cen = 300, cumul = TRUE)
pred_lin_decay = crosspred(cb_lin_decay, mod_lin_decay, at = tl_predvalues_change, cen = 0, cumul = TRUE)


# j decay
persp(x = tl_predvalues, y = lag_seq, l_coefs[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, pred_j_decay$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues_change, y = lag_seq, pred_lin_decay$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, 
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

library(ggeffects)
preds_ra = ggpredict(l_fit_ra[[2]], terms= "ra_t_load [tl_predvalues]", type = "survival")
ggplot(preds_ra, aes(x, predicted)) + geom_line()

preds_ra[[3]]

fit_ra = survfit(Surv(Start, Stop, Event) ~ ra_t_load, data = l_survival_sim_basemethods[[2]])
class(fit_ra)

library(survminer)
ggsurvplot(
  fit_ra,
  size = 0.5,
  # lines become thicker in the pdf version
  conf.int = FALSE,
  risk.table = FALSE,
  tables.col = NULL,
  palette = nih_distinct[4],
  ggtheme = theme_line(), # from ostrc-package
  censor = FALSE,
  ylab = "",
  xlab = "Days",
  title = "Probability of Sports Injury",
  tables.y.text = FALSE,
  fontsize = 10,
  axes.offset = FALSE,
  break.x.by = 2
)
?ggsurvplot

# lists of functions
flist = list(fj, fj, fj, flin)
wlist = list(wconst, wdecay, wexponential_decay, wdirection_flip)
wlist_change = list(wconst, wdecay, wexponential_decay)
