
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
tl_observed_change = 2 - tl_observed
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


dag0 = l_coefs_change[[1]][,1]
plot(x = tl_predvalues_change, y = dag0)

#################################################Simulating exposure histories###################################
# simulate exposure histories
# for nsub number of subjects
# and t_max days under study
set.seed(1234)
l_tl_hist = vector(mode = "list", length = nsub)
l_tl_hist = l_tl_hist %>% map(~sample(tl_valid, t_max))
# add IDs and calc the change in training load from previous day
d_sim_tl_hist = l_tl_hist %>% 
  map(. %>% enframe(name = NULL)) %>% 
  bind_rows() %>% 
  rename(t_load = value) %>% 
  mutate(id = rep(1:nsub, each = t_max),
         day = rep(1:t_max, nsub),
         t_load_change = lead(t_load)-t_load) 

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

# Test: using gasparrini method = same results
# test_sim = d_sim_tl_hist %>% filter(id == 1)
# testsim_exphist = exphist(test_sim %>% pull(t_load), lag = c(lag_min, lag_max))
# l_exposure_seqs = split(testsim_exphist, seq(nrow(testsim_exphist)))
# decay_coefs = do.call(wdecay,list(0:28))
# exposure_coefs = l_exposure_seqs %>% map(., flin)
# combfun_effect = exposure_coefs %>% map(., ~.*decay_coefs)
# cumeffect = combfun_effect %>% map(., sum)
# cumeffect_mat = unlist(cumeffect)
# as.numeric(cumeffect_mat)
# calc_cumeff(test_sim$t_load, fvar = flin, flag = wdecay)

# function for calculating the cumulative effect for each individual in a dataset
# with the choice of var-function (fvar) and lag-function (flag)
# and specification of whether it is the amount of training load "amount"
# or change in training load "change" the calculation is performed on.
calc_cumeffs_all = function(l_sim_tl_hist, t_load_type = "amount", fvar, flag){
  
  if(t_load_type == "amount"){
    l_sim_tl_hist$data = l_sim_tl_hist$data %>% map(~calc_cumeff(.$t_load, fvar, flag))  
    } else if(t_load_type == "change"){
    l_sim_tl_hist$data = l_sim_tl_hist$data %>% map(~calc_cumeff(.$t_load_change, fvar, flag))
  }
  d_cumeffs = unnest(l_sim_tl_hist, cols = c(data)) %>% rename(cumeff = data)
  d_cumeffs
}
# test: calc_cumeffs_all(l_load_nested, fj, wdecay)

# lists of functions
flist = list(fj, fj, fj, flin)
wlist = list(wconst, wdecay, wexponential_decay, wdirection_flip)
wlist_change = list(wconst, wdecay, wexponential_decay)

# labels that allow the identification of the relationship in each element in a list
names_rels = c("constant", "decay", "exponential_decay")
name_extra = "flin_direction_flip"

# nest simulated dataset on each individual
l_load_nested = d_sim_tl_hist %>% nest(data = c(day, t_load, t_load_change)) 

# calc cumulative effects for each function in flist matched to each function in wlist
l_cumeffs = map2(.x = flist,
                 .y = wlist, 
                 ~calc_cumeffs_all(l_load_nested, "amount", .x, .y))
# add labels and collapse to dataset
l_cumeffs = map2(.x = l_cumeffs,
                 .y = c(names_rels, name_extra), 
                 ~.x %>% mutate(relationship = paste0(.y)))

# do the same for change. the difference is that flin is used for all functions in wlist_change.
l_cumeffs_change = wlist_change %>% 
                   map(~calc_cumeffs_all(l_load_nested, "change", flin, .)) %>% 
                   map2(.x = .,
                        .y = names_rels, 
                        ~.x %>% mutate(relationship = paste0(.y)))

# calculate cumulative effect of training load for each individual,
# as a cross-basis of a function f or exposure amount, and function w, lag-time
# l_load_nested = d_sim_tl_hist %>% nest(data = c(day, t_load, t_load_change)) 
# l_load_nested$data = l_load_nested$data %>% map(~calc_cumeff(.$t_load, fj, wdecay))
# d_load_cumeffs_j_decay = unnest(l_load_nested, cols = c(data)) %>% rename(cumeff = data)

# to have 2 examples
# l_load_nested = d_sim_tl_hist %>% nest(data = t_load) 
# l_load_nested$data = l_load_nested$data %>% map(~calc_cumeff(.$t_load, flin, wdecay))
# d_load_cumeffs_lin_decay = unnest(l_load_nested, cols = c(data)) %>% rename(cumeff = data)

########################################Simulate injuries based on cumulative effect of training load#####################
l_cumeffs_mats = l_cumeffs %>% map(. %>% mutate(day = rep(1:t_max, nsub)) %>% 
                    pivot_wider(names_from = id, values_from = cumeff) %>% 
                    select(-day, -relationship) %>% as.matrix)

l_cumeffs_mats_change = l_cumeffs_change %>% map(. %>% mutate(day = lead(rep(1:t_max, nsub))) %>% 
                                                   filter(!is.na(cumeff)) %>% 
                           pivot_wider(names_from = id, values_from = cumeff) %>% 
                           select(-day, -relationship) %>% as.matrix)

set.seed(1234)
l_survival_sim = l_cumeffs_mats %>% map(., ~permalgorithm(nsub, t_max, Xmat = .,
                                   censorRandom = runif(nsub, 1, t_max*2), betas=1))

l_survival_sim_change = l_cumeffs_mats_change %>% map(., ~permalgorithm(nsub, t_max, Xmat = .,
                                  censorRandom = runif(nsub, 1, (t_max)*2), betas=1))


# arrange in a matrix which will be used later
# cumeff_j_decay = l_cumeffs[[2]]
# cumeff_lin_decay = l_cumeffs_change[[2]]
# cumeff_j_decay_i = 
#   d_load_cumeffs_j_decay %>% mutate(day = rep(1:t_max, nsub)) %>% 
#   pivot_wider(names_from = id, values_from = cumeff) %>% select(-day) %>% as.matrix
# 
# cumeff_lin_decay_i = 
#   d_load_cumeffs_lin_decay %>% mutate(day = rep(1:t_max, nsub)) %>% 
#   pivot_wider(names_from = id, values_from = cumeff) %>% select(-day) %>% as.matrix

# simulate survival times and events for each participant
# set the beta-coefficient to be 1, as in 1 times the cumeffect provided in the matrix
# set censoring probability to 0.10
# use cumulative effect as the time-dependent variable
# set.seed(1234)
# d_sim_surv_j_decay = permalgorithm(nsub, t_max, Xmat = cumeff_j_decay_i,
#                       censorRandom = runif(nsub, 1, t_max*2), betas=1)
# 
# d_sim_surv_lin_decay = permalgorithm(nsub, t_max, Xmat = cumeff_lin_decay_i,
#                                    censorRandom = runif(nsub, 1, t_max*2),betas=1)
# 
# d_sim_surv_j_decay %>% summarise(n_events = sum(Event == 1))
# d_sim_surv_lin_decay %>% summarise(n_events = sum(Event == 1))

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
  d_sim_tl_hist  %>% select(-t_load_change) %>% 
  pivot_wider(names_from = day, values_from = t_load) %>% select(-id) %>% as.matrix

d_sim_tl_hist_spread_day_change = 
  d_sim_tl_hist %>% select(-t_load) %>%
  pivot_wider(names_from = day, values_from = t_load_change) %>% select(-id) %>% as.matrix

# list of q-matrices
l_q_matrices = l_survival_sim %>% map(., ~from_sim_surv_to_q(., d_sim_tl_hist_spread_day))
l_q_matrices_change = l_survival_sim_change %>% map(., ~from_sim_surv_to_q(., d_sim_tl_hist_spread_day_change))

#-------------------------------------------Sensitivity analysis to determine best model
# sensitivity analysis
# what combination of lag and var functions
# have the best AIC for each of the Q matrices?
knots = 3:6
arglist_splines = knots %>% map(., ~list(fun="ns", knots = .))
arglist_poly = list(list(fun="poly",degree=2))
arglist_lin = list(list(fun = "lin"))

arglist1 = append(arglist_splines, arglist_poly)
arglist = append(arglist1, arglist_lin)

# instead of making a list of lists, copy-paste
# easier to keep head straight on permutations
# this needs careful coding
crossbasis_j_constant = arglist %>% 
  map(., ~crossbasis(l_q_matrices[[1]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="ns", knots = 3, intercept = FALSE),
                     arglag = .))

crossbasis_j_decay = arglist %>% 
  map(., ~crossbasis(l_q_matrices[[2]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="ns", knots = 3, intercept = FALSE),
                     arglag = .))

crossbasis_j_exponential_decay = arglist %>% 
  map(., ~crossbasis(l_q_matrices[[3]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="ns", knots = 3, intercept = FALSE),
                     arglag = .))

crossbasis_flin_direction_flip = arglist %>% 
  map(., ~crossbasis(l_q_matrices[[4]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="lin", intercept = FALSE),
                     arglag = .))

# need data on counting process form
l_counting_survival_sim = l_survival_sim %>% map(., ~counting_process_form(.))

names_dlnm_funs = c("splines_3", "splines_4", "splines_5", "splines_6", "poly2", "lin")
names_truerel = c(names_rels, name_extra)
aic_vec = c()
for(i in 1:length(arglist)){
  surfit_j_constant = coxph(Surv(enter, exit, event) ~ crossbasis_j_constant[[i]], l_counting_survival_sim[[1]], y = FALSE, ties = "efron")
  surfit_j_decay = coxph(Surv(enter, exit, event) ~ crossbasis_j_decay[[i]], l_counting_survival_sim[[2]], y = FALSE, ties = "efron")
  survfit_j_exponential_decay = coxph(Surv(enter, exit, event) ~ crossbasis_j_exponential_decay[[i]], l_counting_survival_sim[[3]], y = FALSE, ties = "efron")
  survfit_lin_direction_flip = coxph(Surv(enter, exit, event) ~ crossbasis_flin_direction_flip[[i]], l_counting_survival_sim[[4]], y = FALSE, ties = "efron")
  aic = c(AIC(surfit_j_constant), AIC(surfit_j_decay), AIC(survfit_j_exponential_decay), AIC(survfit_lin_direction_flip))
  aic_vec = append(aic_vec, aic)
}

d_aic_amount_load = enframe(aic_vec, name = NULL) %>% 
                    mutate(fun_names = rep(names_dlnm_funs, each = 4), 
                           true_rels = rep(names_truerel, length(arglist))) %>% rename(aic = value)
best_aics_tl_amount = d_aic_amount_load %>% group_by(true_rels) %>% filter(aic == min(aic))

#--------------------do the same for change in load

crossbasis_lin_constant = arglist %>% 
  map(., ~crossbasis(l_q_matrices_change[[1]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="lin"),
                     arglag = .))

crossbasis_lin_decay = arglist %>% 
  map(., ~crossbasis(l_q_matrices_change[[2]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="lin"),
                     arglag = .))

crossbasis_lin_exponential_decay = arglist %>% 
  map(., ~crossbasis(l_q_matrices_change[[3]],
                     lag=c(lag_min, lag_max),
                     argvar = list(fun="lin"),
                     arglag = .))

# need data on counting process form
l_counting_survival_sim_change = l_survival_sim_change %>% map(., ~counting_process_form(.))
aic_vec_change = c()
for(i in 1:length(arglist)){
  survfit_lin_constant = coxph(Surv(enter, exit, event) ~ crossbasis_lin_constant[[i]], l_counting_survival_sim_change[[1]], y = FALSE, ties = "efron")
  survfit_lin_decay = coxph(Surv(enter, exit, event) ~ crossbasis_lin_decay[[i]], l_counting_survival_sim_change[[2]], y = FALSE, ties = "efron")
  survfit_lin_exponential_decay = coxph(Surv(enter, exit, event) ~ crossbasis_lin_exponential_decay[[i]], l_counting_survival_sim_change[[3]], y = FALSE, ties = "efron")
  aic_change = c(AIC(survfit_lin_constant), AIC(survfit_lin_decay), AIC(survfit_lin_exponential_decay))
  aic_vec_change = append(aic_vec_change, aic_change)
}

d_aic_change_load = enframe(aic_vec_change, name = NULL) %>% 
                    mutate(fun_names = rep(names_dlnm_funs, each = 3), 
                    true_rels = rep(names_rels, length(arglist))) %>% rename(aic = value)
best_aics_tl_change = d_aic_change_load %>% group_by(true_rels) %>% filter(aic == min(aic))



####################################### Running the different models of training load ####################################

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
  nested_list = d_sim_hist %>% nest(data = c(t_load, day, t_load_change))
  nested_list$data = nested_list$data %>% map(., ~FUN(.$t_load))
  l_unnest = unnest(nested_list, cols = c(data)) %>% mutate(day = rep(day_start:t_max, nsub)) 
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


# RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
modellist = list()
modellist[[1]] = coxph(Surv(enter, exit, event) ~ cb_j_decay, dataspl_j_decay, y = FALSE, ties = "efron")
modellist[[2]] = coxph(Surv(enter, exit, event) ~ cb_lin_decay, dataspl_lin_decay, y = FALSE, ties = "efron")
mod_j_decay = modellist[[1]]
mod_lin_decay = modellist[[2]]
AIC(mod_j_decay)
AIC(mod_lin_decay)





# DEFINE THE CROSS-BASIS
# you can make a list of functions for argvar and arglag, with splines etc.
cb_j_decay = crossbasis(q_j_decay,
                        lag=c(lag_min, lag_max),
                        argvar = list(fun="ns", knots = 3, intercept = FALSE),
                        arglag = list(fun="lin", intercept=TRUE))

cb_lin_decay = crossbasis(q_lin_decay, 
                          lag=c(lag_min, lag_max), 
                          argvar = list("lin", intercept = FALSE), 
                          arglag = list(fun="lin", intercept = TRUE))


#------j decay
  # We need a Q-matrix to define the crossbasis
  # we will arrange the data in a counting process form
  # first, restrict data to event/censor times
  d_follow_up_times_j = d_sim_surv_j_decay %>% distinct(Id, Fup)
  d_surv_lim_j_decay = map2(.x = d_follow_up_times_j$Id,
       .y = d_follow_up_times_j$Fup,
       ~d_sim_surv_j_decay %>% filter(Id == .x) %>% slice(.y)) %>% 
    bind_rows() %>% select(event = Event, exit = Stop, id = Id)
  
  # extracting timepoints in which an event happened
  ftime = d_surv_lim_j_decay %>% filter(event == 1) %>% distinct(exit) %>% arrange(exit) %>% pull()
  # arrange the survival data so that, for each individual, we have an interval of enter and exit times
  # for each of the exit times above, with the information of whether or not they were injured at that time
  # meaning we will have the same time intervals per participant
  dataspl_j_decay = survSplit(Surv(exit, event)~., d_surv_lim_j_decay, cut = ftime, start="enter") %>% arrange(id)
  
  # for each individual, for each of these exit times, we will extract the exposure history 
  # for the given lag-time which we are interested in
  # This is called the Q-matrix. The Q-matrix should be nrow(dataspl) X 0:lag_max dimensions.
  q_j_decay = 
      map2(.x = dataspl_j_decay$id, 
           .y = dataspl_j_decay$exit, 
           ~exphist(d_sim_tl_hist_spread_day[.x,], .y, c(lag_min, lag_max))) %>% 
    do.call("rbind", .)
  
#------lin decay
  
  # We need a Q-matrix to define the crossbasis
  # we will arrange the data in a counting process form
  # first, restrict data to event/censor times
  d_follow_up_times_lin = d_sim_surv_lin_decay %>% distinct(Id, Fup)
  d_surv_lim_lin_decay = map2(.x = d_follow_up_times_lin$Id,
                            .y = d_follow_up_times_lin$Fup,
                            ~d_sim_surv_lin_decay %>% filter(Id == .x) %>% slice(.y)) %>% 
                             bind_rows() %>% select(event = Event, exit = Stop, id = Id)

  # extracting timepoints in which an event happened
  ftime = d_surv_lim_lin_decay %>% filter(event == 1) %>% distinct(exit) %>% arrange(exit) %>% pull()
  # arrange the survival data so that, for each individual, we have an interval of enter and exit times
  # for each of the exit times above, with the information of whether or not they were injured at that time
  # meaning we will have the same time intervals per participant
  dataspl_lin_decay = survSplit(Surv(exit, event)~., d_surv_lim_lin_decay, cut = ftime, start="enter") %>% arrange(id)
  
  # for each individual, for each of these exit times, we will extract the exposure history 
  # for the given lag-time which we are interested in
  # This is called the Q-matrix. The Q-matrix should be nrow(dataspl) X 0:lag_max dimensions.
  # testsim = d_sim_tl_hist %>% filter(id == 1) %>% pull(t_load)
  # slide(testsim, ~exphist(., times = lag_max, lag = c(0,27)), .before = lag_max-1, step = 1, .complete = FALSE)

  # CREATE THE MATRIX Q OF LAGGED EXPOSURES USING THE FUNCTION exphist()
  # LAGGED EXPOSURES IS BASED ON TIME SINCE STUDY FIRST MEASURE (STUDY START)
  q_lin_decay = map2(.x = dataspl_lin_decay$id, 
                     .y = dataspl_lin_decay$exit, 
                        ~exphist(d_sim_tl_hist_spread_day[.x,], .y, c(lag_min, lag_max))) %>% 
                     do.call("rbind", .) 

# DEFINE THE CROSS-BASIS
# you can make a list of functions for argvar and arglag, with splines etc.
cb_j_decay = crossbasis(q_j_decay,
                lag=c(lag_min, lag_max),
                argvar = list(fun="ns", knots = 3, intercept = FALSE),
                arglag = list(fun="lin", intercept=TRUE))

cb_lin_decay = crossbasis(q_lin_decay, 
                lag=c(lag_min, lag_max), 
                argvar = list("lin", intercept = FALSE), 
                arglag = list(fun="lin", intercept = TRUE))

#list(knots = logknots(20,nk=3))

# RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
modellist = list()
modellist[[1]] = coxph(Surv(enter, exit, event) ~ cb_j_decay, dataspl_j_decay, y = FALSE, ties = "efron")
modellist[[2]] = coxph(Surv(enter, exit, event) ~ cb_lin_decay, dataspl_lin_decay, y = FALSE, ties = "efron")
mod_j_decay = modellist[[1]]
mod_lin_decay = modellist[[2]]
AIC(mod_j_decay)
AIC(mod_lin_decay)

#--------------------------------------------Figures

################################################################################
# PREDICTION ALONG TIME

# OBTAIN THE PREDICTED RISK FOR A SEQUENCE OF TL LEVELS
pred_j_decay = crosspred(cb_j_decay, mod_j_decay, at = tl_predvalues, cen = 300, cumul = TRUE)
pred_lin_decay = crosspred(cb_lin_decay, mod_lin_decay, at = tl_predvalues, cen = 0, cumul = TRUE)

# 3D GRAPHS OF PREDICTED VALUES FOR ASSESSING MODEL FIT

# j decay
persp(x = tl_predvalues, y = lag_seq, l_coefs[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, pred_j_decay$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

# lin decay

persp(x = tl_predvalues, y = lag_seq, l_coefs[[1]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, pred_lin_decay$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, 
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")


# plot the lag-response curves for specific and incremental cumulative effects
plot(pred_j_decay, "slices", var=175, col=2, ylab="HR", ci.arg=list(density=15,lwd=2),
     main="Lag-response curve for a 111-unit increase in sRPE")
plot(pred_j_decay, "slices", var=175, col=2, cumul=TRUE, ylab="Cumulative HR",
     main="Lag-response curve of incremental cumulative effects")

# lag-response curve for dose 60
plot(pred_lin_decay, var=1000, ylab="HR for exposure 175", xlab="Lag (days)", xlim=c(0, 28))
# exposure-response curve for lag 10
plot(pred_lin_decay, lag=0, ylab="HR at lag 10", xlab="sRPE", ylim=c(0, 6), xlim=c(0, 1000))

# ARGUMENTS FOR 3D PLOTS
arg3D <- list(x=seq(0,10,0.25),y=0:20*2,ticktype="detailed",theta=230,
              ltheta=200,phi=30,lphi=30,xlab="Exposure",ylab="Lag",zlab="HR",
              shade = 0.75,r=sqrt(3),d=5,cex.axis=0.7,cex.lab=0.8,border=grey(0.3),
              col=grey(0.99))

# alternative ggplot2 3D figures
#devtools::install_github("AckerDWM/gg3D", force = TRUE)
library(gg3D)
d_coef_1 = l_coefs[[1]] %>%
  reshape2::melt() %>% 
  as_tibble() %>% 
  transmute(t_load = as.numeric(Var1), lag = as.numeric(Var2), coef = value)

ggplot(d_coef_1, aes(x = t_load, y = lag, z = coef, color=coef)) +
  axes_3D() +
  stat_wireframe(alpha=.8, theta=230, phi=40) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors=plot3D::jet2.col()) +
  labs_3D(hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 90))


#----------------------------------------function for simulating an exposure history

# function to simulate an exposure history for one participant
# with number of days t
# and training load sampled from a vector of values t_load
# if you would like to have different starting times for different athletes

sim_exp_history = function(t, t_load){
  
  # we assume that athletes are measured at the end of vacation,
  # which is a few days of 0
  # before preseason starts
  # but that the startday is different for different athletes
  start = round(runif(1, 1, 5), 0) # individual start date
  duration =  7 + 7*rpois(1,3) # duration in days
  tl = sample(t_load) # training load exposure
  vec = c(rep(0, start-1), rep(tl, duration))
  
  while (length(vec) <= t){
    intermission = 21 + 7*rpois(1,3) # in days
    duration =  7 + 7*rpois(1,3) # in days
    tl =  sample(tl_observed)
    vec = append(vec, c(rep(0, start), rep(tl, duration)))
  }
  
  return(vec[1:t])
}



