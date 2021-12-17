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
l_load_nested = d_sim_tl_hist %>% group_by(id) %>% nest() 

# calc cumulative effects for each function in flist matched to each function in wlist
l_cumeffs = map2(.x = flist,
                 .y = wlist, 
                 ~calc_cumeffs_all(l_load_nested, t_load_type = "amount", .x, .y))
# add labels and collapse to dataset
l_cumeffs = map2(.x = l_cumeffs,
                 .y = c(names_rels, name_extra), 
                 ~.x %>% mutate(relationship = paste0(.y)))

# do the same for change. the difference is that flin is used for all functions in wlist_change.
l_cumeffs_change = wlist_change %>% 
  map(~calc_cumeffs_all(l_load_nested, t_load_type = "change", flin, .)) %>% 
  map2(.x = .,
       .y = names_rels, 
       ~.x %>% mutate(relationship = paste0(.y)))

########################################Simulate injuries based on cumulative effect of training load#####################
l_cumeffs_mats = l_cumeffs %>% map(. %>% mutate(day = rep(1:t_max, nsub)) %>% 
                                     pivot_wider(names_from = id, values_from = cumeff) %>% 
                                     select(-relationship) %>% as.matrix)

l_cumeffs_mats_change = l_cumeffs_change %>% 
  map(. %>% mutate(day = lead(rep(1:t_max, nsub))) %>% 
        filter(!is.na(cumeff)) %>% 
        pivot_wider(names_from = id, values_from = cumeff) %>% 
        select(-relationship) %>% as.matrix)

set.seed(1234)
l_survival_sim = l_cumeffs_mats %>% map(., ~permalgorithm(nsub, t_max, Xmat = .,
                                                          censorRandom = runif(nsub, 1, t_max*2), betas=1))

l_survival_sim_change = l_cumeffs_mats_change %>% map(., ~permalgorithm(nsub, t_max, Xmat = .,
                                                                        censorRandom = runif(nsub, 1, (t_max)*2), betas=1))

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