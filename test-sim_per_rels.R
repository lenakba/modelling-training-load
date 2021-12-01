
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

# 200 pers, 1 sesong er ca. 150 skader 
nsub = 250 # number of subjects
t_max = 300 # number of days (length of study)

symmetrized_change = function(x, y){
  100*((x - y)/(x + y))
}

# observed values
tl_observed = (d_load %>% filter(srpe <= 1200))$srpe
#tl_observed_change = lead(tl_observed)-tl_observed
tl_observed_change = symmetrized_change(lead(tl_observed), tl_observed)
tl_observed_change = tl_observed_change[-length(tl_observed_change)]
tl_valid = min(tl_observed):max(tl_observed)
# lag set at 4 weeks (28) as is often used in tl studies
lag_min = 0
lag_max = 27
lag_seq = lag_min:lag_max # number of days before current day assumed to affect risk of injury

# vector of tl values used in visualizations of predictions
tl_predvalues = seq(min(tl_observed), max(tl_observed), 10)
tl_predvalues_change = seq(min(tl_observed_change, na.rm = TRUE), max(tl_observed_change, na.rm = TRUE))
tl_predvalues_acwr = seq(0.1, 3.5)
tl_predvalues_weekly_change = seq(-80, 80)

###################################Training load and lag structure functions###########################################
# FUNCTIONS TO COMPUTE THE LOGIT AND INVERSE LOGIT TRANSFORMATIONS
# logit = function(prob) log(prob/(1-prob))
# invlogit = function(linpred) exp(linpred)/(1+exp(linpred))

# functions for simulating the effect of the amount of training load on injury risk (J-shape)
# and the effect of amount of change of training load on injury risk (linear)
# the only exception is for the lag-function wdirection_flip, the linear function will be used instead of J-shape
fj = function(x)case_when(x < 600 ~ ((600-x)/200)^1.5/10,
                          x >= 600 ~ ((x-600)/200)^3/30)
# for absolute change (no longer used)
flin_amount = function(x) 0.0009*x
flin = function(x) 0.009*x

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
flindirection_flip = function(x, lag) flin_amount(x) * wdirection_flip(lag)
fw_funs = list(fjconst, fjdecay, fjexponential_decay, flindirection_flip)

# f*w functions for change in load
flinconst = function(x, lag) flin(x) * wconst(lag)
flindecay = function(x, lag) flin(x) * wdecay(lag)
flinexponential_decay = function(x, lag) flin(x) * wexponential_decay(lag)
fw_funs_change = list(flinconst, flindecay, flinexponential_decay)

# f*w functions for relative change in load
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
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[3]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="RR", shade=0.75, r=sqrt(3), d=5, cex.axis=0.7, cex.lab=0.8,
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)")


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
           day = rep(1:t_max, nsub)) %>% 
    group_by(id) %>% 
    mutate(t_load_change = lag(symmetrized_change(lead(t_load), t_load))) %>% 
    ungroup()  
  d_sim_tl_hist %>% select(id, day, t_load, t_load_change)
}

set.seed(1234)
d_sim_tl_hist = sim_tl_history(nsub, t_max, tl_valid)

#############################Calculate cumulative effect of training load given exposure histories#######################

# function to compute the cumulative effect of training load exposure, requires:
# tl_hist     the tl exposure history from t0 to max(t) for one participant
# fvar        the function for the effect of training load
# flag        the function for the time-lagged effect of training load
calc_cumeff = function(tl_hist, fvar, flag){
  l_effect_seqs = slide(tl_hist, ~pluck(.), .before = lag_max, step = 1, .complete = FALSE)
  l_effect_lags = compact(l_effect_seqs) %>% map(~(NROW(.)-1):0)
  l_flin = l_effect_seqs %>% map(~fvar(.))
  l_wdecay = l_effect_lags %>% map(~flag(.))
  combfun_effect = l_flin %>% map2(.x = ., .y = l_wdecay, ~ .x * .y)
  cumeffect = combfun_effect %>% map(., ~sum(., na.rm = TRUE))
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

d_cumeff = calc_cumeffs_all(d_sim_tl_hist, t_load_type = "change", flin, wdecay)

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
d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "change", fvar = flin, flag = wdecay)

# estimate common number of events
set.seed(1234)
vec_fjwconst = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "amount", fvar = fj, flag = wconst)
  vec_fjwconst[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}
vec_fjwdecay = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "amount", fvar = fj, flag = wdecay)
  vec_fjwdecay[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}
vec_fjwexponential_decay = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "amount", fvar = fj, flag = wexponential_decay)
  vec_fjwexponential_decay[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}

vec_flinwconst = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "change", fvar = flin, flag = wconst)
  vec_flinwconst[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}
vec_flinwdecay = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "change", fvar = flin, flag = wdecay)
  vec_flinwdecay[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}
vec_flinwexponential_decay = rep(0, 100)
for(i in 1:100){
  d_survival_sim = sim_survdata(nsub, t_max, tl_valid, t_load_type = "change", fvar = flin, flag = wexponential_decay)
  vec_flinwexponential_decay[i] = d_survival_sim %>% summarise(sum(Event)) %>% pull()
}

overall_mean_injuries = mean(c(vec_fjwconst, vec_fjwdecay, vec_fjwexponential_decay, vec_flinwconst, vec_flinwdecay, vec_flinwexponential_decay))
d_n_injuries_rels = tibble(vec_fjwconst, vec_fjwdecay, vec_fjwexponential_decay, vec_flinwconst, vec_flinwdecay, vec_flinwexponential_decay) 
mean_injuries = d_n_injuries_rels %>% summarise_all(mean)


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

# arrange the exposure history in wide format in a matrix
# which is neeeded for calculating the q-matrix for the crossbasis
d_sim_tl_hist_spread_day = 
  d_sim_tl_hist %>% select(-t_load) %>% 
  pivot_wider(names_from = day, values_from = t_load_change) %>% select(-id) %>% as.matrix
d_survival_sim_cpform = counting_process_form(d_survival_sim)
q_mat = calc_q_matrix(d_survival_sim_cpform, d_sim_tl_hist_spread_day)
q_mat_same_n = q_mat %>% as_tibble() %>% mutate(day = as.numeric(rownames(q_mat))) %>% 
               filter(day >= lag_max+1) %>% select(-day) %>% as.matrix

####################################### Modify training load with different methods ####################################

# Perform the typical methods for handling training load amount: rolling average and EWMA

# Function for calculating rolling averages on a chooseable number of days
# Based on rollapplyr, not rollmean, as rollmean will only start calculating averages
# at n values, while rollapplyr allows the user to decide preliminary values.
ra = function(x, n_days = lag_max+1, window = TRUE, ...){
  zoo::rollapplyr(x, n_days, mean, partial = FALSE)
}

# function for calculating exponentially waited moving averages
# using similar syntax as the RA-function
# wilder=FALSE (the default) uses an exponential smoothing ratio of 2/(n+1)
# same as in williams et al. 2016
ewma = function(x, n_days = lag_max+1){
  TTR::EMA(x, n = n_days, wilder = FALSE)
}

# robust exponential decreasing index (REDI)
# for details, see: http://dx.doi.org/10.1136/bmjsem-2019-000573
redi = function(x, n_days = lag_max+1, lag = lag_seq, lambda = 0.1){
  lambda = lambda
  lag_effect = exp((-lambda)*lag)
  t_max = length(x)
  vec = n_days:t_max
  for(i in 1:(t_max-(n_days-1))){
    tl = lead(x, i-1)[1:n_days]
    vec[i] = sum(lag_effect*tl)/sum(lag_effect)
  }
  vec 
}

# functions for calculating ra on a sliding window that moves one day at a time.
# this ensures that RA isn't calculated until the first 4 weeks have passed
slide_ra = function(x){
  l = slide(x, ~ra(., lag_max+1), .before = lag_max, step = 1, .complete = TRUE) %>% map(last)
  l = compact(l)
  l = unlist(l)
  l
}

# function to nest the exposure history data by each individual, 
# and run a user-specified function on each of their datasets in the list
function_on_list = function(d_sim_hist, FUN = NULL, day_start){
  nested_list = d_sim_hist %>% group_by(id) %>% nest()
  nested_list$data = nested_list$data %>% map(., ~FUN(.$t_load))
  l_unnest = unnest(nested_list, cols = c(data)) %>% ungroup() %>% 
    filter(!is.na(data)) %>% mutate(day = rep(day_start:t_max, nsub)) 
  l_unnest
}

d_sim_hist_ra = function_on_list(d_sim_tl_hist, ra, lag_max+1) %>% rename(ra_t_load = data)
d_sim_hist_ewma = function_on_list(d_sim_tl_hist, ewma, lag_max+1) %>% rename(ewma_t_load = data)
d_sim_hist_redi = function_on_list(d_sim_tl_hist, redi, lag_max+1) %>% rename(redi_t_load = data)

# calc rolling average and ewma on training load amount
d_survival_sim_cpform_mods = d_survival_sim_cpform %>% filter(exit >= lag_max+1)
d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% left_join(d_sim_hist_ra, by = c("id", "exit" = "day"))
d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% left_join(d_sim_hist_ewma, by = c("id", "exit" = "day"))
d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% left_join(d_sim_hist_redi, by = c("id", "exit" = "day"))
ob_ra = onebasis(d_survival_sim_cpform_mods$ra_t_load, "poly", degree = 2)
ob_ewma = onebasis(d_survival_sim_cpform_mods$ewma_t_load, "poly", degree = 2)
ob_redi = onebasis(d_survival_sim_cpform_mods$redi_t_load, "poly", degree = 2)
cb_dlnm = crossbasis(q_mat_same_n, lag=c(lag_min, lag_max), 
                     argvar = list(fun="poly", degree = 2),
                     arglag = list(fun="ns", knots = 3))

# Perform typical methods for handling change in training load, ACWR and week-to-week change

# calculate 7:28 coupled ACWR (this becomes, in theory, a measure of change)
# use equation in Lolli et al. 2017, most common in football studies Wang et al. 2021.
slide_sum = function(x){
  l = slide(x, ~sum(.), .before = 6, step = 1, .complete =TRUE)
  l = compact(l)
  l = unlist(l)
  l
}

slide_chronic = function(x){
  l = slide(x, ~sum(.), .before = lag_max, step = 1, .complete =TRUE) %>% map(~./4)
  l = compact(l)
  l = unlist(l)
  l
}

# the first acute load can be calculated from day 22 to day 28
# to have an equal number acute and chronic values
d_sim_tl_hist_acwr_acute = d_sim_tl_hist %>% filter(day >= lag_max-5)
d_sim_hist_acute = function_on_list(d_sim_tl_hist_acwr_acute, FUN = slide_sum, lag_max+1) %>% rename(acute_load = data)
d_sim_hist_chronic = function_on_list(d_sim_tl_hist, FUN = slide_chronic, lag_max+1) %>% rename(chronic_load = data)

# calculate the ACWR be acute/chronic
d_sim_hist_acwr = d_sim_hist_acute %>% 
  left_join(d_sim_hist_chronic, by = c("id", "day")) %>% 
  mutate(acwr = acute_load/chronic_load)

# week to week change requires 2 full weeks before calculation
# the sliding window thereafter jumps one day at a time
# but the calculation is "uncoupled"
d_sim_hist_weekly = function_on_list(d_sim_tl_hist, FUN = slide_sum, 7) %>% 
  rename(week_sum = data) %>% 
  mutate(week_sum_lead = lead(week_sum, 7),
         weekly_change = symmetrized_change(week_sum_lead, week_sum))

# the difference can't be measured until a week after the first week
d_sim_hist_weekly = d_sim_hist_weekly %>% 
  select(-week_sum, -week_sum_lead) %>% 
  group_by(id) %>% 
  mutate(day = lead(day, 7)) %>% 
  filter(!is.na(day)) %>% ungroup()

# couple survival data with change in load data
d_survival_sim_cpform_mods = d_survival_sim_cpform %>% filter(exit >= lag_max+1)

d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% 
                             left_join(d_sim_hist_acwr, by = c("id", "exit" = "day"))

d_survival_sim_cpform_mods = d_survival_sim_cpform_mods %>% 
                             left_join(d_sim_hist_weekly, by = c("id", "exit" = "day"))

ob_acwr = onebasis(d_survival_sim_cpform_mods$acwr, "lin")
ob_weekly_change = onebasis(d_survival_sim_cpform_mods$weekly_change, "lin")
cb_dlnm_change = crossbasis(q_mat_same_n, lag=c(lag_min, lag_max), 
                     argvar = list(fun="lin"),
                     arglag = list(fun="ns", knots = 3))

################################################### Fit the models ##################################################
library(rlang)
cox_basis = function(d_sim_surv, basis){
  basis = enexpr(basis)
  eval_bare(expr(coxph(Surv(enter, exit, event) ~ !!basis, d_sim_surv, y = FALSE, ties = "efron")))
}

fit_ra = cox_basis(d_survival_sim_cpform_mods, ob_ra)
fit_ewma = cox_basis(d_survival_sim_cpform_mods, ob_ewma)
fit_redi = cox_basis(d_survival_sim_cpform_mods, ob_redi)
fit_dlnm = cox_basis(d_survival_sim_cpform_mods, cb_dlnm)

fit_acwr = cox_basis(d_survival_sim_cpform_mods, ob_acwr)
fit_weekly_change = cox_basis(d_survival_sim_cpform_mods, ob_weekly_change)
fit_dlnm = cox_basis(d_survival_sim_cpform_mods, cb_dlnm_change)

################################################## Calculate numeric performance measures ############################
  
aic_ra = AIC(fit_ra)
aic_ewma = AIC(fit_ewma)
aic_redi = AIC(fit_redi)
aic_dlnm = AIC(fit_dlnm)  

anova(fit_ra, fit_ewma, fit_dlnm)

aic_acwr = AIC(fit_acwr)
aic_weekly_change = AIC(fit_weekly_change)
aic_dlnm = AIC(fit_dlnm)  

anova(fit_acwr, fit_weekly_change, fit_dlnm)

cp_preds_acwr = crosspred(ob_acwr, fit_acwr, at = tl_predvalues_acwr, cumul = TRUE)
cp_preds_weekly_change = crosspred(ob_weekly_change, fit_weekly_change, at = tl_predvalues_weekly_change, cumul = TRUE)
cp_preds_dlnm = crosspred(cb_dlnm_change, fit_dlnm, at = tl_predvalues_change, cumul = TRUE)

# p -values?
# RMSE for model fit?
sqrt(mean(fit_acwr$residuals^2))
sqrt(mean(fit_weekly_change$residuals^2))
sqrt(mean(fit_dlnm$residuals^2))


## compare true coefs vs. predicted coefs

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

######################################################## GRAPHS OF PREDICTED VALUES FOR ASSESSING MODEL FIT
d_cumulative_preds = bind_cols(t_load = tl_predvalues, 
                               ra = cp_preds_ra$allfit, 
                               ewma = cp_preds_ewma$allfit,
                               dlnm = cp_preds_dlnm$allfit) %>% 
                               pivot_longer(!t_load, names_to = "model", 
                                            values_to = "cumulative_effect")

calc_ci = function(estimate, se, dir = NULL){
  qn = qnorm(0.95)  
  if(dir == "high"){
    ci = estimate + qn*se
  } else if (dir == "low"){
    ci = estimate - qn*se
  }
}
cumeff_preds_ci_high_ra = calc_ci(cumeff_preds_ra, cp_preds_ra$allse, "high")
cumeff_preds_ci_low_ra = calc_ci(cumeff_preds_ra, cp_preds_ra$allse, "low")
cumeff_preds_ci_high_ewma = calc_ci(cumeff_preds_ewma, cp_preds_ewma$allse, "high")
cumeff_preds_ci_low_ewma = calc_ci(cumeff_preds_ewma, cp_preds_ewma$allse, "low")
cumeff_preds_ci_high_dlnm = calc_ci(cumeff_preds_dlnm, cp_preds_dlnm$allse, "high")
cumeff_preds_ci_low_dlnm = calc_ci(cumeff_preds_dlnm, cp_preds_dlnm$allse, "low")

cis_high = bind_cols(t_load = tl_predvalues, 
          ra = cumeff_preds_ci_high_ra, 
          ewma = cumeff_preds_ci_high_ewma,
          dlnm = cumeff_preds_ci_high_dlnm) %>% 
  pivot_longer(!t_load, names_to = "model", 
               values_to = "ci_high")

cis_low = bind_cols(t_load = tl_predvalues, 
                     ra = cumeff_preds_ci_low_ra, 
                     ewma = cumeff_preds_ci_low_ewma,
                     dlnm = cumeff_preds_ci_low_dlnm) %>% 
  pivot_longer(!t_load, names_to = "model", 
               values_to = "ci_low")

d_cumul_preds = d_cumulative_preds %>% left_join(cis_high, by = c("t_load", "model")) %>% left_join(cis_low, by = c("t_load", "model")) 
# true relationship 
d_true_coefs = truecumcoefs %>% enframe() %>% rename(t_load = name, truerel = value) %>% mutate(t_load = as.numeric(t_load))

text_size = 14
ggplot(d_cumul_preds, aes(x = t_load)) +
  facet_wrap(~model, ncol = 3, scales = "free") +
  geom_line(data = d_true_coefs, aes(y = truerel, color = "True relationship"), size = 0.5) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) + 
  geom_line(aes(y = cumulative_effect, color = "Estimation"), size = 0.5) +
  scale_y_continuous(limit = c(-15, 25), breaks = scales::breaks_width(5))  +
  scale_color_manual(values = c(nih_distinct[1], "Black")) + 
  ylab("Cumulative Hazard") +
  xlab("sRPE (AU) on Day 0") +
  theme_line(text_size, legend = TRUE) +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = nih_distinct[4]),
        strip.background = element_blank(),
        strip.text.x = element_text(size = text_size, family="Trebuchet MS", colour="black", face = "bold", hjust = -0.01),
        axis.ticks = element_line(color = nih_distinct[4]),
        legend.position = "bottom")

# time x amount effect (DLNM only)

# j decay
persp(x = tl_predvalues, y = lag_seq, exp(true_effect), ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, cp_preds_dlnm$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

#------------------Figures for change in load


cp_preds_acwr = crosspred(ob_acwr, fit_acwr, at = tl_predvalues_acwr, cumul = TRUE)
cp_preds_weekly_change = crosspred(ob_weekly_change, fit_weekly_change, at = tl_predvalues_weekly_change, cumul = TRUE)
cp_preds_dlnm = crosspred(cb_dlnm_change, fit_dlnm, at = tl_predvalues_change, cumul = TRUE)

calc_ci = function(estimate, se, dir = NULL){
  qn = qnorm(0.95)  
  if(dir == "high"){
    ci = estimate + qn*se
  } else if (dir == "low"){
    ci = estimate - qn*se
  }
}
cumeff_preds_ci_high_acwr = calc_ci(cp_preds_acwr$allfit, cp_preds_acwr$allse, "high")
cumeff_preds_ci_low_acwr = calc_ci(cp_preds_acwr$allfit, cp_preds_acwr$allse, "low")
cumeff_preds_ci_high_weekly_change = calc_ci(cp_preds_weekly_change$allfit, cp_preds_weekly_change$allse, "high")
cumeff_preds_ci_low_weekly_change = calc_ci(cp_preds_weekly_change$allfit, cp_preds_weekly_change$allse, "low")
cumeff_preds_ci_high_dlnm = calc_ci(cp_preds_dlnm$allfit, cp_preds_dlnm$allse, "high")
cumeff_preds_ci_low_dlnm = calc_ci(cp_preds_dlnm$allfit, cp_preds_dlnm$allse, "low")

d_cumulative_preds_acwr = bind_cols(t_load = tl_predvalues_acwr, 
                                    cumulative_effect = cp_preds_acwr$allfit) %>% 
                                    mutate(method = "ACWR",
                                           ci_high = cumeff_preds_ci_high_acwr,
                                           ci_low = cumeff_preds_ci_low_acwr) 

d_cumulative_preds_weekly_change = bind_cols(t_load = tl_predvalues_weekly_change, 
                                             cumulative_effect = cp_preds_weekly_change$allfit) %>% 
                                             mutate(method = "Week-to-week %-change",
                                                    ci_high = cumeff_preds_ci_high_weekly_change,
                                                    ci_low = cumeff_preds_ci_low_weekly_change) 

d_cumulative_preds = bind_cols(t_load = tl_predvalues_change, 
                               cumulative_effect = cp_preds_dlnm$allfit) %>% 
                               mutate(method = "DLNM %-change",
                                      ci_high = cumeff_preds_ci_high_dlnm,
                                      ci_low = cumeff_preds_ci_low_dlnm)

d_cumulative_preds = bind_rows(d_cumulative_preds_acwr, d_cumulative_preds_weekly_change, d_cumulative_preds)

d_cumulative_preds %>% View()

text_size = 14
ggplot(d_cumulative_preds, aes(x = t_load)) +
  facet_wrap(~method, ncol = 3, scales = "free") +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) + 
  geom_line(aes(y = cumulative_effect, color = "Estimation"), size = 0.5) +
  scale_y_continuous(limit = c(-15, 20), breaks = scales::breaks_width(5))  +
  scale_color_manual(values = c(nih_distinct[1], "Black")) + 
  ylab("Cumulative Hazard") +
  xlab("Relative change in sRPE (AU) on Day 0") +
  theme_line(text_size, legend = TRUE) +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = nih_distinct[4]),
        strip.background = element_blank(),
        strip.text.x = element_text(size = text_size, family="Trebuchet MS", colour="black", face = "bold", hjust = -0.01),
        axis.ticks = element_line(color = nih_distinct[4]),
        legend.position = "bottom")

d_cumul_preds = d_cumulative_preds %>% left_join(cis_high, by = c("t_load", "model")) %>% left_join(cis_low, by = c("t_load", "model")) 
# true relationship 
true_effect = calc_coefs(tl_predvalues_change, lag_seq, flindecay)
truecumcoefs = rowSums(true_effect)
d_true_coefs = truecumcoefs %>% enframe() %>% rename(t_load = name, truerel = value) %>% mutate(t_load = as.numeric(t_load))



text_size = 14
ggplot(d_cumulative_preds, aes(x = t_load)) +
  facet_wrap(~method, ncol = 3, scales = "free") +
  geom_line(data = d_true_coefs, aes(y = truerel, color = "True relationship"), size = 0.5) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) + 
  geom_line(aes(y = cumulative_effect, color = "Estimation"), size = 0.5) +
  scale_y_continuous(limit = c(-15, 20), breaks = scales::breaks_width(5))  +
  scale_color_manual(values = c(nih_distinct[1], "Black")) + 
  ylab("Cumulative Hazard") +
  xlab("Relative change in sRPE (AU) on Day 0") +
  theme_line(text_size, legend = TRUE) +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = nih_distinct[4]),
        strip.background = element_blank(),
        strip.text.x = element_text(size = text_size, family="Trebuchet MS", colour="black", face = "bold", hjust = -0.01),
        axis.ticks = element_line(color = nih_distinct[4]),
        legend.position = "bottom")

######################################################## 3D GRAPHS OF PREDICTED VALUES FOR ASSESSING MODEL FIT


cp_preds_dlnm = crosspred(cb_dlnm, fit_dlnm, at = tl_predvalues, cen = 600, cumul = TRUE)

# j decay
persp(x = tl_predvalues, y = lag_seq, l_coefs[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, cp_preds_dlnm$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")


