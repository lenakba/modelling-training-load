
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

# FUNCTIONS TO COMPUTE THE LOGIT AND INVERSE LOGIT TRANSFORMATIONS
# logit = function(prob) log(prob/(1-prob))
# invlogit = function(linpred) exp(linpred)/(1+exp(linpred))

# functions for simulating the effect of the amount of training load on injury risk (J-shape)
# and the effect of amount of change of training load on injury risk (linear)
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


# making  lists with the functions
flist = list(flin, flin, fj, fj)
wlist = list(wconst, wconst, wdecay, wdecay)


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

# arrange in a matrix which will be used later
d_sim_tl_hist_spread_day = 
  d_sim_tl_hist  %>% select(-t_load_change) %>% 
  pivot_wider(names_from = day, values_from = t_load) %>% select(-id) %>% as.matrix

d_sim_tl_hist_spread_day_change = 
  d_sim_tl_hist %>% select(-t_load) %>%
  pivot_wider(names_from = day, values_from = t_load_change) %>% select(-id) %>% as.matrix


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

# calculate cumulative effect of training load for each individual,
# as a cross-basis of a function f or exposure amount, and function w, lag-time
l_load_nested = d_sim_tl_hist %>% nest(data = t_load) 

# can run a for-loop for each element in a list of functions
# for the final simulation
# but we will run 1 at a time during testing of the code.
# for(i in 1:nrow(combsim)){
# l = l_load_nested$data %>% map(~fcumeff_tot(.$t_load, lag = c(lag_min, lag_max), combsim$x_funs[i], combsim$lag_funs[i]))
# }
l_load_nested$data = l_load_nested$data %>% map(~calc_cumeff(.$t_load, fj, wdecay))
d_load_cumeffs_j_decay = unnest(l_load_nested, cols = c(data)) %>% rename(cumeff = data)

# to have 2 examples
l_load_nested = d_sim_tl_hist %>% nest(data = t_load) 
l_load_nested$data = l_load_nested$data %>% map(~calc_cumeff(.$t_load, flin, wdecay))
d_load_cumeffs_lin_decay = unnest(l_load_nested, cols = c(data)) %>% rename(cumeff = data)

# # permalgorithm requires a matrix of covariates
# cumeff_j_decay = as.matrix(d_load_cumeffs_j_decay$cumeff)
# cumeff_lin_decay = as.matrix(d_load_cumeffs_lin_decay$cumeff)

# arrange in a matrix which will be used later
cumeff_j_decay_i = 
  d_load_cumeffs_j_decay %>% mutate(day = rep(1:t_max, nsub)) %>% 
  pivot_wider(names_from = id, values_from = cumeff) %>% select(-day) %>% as.matrix

cumeff_lin_decay_i = 
  d_load_cumeffs_lin_decay %>% mutate(day = rep(1:t_max, nsub)) %>% 
  pivot_wider(names_from = id, values_from = cumeff) %>% select(-day) %>% as.matrix

# simulate survival times and events for each participant
# set the beta-coefficient to be 1, as in 1 times the cumeffect provided in the matrix
# set censoring probability to 0.10
# use cumulative effect as the time-dependent variable
set.seed(1234)
d_sim_surv_j_decay = permalgorithm(nsub, t_max, Xmat = cumeff_j_decay_i,
                      censorRandom = runif(nsub, 1, t_max*2), betas=1)

d_sim_surv_lin_decay = permalgorithm(nsub, t_max, Xmat = cumeff_lin_decay_i,
                                   censorRandom = runif(nsub, 1, t_max*2),betas=1)

d_sim_surv_j_decay %>% summarise(n_events = sum(Event == 1))
d_sim_surv_lin_decay %>% summarise(n_events = sum(Event == 1))


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



