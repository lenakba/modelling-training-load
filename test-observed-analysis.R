# Packages
library(tidyverse) # for everything
library(mice) # for pooling an imputed model
library(dlnm) # distributed lag nonlinear models
library(survival)

# questions:
# 1 Only load was imputed. I did not impute missing injuries. 
# 2 Time period for injuries is from 1-4 days after training load ACWR. 
# If there are more than one injuries, they are still counted as 1 injury (since we are using logit link). Right?

# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors:
options(scipen = 30, 
        stringsAsFactors = FALSE)

symmetrized_change = function(x, y){
  100*((x - y)/(x + y))
}

# loading imputed datasets
folder_handball = paste0("O:\\Prosjekter\\Bache-Mathiesen-Biostatistikk\\Data\\handball\\") # location of handball data
l_handball = readRDS(paste0(folder_handball,"handball_imputed.RDS"))
l_handball = l_handball %>% map(. %>% dplyr::select(., p_id, date_training, load, age, sex, injury))
l_index = l_handball %>% map(. %>% distinct(p_id) %>% rownames_to_column)
l_handball = map2(.x = l_handball,
     .y = l_index,
     ~left_join(.x, .y, by = "p_id")) %>% map(. %>% select(-p_id) %>% rename(p_id = rowname))

d_confounders = l_handball[[1]] %>% distinct(p_id, age, sex) %>% mutate(sex = factor(sex))

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

# add number of days
l_handball = l_handball %>% map(. %>% arrange(p_id, date_training) %>% group_by(p_id) %>% mutate(day = 1:n()))
l_handball = l_handball %>% map(. %>% group_by(p_id) %>% mutate(Fup = ifelse(injury == 1, day, NA)) %>% fill(Fup, .direction = "up") %>% ungroup())

########################################## with frailty for recurrent events
l_surv = l_handball %>% map(. %>% group_by(p_id) %>% 
                                        rename(Stop = day, Id = p_id, Event = injury) %>% 
                                        mutate(Start = lag(Stop),
                                               Start = ifelse(is.na(Start), 0, Start)) %>% ungroup())

# rearrange to counting process form to calculate the Q matrix
l_surv_cpform = l_surv %>% map(~counting_process_form(.) %>% mutate(id = as.numeric(id)))


# arrange the exposure history in wide format in a matrix
l_tl_hist = l_surv %>% map(. %>% group_by(Id) %>% ungroup() %>% select(Id, load, Stop))

# arrange the exposure history in wide format in a matrix
l_tl_hist_spread_day = 
  l_tl_hist %>% map(. %>% pivot_wider(names_from = Stop, values_from = load) %>% select(-Id) %>% as.matrix)

# find the min and max lag
# lag set at 4 weeks (28) as is often used in tl studies
# since the first day is day 0, the 28th day is day 27
lag_min = 0
lag_max = 27

# calc Q matrices
l_q_mat = map2(.x = l_surv_cpform,
               .y = l_tl_hist_spread_day, 
               ~calc_q_matrix(.x, .y))

# make the crossbasis
l_cb_dlnm = l_q_mat %>% map(~crossbasis(., lag=c(lag_min, lag_max), 
                                        argvar = list(fun="ns", knots = 3),
                                        arglag = list(fun="ns", knots = 3)))

# add confounders back to datasets
#l_surv_cpform = l_surv_cpform %>% map(. %>% mutate(id = as.character(id)) %>% left_join(d_confounders, by = c("id" = "p_id")))

# if we do not use the whole data, we only have 12 injuries
# fit DLNM
l_fit_dlnm = map2(.x = l_surv_cpform,
                  .y = l_cb_dlnm,
                  ~coxph(Surv(enter, exit, event) ~ .y, .x, y = FALSE, ties = "efron"))

l_fit_dlnm %>% pool()


library(parfm)
# b) We then fit a model with Weibull baseline and gamma frailty:
fit_frailty = parfm(Surv(time,status)~factor(trt),cluster="id", frailty="gamma", data=eye)
exp(-0.960)






########################################## Time to first injury below

# calc follow up time for those with and without an injury event
l_daystoinjury = l_handball_nomissing %>% map(. %>% filter(injury == 1) %>% group_by(p_id) %>% distinct(injury, .keep_all = TRUE) %>% rename(Fup = day))
l_daystocensor = l_handball_nomissing %>% map(. %>% group_by(p_id) %>% mutate(n_inj = sum(injury)) %>% filter(n_inj == 0) %>% mutate(Fup = max(day)) %>% distinct(Fup, .keep_all = TRUE) %>% select(-n_inj, -day))

l_fup = map2(.x = l_daystoinjury,
     .y = l_daystocensor,
     ~bind_rows(.x, .y))

l_fup = l_fup %>% map(. %>% select(p_id, Fup))

# format data with start and stop times
l_surv_lim = map2(.x = l_handball_nomissing,
              .y = l_fup,
              ~left_join(.x, .y, by = "p_id"))

l_surv = l_surv_lim %>% map(. %>% group_by(p_id) %>% 
                          filter(day <= Fup) %>% 
                          rename(Stop = day, Id = p_id, Event = injury) %>% 
                          mutate(Start = lag(Stop),
                          Start = ifelse(is.na(Start), 0, Start)) %>% ungroup())

# rearrange to counting process form to calculate the Q matrix
l_surv_cpform = l_surv %>% map(~counting_process_form(.) %>% mutate(id = as.numeric(id)))

# arrange the exposure history  in wide format in a matrix
l_tl_hist = l_surv_lim %>% map(. %>% group_by(p_id) %>% 
                      filter(day <= Fup) %>% ungroup() %>% select(p_id, load, day))

# arrange the exposure history in wide format in a matrix
l_tl_hist_spread_day = 
  l_tl_hist %>% map(. %>% pivot_wider(names_from = day, values_from = load) %>% select(-p_id) %>% as.matrix)

# find the min and max lag
# lag set at 4 weeks (28) as is often used in tl studies
# since the first day is day 0, the 28th day is day 27
lag_min = 0
lag_max = 27

# calc Q matrices
l_q_mat = map2(.x = l_surv_cpform,
             .y = l_tl_hist_spread_day, 
             ~calc_q_matrix(.x, .y))

# make the crossbasis
l_cb_dlnm = l_q_mat %>% map(~crossbasis(., lag=c(lag_min, lag_max), 
                     argvar = list(fun="poly", degree = 2),
                     arglag = list(fun="ns", knots = 3)))

# if we do not use the whole data, we only have 12 injuries
# fit DLNM
l_fit_dlnm = map2(.x = l_surv_cpform,
                  .y = l_cb_dlnm,
                  ~coxph(Surv(enter, exit, event) ~ .y, .x, y = FALSE, ties = "efron"))

l_fit_dlnm %>% pool()



# Function for calculating rolling averages on a chooseable number of days
# Based on rollapplyr, not rollmean, as rollmean will only start calculating averages
# at n values, while rollapplyr allows the user to decide preliminary values.
ra = function(x, n_days = lag_max+1, window = TRUE, ...){
  zoo::rollapplyr(x, n_days, mean, partial = TRUE)
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

# function to nest the exposure history data by each individual, 
# and run a user-specified function on each of their datasets in the list
function_on_list = function(d_sim_hist, FUN = NULL, day_start){
  nested_list = d_sim_hist %>% group_by(p_id) %>% nest()
  nested_list$data = nested_list$data %>% map(., ~FUN(.$load))
  l_unnest = unnest(nested_list, cols = data) %>% ungroup() %>% 
    filter(!is.na(data))
  l_unnest
}

# rolling average
l_hist_ra = l_tl_hist %>% map(. %>% function_on_list(., ra, lag_max+1) %>% rename(ra_t_load = data, id = p_id) %>% group_by(id) %>% mutate(day = 1:n()) %>% ungroup)
l_surv_cpform = l_surv_cpform %>% map(. %>% mutate(id = as.character(id)))
l_surv_cpform_mods = map2(.x = l_hist_ra,
     .y = l_surv_cpform,
     ~left_join(.y, .x, by = c("id", "exit" = "day")))

l_ob_ra = l_surv_cpform_mods %>% map(~onebasis(.$ra_t_load, "poly", degree = 2))                                    
l_fit_ra = map2(.x = l_surv_cpform_mods,
                  .y = l_ob_ra,
                  ~coxph(Surv(enter, exit, event) ~ .y, .x, y = FALSE, ties = "efron"))

# redi
lag_seq = lag_min:lag_max
d_sim_hist_redi = l_tl_hist %>% map(. %>% function_on_list(., redi, lag_max+1) %>% rename(redi_t_load = data))


l_tl_hist[[1]]

function_on_list(l_tl_hist[[1]], redi, lag_max+1)


d_sim_hist_ewma = l_tl_hist %>% map(. %>% function_on_list(., ewma, lag_max+1) %>% rename(ewma_t_load = data))


d_survival_sim_cpform_mods_a = d_survival_sim_cpform_mods_a %>% left_join(d_sim_hist_ra, by = c("id", "exit" = "day"))
d_survival_sim_cpform_mods_a = d_survival_sim_cpform_mods_a %>% left_join(d_sim_hist_ewma, by = c("id", "exit" = "day"))
d_survival_sim_cpform_mods_a = d_survival_sim_cpform_mods_a %>% left_join(d_sim_hist_redi, by = c("id", "exit" = "day"))

ob_ewma = onebasis(d_survival_sim_cpform_mods_a$ewma_t_load, "poly", degree = 2)
ob_redi = onebasis(d_survival_sim_cpform_mods_a$redi_t_load, "poly", degree = 2)


# make sure sex is treated as categorical
# we also create it as a dummy variable 0 = female, 1 = male
l_handball = l_handball %>% map(. %>% mutate(sex_fac = factor(sex)))

d_survival_sim_cpform = counting_process_form(d_survival_sim)
q_mat = calc_q_matrix(d_survival_sim_cpform, d_sim_tl_hist_spread_day)




l_fit_dlnm %>% map(~AIC(.))
l_fit_ra %>% map(~AIC(.))


####################--------------------------------Frailty?--------------------------------------------

l_surv = l_handball_nomissing %>% map(. %>% group_by(p_id) %>% 
                               rename(Stop = day, Id = p_id, Event = injury) %>% 
                               mutate(Start = lag(Stop),
                                      Start = ifelse(is.na(Start), 0, Start),
                                      Fup = max(Stop)) %>% ungroup())

# rearrange to counting process form to calculate the Q matrix
l_surv_cpform = l_surv %>% map(~counting_process_form(.) %>% mutate(id = as.numeric(id)))

# arrange the exposure history  in wide format in a matrix
l_tl_hist = l_surv %>% map(. %>% group_by(Id) %>% ungroup() %>% select(Id, load, Stop))

# arrange the exposure history in wide format in a matrix
l_tl_hist_spread_day = 
  l_tl_hist %>% map(. %>% pivot_wider(names_from = Stop, values_from = load) %>% select(-Id) %>% as.matrix)

# find the min and max lag
# lag set at 4 weeks (28) as is often used in tl studies
# since the first day is day 0, the 28th day is day 27
lag_min = 0
lag_max = 27

# calc Q matrices
l_q_mat = map2(.x = l_surv_cpform,
               .y = l_tl_hist_spread_day, 
               ~calc_q_matrix(.x, .y))

# make the crossbasis
l_cb_dlnm = l_q_mat %>% map(~crossbasis(., lag=c(lag_min, lag_max), 
                                        argvar = list(fun="poly", degree = 2),
                                        arglag = list(fun="ns", knots = 3)))

# if we do not use the whole data, we only have 12 injuries
# fit DLNM
l_fit_dlnm = map2(.x = l_surv_cpform,
                  .y = l_cb_dlnm,
                  ~coxph(Surv(enter, exit, event) ~ .y, .x, y = FALSE, ties = "efron"))

l_fit_dlnm %>% pool()


library(parfm)
# b) We then fit a model with Weibull baseline and gamma frailty:
fit_frailty = parfm(Surv(time,status)~factor(trt),cluster="id", frailty="gamma", data=eye)
exp(-0.960)

# fit models
fit_ra = coxph(Surv(enter, exit, event) ~ ob_ra, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
fit_ewma = coxph(Surv(enter, exit, event) ~ ob_ewma, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")
fit_redi = coxph(Surv(enter, exit, event) ~ ob_redi, d_survival_sim_cpform_mods, y = FALSE, ties = "efron")

list_fits = list(fit_ra, fit_ewma, fit_redi, fit_dlnm)

# run predictions
cp_preds_ra = crosspred(ob_ra, fit_ra, at = predvalues, cen = 600, cumul = TRUE)
cp_preds_ewma = crosspred(ob_ewma, fit_ewma, at = predvalues, cen = 600, cumul = TRUE)
cp_preds_redi = crosspred(ob_redi, fit_redi, at = predvalues, cen = 600, cumul = TRUE)
cp_preds_dlnm = crosspred(cb_dlnm, fit_dlnm, at = predvalues, cen = 600, cumul = TRUE)