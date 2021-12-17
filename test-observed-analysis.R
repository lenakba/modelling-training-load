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

# find the min and max lag
# lag set at 4 weeks (28) as is often used in tl studies
# since the first day is day 0, the 28th day is day 27
lag_min = 0
lag_max = 27

# loading imputed datasets
folder_handball = paste0("O:\\Prosjekter\\Bache-Mathiesen-Biostatistikk\\Data\\handball\\") # location of handball data
l_handball = readRDS(paste0(folder_handball,"handball_imputed.RDS"))
l_handball = l_handball %>% map(. %>% dplyr::select(., p_id, date_training, load, age, sex, injury))
l_index = l_handball %>% map(. %>% distinct(p_id) %>% rownames_to_column)
l_handball = map2(.x = l_handball,
     .y = l_index,
     ~left_join(.x, .y, by = "p_id")) %>% map(. %>% select(-p_id) %>% rename(p_id = rowname))

d_confounders = l_handball[[1]] %>% distinct(p_id, .keep_all = TRUE) %>% select(p_id, sex, age) %>% mutate(sex = factor(sex))


l_handball[[1]] %>% count(injury)

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
l_handball = l_handball %>% map(. %>% arrange(p_id, date_training) %>% 
                                  group_by(p_id) %>% 
                                  mutate(day = 1:n()))
l_handball = l_handball %>% map(. %>% group_by(p_id) %>%
                                  mutate(Fup = ifelse(injury == 1, day, NA)) %>% 
                                  fill(Fup, .direction = "up") %>% 
                                  ungroup())

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



# calc Q matrices
l_q_mat = map2(.x = l_surv_cpform,
               .y = l_tl_hist_spread_day, 
               ~calc_q_matrix(.x, .y))

# add confounders back to datasets
l_surv_cpform = l_surv_cpform %>% map(. %>% mutate(id = as.character(id)) %>% as_tibble() %>%
                                        left_join(d_confounders, by = c("id" = "p_id")))

# make the crossbasis
l_cb_dlnm = l_q_mat %>% map(~crossbasis(., lag=c(lag_min, lag_max), 
                                        argvar = list(fun="ns", knots = 3),
                                        arglag = list(fun="lin")))
# fit DLNM
l_fit_dlnm = map2(.x = l_surv_cpform,
                  .y = l_cb_dlnm,
                  ~coxph(Surv(enter, exit, event) ~ .y + sex + age + 
                           frailty.gaussian(id), .x, y = FALSE, ties = "efron"))

d_pooled = summary(l_fit_dlnm %>% pool(), conf.int = TRUE, exponentiate = TRUE) %>% as_tibble() %>% mutate_if(is.numeric, ~round(.,3))

write_excel_csv(d_pooled, "sup_tableS4_modelcoefs.csv", delim = ";", na = "")

#################################################### Figures for exploring effects

library(lmisc)
# shared figure options
text_size = 14
ostrc_theme =  theme(panel.border = element_blank(), 
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.line = element_line(color = nih_distinct[4]),
                     strip.background = element_blank(),
                     strip.text.x = element_text(size = text_size, family="Trebuchet MS", colour="black", face = "bold", hjust = -0.01),
                     axis.ticks = element_line(color = nih_distinct[4]),
                     legend.position = "bottom")

# vector of tl values used in visualizations of predictions
predvalues = seq(min(l_handball[[1]]$load), max(l_handball[[1]]$load), 100)
lag_seq = lag_min:lag_max 

# predict hazards
l_cp_preds_dlnm = 
  map2(.x = l_fit_dlnm,
       .y = l_cb_dlnm,
       ~crosspred(.y, .x, at = predvalues, cen = 0, cumul = TRUE))
glimpse(l_cp_preds_dlnm)
# function for plucking the right matrix out of the crosspred list within the list of crosspred lists
pluck_mat = function(x, pos){pluck(l_cp_preds_dlnm, x, pos)}
# 13 is matRRfit
allRRfit = 16
d_preds_cumul1 = pluck_mat(1, allRRfit)
d_preds_cumul2 = pluck_mat(2, allRRfit)
d_preds_cumul3 = pluck_mat(3, allRRfit)
d_preds_cumul4 = pluck_mat(4, allRRfit)
d_preds_cumul5 = pluck_mat(5, allRRfit)
l_cumulRRfit = list(d_preds_cumul1, d_preds_cumul2, d_preds_cumul3, d_preds_cumul4, d_preds_cumul5)
# average across preds
mat_cumulRRfit = reduce(l_cumulRRfit, `+`) / length(l_cumulRRfit)

# conflow
allRRfit_low = 17
d_preds_cumullow1 = pluck_mat(1, allRRfit_low)
d_preds_cumullow2 = pluck_mat(2, allRRfit_low)
d_preds_cumullow3 = pluck_mat(3, allRRfit_low)
d_preds_cumullow4 = pluck_mat(4, allRRfit_low)
d_preds_cumullow5 = pluck_mat(5, allRRfit_low)
l_cumulRRfit_low = list(d_preds_cumullow1, d_preds_cumullow2, d_preds_cumullow3, d_preds_cumullow4, d_preds_cumullow5)
# average across preds
mat_cumulRRfit_low = reduce(l_cumulRRfit_low, `+`) / length(l_cumulRRfit_low)

# confhigh
allRRfit_high = 18
d_preds_cumulhigh1 = pluck_mat(1, allRRfit_high)
d_preds_cumulhigh2 = pluck_mat(2, allRRfit_high)
d_preds_cumulhigh3 = pluck_mat(3, allRRfit_high)
d_preds_cumulhigh4 = pluck_mat(4, allRRfit_high)
d_preds_cumulhigh5 = pluck_mat(5, allRRfit_high)
l_cumulRRfit_cumulhigh = list(d_preds_cumulhigh1, d_preds_cumulhigh2, d_preds_cumulhigh3, d_preds_cumulhigh4, d_preds_cumulhigh5)
# average across preds
mat_cumulRRfit_high = reduce(l_cumulRRfit_cumulhigh, `+`) / length(l_cumulRRfit_cumulhigh)

d_cumul = as_tibble(mat_cumulRRfit) %>% 
          mutate(srpe = predvalues, ci_low = mat_cumulRRfit_low, ci_high = mat_cumulRRfit_high)
plot_cumul = ggplot(d_cumul, aes(x = srpe, y = value, group = 1)) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1) +
  geom_line(size = 0.75, color = nih_distinct[4]) +
  theme_base(text_size) +
  ostrc_theme +
  xlab("sRPE (AU)") +
  ylab("Cumulative HR on Day 0")  +
  coord_cartesian(ylim = c(0.2, 6), expand = FALSE)

# 13 is matRRfit
matRRfit = 13
d_preds1 = pluck_mat(1, matRRfit)
d_preds2 = pluck_mat(2, matRRfit)
d_preds3 = pluck_mat(3, matRRfit)
d_preds4 = pluck_mat(4, matRRfit)
d_preds5 = pluck_mat(5, matRRfit)
l_matRRfit = list(d_preds1, d_preds2, d_preds3, d_preds4, d_preds5)
# average across preds
mat_matRRfit = reduce(l_matRRfit, `+`) / length(l_matRRfit)

# conflow
matRRfit_low = 14
d_preds_low1 = pluck_mat(1, matRRfit_low)
d_preds_low2 = pluck_mat(2, matRRfit_low)
d_preds_low3 = pluck_mat(3, matRRfit_low)
d_preds_low4 = pluck_mat(4, matRRfit_low)
d_preds_low5 = pluck_mat(5, matRRfit_low)
l_matRRfit_low = list(d_preds_low1, d_preds_low2, d_preds_low3, d_preds_low4, d_preds_low5)
# average across preds
mat_matRRfit_low = reduce(l_matRRfit_low, `+`) / length(l_matRRfit_low)

# confhigh
matRRfit_high = 15
d_preds_high1 = pluck_mat(1, matRRfit_high)
d_preds_high2 = pluck_mat(2, matRRfit_high)
d_preds_high3 = pluck_mat(3, matRRfit_high)
d_preds_high4 = pluck_mat(4, matRRfit_high)
d_preds_high5 = pluck_mat(5, matRRfit_high)
l_matRRfit_high = list(d_preds_high1, d_preds_high2, d_preds_high3, d_preds_high4, d_preds_high5)
# average across preds
mat_matRRfit_high = reduce(l_matRRfit_high, `+`) / length(l_matRRfit_high)

# lag-response curve for sRPE 3000
spre_fixed = "100"
rownumber = which(rownames(mat_matRRfit)==spre_fixed)
d_preds_per_lag = as_tibble(mat_matRRfit[rownumber,]) %>% 
  rename(coef = value) %>% 
  mutate(lag = 0:27,
         ci_low = mat_matRRfit_low[rownumber,],
         ci_high = mat_matRRfit_high[rownumber,])

png("figure5_part2_3d.png", units = "in", width = 10, height = 5, res = 600)
# figure for centering at average sRPE, 458
persp(x = predvalues, y = lag_seq, mat_matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE", main = "D 3D plane of effects", zlim = c(0.6, 1.6))
dev.off()


cairo_pdf("figure5_part2_3d.pdf", width = 12, height = 7)
# figure for centering at average sRPE, 458
persp(x = predvalues, y = lag_seq, mat_matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.2, cex.lab=1.2,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE", main = "D 3D plane of effects", zlim = c(0.6, 1.6))
dev.off()


########## create figures for lag effects only

# exposure-response curve for lag 0
lag_fixed = "lag0"
colnumber = which(colnames(mat_matRRfit) == lag_fixed)
d_preds_per_srpe = as_tibble(mat_matRRfit[,colnumber]) %>% 
  rename(coef = value) %>% 
  mutate(srpe = predvalues,
         ci_low = mat_matRRfit_low[,colnumber],
         ci_high = mat_matRRfit_high[,colnumber])

plot_dlnm_2d1 = ggplot(d_preds_per_srpe, aes(x = srpe, y = coef, group = 1)) +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) +
  geom_line(size = 0.75, color = nih_distinct[4]) +
  theme_base(text_size) +
  ostrc_theme +
  xlab("sRPE (AU)") +
  ylab("HR on Day 0") +
  coord_cartesian(ylim = c(0.5, 1.5))

lag_fixed = "lag27"
colnumber = which(colnames(mat_matRRfit) == lag_fixed)
d_preds_per_srpe = as_tibble(mat_matRRfit[,colnumber]) %>% 
  rename(coef = value) %>% 
  mutate(srpe = predvalues,
         ci_low = mat_matRRfit_low[,colnumber],
         ci_high = mat_matRRfit_high[,colnumber])

plot_dlnm_2d2 = ggplot(d_preds_per_srpe, aes(x = srpe, y = coef, group = 1)) +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) +
  geom_line(size = 0.75, color = nih_distinct[4]) +
  theme_base(text_size) +
  ostrc_theme +
  xlab("sRPE (AU)") +
  ylab("HR on Day 27") +
  coord_cartesian(ylim = c(0.5, 1.5))

png("figure5_part1_2d.png", units = "in", width = 10, height = 8, res = 600)
ggpubr::ggarrange(plot_dlnm_2d1, plot_dlnm_2d2, plot_cumul, ncol = 2, nrow = 2, labels = c("A Risk on current day", "B Risk on 28th lag day", "C Cumulative effect"))
dev.off()

cairo_pdf("figure5_part1_2d.pdf", width = 10, height = 8)
ggpubr::ggarrange(plot_dlnm_2d1, plot_dlnm_2d2, plot_cumul, ncol = 2, nrow = 2, labels = c("A Risk on current day", "B Risk on 27th day", "C Cumulative effect"))
dev.off()