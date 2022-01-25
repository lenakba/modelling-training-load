
library(tidyverse) # for datawrangling
library(lmisc) # for OSTRC colors etc.

# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors:
options(scipen = 40, 
        stringsAsFactors = FALSE)

# assume the training load values are in the same folder
# as the r script, and that work repository is in the source file location
d_load = read_delim("norwegian_premier_league_football_td_vec.csv", delim = ";")

# fetching functions for estimating the true relationship
# we assume working directory is the same location as this script
source("functions-relationships.R", encoding = "UTF-8")

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
tl_predvalues = seq(min(tl_observed), max(tl_observed), 25)
tl_predvalues_change = seq(min(tl_observed_change, na.rm = TRUE), max(tl_observed_change, na.rm = TRUE), 10)

# list of the exposure-response (f) and lag-response (w) combination functions
fw_funs = list(fjconst, fjdecay, fjexponential_decay, flindirection_flip)
fw_funs_change = list(flinconst, flindecay, flinexponential_decay)

# calculate true coefs
l_coefs = fw_funs %>% map(., ~calc_coefs(tl_predvalues, lag_seq, .)) %>% map(., ~exp(.))
l_coefs_change = fw_funs_change %>% map(., ~calc_coefs(tl_predvalues_change, lag_seq, .)) %>% map(., ~exp(.))

# figures of f*w functions for amount of training load
png("figure1_relationships_amount.png", units = "in", width = 10, height = 10, res = 600)
par(mfrow=c(2,2), mar=c(0.8,0.1,1,0.1), mgp = c(4, 1, 0))
persp(x = tl_predvalues, y = lag_seq, l_coefs[[1]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE (AU)", main = "A Constant")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE(AU)", main = "B Decay")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[3]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE (AU)", main = "C Exponential Decay")

persp(x = tl_predvalues, y = lag_seq, l_coefs[[4]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE (AU)", main = "D Direct, then inverse")
dev.off()

png("figure2_relationships_change.png", units = "in", width = 10, height = 10, res = 600)
par(mfrow=c(2,2), mar=c(0.8,0.1,1,0.1), mgp = c(4, 1, 0))
# figures for change in load
persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[1]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)", main = "A Constant")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[2]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)", main = "B Decay")

persp(x = tl_predvalues_change, y = lag_seq, l_coefs_change[[3]], ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5, cex.axis=1.0, cex.lab=1.0,
      border=grey(0.2), col = nih_distinct[1], xlab = "%ΔsRPE (AU)", main = "C Exponential Decay")
dev.off()


##################################################### figures of just the f(x) functions

library(devEMF)
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

x = c(tl_predvalues, tl_predvalues, tl_predvalues_change)
y = c(exp(fj(tl_predvalues)), exp(flin_amount(tl_predvalues)), exp(flin(tl_predvalues_change)))
d_fx_coefs = tibble(tl = x, coef = y) %>% 
             mutate(index = 1:n(),
                    group = case_when(index <= length(tl_predvalues) ~ "A sRPE (AU) J shaped", 
                                      index > length(tl_predvalues) & index <= length(tl_predvalues)*2 ~ "B sRPE (AU) Linear", 
                                       TRUE ~ "C %ΔsRPE (AU) Linear"),
                    group = fct_inorder(group))



fplot = ggplot(d_fx_coefs, aes(x = tl, y = coef)) +
  facet_wrap(~group, scales = "free", ncol = 2) +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1.2, color = nih_distinct[1]) +
  geom_line(color = nih_distinct[2], size = 0.8) +
  theme_line(text_size) +
  ylab("Hazard Ratio") +
  xlab("Training load") +
  ostrc_theme +
  scale_y_continuous(limits = c(0.5, 3.0), breaks = scales::breaks_width(0.5, 0))

emf("supfig_f_funksjon_coefs.emf", width = 8, height = 7)
fplot
dev.off()

##################################################### figures of just the w(l) functions

# make list of lag functions
wfun_list = list(wconst, wdecay, wexponential_decay, wdirection_flip)
names_lag = list("A Constant", "B Decay", "C Exponential Decay", "D Direct, then inverse")

# calculate true coefs
predvalue_list = list(lag_seq, lag_seq, lag_seq, lag_seq)
l_coefs_lags = predvalue_list %>% map2(.x = .,
                                       .y = wfun_list,
                                        ~exp(.y(.x)))
# add labels etc.
l_coefs_lags = l_coefs_lags %>% 
  map(~enframe(., name = NULL, value = "coef")) %>% 
  map(. %>% mutate(lag = lag_seq)) %>% 
  map2(.x = .,
       .y = names_lag,
       ~mutate(., group = .y))

# collapse to dataset
d_coefs_lags = bind_rows(l_coefs_lags)

wplot = ggplot(d_coefs_lags, aes(x = lag, y = coef)) +
  facet_wrap(~group, ncol = 2, scales = "free") +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1.2, color = nih_distinct[1]) +
  geom_line(color = nih_distinct[2], size = 0.8) +
  theme_line(text_size) +
  ylab("Hazard Ratio") +
  xlab("Lag (Days)") +
  scale_x_continuous(limits = c(0, 27), breaks = scales::breaks_width(3, 0)) +
  ostrc_theme +
  coord_cartesian(ylim = c(0, 3.0))



emf("supfig_w_funksjon_coefs.emf", width = 10, height = 6)
wplot
dev.off()
