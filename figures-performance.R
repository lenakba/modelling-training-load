
# so we don't have to deal with scientific notations
# and strings aren't automatically read as factors
options(scipen = 17, 
        stringsAsFactors = FALSE)

# loading packages
library(tidyverse) # for datawrangling
library(lmisc) # OSTRC colors

# read data
d_res_amount = read_delim("simulation_results_amount.csv", delim = ";")
d_res_change = read_delim("simulation_results_change.csv", delim = ";") 

# change names for figures
d_res_amount = d_res_amount %>% mutate(method = str_to_upper(method),
                                       method = ifelse(method == "RA", "Rolling Average", method),
                                method_fac = fct_inorder(method))

d_res_change = d_res_change %>% mutate(method = str_to_upper(method),
                        method = ifelse(method == "WEEKLY_CHANGE", "Week-to-week %Δ", method),
                        method_fac = fct_inorder(method))

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

############################################### Predicted values vs. true coefs (amount only) ###############################################

# pred vs. true coefs, 3 main relationships
d_amount_1rep = d_res_amount %>% filter(rep == 5)
d_amount_3rels = d_amount_1rep %>% filter(relationship != "Linear Direction Change")
d_directionflip = d_amount_1rep %>% filter(relationship == "Linear Direction Change")
  
preds_plot = function(d, x, xlab, compare = TRUE){
  x = enquo(x)
  
  plot = ggplot(d, aes(x = !!x))
    if(compare){
      plot = plot + geom_line(data = d, aes(y = true_cumul_coefs, color = "True relationship"), size = 0.5)   
    }
  plot = plot + geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) + 
    geom_line(aes(y = cumul, color = "Estimation"), size = 0.5) +
    scale_color_manual(values = c(nih_distinct[1], "Black")) + 
    ylab("Cumulative Hazard") +
    xlab(xlab) +
    theme_line(text_size, legend = TRUE) +
    ostrc_theme 
  plot
}

preds_plot(d_amount_3rels, t_load, "sRPE (AU) on Day 0")  +
  facet_wrap(c("relationship", "method_fac"), ncol = 4, scales = "free") +
  scale_y_continuous(limit = c(-15, 20), breaks = scales::breaks_width(5)) 

preds_plot(d_directionflip, t_load, "sRPE (AU) on Day 0") + 
  facet_wrap(~method_fac, ncol = 2, scales = "free") 

# for labels?
# rel_labs <- c("A Constant", "B Constant", "C Constant", "D Decay", "E Decay", "F Decay", "G Exponential Decay", "H Exponential Decay", "I Exponential Decay")
# names(rel_labs) <- rep(c("Constant", "Decay", "Exponential Decay"), each = 3)

############################################### Predicted values (change only) ###############################################
d_change_1rep = d_res_change %>% filter(rep == 1)
d_change_long = d_change_1rep %>% pivot_longer(cols = c(t_load_acwr, t_load_weekly_change, t_load_change), values_to = "t_load_relative")
preds_plot(d_change_long, t_load_relative, "Relative change in sRPE (AU) on Day 0", compare = FALSE) +
  facet_wrap(c("relationship", "method_fac"), ncol = 3, scales = "free") +
  theme(legend.position = "none")

############################################## Dotplot ranking figures (RMSE, AIC) ###########################################################################

perf_internal = c("aic", "rmse_residuals", "aw")
perf_external = c("rmse", "coverage")
d_perf_params_amount = d_res_amount %>% group_by(relationship, method) %>% 
                       summarise_at(vars(all_of(perf_internal), all_of(perf_external), starts_with("mcse")), mean)
d_perf_params_change = d_res_change %>% group_by(relationship, method) %>% 
                       summarise_at(vars(all_of(perf_internal)), mean)

d_perf_amount_3rels = d_perf_params_amount %>% filter(relationship != "Linear Direction Change")
d_perf_directionflip = d_perf_params_amount %>% filter(relationship == "Linear Direction Change")

rmse_plot = function(d, x, xlab){
  x = enquo(x)

  d = d %>% mutate(relationship = case_when(relationship == "Constant" ~ "A Constant",
                                            relationship == "Decay" ~ "B Decay",
                                            relationship == "Exponential Decay" ~ "C Exponential Decay"))
  
  ggplot(d, aes(x = !!x, y = method)) +
    geom_point(size = 4, color = nih_distinct[4]) + 
    theme_dot(text_size) + 
    xlab(xlab) +
    ylab(NULL) + 
    theme(axis.line = element_line(color = nih_distinct[4]),
          axis.ticks = element_line(color = nih_distinct[4]),
          panel.border = element_blank(), 
          panel.background = element_blank(),
          title = element_text(face = "bold", family = "Trebuchet MS"),
          strip.text.x = element_text(size = text_size, family="Trebuchet MS", colour="black", face = "bold", hjust = -0.01)) 
}

rmse_plot(d_perf_amount_3rels, rmse, xlab = "External RMSE") + facet_wrap(~relationship, scales = "free")
rmse_plot(d_perf_directionflip, rmse, xlab = "External RMSE")
rmse_plot(d_perf_params_change, rmse_residuals, xlab = "Internal RMSE") + facet_wrap(~relationship, scales = "free")
rmse_plot(d_perf_params_change, aic, xlab = "AIC") + facet_wrap(~relationship, scales = "free")

############################################## DLNM effect at each lag level #################################################################

# lag set at 4 weeks (28) as is often used in tl studies
lag_min = 0
lag_max = 27
lag_seq = lag_min:lag_max # number of days before current day assumed to affect risk of injury
tl_predvalues = d_res_amount %>% filter(rep == 1, method == "EWMA", relationship == "Constant") %>% pull(t_load)
tl_predvalues_change = d_res_lin_constant %>% filter(rep == 1, method == "dlnm") %>% pull(t_load_change)

# time x amount effect (DLNM only)

const_coefs = d_res_amount %>% filter(relationship == "Constant") %>% pull(true_cumul_coefs)

# j decay
persp(x = tl_predvalues, y = lag_seq, exp(const_coefs), ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")

persp(x = tl_predvalues, y = lag_seq, cp_preds_dlnm$matRRfit, ticktype="detailed", theta=230, ltheta=150, phi=40, lphi=30,
      ylab="Lag (Days)", zlab="HR", shade=0.75, r=sqrt(3), d=5,
      border=grey(0.2), col = nih_distinct[1], xlab = "sRPE")
