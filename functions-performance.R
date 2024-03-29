
###################################Performance measure functions###########################################
# this script has a number of functions used in other scripts

# calculate confidence intervals
calc_ci = function(estimate, se, ci_direction){
  qn = qnorm(0.95)  
  if(ci_direction == "high"){
    ci = estimate + qn*se
  } else if (ci_direction == "low"){
    ci = estimate - qn*se   
  }
  ci
}

# calculate coverage of 95% confidence intervals
coverage = function(ci_high, ci_low, target){
  coverage = target <= ci_high & target >= ci_low
  numerator = sum(coverage == TRUE)
  denominator = length(coverage)
  coverage_prop = numerator/denominator
  coverage_prop
}

# calculate the average width of confidence intervals
average_width = function(ci_high, ci_low){
  aw = mean(ci_high-ci_low)
  aw
}

# the absolute bias
raw_bias = function(estimate, target){
  mean(abs(estimate - target))
}

# Calculate the root-mean-squared-error compared to the true cumulative effects
rmse = function(estimate, target){
  sqrt(mean((estimate - target)^2)) 
}

# calculate the root-mean-squared-error of predicted vs. real values (for change in load)
rmse_residuals = function(residuals){
  sqrt(mean(residuals^2)) 
}


# Monte Carlo Standard Error functions for the different performance measures
# monte carlo standard error also requires the number of simulations (runs, permutations)
mcse_bias = function(estimate, target, nsim){
  
  d_se = bind_cols(data.frame(estimate), data.frame(target)) 
  d_est = data.frame(numeric(nrow(d_se)))
  colnames(d_est) = "bias_j"
  for(i in 1:nrow(d_se)){
    d_temp = d_se[-i,]
    rb = raw_bias(d_temp$estimate, d_temp$target)
    d_est[i,1] = rb
  }
  
  bias_j = d_est$bias_j
  main_rb = raw_bias(estimate, target)
  mcse = sqrt(sum((bias_j-main_rb)^2)/(nsim*(nsim-1)))
  mcse
}

mcse_rmse = function(estimate, target, nsim){
  
  d_se = bind_cols(data.frame(estimate), data.frame(target)) 
  d_est = data.frame(numeric(nrow(d_se)))
  colnames(d_est) = "rmse_j"
  for(i in 1:nrow(d_se)){
    d_temp = d_se[-i,]
    rmse = rmse(d_temp$estimate, d_temp$target)
    d_est[i,1] = rmse
  }
  
  rmse_j = d_est$rmse_j
  main_rmse = rmse(estimate, target)
  mcse = sqrt(sum((rmse_j-main_rmse)^2)/(nsim*(nsim-1)))
  mcse
}

# the denominator is the number of CI values
mcse_coverage = function(ci_low, ci_high, target, denominator, nsim){
  is_covered = ifelse((ci_low < target) & (target < ci_high), 1, 0)
  cr = 100*(sum(is_covered == 1, na.rm = TRUE)/denominator)
  mcse = sqrt(abs(((95-cr)*(5-cr)))/nsim)
  mcse
}
