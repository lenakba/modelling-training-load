
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
average_width = function(coef_high, coef_low){
  aw = mean(coef_high-coef_low)
  aw
}

# Calculate the root-mean-squared-error compared to the true cumulative effects
rmse = function(estimate, target){
  sqrt(mean((estimate - target)^2)) 
}