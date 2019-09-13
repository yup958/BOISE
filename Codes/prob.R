## Function 1 for Density Calculation
## Input: x is a 0-1 vector of dim 366, phi is pars of dim 367 (alpha + beta[366])
## Output: Probability of x with parameter phi
prob = function(x, phi){
  alpha = phi[1]
  beta = phi[2: length(phi)]
  p = exp(alpha + beta) / (1 + exp(alpha + beta))
  loglike = sum(x * (alpha + beta)) + sum(log(1 - p))
  return(exp(loglike))
}
