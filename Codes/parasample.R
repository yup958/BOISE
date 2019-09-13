## Function 2 for Auxillary Parameters Sampling
## Input: Number of Auxillarys m; Number of Pars n(Prior information if needed in the future)
## Output: A n * m matrix with each row a i.i.d sample of pars
parasample <- function(aux, dimen, mu1, sigma1, mu2, sigma2){
  alpha = rnorm(aux, mean = mu1, sd = sigma1)
  beta = mvtnorm::rmvnorm(aux, mean = rep(mu2, dimen - 1), sigma = diag(sigma2,dimen - 1, dimen -1))
  pars = rep(0, aux * dimen)
  dim(pars) = c(aux, dimen)
  pars[ ,1] = alpha
  pars[ ,2:dimen] = beta
  return(pars)
}
