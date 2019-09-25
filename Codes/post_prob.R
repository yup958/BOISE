## Func 2 for BOISE: given xA, certain block of original data xb, calculate posterior probability p(xA|x0, class = b).
## Input: original data block xb, number of targets Nb, informer data xA, informer index A
## Output: posterior probability calculated using MCMC. As p(xA|xb) = \int p(xA|theta)*p(theta|xb) dtheta

post_prob <- function(xb,Nb,xA,A,nA){
  source("prob.R")
  if(Nb == 1){
    tmp <- list(J = Nb, K = length(xb), x = xb)
    fit <- stan(file = "Rasch_Vector.stan",data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
  } else{
    tmp <- list(J = Nb, K = dim(xb)[2], x = xb)
    fit <- stan(file = "Rasch_Matrix.stan",data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
  }
  la <- extract(fit, permuted = TRUE) # return a list of arrays 
  alphahat <- la$alpha
  if(nA >= 2){
    btahat <- la$bta[ ,A]
    postprob = 0
    for (i in 1:length(alphahat)) {
      p = exp(alphahat[i] + btahat[i,]) / (1 + exp(alphahat[i] + btahat[i,]))
      loglike = sum(xA * (alphahat[i] + btahat[i,])) + sum(log(1 - p))
      postprob = postprob + exp(loglike)
    }
  } else{
    btahat <- la$bta[ ,A]
    postprob = 0
    for (i in 1:length(alphahat)) {
      p = exp(alphahat[i] + btahat[i]) / (1 + exp(alphahat[i] + btahat[i]))
      loglike = sum(xA * (alphahat[i] + btahat[i])) + sum(log(1 - p))
      postprob = postprob + exp(loglike)
    }
  }
  postprob = postprob / length(alphahat)
  return(postprob)
}

#Test
# xb <- dat[1,]
# Nb <- 1
# xA <- dat[100, c(1,10,100,200)]
# A <- c(1,10,100,200)
# p <- post_prob(xb,Nb,xA,A)
