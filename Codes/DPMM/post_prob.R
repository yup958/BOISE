## Func 2 for BOISE: given xA, certain block of original data xb, calculate posterior probability p(xA|x0, class = b).
## Input: original data block xb, number of targets Nb, informer data xA, informer index A
## Output: posterior probability calculated using MCMC. As p(xA|xb) = \int p(xA|theta)*p(theta|xb) dtheta

post_prob <- function(xb,Nb,xA,A,nA){
  if(Nb == 1){
    tmp <- list(J = Nb, K = length(xb), x = xb)
    fit <- stan(file = "PEL_Rasch_Vector.stan",data = tmp, control = list(max_treedepth = 15), iter = 4000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- la$alpha
  } else{
    tmp <- list(J = Nb, K = dim(xb)[2], x = xb)
    fit <- stan(file = "PEL_Rasch_Matrix.stan",data = tmp, control = list(max_treedepth = 15), iter = 4000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- apply(la$alpha,1,mean)
  }
  
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
# xb <- dat[c(1,2),]
# Nb <- 2
# xA <- dat[100, 1]
# A <- 1
# nA = 1
# p <- post_prob(xb,Nb,xA,A,nA)
