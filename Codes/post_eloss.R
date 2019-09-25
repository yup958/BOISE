## Func 3 for BOISE: Given an informer result xA, certain block of original data xb,
## informer set A with size nA, top set size nT, return a list of {1. Top set index. 2. Posterior Expected Loss for xA, xb}
## Input: xb, Nb, xA, nA, A,nT
## Output: list of Top set and PEL_b

post_eloss <- function(xb, Nb, xA, nA, A,nT = 2,iter = 2){
  source("post_prob.R")
  source("post_theta.R")
  if(nT <= 1){
    print("Please require at least 2 top element.")
    return(0)
  }
  postls <- list(loss = 100, top = rep(0,nT))
  theta = post_theta(xb, Nb, xA, A, nA, thres = 10,iter = iter)
  postls$top = order(theta, decreasing = T)[1:nT]
  postls$loss = post_prob(xb, Nb,xA,A,nA) * sum(1 - theta[postls$top])
  return(postls)
}

##Test
# xb <- dat[1:3,]
# Nb <- 3
# xA <- dat[100, 100]
# A <- 100
# nA = 1
# nT = 2
# postls <- post_eloss(xb,Nb,xA,nA,A,nT, iter = 2)
