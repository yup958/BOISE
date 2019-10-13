## Func 1 for BOISE: Posterior estimation for theta (alpha + bta)
## Input: Some block of original data xb, number of targets in xb, Nb, (possible) informer set xA, informer index set A
## If xb has more than thres = 5 targets, may need re-clustering use DPMM
## Output: Posterior mean of theta

post_theta <- function(xb, Nb, xA, A, nA){
  source("single_post.R")
  source("multi_post.R")
  if(nA == 1){
    theta = single_post(xb, Nb, xA, A)
  } else{
    theta = multi_post(xb,Nb,xA,A)
  }
   return(theta)
}


# #Test
# xb <- dat[1:10,]
# Nb <- 10
# xA <- dat[100, c(1,10,100,200)]
# A <- c(1,10,100,200)
# nA = 4
# theta <- post_theta(xb,Nb,xA,A,nA)
