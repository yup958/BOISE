## Func 1 for BOISE: Posterior estimation for theta (alpha + bta)
## Input: Some block of original data xb, number of targets in xb, Nb, (possible) informer set xA, informer index set A
## If xb has more than thres = 5 targets, may need re-clustering use DPMM
## Output: Posterior mean of theta

post_theta <- function(xb, Nb, xA, A, nA, thres = 10,iter = 2){
  source("dpmm.R")
  source("single_post.R")
  source("multi_post.R")
  if(nA == 1){
    theta = single_post(xb, Nb, xA, A)
  } else{
    theta = multi_post(xb,Nb,xA,A)
  }
  ## Only 1 target in this block
  # if(Nb <= thres){          
  #   if(nA == 1){
  #     theta = single_post(xb, Nb, xA, A)
  #   } else{
  #     theta = multi_post(xb,Nb,xA,A)
  #   }
  # } else{
  #   subcl <- dpmm(data = xb, a = Nb %/% 10 + 1, aux = max(2, Nb %/% 20), iter = iter)
  #   Clas = subcl$K
  #   N = subcl$N
  #   C = subcl$C
  #   theta <- matrix(0, Clas, dim(xb)[2])
  #   for (i in 1:Clas) {
  #     if(nA == 1){
  #       theta[i, ] <- single_post(xb[which(C == i),], N[i], xA, A)
  #     } else{
  #       theta[i, ] <- multi_post(xb[which(C == i),], N[i], xA, A)
  #     }
  #   }
  #   theta <- apply(theta,2, mean)
  # }
   return(theta)
}


# #Test
# xb <- dat[1:10,]
# Nb <- 10
# xA <- dat[100, c(1,10,100,200)]
# A <- c(1,10,100,200)
# nA = 4
# theta <- post_theta(xb,Nb,xA,A,nA)
