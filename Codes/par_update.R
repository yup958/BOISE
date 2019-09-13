## Function 3: Update Posterior Parameter w.r.t Clusters
## Input: A list of {Number of clusters K, # of targets within each cluster N[K], 
##        Parameters w.r.t each cluster phi[K], Cluster assignment for each target C[n]}, data X
## Output: A list with only parameters updated
## We use Stan in a loop. See what happens

par_update <- function(cl = list("K" = 2, "N" = rep(1,2), "phi" = matrix(0,2,2), "C" = c(1,2)), data){
  n = dim(data)[1]
  m = dim(data)[2]
  for (i in 1:cl$K) {
    if(cl$N[i] == 0){
      next
    }
    subdat = data[which(cl$C == i), ]
    tmp <- list(J = cl$N[i], K = m, x = subdat)
    if(cl$N[i] == 1){
      fit <- stan(file = 'Rasch_Vector.stan', data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    } else{
      fit <- stan(file = 'Rasch_Matrix.stan', data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    }
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- mean(la$alpha)
    btahat <- apply(la$bta, 2, mean)
    cl$phi[i, 1] = alphahat
    cl$phi[i, 2:(m+1)] = btahat
  }
  return(cl)
}

# ## Test
# setwd("~/Drug_Discovery/BOISE")
# library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# 
# phi = matrix(0,3,367)
# cl = list("K" = 3, "N" = c(2,0,4), "phi" = phi, "C" = c(3,3,1,3,1,3))
# dat = data[1:6, ]
# new_cl = par_update(cl, dat)

