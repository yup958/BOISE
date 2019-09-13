## Function 5 Rearrange Cluster Assignment
## Input: A list of {Number of clusters K, # of targets within each cluster N[K], 
##        Parameters w.r.t each cluster phi[K], Cluster assignment for each target C[n]}
##        that has been Innitialized and updated parameters.
##        
## Output: A same type list with label already updated (i.e. no empty clusters and clusters are ranked)

Rearrange <- function(cl = list("K" = 2, "N" = rep(1,2), "phi" = matrix(0,2,2), "C" = c(1,2)), 
                      a = 0, aux = 5, mu1 = 0, sigma1=1, mu2=0, sigma2=1, data){
  source("prob.R")
  source("parasample.R")
  source("label_update.R")
  if(aux <= 1){
    print("Please Allow least 2 Auxillary Parameters")
    return(1)
  }
  n = dim(data)[1]
  m = dim(data)[2]
  new_K = cl$K
  new_N = rep(0,dim(data)[1])
  new_N[1:cl$K] = cl$N
  new_C = cl$C
  new_phi = matrix(0, n, m + 1)
  new_phi[1:cl$K, ] = cl$phi
  
  for (j in 1:n) {
    aux_pars <- parasample(aux = aux, dimen = m + 1, mu1,sigma1,mu2,sigma2)
    proba = rep(0, new_K + aux)
    for (i in 1:new_K) {
      proba[i] = new_N[i] * prob(data[j, ], phi = new_phi[i, ]) / (n - 1 + a)
    }
    for (i in (new_K + 1):(new_K + aux)) {
      proba[i] = prob(data[j, ], phi = aux_pars[(i - new_K), ]) * a / (aux * (n - 1 + a))
    }
    proba = proba / sum(proba)
    tmp = which(rmultinom(1,1,proba) == 1)
    if(tmp <= new_K){
      new_N[tmp] = new_N[tmp] + 1
      new_N[new_C[j]] = new_N[new_C[j]] - 1
      new_C[j] = tmp
    } else{
      new_N[new_C[j]] = new_N[new_C[j]] - 1
      new_K = new_K + 1
      new_C[j] = new_K
      new_N[new_K] = 1
      new_phi[new_K, ] = aux_pars[(tmp - new_K + 1), ]
    }
  }
  
  new_cl = list("K" = new_K, "N" = new_N[1:new_K], "phi" = new_phi[1:new_K, ], "C" = new_C)
  new_cl = label_update(new_cl)
  return(new_cl)
}

# #test
# source("Initialization.R")
# dat <- data[1:5, ]
# cl = Initialization(dat, a = 2, aux = 2)
# new_cl = Rearrange(cl, a= 2, aux = 3, data = dat)
