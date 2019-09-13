## Function 4 Initializing Chinese Restaurant Process
## Input: data x, concentration parameter a, base distribution for alpha (mu1, sigma1) for beta (mu2, sigma2)
## Output: A list of {Number of clusters K, # of targets within each cluster N[K], 
##        Parameters w.r.t each cluster phi[K], Cluster assignment for each target C[n]}

Initialization <- function(data, a = 0, aux = 5, mu1 = 0, sigma1=1, mu2=0, sigma2=1){
  source("parasample.R")
  source("prob.R")
  if(aux <= 1){
    print("Please Allow least 2 Auxillary Parameters")
    return(1)
  }
  aux_pars <- parasample(aux = aux, dimen = dim(data)[2] + 1, mu1,sigma1,mu2,sigma2)
  K = 1
  N = rep(0,dim(data)[1])
  C = rep(0,dim(data)[1])
  phi = matrix(0, dim(data)[1], dim(data)[2] + 1)
  
  N[1] = 1
  C[1] = 1
  proba = rep(0, aux)
  for (i in 1:aux) {
    proba[i] = prob(data[1, ], aux_pars[i, ])
  }
  proba = proba / sum(proba)
  phi[1, ] = aux_pars[which(rmultinom(1,1,proba) == 1), ]
  for (j in 2:dim(data)[1]) {
    aux_pars <- parasample(aux = aux, dimen = dim(data)[2] + 1, mu1,sigma1,mu2,sigma2)
    proba = rep(0, K + aux)
    for (i in 1:K) {
      proba[i] = N[i] * prob(data[j, ], phi = phi[i, ]) / (j - 1 + a)
    }
    for (i in (K + 1):(K + aux)) {
      proba[i] = prob(data[j, ], phi = aux_pars[(i - K), ]) * a / (aux * (j - 1 + a))
    }
    proba = proba / sum(proba)
    tmp = which(rmultinom(1,1,proba) == 1)
    if(tmp > K){
      K = K + 1
      N[K] = 1
      C[j] = K
      phi[K, ] = aux_pars[(tmp - K + 1), ] 
    } else{
      N[tmp] = N[tmp] + 1
      C[j] = tmp
    }
  }
  cl = list("K" = K, "N" = N[1:K], "phi" = phi[1:K, ],"C" = C)
  return(cl)
}


## test
# dat <- data[1:5, ]
# cl = Initialization(dat, a = 2, aux = 2)
# source("label_update.R")
# source("par_update.R")
# cl = label_update(cl)
