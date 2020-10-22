## Calculating entropy for Accelerated BOISE
## Input: P, a matrix summarize information of each clustering;
##        inform, candidate informer set;
##        XX, posterior samples of x_i*;
##        priors a and b, prior mass alpha.
## Output: Entropy under this clustering assignment and candidate informer set.


Entropy <- function(P, inform, XX, alpha, a, b){
  K = nrow(P)
  m = ncol(P) - 1
  size = nrow(XX)
  XA = as.matrix(XX[, inform])
  XA = apply(XA, 1, function(x){
    return(paste(x, collapse = ""))
  })
  atom = unique(XA)
  entropy = sapply(atom, function(s){
    Ps = length(which(XA == s)) / size
    xA = as.numeric(unlist(strsplit(s, split = "")))
    pc = rep(0, K+1)
    pc[K+1] = log(alpha) + sum(xA*(log(a/(a+b))[inform])) + sum((1-xA)*(log(b/(a+b))[inform]))
    pc[1:K] = apply(P,1,function(x){
      return(log(x[m+1]) + sum(xA*(log(x[1:m])[inform])) + sum((1-xA)*(log(1-x[1:m])[inform])))
    })
    pc = exp(pc) / (sum(exp(pc)))
    return(Ps * sum(-pc*log(pc)))
  })
  entropy = sum(entropy)
  return(entropy)
}
