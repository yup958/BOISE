## Simplified verision of Entropy function
## Input: cluster assignment cl, informer set inform, original data matrix dat
## Output: Entropy under this cluster assignment and informer set.


Entropy2 <- function(cl, inform, dat){
  n = nrow(dat)
  K = cl$K
  subdat = as.matrix(dat[,inform])
  subdat = apply(subdat, 1, function(x){
    return(paste(x, collapse = ""))
  })
  atom = unique(subdat)
  entropy = sapply(atom, function(s){
    target = which(subdat == s)
    Ps = length(target) / n
    P = sapply(1:K, function(k){
      classk = which(cl$C == k)
      return(length(intersect(classk,target)))
    })
    P = P / sum(P)
    P = P[which(P>0)]
    return(Ps * sum(-P*log(P)))
  })
  entropy = sum(entropy)
  return(entropy)
}
