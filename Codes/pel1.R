## Func 4 for BOISE: Given informer index A, and size nA, given DPMM result cl, calculate PEL1
## Input: A, nA, nT, cl, data
## Output: PEL1 of A
## Use package R.utils

pel1 <- function(A, nA, nT = 2, cl, data, iter = 2){
  if(cl$K == 1){
    print("Please run DPMM before calculating PEL")
    return(1)
  }
  source("post_eloss.R")
  B = cl$K
  pel = 0
  if(nA == 1){
    for (b in 1:B) {
      xb = data[which(cl$C == b), ]
      Nb = cl$N[b]
      tmp = post_eloss(xb, Nb, 0, nA, A, nT, iter)
      pel = pel + tmp$loss
      tmp = post_eloss(xb, Nb, 1, nA, A, nT, iter)
      pel = pel + tmp$loss
    }
    pel = pel / B
  } else{
    X = matrix(0, 2^nA, nA)
    for (i in 1:(2^nA-1)) {
      s = R.utils::intToBin(i)
      n = nchar(s)
      for (j in (nA-n+1):nA) {
        X[i+1,j] = as.integer(substr(s,j-nA+n,j-nA+n))
      }
    }
    for (b in 1:B) {
      xb = data[which(cl$C == b), ]
      Nb = cl$N[b]
      Loss = apply(X, 1, function(x){
        tmp = post_eloss(xb, Nb, x, nA, A, nT,iter)
        return(tmp$loss)
      })
      pel = pel + sum(Loss)
    }
    pel = pel / B
  }
  return(pel)  
}

##Test
# A = c(10,100)
# nA = 2
# nT = 3
# cl = list(K = 2, N = c(11, 4), C = c(rep(1,11), rep(2,4)))
# data = dat[1:15, ]
# pel1(A,nA,nT,cl,data,iter = 2)
