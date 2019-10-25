### Innitialization based on Chinese Restaurant Process for conjugate prior DPMM
## Input: number of targets to be clustered n, divergency parameter alpha
## Output: A list of {K, N, C}. K = number of clusters, N = #of targets in each cluster, C = cluster assignment

Initial_beta <- function(n, alpha = 2){
  N = rep(0, 2 * n)
  C = rep(0, n)
  K = 1
  
  N[1] = 1
  C[1] = 1
  for (i in 2:n) {
    p = rep(0, K + 1)
    p[K + 1] = alpha / (i - 1 + alpha)
    for (k in 1:K) {
      p[k] = N[k] / (i - 1 + alpha)
    }
    class = which(rmultinom(1, 1, p) == 1)
    C[i] = class
    N[class] = N[class] + 1
    if(class == K + 1){K = K + 1}
  }
  
  cl = list("K" = K, "N" = N, "C" = C)
  return(cl)
}


# #Test
# cl = Initial_beta(200,10)
