## Function 4: Update labels and eliminate clusters with 0 items..
## Input: A list of {Number of clusters K, # of targets within each cluster N[K], 
##        Parameters w.r.t each cluster phi[K], Cluster assignment for each target C[n]}
## Output: A list with K,N,C and phi updated correspondingly

label_update <- function(cl = list("K" = 2, "N" = rep(1,2), "phi" = matrix(0,2,2), "C" = c(1,2))){
  ord = order(cl$N,decreasing = T)
  cl$K = length(which(cl$N > 0))
  new_N = rep(0, cl$K)
  new_phi = matrix(0, cl$K, dim(cl$phi)[2])
  new_C = cl$C
  for (i in 1:cl$K) {
    tmp = ord[i]
    new_N[i] = cl$N[tmp]
    new_phi[i, ] = cl$phi[tmp, ]
    new_C[which(cl$C == tmp)] = i
  }
  cl$N = new_N
  cl$phi = new_phi
  cl$C = new_C
  return(cl)
}
