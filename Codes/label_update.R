## Function 4: Update labels and eliminate clusters with 0 items..
## Input: A list of {Number of clusters K, # of targets within each cluster N[K], 
##        Parameters w.r.t each cluster phi[K], Cluster assignment for each target C[n]}
## Output: A list with K,N,C and phi updated correspondingly

label_update <- function(cl = list("K" = 2, "N" = rep(1,2), "phi" = matrix(0,2,2), "C" = c(1,2)),col){
  cl$K = length(which(cl$N > 0))
  if(cl$K == 1){
    ## If there are only 1 cluster, we add some uncertainty by forcing one customer leave the table.
    cl$K = 2
    cl$N[which(cl$N > 0)] = cl$N[which(cl$N > 0)] - 1
    cl$N = c(cl$N,1)
    cl$C[1] = length(cl$N)
    cl$phi = matrix(0, length(cl$N), col)
  }
  ord = order(cl$N,decreasing = T)
  new_N = rep(0, cl$K)
  new_phi = matrix(0, cl$K, col)
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
