### function to return scores of each cpd when having the intermediate data
### test should be a vector

evaluate_interm <- function(cl, inform, train, test, m0_selection, block, row_sample_size){
  grp = max(cl)
  Scores = rep(0, ncol(train))
  
  ## fill in the Scores vector by iterating through the column clusters
  for (i in 1:grp) {
    sub_cols = which(cl == i) # column idx in ith column cluster
    sub_train = train[ , sub_cols]
    m = nrow(sub_train)
    n = ncol(sub_train)
    a = rep(mean(sub_train, na.rm = T), ncol(sub_train))
    b = 1 - a
    ### revised for chemical clustering!!! V2 -> V1, 3 -> 2
    m0 = 1
    if(i %in% m0_selection$V2){
      m0 = m0_selection[which(m0_selection$V2 == i), 3] # retrieve m0 prior
    }
    
    ## recover cl_sample class from block list (just save the CC matrix)
    row_cl_sample = list(KK=apply(block[[i]], 1, max), NN = matrix(0, row_sample_size, nrow(sub_train)),
                         CC=block[[i]])
    for (j in 1:row_sample_size) {
      row_cl_sample$NN[j, 1:row_cl_sample$KK[j]] = as.numeric(table(block[[i]][j,]))
    }
    
    ## traditional BOISE 
    P = clust_sum(row_cl_sample, sub_train, row_sample_size, a, b)
    ### sub_xA and sub_inform of all selected informers in sub-train matrix
    sub_inform = c() 
    sub_xA = c()
    for (informer in inform) {
      if(informer %in% sub_cols & !is.na(test[informer])){
        sub_inform = c(sub_inform, which(sub_cols == informer))
        sub_xA = c(sub_xA, test[informer])
      }
    }
    if(is.null(sub_inform)){
      ### if no informer is in this subgroup of columns
      post_thetas = matrix(0, row_sample_size, n)
      for (k in 1:row_sample_size){
        w = c(P[[k]][, n+1], m0)
        w = as.matrix(w / sum(w))
        post_theta = rbind(P[[k]][, 1:n], a/(a+b))
        post_thetas[k, ] = t(w) %*% post_theta
      }
      tmp_score = apply(post_thetas, 2, mean)
    } else{
      post_probs = matrix(0, 1, row_sample_size)
      post_thetas = matrix(0, row_sample_size, n)
      for (k in 1:row_sample_size){
        postls = pel2_beta(P = P[[k]], x0 = sub_train, xA=sub_xA, A=sub_inform, alpha=a, beta=b, m0=m0)
        post_probs[k] = postls$post_prob
        post_thetas[k, ] = postls$post_theta
      }
      post_probs = post_probs / (sum(post_probs))
      tmp_score = post_probs %*% post_thetas
    }
    Scores[sub_cols] = tmp_score
  }
  return(Scores)
}
