### New codes to score a candidate informer set for a given column clustering. 
### Analog to pel1_beta.R

score_inform <- function(cl, inform, train, m0s, block,
                         row_sample_size, interm_size, nT,
                         simplified = TRUE, thres = 15000){
  grp = max(cl$C)
  final_Scores = list()
  
  ## fill in the Scores vector by iterating through the column clusters
  for (i in 1:grp) {
    sub_cols = which(cl$C == i) # column idx in ith column cluster
    sub_train = train[ , sub_cols]
    m = nrow(sub_train)
    n = ncol(sub_train)
    a = rep(mean(sub_train, na.rm = T), ncol(sub_train))
    b = 1 - a
    m0 = m0s[i]
    
    ## recover cl_sample class from block list (just save the CC matrix)
    row_cl_sample = list(KK=apply(block[[i]], 1, max), NN = matrix(0, row_sample_size, nrow(sub_train)),
                         CC=block[[i]])
    for (j in 1:row_sample_size) {
      row_cl_sample$NN[j, 1:row_cl_sample$KK[j]] = as.numeric(table(block[[i]][j,]))
    }
    
    ## traditional BOISE 
    P = clust_sum(row_cl_sample, sub_train, row_sample_size, a, b)
    ### sub_inform is the index of all selected informers in sub-train matrix
    sub_inform = c() 
    for (informer in inform) {
      if(informer %in% sub_cols){
        sub_inform = c(sub_inform, which(sub_cols == informer))
      }
    }
    ## new scores list
    Scores = score_cpd(P, sub_train, sub_inform, nT, m0, a, b, row_sample_size, interm_size)
    ## update final_Scores list
    if(length(final_Scores) == 0){
      final_Scores = Scores
    } else if(length(Scores) == 1){
      new_xa = names(Scores)[1]
      final_Scores = lapply(final_Scores, function(x){
        new_sc = c(x$sc, Scores[[new_xa]]$sc)
        new_sc = sort(new_sc, decreasing = T)[1:min(nT, length(new_sc))]
        return(list(lg_wt = x$lg_wt, sc = new_sc))
      })
    } else{
      prev_lg_wts = unname(unlist(lapply(final_Scores, function(x){return(x$lg_wt)})))
      new_lg_wts = unname(unlist(lapply(Scores, function(x){return(x$lg_wt)})))
      sorted_lg_wts = sort(c(prev_lg_wts, new_lg_wts), decreasing = T)
      for (lg_wt_thres in sorted_lg_wts) {
        count = length(which(prev_lg_wts >= lg_wt_thres)) * length(which(new_lg_wts >= lg_wt_thres))
        if(count >= thres)
          break
      }
      final_Scores = Filter(function(x) x$lg_wt >= lg_wt_thres, final_Scores)
      Scores = Filter(function(x) x$lg_wt >= lg_wt_thres, Scores)
      tmp_Scores = list()
      for (xa in names(final_Scores)) {
        for (new_xa in names(Scores)) {
          key = paste(xa, new_xa, sep = '')
          new_sc = c(final_Scores[[xa]]$sc, Scores[[new_xa]]$sc)
          new_sc = sort(new_sc, decreasing = T)[1:min(nT, length(new_sc))]
          new_lg_wt = final_Scores[[xa]]$lg_wt + Scores[[new_xa]]$lg_wt
          tmp_Scores[[key]] = list(lg_wt = new_lg_wt, sc = new_sc)
        }
      }
      if(simplified & length(tmp_Scores) > thres){
        lg_wt_low = unlist(lapply(tmp_Scores, function(x){return(unname(x$lg_wt))})) # list of lg_wts
        lg_wt_low = sort(lg_wt_low, decreasing = T)[thres]
        tmp_Scores = Filter(function(x) x$lg_wt >= lg_wt_low, tmp_Scores)
      }
      final_Scores = tmp_Scores
    }
  }
  return(final_Scores)
}

