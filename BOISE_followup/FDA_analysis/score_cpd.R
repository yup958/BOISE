### New codes to score compounds in a submatrix. Analog to pel2_beta.R
### input a summary list P, sub_train and sub_inform for a subgroup of columns
### return a list of Scores = list(xA = list(lg_wt, sc)), where xA is possible outcome on
### sub_inform.

score_cpd <- function(P, sub_train, sub_inform, nT,  m0, a, b, row_sample_size, interm_size){
  m = nrow(sub_train)
  n = ncol(sub_train)
  Scores = list()
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
    Scores[[' ']] = list(lg_wt = 0, 
                         sc = sort(tmp_score, decreasing = T)[1:min(nT, ncol(sub_train))])
  } else{
    ### if there are informers in the 
    nA = length(sub_inform)
    XA = matrix(0, row_sample_size * interm_size, nA)
    
    ## Sample for intermediate x_A with sample size "interm_size"
    cur_position = 1
    for (k in 1:row_sample_size) {
      K = nrow(P[[k]])
      p = rep(0, K + 1)
      p[K + 1] = m0 / (m + m0)
      p[1:K] = P[[k]][ ,n+1] / (m + m0)
      for (i in 1:interm_size) {
        classi = which(rmultinom(1, 1, p) == 1)
        if(classi == K+1){
          post_theta = a[sub_inform] / (a[sub_inform] + b[sub_inform])
        } else{
          post_theta = P[[k]][classi,sub_inform]
        }
        interm_xA = as.numeric(rbinom(nA, 1, p = post_theta))
        XA[cur_position, ] = interm_xA
        cur_position = cur_position + 1
      }
    }
    
    ## Score computation
    XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
    tab = table(XA)
    YA = names(tab)
    l = length(tab)
    ### compute post_prob and post_thetas
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      post_probs = matrix(0, 1, row_sample_size)
      post_thetas = matrix(0, row_sample_size, n)
      for (k in 1:row_sample_size){
        postls = pel2_beta(P = P[[k]], x0 = sub_train, xA=xA, A=sub_inform, alpha=a, beta=b, m0=m0)
        post_probs[k] = postls$post_prob
        post_thetas[k, ] = postls$post_theta
      }
      # p(xA | x0)
      post_prob = mean(post_probs)
      # E(theta | x0, xA)
      post_probs = post_probs / (sum(post_probs))
      tmp_score = post_probs %*% post_thetas
      Scores[[YA[i]]] = list(lg_wt = log(post_prob), 
                             sc = sort(tmp_score, decreasing = T)[1:min(nT, ncol(sub_train))])
    }
  }
  return(Scores)
}
