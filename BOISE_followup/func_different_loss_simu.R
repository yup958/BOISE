### for L_0 loss: sum of (1-\theta_ij)
pel2_l0 <-
  function (P, x0, xA, A, nT = 10, alpha, beta, m0){
    K = nrow(P)
    m = ncol(x0)
    n = nrow(x0)
    nA = length(A)
    lp = rep(0, K + 1) #posterior probabilities of log p(xA|ck,C,x0)
    theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
    
    theta[K+1,] = alpha / (alpha + beta)
    #theta[K+1,A] = xA
    
    # p(xA|i in ck, C, x0)
    lp[K + 1] = sum(xA * log(alpha[A] / (alpha[A] + beta[A]))) + 
      sum((1 - xA) * log(beta[A] / (alpha[A] + beta[A])))
    if (K==1){
      lp[K] = sum(xA * log(P[A])) + sum((1-xA) * log(1 - P[A]))
    } else{
      lp[1:K] = apply(P[,1:m], 1, function(q){
        return(sum(xA * log(q[A])) + sum((1-xA) * log(1 - q[A])))
      })
    }
    
    theta[1:K,] = P[,1:m]
    if(length(A) == 1){
      theta[1:K, A] = (P[,A] * (P[,m+1] + rep(alpha[1]+beta[1], K)) + xA) /
        (P[,m+1] + rep(alpha[1]+beta[1], K) + 1)
    } else{
      theta[1:K,A] = t(as.matrix(apply(P[,c(A,m+1)], 1, function(q){
        l = length(A)
        counts = q[1:l] * (q[l+1]+alpha[1]+beta[1])
        counts = counts + xA
        return(counts / (q[l+1]+alpha[1]+beta[1]+1))
      })))
    }
    #theta[1:K, A] = t(matrix(rep(xA, K), nrow = nA, ncol = K))
    
    # P(xA|C,x0)
    weight = c(P[,m+1], m0)
    log_weight = log(weight / (sum(weight)))
    log_weight = log_weight + lp
    c = max(log_weight)
    log_weight = log_weight - c
    log_post_prob = c + log(sum(exp(log_weight)))
    
    ##E(theta|C,x0,xA)
    lw = log(c(P[, m + 1], m0))
    lw = lw + lp
    lw = lw - max(lw)
    w = as.matrix(exp(lw) / sum(exp(lw)))
    post_theta = t(w) %*% theta
    postls <- list(post_theta = post_theta, log_post_prob = log_post_prob)
    return(postls)
  }
pel1_l0 <-
  function(cl_sample, P, sample_size, interm_size, A, nT = 10, alpha, beta, x0, m0 = 2){
    ### Compute PEL1
    nA = length(A)
    XA = matrix(0, sample_size * interm_size, nA)
    m = nrow(x0)
    n = ncol(x0)
    ## Sample for intermediate x_A with sample size "interm_size"
    cur_position = 1
    for (k in 1:sample_size) {
      K = cl_sample$KK[k]
      p = rep(0, K + 1)
      p[K + 1] = m0 / (m + m0)
      p[1:K] = P[[k]][ ,n+1] / (m + m0)
      for (i in 1:interm_size) {
        classi = which(rmultinom(1, 1, p) == 1)
        if(classi == K+1){
          post_theta = a[A] / (a[A] + b[A])
        } else{
          post_theta = P[[k]][classi,A]
        }
        interm_xA = as.numeric(rbinom(nA, 1, p = post_theta))
        XA[cur_position, ] = interm_xA
        cur_position = cur_position + 1
      }
    }
    XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
    tab = table(XA)
    YA = names(tab)
    l = length(tab)
    lg_wts = rep(0, l)
    pel2s = rep(0, l)
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, n)
      for (k in 1:sample_size){
        postls = pel2_l0(P = P[[k]], x0 = x0, xA=xA, A=A, alpha=alpha, beta=beta, m0=m0)
        log_post_probs[k] = postls$log_post_prob
        post_thetas[k, ] = postls$post_theta
      }
      # p(xA | x0)
      c = max(log_post_probs)
      log_post_probs = log_post_probs - c
      log_post_prob = -log(sample_size) + c + log(sum(exp(log_post_probs))) ## log(mean of post_probs)
      lg_wts[i] = log_post_prob
      # E(theta | x0, xA)
      post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
      Score = post_probs %*% post_thetas
      pel2s[i] = sum(sort(1-Score)[1:nT])
    }
    ## calculate PEL1
    lg_wts = lg_wts - max(lg_wts)
    wts = exp(lg_wts) / sum(exp(lg_wts))
    pel1 = sum(wts * pel2s)
    return(pel1)
  }
Boise_l0 <-
  function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0,
           mcParallel = TRUE){
    ### Input check
    m = nrow(x0)
    n = ncol(x0)
    
    P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
    ## BOISE selection based on pel1
    step = 1
    if(!mcParallel){
      pel1 = unlist(lapply(1:ncol(x0), function(x){
        return(pel1_l0(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))}))
    } else{
      pel1 = unlist(mclapply(1:ncol(x0), function(x){
        return(pel1_l0(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))},
        mc.cores = detectCores()))
    }  
    tmp = order(pel1)[1]
    inform = tmp
    candidate = order(pel1)
    while (step < nA) {
      step = step +1
      candidate = candidate[-which(candidate == tmp)]
      pel1 = rep(0,length(candidate))
      if(!mcParallel){
        pel1 = unlist(lapply(candidate, function(x){
          return(pel1_l0(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))}))
      } else{
        pel1 = unlist(mclapply(candidate, function(x){
          return(pel1_l0(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
          mc.cores = detectCores()))
      }
      tmp = candidate[order(pel1)[1]]
      inform = c(inform, tmp)
    }
    return(inform)
  }

### for L_1 loss: sum of (1-x_ij)
pel2 <-
  function (P, x0, xA, A, nT = 10, alpha, beta, m0){
    K = nrow(P)
    m = ncol(x0)
    n = nrow(x0)
    nA = length(A)
    lp = rep(0, K + 1) #posterior probabilities of log p(xA|ck,C,x0)
    theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
    
    theta[K+1,] = alpha / (alpha + beta)
    theta[K+1,A] = xA
    
    # p(xA|i in ck, C, x0)
    lp[K + 1] = sum(xA * log(alpha[A] / (alpha[A] + beta[A]))) + 
      sum((1 - xA) * log(beta[A] / (alpha[A] + beta[A])))
    if (K==1){
      lp[K] = sum(xA * log(P[A])) + sum((1-xA) * log(1 - P[A]))
    } else{
      lp[1:K] = apply(P[,1:m], 1, function(q){
        return(sum(xA * log(q[A])) + sum((1-xA) * log(1 - q[A])))
      })
    }
    
    theta[1:K,] = P[,1:m]
    if(length(A) == 1){
      theta[1:K, A] = (P[,A] * (P[,m+1] + rep(alpha[1]+beta[1], K)) + xA) /
        (P[,m+1] + rep(alpha[1]+beta[1], K) + 1)
    } else{
      theta[1:K,A] = t(as.matrix(apply(P[,c(A,m+1)], 1, function(q){
        l = length(A)
        counts = q[1:l] * (q[l+1]+alpha[1]+beta[1])
        counts = counts + xA
        return(counts / (q[l+1]+alpha[1]+beta[1]+1))
      })))
    }
    theta[1:K, A] = t(matrix(rep(xA, K), nrow = nA, ncol = K))
    
    # P(xA|C,x0)
    weight = c(P[,m+1], m0)
    log_weight = log(weight / (sum(weight)))
    log_weight = log_weight + lp
    c = max(log_weight)
    log_weight = log_weight - c
    log_post_prob = c + log(sum(exp(log_weight)))
    
    ##E(theta|C,x0,xA)
    lw = log(c(P[, m + 1], m0))
    lw = lw + lp
    lw = lw - max(lw)
    w = as.matrix(exp(lw) / sum(exp(lw)))
    post_theta = t(w) %*% theta
    postls <- list(post_theta = post_theta, log_post_prob = log_post_prob)
    return(postls)
  }
pel1_l1 <-
  function(cl_sample, P, sample_size, interm_size, A, nT = 10, alpha, beta, x0, m0 = 2){
    #source("npel2.R")
    ### Compute PEL1
    nA = length(A)
    XA = matrix(0, sample_size * interm_size, nA)
    m = nrow(x0)
    n = ncol(x0)
    ## Sample for intermediate x_A with sample size "interm_size"
    cur_position = 1
    for (k in 1:sample_size) {
      K = cl_sample$KK[k]
      p = rep(0, K + 1)
      p[K + 1] = m0 / (m + m0)
      p[1:K] = P[[k]][ ,n+1] / (m + m0)
      for (i in 1:interm_size) {
        classi = which(rmultinom(1, 1, p) == 1)
        if(classi == K+1){
          post_theta = a[A] / (a[A] + b[A])
        } else{
          post_theta = P[[k]][classi,A]
        }
        interm_xA = as.numeric(rbinom(nA, 1, p = post_theta))
        XA[cur_position, ] = interm_xA
        cur_position = cur_position + 1
      }
    }
    XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
    tab = table(XA)
    YA = names(tab)
    l = length(tab)
    lg_wts = rep(0, l)
    pel2s = rep(0, l)
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, n)
      for (k in 1:sample_size){
        postls = pel2(P = P[[k]], x0 = x0, xA=xA, A=A, alpha=alpha, beta=beta, m0=m0)
        log_post_probs[k] = postls$log_post_prob
        post_thetas[k, ] = postls$post_theta
      }
      # p(xA | x0)
      c = max(log_post_probs)
      log_post_probs = log_post_probs - c
      log_post_prob = -log(sample_size) + c + log(sum(exp(log_post_probs))) ## log(mean of post_probs)
      lg_wts[i] = log_post_prob
      # E(theta | x0, xA)
      post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
      Score = post_probs %*% post_thetas
      pel2s[i] = sum(sort(1-Score)[1:nT])
    }
    ## calculate PEL1
    lg_wts = lg_wts - max(lg_wts)
    wts = exp(lg_wts) / sum(exp(lg_wts))
    pel1 = sum(wts * pel2s)
    return(pel1)
  }

Boise_l1 <-
  function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0,
           mcParallel = TRUE){
    
    ### Input check
    m = nrow(x0)
    n = ncol(x0)
    
    P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
    ## BOISE selection based on pel1
    step = 1
    if(!mcParallel){
      pel1 = unlist(lapply(1:ncol(x0), function(x){
        return(pel1_l1(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))}))
    } else{
      pel1 = unlist(mclapply(1:ncol(x0), function(x){
        return(pel1_l1(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))},
        mc.cores = detectCores()))
    }  
    tmp = order(pel1)[1]
    inform = tmp
    candidate = order(pel1)
    while (step < nA) {
      step = step +1
      candidate = candidate[-which(candidate == tmp)]
      pel1 = rep(0,length(candidate))
      if(!mcParallel){
        pel1 = unlist(lapply(candidate, function(x){
          return(pel1_l1(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))}))
      } else{
        pel1 = unlist(mclapply(candidate, function(x){
          return(pel1_l1(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
          mc.cores = detectCores()))
      }
      tmp = candidate[order(pel1)[1]]
      inform = c(inform, tmp)
    }
    return(inform)
  }

### for L_2 loss: ranking loss
pel1_l2 <-
  function(cl_sample, P, sample_size, interm_size, A, nT = 10, alpha, beta, x0, m0 = 2){
    #source("npel2.R")
    ### Compute PEL1
    nA = length(A)
    XA = matrix(0, sample_size * interm_size, nA)
    m = nrow(x0)
    n = ncol(x0)
    ## Sample for intermediate x_A with sample size "interm_size"
    cur_position = 1
    for (k in 1:sample_size) {
      K = cl_sample$KK[k]
      p = rep(0, K + 1)
      p[K + 1] = m0 / (m + m0)
      p[1:K] = P[[k]][ ,n+1] / (m + m0)
      for (i in 1:interm_size) {
        classi = which(rmultinom(1, 1, p) == 1)
        if(classi == K+1){
          post_theta = a[A] / (a[A] + b[A])
        } else{
          post_theta = P[[k]][classi,A]
        }
        interm_xA = as.numeric(rbinom(nA, 1, p = post_theta))
        XA[cur_position, ] = interm_xA
        cur_position = cur_position + 1
      }
    }
    XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
    tab = table(XA)
    YA = names(tab)
    l = length(tab)
    lg_wts = rep(0, l)
    pel2s = rep(0, l)
    # unit_vec = matrix(1, ncol(train),1)
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, n)
      for (k in 1:sample_size){
        postls = pel2(P = P[[k]], x0 = x0, xA=xA, A=A, alpha=alpha, beta=beta, m0=m0)
        log_post_probs[k] = postls$log_post_prob
        post_thetas[k, ] = postls$post_theta
      }
      # p(xA | x0)
      c = max(log_post_probs)
      log_post_probs = log_post_probs - c
      log_post_prob = -log(sample_size) + c + log(sum(exp(log_post_probs))) ## log(mean of post_probs)
      lg_wts[i] = log_post_prob
      # E(theta | x0, xA)
      post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
      Score = post_probs %*% post_thetas
      Score = sort(Score, decreasing = T)
      weights = (1:n) - 1
      pel2s[i] = sum(Score * weights)
    }
    ## calculate PEL1
    lg_wts = lg_wts - max(lg_wts)
    wts = exp(lg_wts) / sum(exp(lg_wts))
    pel1 = sum(wts * pel2s)
    return(pel1)
  }

Boise_l2 <-
  function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0,
           mcParallel = TRUE){
    ### Input check
    m = nrow(x0)
    n = ncol(x0)
    
    P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
    ## BOISE selection based on pel1
    step = 1
    if(!mcParallel){
      pel1 = unlist(lapply(1:ncol(x0), function(x){
        return(pel1_l2(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))}))
    } else{
      pel1 = unlist(mclapply(1:ncol(x0), function(x){
        return(pel1_l2(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))},
        mc.cores = detectCores()))
    }  
    tmp = order(pel1)[1]
    inform = tmp
    candidate = order(pel1)
    while (step < nA) {
      step = step +1
      candidate = candidate[-which(candidate == tmp)]
      pel1 = rep(0,length(candidate))
      if(!mcParallel){
        pel1 = unlist(lapply(candidate, function(x){
          return(pel1_l2(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))}))
      } else{
        pel1 = unlist(mclapply(candidate, function(x){
          return(pel1_l2(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
          mc.cores = detectCores()))
      }
      tmp = candidate[order(pel1)[1]]
      inform = c(inform, tmp)
    }
    return(inform)
  }

## Evaluating function
Evaluate_l0 <-
  function(cl_sample, inform, measure, percent,
           test, train, nT, sample_size, alpha, beta, m0){
    P = clust_sum(cl_sample,train,sample_size, alpha, beta)
    Score = rep(0, ncol(train))
    xA = test[inform]
    nA=length(inform)
    m = ncol(train)
    log_post_probs = matrix(0, 1, sample_size)
    post_thetas = matrix(0, sample_size, m)
    for (k in 1:sample_size){
      postls = pel2_l0(P[[k]], x0=train, xA, A=inform, nT, alpha, beta, m0)
      log_post_probs[k] = postls$log_post_prob
      post_thetas[k, ] = postls$post_theta
    }
    log_post_probs = log_post_probs - max(log_post_probs)
    post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
    Score = post_probs %*% post_thetas
    test = as.vector(test)
    Score = as.vector(Score)
    if(measure == "nef"){
      nTop = round(ncol(train) * percent, 0)
      top = order(Score,decreasing = T)[1:nTop]
      pred_hit = sum(test[top])
      hit = sum(test)
      maxhit = min(hit,nTop)
      result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
    } else if(measure == "rocauc"){
      rocobj = pROC::roc(test,Score)
      result = rocobj$auc
    } else if(measure == "prauc"){
      fg = Score[test == 1]
      bg = Score[test == 0]
      pr = pr.curve(fg, bg)
      result = pr$auc.integral
    } else if(measure %in% c("mat", "f")){
      pred.obj = ROCR::prediction(Score, test)
      perform.obj = ROCR::performance(pred.obj, measure)
      result = max(unlist(perform.obj@y.values),na.rm = T)
    } else{
      print("Criteria is not supported.")
      result = 0
    }
    return(result)
  }
Evaluate_l1 <-
  function(cl_sample, inform, measure, percent,
           test, train, nT, sample_size, alpha, beta, m0){
    P = clust_sum(cl_sample,train,sample_size, alpha, beta)
    Score = rep(0, ncol(train))
    xA = test[inform]
    nA=length(inform)
    m = ncol(train)
    log_post_probs = matrix(0, 1, sample_size)
    post_thetas = matrix(0, sample_size, m)
    for (k in 1:sample_size){
      postls = pel2(P[[k]], x0=train, xA, A=inform, nT, alpha, beta, m0)
      log_post_probs[k] = postls$log_post_prob
      post_thetas[k, ] = postls$post_theta
    }
    log_post_probs = log_post_probs - max(log_post_probs)
    post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
    Score = post_probs %*% post_thetas
    test = as.vector(test)
    Score = as.vector(Score)
    if(measure == "nef"){
      nTop = round(ncol(train) * percent, 0)
      top = order(Score,decreasing = T)[1:nTop]
      pred_hit = sum(test[top])
      hit = sum(test)
      maxhit = min(hit,nTop)
      result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
    } else if(measure == "rocauc"){
      rocobj = pROC::roc(test,Score)
      result = rocobj$auc
    } else if(measure == "prauc"){
      fg = Score[test == 1]
      bg = Score[test == 0]
      pr = pr.curve(fg, bg)
      result = pr$auc.integral
    } else if(measure %in% c("mat", "f")){
      pred.obj = ROCR::prediction(Score, test)
      perform.obj = ROCR::performance(pred.obj, measure)
      result = max(unlist(perform.obj@y.values),na.rm = T)
    } else{
      print("Criteria is not supported.")
      result = 0
    }
    return(result)
  }

