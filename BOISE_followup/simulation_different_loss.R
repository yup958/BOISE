setwd('~/RAwork/BOISE/BOISE_followup')
set.seed(817)
library(BOISE)
### Define a function to find m0 that minimize prior and posterior difference
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

m0_Find <- function(M, lower, upper){
  a = rep(mean(M), ncol(M))
  b = 1 - a
  n = nrow(M)
  min_diff = 1000
  best_m0 = 1
  
  for(m0 in lower:upper){
    prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
    cl_sample = dpmm_beta(M, a, b, m0, burn_in = 500, sample_size=10, thinning = 5)
    diff = abs(mean(cl_sample$KK) - prior_m0_expect)
    if(diff < min_diff){
      min_diff = diff
      best_m0 = m0
    }
  }
  return(best_m0)
}

M_generate <- function(P){
  M = apply(P, 2, function(p){
    return(rbinom(length(p), 1, p))
  })
  M_sum = apply(M, 1, sum)
  while (length(which(M_sum == 0)) > 0) {
    M = apply(P, 2, function(p){
      return(rbinom(length(p), 1, p))
    })
    M_sum = apply(M, 1, sum)
  }
  return(M)
}

new_pel1_beta <-
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
    unit_vec = matrix(1, ncol(x0),1)
    lg_wts = rep(0, l)
    pel2s = rep(0, l)
    unit_vec = matrix(1, ncol(train),1)
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, n)
      for (k in 1:sample_size){
        postls = pel2_beta(P = P[[k]], x0 = x0, xA=xA, A=A, alpha=alpha, beta=beta, m0=m0)
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
      # pel2s[i] = sum(sort(1-Score)[1:nT])
      # Score = matrix(Score, ncol=1)
      # pel2s[i] = -sum(abs(Score %*% t(unit_vec) - unit_vec %*% t(Score)))
      Score = sort(Score, decreasing = T)
      nT = which(Score <= mean(Score))[1] - 1
      pel2s[i] = (n-nT) * sum(Score[1:nT]) - nT * sum(Score[(nT+1):n])
    }
    ## calculate PEL1
    lg_wts = lg_wts - max(lg_wts)
    wts = exp(lg_wts) / sum(exp(lg_wts))
    pel1 = sum(wts * pel2s)
    return(-pel1)
  }

new_Boise <-
  function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0,
           mcParallel = TRUE){
    #source("clust_sum.R")
    #source("npel1.R")
    if (!require('parallel')) {
      install.packages("parallel")
      library(parallel)
    }
    
    ### Input check
    m = nrow(x0)
    n = ncol(x0)
    
    P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
    ## BOISE selection based on pel1
    step = 1
    if(!mcParallel){
      pel1 = unlist(lapply(1:ncol(x0), function(x){
        return(new_pel1_beta(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))}))
    } else{
      pel1 = unlist(mclapply(1:ncol(x0), function(x){
        return(new_pel1_beta(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))},
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
          return(new_pel1_beta(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))}))
      } else{
        pel1 = unlist(mclapply(candidate, function(x){
          return(new_pel1_beta(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
          mc.cores = detectCores()))
      }
      tmp = candidate[order(pel1)[1]]
      inform = c(inform, tmp)
    }
    return(inform)
  }
pel1_beta <-
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
    unit_vec = matrix(1, ncol(x0),1)
    lg_wts = rep(0, l)
    pel2s = rep(0, l)
    unit_vec = matrix(1, ncol(train),1)
    for (i in 1:l) {
      xA = as.numeric(strsplit(YA[i], "")[[1]])
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, n)
      for (k in 1:sample_size){
        postls = pel2_beta(P = P[[k]], x0 = x0, xA=xA, A=A, alpha=alpha, beta=beta, m0=m0)
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

Boise <-
  function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0,
           mcParallel = TRUE){
    #source("clust_sum.R")
    #source("npel1.R")
    if (!require('parallel')) {
      install.packages("parallel")
      library(parallel)
    }
    
    ### Input check
    m = nrow(x0)
    n = ncol(x0)
    
    P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
    ## BOISE selection based on pel1
    step = 1
    if(!mcParallel){
      pel1 = unlist(lapply(1:ncol(x0), function(x){
        return(pel1_beta(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))}))
    } else{
      pel1 = unlist(mclapply(1:ncol(x0), function(x){
        return(pel1_beta(cl_sample, P, sample_size, interm_size, A = x, nT,alpha, beta, x0, m0))},
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
          return(pel1_beta(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))}))
      } else{
        pel1 = unlist(mclapply(candidate, function(x){
          return(pel1_beta(cl_sample, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
          mc.cores = detectCores()))
      }
      tmp = candidate[order(pel1)[1]]
      inform = c(inform, tmp)
    }
    return(inform)
  }
#####################################################################
#########                  Simulations                    ###########
#####################################################################
N = 10
sample_size = 20
interm_size = 1000
# old_roc_results_mat = matrix(0, nrow = N, ncol = 30)
# old_nef_results_mat = matrix(0, nrow = N, ncol = 30)
old_pr_results_mat = matrix(0, nrow = N, ncol = 30)
# new_roc_results_mat = matrix(0, nrow = N, ncol = 30)
# new_nef_results_mat = matrix(0, nrow = N, ncol = 30)
new_pr_results_mat = matrix(0, nrow = N, ncol = 30)
for (i in 1:N) {
  print(i)
  alpha = 0.1
  beta = 0.9
  # # Full clustered
  # P_full = matrix(0,30,50)
  # probs = rbeta(150, alpha, beta)
  # dim(probs) = c(3, 50)
  # for (r in 1:30) {
  #   P_full[r, ] = probs[((r-1) %/% 10 + 1), ]
  # }
  # M = M_generate(P_full)
  
  # Half clustered
  P_half = matrix(0,30,50)
  left_probs = rbeta(75, alpha, beta)
  dim(left_probs) = c(3, 25)
  for (r in 1:30) {
    P_half[r, 1:25] = left_probs[((r-1) %/% 10 + 1), ]
    P_half[r, 26:50] = rbeta(25, alpha, beta)
  }
  M = M_generate(P_half)
  
  # # None clustered
  # P_none = matrix(0,30,50)
  # for (r in 1:30) {
  #   P_none[r, ] = rbeta(50, alpha, beta)
  # }
  # M = M_generate(P_none)
  
  m0 = m0_Find(M, lower = 2, upper = 10)
  old_roc_results = rep(0, nrow(M))
  old_nef_results = rep(0, nrow(M))
  old_pr_results = rep(0, nrow(M))
  new_roc_results = rep(0, nrow(M))
  new_nef_results = rep(0, nrow(M))
  new_pr_results = rep(0, nrow(M))
  for (j in 1:nrow(M)) {
    train = M[-j, ]
    test = M[j, ]
    a = rep(mean(train), ncol(train))
    b = 1 - a
    cl_sample = dpmm_beta(train, a, b, m0, burn_in = 500, sample_size=20, thinning = 5)
    ## new informer search
    new_inform = new_Boise(cl_sample, sample_size, interm_size, nA=5, nT=10, a, b, train, m0)
    new_roc_results[j] = Evaluate(cl_sample, new_inform, "rocauc", percent=0.1, test, train, nT,
                                  sample_size=20, a, b, m0)
    new_nef_results[j] = Evaluate(cl_sample, new_inform, "nef", percent=0.1,test, train, nT,
                                  sample_size=20, a, b, m0)
    new_pr_results[j] = Evaluate(cl_sample, new_inform, "prauc", percent=0.1, test, train, nT,
                                 sample_size=20, a, b, m0)
    ## old informer search
    inform = Boise(cl_sample, sample_size, interm_size, nA=5, nT=10, a, b, train, m0)
    old_roc_results[j] = Evaluate(cl_sample, inform, "rocauc", percent=0.1, test, train, nT,
                              sample_size=20, a, b, m0)
    old_nef_results[j] = Evaluate(cl_sample, inform, "nef", percent=0.1,test, train, nT,
                              sample_size=20, a, b, m0)
    old_pr_results[j] = Evaluate(cl_sample, inform, "prauc", percent=0.1, test, train, nT,
                              sample_size=20, a, b, m0)
  }
  # old_roc_results_mat[i, ] = old_roc_results
  # old_nef_results_mat[i, ] = old_nef_results
  old_pr_results_mat[i, ] = old_pr_results
  # new_nef_results_mat[i, ] = new_nef_results
  # new_roc_results_mat[i, ] = new_roc_results
  new_pr_results_mat[i, ] = new_pr_results
}
new_roc_results_unclust = as.vector(roc_results_mat)
new_nef_results_unclust = as.vector(nef_results_mat)
old_roc_results_unclust = as.vector(roc_results_mat)
old_nef_results_unclust = as.vector(nef_results_mat)
load("tmp_simu_result.RData")
# old_results = data.frame(old_roc_full = old_roc_results, old_roc_half = old_roc_results_half,
#                      old_roc_unclust = old_roc_results_unclust, old_nef_full = old_nef_results,
#                      old_nef_half = old_nef_results_half, old_nef_unclust = old_nef_results_unclust)
# new_results = data.frame(new_roc_full = new_roc_results, new_roc_half = new_roc_results_half,
#                          new_roc_unclust = new_roc_results_unclust, new_nef_full = new_nef_results,
#                          new_nef_half = new_nef_results_half, new_nef_unclust = new_nef_results_unclust)
old_results$old_pr_none = as.vector(old_pr_results_mat)
new_results$new_pr_none = as.vector(new_pr_results_mat)
save(list = c("old_results", "new_results"), 
     file = "tmp_simu_result.RData")

m = 30
k = 10
title_size = 18
axis_title_size = 16
axis_text_size = 16
load("tmp_simu_result.RData")
temp = data.frame("ROCAUC" = rep(0, 6*k*m), "NEF10" = rep(0, 6*k*m), "IBR" = c(rep("Boise_new", 3*k*m), rep("Boise_old", 3*k*m)),
                  "Setting" = rep(c(rep("Fully-clustered", k*m), rep("Half-clustered", k*m), 
                                    rep("Unclustered", k*m)), 2))
temp$ROCAUC = c(new_results$new_roc_full, new_results$new_roc_half, new_results$new_roc_unclust,
                old_results$old_roc_full, old_results$old_roc_half, old_results$old_roc_unclust)
temp$NEF10 = c(new_results$new_nef_full, new_results$new_nef_half, new_results$new_nef_unclust,
               old_results$old_nef_full, old_results$old_nef_half, old_results$old_nef_unclust)
pd <- position_dodge(0.5)
p1 <- ggplot(temp, aes(x = Setting, y = ROCAUC, fill = IBR)) +
  geom_boxplot(width = 0.5)+
  scale_fill_manual(values = c("#D55E00", "#CC79A7",
                               "#000000"))+
  stat_summary(fun = "mean", geom="point", size=1, color="white", position = pd) +
  labs(title="ROCAUC Comparison", x="Clustering Setting") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(size = axis_text_size))
p2 <- ggplot(temp, aes(x = Setting, y = NEF10, fill = IBR)) +
  geom_boxplot(width = 0.5)+
  scale_fill_manual(values = c("#D55E00", "#CC79A7",
                               "#000000"))+
  stat_summary(fun = "mean", geom="point", size=1, color="white", position = pd) +
  labs(title="NEF10 Comparison", x="Clustering Setting") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(size = axis_text_size))

myGrobs <- list(p1,p2)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 1,ncol = 2)

p1 <- ggplot()+
  geom_point(aes(x = old_results$old_pr_none, y = new_results$new_pr_none))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of PRAUC old vs. new loss', x = 'old Boise PRAUC', y = 'new Boise PRAUC')
p1
length(which(new_results$new_roc_half > old_results$old_roc_half))

