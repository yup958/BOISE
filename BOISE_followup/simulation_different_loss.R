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


#####################################################################
#########                  Simulations                    ###########
#####################################################################
N = 10
sample_size = 20
interm_size = 1000
l0_roc_results_mat = matrix(0, nrow = N, ncol = 30)
l0_nef_results_mat = matrix(0, nrow = N, ncol = 30)
l1_roc_results_mat = matrix(0, nrow = N, ncol = 30)
l1_nef_results_mat = matrix(0, nrow = N, ncol = 30)
l2_roc_results_mat = matrix(0, nrow = N, ncol = 30)
l2_nef_results_mat = matrix(0, nrow = N, ncol = 30)
for (i in 1:N) {
  print(i)
  alpha = 0.1
  beta = 0.9
  # # Full clustered
  # P_full = matrix(0,30,60)
  # probs = rbeta(180, alpha, beta)
  # dim(probs) = c(3, 60)
  # for (r in 1:30) {
  #   P_full[r, ] = probs[((r-1) %/% 10 + 1), ]
  # }
  # M = M_generate(P_full)
  
  # # Half clustered
  # P_half = matrix(0,30,60)
  # left_probs = rbeta(90, alpha, beta)
  # dim(left_probs) = c(3, 30)
  # for (r in 1:30) {
  #   P_half[r, 1:30] = left_probs[((r-1) %/% 10 + 1), ]
  #   P_half[r, 31:60] = rbeta(30, alpha, beta)
  # }
  # M = M_generate(P_half)
  
  # None clustered
  P_none = matrix(0,30,60)
  for (r in 1:30) {
    P_none[r, ] = rbeta(60, alpha, beta)
  }
  M = M_generate(P_none)
  
  m0 = m0_Find(M, lower = 1, upper = 10)
  roc_results_l0 = rep(0, nrow(M))
  roc_results_l1 = rep(0, nrow(M))
  roc_results_l2 = rep(0, nrow(M))
  nef_results_l0 = rep(0, nrow(M))
  nef_results_l1 = rep(0, nrow(M))
  nef_results_l2 = rep(0, nrow(M))
  for (j in 1:nrow(M)) {
    train = M[-j, ]
    test = M[j, ]
    a = rep(mean(train), ncol(train))
    b = 1 - a
    cl_sample = dpmm_beta(train, a, b, m0, burn_in = 500, sample_size, thinning = 5)
    ## L_2 informer search
    inform_l2 = Boise_l2(cl_sample, sample_size, interm_size, nA=3, nT=6, a, b, train, m0)
    roc_results_l2[j] = Evaluate_l1(cl_sample, inform_l2, "rocauc", percent=0.1, test, train, nT,
                                    sample_size, a, b, m0)
    nef_results_l2[j] = Evaluate_l1(cl_sample, inform_l2, "nef", percent=0.1,test, train, nT,
                                    sample_size, a, b, m0)
    ## L_1 informer search
    inform_l1 = Boise_l1(cl_sample, sample_size, interm_size, nA=3, nT=6, a, b, train, m0)
    roc_results_l1[j] = Evaluate_l1(cl_sample, inform_l1, "rocauc", percent=0.1, test, train, nT,
                                    sample_size, a, b, m0)
    nef_results_l1[j] = Evaluate_l1(cl_sample, inform_l1, "nef", percent=0.1,test, train, nT,
                                    sample_size, a, b, m0)
    ## L_0 informer search
    inform_l0 = Boise_l0(cl_sample, sample_size, interm_size, nA=3, nT=6, a, b, train, m0)
    roc_results_l0[j] = Evaluate_l0(cl_sample, inform_l0, "rocauc", percent=0.1, test, train, nT,
                                  sample_size, a, b, m0)
    nef_results_l0[j] = Evaluate_l0(cl_sample, inform_l0, "nef", percent=0.1,test, train, nT,
                                  sample_size, a, b, m0)
  }
  l0_roc_results_mat[i, ] = roc_results_l0
  l1_roc_results_mat[i, ] = roc_results_l1
  l2_roc_results_mat[i, ] = roc_results_l2
  l0_nef_results_mat[i, ] = nef_results_l0
  l1_nef_results_mat[i, ] = nef_results_l1
  l2_nef_results_mat[i, ] = nef_results_l2
}
load("simu_result_l012.RData")
# l0_results = data.frame(roc_full = as.vector(l0_roc_results_mat), nef_full = as.vector(l0_nef_results_mat))
# l1_results = data.frame(roc_full = as.vector(l1_roc_results_mat), nef_full = as.vector(l1_nef_results_mat))
# l2_results = data.frame(roc_full = as.vector(l2_roc_results_mat), nef_full = as.vector(l2_nef_results_mat))
l0_results$roc_none = as.vector(l0_roc_results_mat)
l0_results$nef_none = as.vector(l0_nef_results_mat)
l1_results$roc_none = as.vector(l1_roc_results_mat)
l1_results$nef_none = as.vector(l1_nef_results_mat)
l2_results$roc_none = as.vector(l2_roc_results_mat)
l2_results$nef_none = as.vector(l2_nef_results_mat)

save(list = c("l0_results", "l1_results", "l2_results"), 
     file = "simu_result_l012.RData")
load("tmp_simu_result.RData")
# old_results = data.frame(old_roc_full = old_roc_results, old_roc_half = old_roc_results_half,
#                      old_roc_unclust = old_roc_results_unclust, old_nef_full = old_nef_results,
#                      old_nef_half = old_nef_results_half, old_nef_unclust = old_nef_results_unclust)
# new_results = data.frame(new_roc_full = new_roc_results, new_roc_half = new_roc_results_half,
#                          new_roc_unclust = new_roc_results_unclust, new_nef_full = new_nef_results,
#                          new_nef_half = new_nef_results_half, new_nef_unclust = new_nef_results_unclust)
save(list = c("old_results", "new_results"), 
     file = "tmp_simu_result.RData")

m = 30
k = 10
title_size = 14
axis_title_size = 14
axis_text_size = 13
load("simu_result_l012.RData")
temp = data.frame("ROCAUC" = rep(0, 9*k*m), "NEF10" = rep(0, 9*k*m), 
                  "IBR" = c(rep("Boise_L0", 3*k*m), rep("Boise_L1", 3*k*m),rep("Boise_L2", 3*k*m)),
                  "Setting" = rep(c(rep("Fully-clustered", k*m), rep("Half-clustered", k*m), 
                                    rep("Unclustered", k*m)), 3))
temp$ROCAUC = c(l0_results$roc_full, l0_results$roc_half, l0_results$roc_none,
                l1_results$roc_full, l1_results$roc_half, l1_results$roc_none,
                l2_results$roc_full, l2_results$roc_half, l2_results$roc_none)
temp$NEF10 = c(l0_results$nef_full, l0_results$nef_half, l0_results$nef_none,
               l1_results$nef_full, l1_results$nef_half, l1_results$nef_none,
               l2_results$nef_full, l2_results$nef_half, l2_results$nef_none)
pd <- position_dodge(0.5)
p1 <- ggplot(temp, aes(x = Setting, y = ROCAUC, fill = IBR)) +
  geom_boxplot(width = 0.5)+
  scale_fill_manual(values = c("#F0E442","#D55E00", "#CC79A7",
                               "#000000"))+
  stat_summary(fun = "mean", geom="point", size=1, color="white", position = pd) +
  labs(title="ROCAUC Comparison", x="Clustering Setting") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(size = axis_text_size),
        legend.position = c(0.85,0.9))
p2 <- ggplot(temp, aes(x = Setting, y = NEF10, fill = IBR)) +
  geom_boxplot(width = 0.5)+
  scale_fill_manual(values = c("#F0E442","#D55E00", "#CC79A7",
                               "#000000"))+
  stat_summary(fun = "mean", geom="point", size=1, color="white", position = pd) +
  labs(title="NEF10 Comparison", x="Clustering Setting") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = axis_title_size),
        legend.text = element_text(size = axis_text_size),
        legend.position = c(0.85,0.9))

myGrobs <- list(p1,p2)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 1,ncol = 2)

p1 <- ggplot()+
  geom_point(aes(x = old_results$old_pr_none, y = new_results$new_pr_none))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of PRAUC old vs. new loss', x = 'old Boise PRAUC', y = 'new Boise PRAUC')
p1
length(which(new_results$new_roc_half > old_results$old_roc_half))

