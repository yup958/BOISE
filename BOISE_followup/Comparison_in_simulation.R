### Import updated functions that take penalization into account (as in revision paper)
setwd('~/RAwork/BOISE/BOISE_followup/')
source_files = c('Boise.R', 'clust_sum.R', 'dpmm_beta.R', 'Evaluate.R', 'Initial_beta.R',
                 'pel1_beta.R', 'pel2_beta.R', 'Update_beta.R', 'post_theta_j.R')
for(f in source_files){
  source(f)
}
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

package_cvboise <- function(M, m0, nA, nT){
  roc_results = rep(0, nrow(M))
  nef_results = rep(0, nrow(M))
  selective_nef = rep(0, nrow(M))
  for (i in 1:nrow(M)) {
    train = M[-i, ]
    test = M[i, ]
    a = rep(mean(train), ncol(train))
    b = 1 - a
    cl_sample = BOISE::dpmm_beta(train, a, b, m0, burn_in = 500, sample_size=20, thinning = 5)
    informer = BOISE::Boise(cl_sample, sample_size=20, interm_size = 1000, nA, nT, a, b,
                     train, m0, mcParallel = T)
    roc_results[i] = BOISE::Evaluate(cl_sample, informer, "rocauc", percent=0.1, test, train, nT,
                              sample_size=20, a, b, m0)
    nef_results[i] = BOISE::Evaluate(cl_sample, informer, "nef", percent=0.1,test, train, nT,
                              sample_size=20, a, b, m0)
    selective_nef[i] = BOISE::Evaluate(cl_sample, informer, "selectivity_nef", percent=0.1,test, train, nT,
                                sample_size=20, a, b, m0)
  }
  return(c(roc_results, nef_results, selective_nef))
}

cvboise <- function(M, m0, nA, nT){
  roc_results = rep(0, nrow(M))
  nef_results = rep(0, nrow(M))
  selective_nef = rep(0, nrow(M))
  for (i in 1:nrow(M)) {
    train = M[-i, ]
    test = M[i, ]
    a = rep(mean(train), ncol(train))
    b = 1 - a
    cl_sample = dpmm_beta(train, a, b, m0, burn_in = 500, sample_size=20, thinning = 5)
    informer = Boise(cl_sample, sample_size=20, interm_size = 1000, nA, nT, a, b,
                     train, m0, mcParallel = T)
    roc_results[i] = Evaluate(cl_sample, informer, "rocauc", percent=0.1, test, train, nT,
                              sample_size=20, a, b, m0)
    nef_results[i] = Evaluate(cl_sample, informer, "nef", percent=0.1,test, train, nT,
                              sample_size=20, a, b, m0)
    selective_nef[i] = Evaluate(cl_sample, informer, "selectivity_nef", percent=0.1,test, train, nT,
                              sample_size=20, a, b, m0)
  }
  return(c(roc_results, nef_results, selective_nef))
}

### Simulation generated from CRP. 
set.seed(11)
N = 10
results_mat = matrix(0, nrow = N, ncol = 90)
results_m0_lo = list(old = results_mat, new = results_mat)
results_m0_mid = list(old = results_mat, new = results_mat)
results_m0_hi = list(old = results_mat, new = results_mat)

for (i in 1:N) {
  print(i)
  alpha = 0.1
  beta = 0.9
  
  # m0 = 3 generating matrices
  cl = BOISE::Initial_beta(x0 = matrix(0, 30, 50), m0 = 3)
  P_m0_lo = matrix(0,30,50)
  probs = rbeta(50 * cl$K, alpha, beta)
  dim(probs) = c(cl$K, 50)
  for (r in 1:30) {
    P_m0_lo[r, ] = probs[cl$C[r], ]
  }
  M = M_generate(P_m0_lo)
  m0 = m0_Find(M, lower = 2, upper = 10)
  results_m0_lo$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_m0_lo$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
  
  # m0 = 6 generating matrices
  cl = BOISE::Initial_beta(x0 = matrix(0, 30, 50), m0 = 6)
  P_m0_mid = matrix(0,30,50)
  probs = rbeta(50 * cl$K, alpha, beta)
  dim(probs) = c(cl$K, 50)
  for (r in 1:30) {
    P_m0_mid[r, ] = probs[cl$C[r], ]
  }
  M = M_generate(P_m0_mid)
  m0 = m0_Find(M, lower = 2, upper = 10)
  results_m0_mid$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_m0_mid$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
  
  # m0 = 10 generating matrices
  cl = BOISE::Initial_beta(x0 = matrix(0, 30, 50), m0 = 9)
  P_m0_hi = matrix(0,30,50)
  probs = rbeta(50 * cl$K, alpha, beta)
  dim(probs) = c(cl$K, 50)
  for (r in 1:30) {
    P_m0_hi[r, ] = probs[cl$C[r], ]
  }
  M = M_generate(P_m0_hi)
  m0 = m0_Find(M, lower = 1, upper = 10)
  results_m0_hi$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_m0_hi$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
}
save(list = c("results_m0_lo","results_m0_mid","results_m0_hi"), 
     file = "penalized_loss_compar_synthetic_with_CRP_result.RData")


### Simulation as in revised BOISE paper, fully / half / none clustering
set.seed(11)
N = 10
results_mat = matrix(0, nrow = N, ncol = 90)
results_full_cluster = list(old = results_mat, new = results_mat)
results_half_cluster = list(old = results_mat, new = results_mat)
results_none_cluster = list(old = results_mat, new = results_mat)
for (i in 1:N) {
  print(i)
  alpha = 0.1
  beta = 0.9
  # Full clustered
  P_full = matrix(0,30,50)
  probs = rbeta(150, alpha, beta)
  dim(probs) = c(3, 50)
  for (r in 1:30) {
    P_full[r, ] = probs[((r-1) %/% 10 + 1), ]
  }
  M = M_generate(P_full)
  m0 = m0_Find(M, lower = 2, upper = 10)
  results_full_cluster$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_full_cluster$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
  # Half clustered
  P_half = matrix(0,30,50)
  left_probs = rbeta(75, alpha, beta)
  dim(left_probs) = c(3, 25)
  for (r in 1:30) {
    P_half[r, 1:25] = left_probs[((r-1) %/% 10 + 1), ]
    P_half[r, 26:50] = rbeta(25, alpha, beta)
  }
  M = M_generate(P_half)
  m0 = m0_Find(M, lower = 2, upper = 10)
  results_half_cluster$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_half_cluster$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
  # None clustered
  P_none = matrix(0,30,50)
  for (r in 1:30) {
    P_none[r, ] = rbeta(50, alpha, beta)
  }
  M = M_generate(P_none)
  m0 = m0_Find(M, lower = 2, upper = 10)
  results_none_cluster$old[i,] = package_cvboise(M, m0 = m0, nA = 5, nT = 10)
  results_none_cluster$new[i,] = cvboise(M, m0 = m0, nA = 5, nT = 10)
}
save(list = c("results_full_cluster","results_half_cluster","results_none_cluster"), 
     file = "penalized_loss_compar_result.RData")

### Plot comparison
m = 30
k = 10
title_size = 18
axis_title_size = 16
axis_text_size = 16
temp = data.frame("ROCAUC" = rep(0, 2*k*m), "Penalization" = c(rep("No", k*m),
                                                     rep("Yes",k*m)))
temp$ROCAUC = c(as.vector(results_m0_lo$old[1:k, 1:m]), as.vector(results_m0_lo$new[1:k, 1:m]))
temp$Penalization = as.factor(temp$Penalization)
temp$Penalization <- factor(temp$Penalization, levels =c("No","Yes"))
p1 <- ggplot(temp, aes(x = Penalization, y = ROCAUC, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="ROCAUC with m0 = 3") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$ROCAUC = c(as.vector(results_m0_mid$old[1:k, 1:m]), as.vector(results_m0_mid$new[1:k, 1:m]))
p2 <- ggplot(temp, aes(x = Penalization, y = ROCAUC, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="ROCAUC with m0 = 6") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$ROCAUC = c(as.vector(results_m0_hi$old[1:k, 1:m]), as.vector(results_m0_hi$new[1:k, 1:m]))
p3 <- ggplot(temp, aes(x = Penalization, y = ROCAUC, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="ROCAUC with m0 = 9") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp = data.frame("NEF10" = rep(0, 2*k*m), "Penalization" = c(rep("No", k*m),
                                                              rep("Yes",k*m)))
temp$NEF10 = c(as.vector(results_m0_lo$old[1:k, (m+1):(2*m)]), as.vector(results_m0_lo$new[1:k, (m+1):(2*m)]))
temp$Penalization = as.factor(temp$Penalization)
temp$Penalization <- factor(temp$Penalization, levels =c("No","Yes"))
p4 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="NEF10 with m0 = 3") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$NEF10 = c(as.vector(results_m0_mid$old[1:k, (m+1):(2*m)]), as.vector(results_m0_mid$new[1:k, (m+1):(2*m)]))
p5 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="NEF10 with m0 = 6") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$NEF10 = c(as.vector(results_m0_hi$old[1:k, (m+1):(2*m)]), as.vector(results_m0_hi$new[1:k, (m+1):(2*m)]))
p6 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="NEF10 with m0 = 6") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

myGrobs <- list(p4,p5,p6,p1,p2,p3)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 3 )


temp$NEF10 = c(as.vector(results_m0_lo$old[1:k, (2*m+1):(3*m)]), as.vector(results_m0_lo$new[1:k, (2*m+1):(3*m)]))
temp$Penalization = as.factor(temp$Penalization)
temp$Penalization <- factor(temp$Penalization, levels =c("No","Yes"))
p7 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="Synthetic Data with m0 = 3", y = "Selectivity NEF10") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$NEF10 = c(as.vector(results_m0_mid$old[1:k, (2*m+1):(3*m)]), as.vector(results_m0_mid$new[1:k, (m+1):(2*m)]))
p8 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="Synthetic Data with m0 = 6", y = "Selectivity NEF10") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))

temp$NEF10 = c(as.vector(results_m0_hi$old[1:k, (2*m+1):(3*m)]), as.vector(results_m0_hi$new[1:k, (m+1):(2*m)]))
p9 <- ggplot(temp, aes(x = Penalization, y = NEF10, fill = Penalization)) +  
  # geom_violin(scale = "width", trim = T, adjust = 0.5 )+
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c("#009E73","#D55E00",
                               "#000000"))+
  theme(legend.position="none") + 
  stat_summary(fun.y=median, geom="point", size=1, color="white") +
  labs(title="Synthetic Data with m0 = 9", y = "Selectivity NEF10") +
  theme(plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size))
myGrobs <- list(p7,p8,p9)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 1,ncol = 3 )
