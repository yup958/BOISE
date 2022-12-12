setwd("~/RAwork/BOISE/BOISE_followup/")
library(BOISE)
data(pkis1)
x0 <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
counts = apply(x0, 2, sum)
x0 = x0[,which(counts > 0)]
alpha = rep(mean(x0,na.rm = T),ncol(x0))
beta = 1-alpha

### inside the dpmm_beta function
m = nrow(x0)
n = ncol(x0)
m0 = 15
sample_size = 100
burn_in = 1000
thinning = 50
## MCMC for random clustering
KK = rep(0,sample_size)
NN = matrix(0, sample_size, 2*m)
CC = matrix(0,sample_size, m)

### Initialization: for different starting point
## for random starts
cl_sample = dpmm_beta(x0, alpha, beta, m0, burn_in, sample_size, thinning)
## for all one cluster
N = rep(0, 2 * m)
N[1] = m
C = rep(1, m)
K = 1
cl = list("K" = K, "N" = N, "C" = C)
## for all separate clusters
N = rep(0, 2 * m)
N[1:m] = rep(1, m)
C = 1:m
K = m
cl = list("K" = K, "N" = N, "C" = C)

### burn-in stage
for (i in 1:burn_in) {
  cl = Update_beta(cl, x0, alpha, beta, m0)
}
### Sample clustering assignments
for (i in 1:sample_size) {
  tmp = 0
  while (tmp < thinning) {
    cl = Update_beta(cl, x0, alpha, beta, m0)
    tmp = tmp + 1
  }
  KK[i] = cl$K
  NN[i, ] = cl$N
  CC[i, ] = cl$C
}

cl_sample_whole = list(KK = KK, NN = NN, CC = CC)
cl_sample_single = list(KK = KK, NN = NN, CC = CC)

### Plot
adj_mat_whole = matrix(0, nrow(x0), nrow(x0))
for (i in 1:nrow(x0)) {
  for (j in 1:nrow(x0)) {
    adj_mat_whole[i,j] = sum(cl_sample_whole$CC[,i]==cl_sample_whole$CC[,j])
  }
}
adj_mat_whole = adj_mat_whole / sample_size

adj_mat_single = matrix(0, nrow(x0), nrow(x0))
for (i in 1:nrow(x0)) {
  for (j in 1:nrow(x0)) {
    adj_mat_single[i,j] = sum(cl_sample_single$CC[,i]==cl_sample_single$CC[,j])
  }
}
adj_mat_single = adj_mat_single / sample_size
par(mfrow = c(1,2))
image(adj_mat_whole, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 220, 20), side=2, line=0.5, at=seq(0, 220, 20) / 224, las=1, cex=1)
mtext(text=seq(0, 220, 20), side=1, line=0.5, at=seq(0, 220, 20) / 224, las=1, cex=1)

image(adj_mat_single, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 220, 20), side=2, line=0.5, at=seq(0, 220, 20) / 224, las=1, cex=1)
mtext(text=seq(0, 220, 20), side=1, line=0.5, at=seq(0, 220, 20) / 224, las=1, cex=1)

