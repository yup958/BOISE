#load data
set.seed(817)
source("/ua/pyu58/codes/clust_sum.R")
source("/ua/pyu58/codes/pel2_beta.R")

## candidate inform
args <- commandArgs()
id = as.numeric(args[3])
A = as.numeric(args[4])

# load original RData file
load(paste('/ua/pyu58/data/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))

a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
sample_size = length(cl_sample$KK)#500
P = clust_sum(cl_sample, train, sample_size, a, b)

inform = A
m = nrow(train)
n = ncol(train)
nA = length(inform)
nT = as.integer(0.1 * n)
## PEL 1 calculation, compatible with parallel computing
lg_wts = rep(0, 2)
pel2s = rep(0, 2)
XA = c(0,1)
for (i in 1:2) {
  xA = XA[i]
  log_post_probs = matrix(0, 1, sample_size)
  post_thetas = matrix(0, sample_size, n)
  for (k in 1:sample_size){
    postls = pel2_beta(P = P[[k]], x0 = train, xA=xA, A=inform, alpha=a, beta=b, m0=m0)
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

## Write result
result = data.frame('id' = id, 'A' = A,  "PEL" = pel1) 
if(A < 10){
  write.table(result,file = paste("testid_", as.character(id),  "_00", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else if(A < 100){
  write.table(result,file = paste("testid_", as.character(id), "_0", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else{
  write.table(result,file = paste("testid_", as.character(id), "_", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
}

