#load data
set.seed(817)
source("/ua/pyu58/codes/clust_sum.R")
source("/ua/pyu58/codes/pel2_beta.R")

## candidate inform
args <- commandArgs()
id = as.numeric(args[3])
A = as.numeric(args[4])

# load original RData file
load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))

a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
sample_size = length(cl_sample$KK)
P = clust_sum(cl_sample, train, sample_size, a, b)

informs = read.table('/z/Comp/boise/results/orig_new_informer_1.txt')
pre_inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
#pre_inform = c()
inform = c(pre_inform, A)
m = nrow(train)
n = ncol(train)
nA = length(inform)
nT = as.integer(0.1 * n)
interm_size = 1000
XA = matrix(0, sample_size * interm_size, nA)

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
      post_theta = a[inform] / (a[inform] + b[inform])
    } else{
      post_theta = P[[k]][classi,inform]
    }
    interm_xA = as.numeric(rbinom(nA, 1, p = post_theta))
    XA[cur_position, ] = interm_xA
    cur_position = cur_position + 1
  }
}

## PEL 1 calculation, compatible with parallel computing
XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
tab = table(XA)
YA = names(tab)
l = length(tab)
rm(XA)
### PEL 2 calculation
lg_wts = rep(0, l)
pel2s = rep(0, l)
unit_vec = matrix(1, ncol(train),1)
for (i in 1:l) {
  xA = as.numeric(strsplit(YA[i], "")[[1]])
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
  #pel2s[i] = sum(sort(1-Score)[1:nT])
  Score = sort(Score, decreasing = T)
  weights = (1:n) - 1
  pel2s[i] = sum(Score * weights)
}
## calculate PEL1
lg_wts = lg_wts - max(lg_wts)
wts = exp(lg_wts) / sum(exp(lg_wts))
pel1 = sum(wts * pel2s)
# tmp_pel2 = sapply(1:l, function(i){
#   xA = as.numeric(strsplit(YA[i], "")[[1]])
#   post_probs = matrix(0, 1, sample_size)
#   post_thetas = matrix(0, sample_size, n)
#   for (k in 1:sample_size){
#     postls = pel2_beta(P = P[[k]], x0 = train, xA=xA, A=inform, nT=nT, alpha=a, beta=b, m0=m0)
#     post_probs[k] = postls$post_prob
#     post_thetas[k, ] = postls$post_theta
#   }
#   post_probs = post_probs / (sum(post_probs))
#   Score = post_probs %*% post_thetas
#   return(sum(sort(1-Score)[1:nT]))
# })
# pel1 = (tmp_pel2 * tab) / sample_size
# pel1 = sum(pel1) / interm_size

## Write result
result = data.frame('id' = id, 'A' = A,  "PEL" = pel1) 
if(A < 10){
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA), "_00", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else if(A < 100){
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA),"_0", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else{
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA),"_", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
}

