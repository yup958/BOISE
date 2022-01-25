load('/ua/pyu58/results/col_clust_res_5.RData')
source("/ua/pyu58/codes/dpmm_beta.R")
source("/ua/pyu58/codes/Initial_beta.R")
source("/ua/pyu58/codes/Update_beta.R")
set.seed(2)
m0_selections = read.table('/ua/pyu58/submit_files/Prior_mass_FDA_separate.txt')

train = t(train)
set.seed(2)

args <- commandArgs()
k = as.numeric(args[3])
sample_size = 100
m0_selections = m0_selections[which(m0_selections$V1==k),]
block = list()

## merge all column clusters with size <= 2
for (grp in 1:cl_sample$KK[k]){
  if(cl_sample$NN[k, grp] <= 2)
    break
}
cl_sample$KK[k] = grp
cl_sample$CC[k, which(cl_sample$CC[k,] > grp)] = grp

for (grp in 1:cl_sample$KK[k]) {
  col_grp = which(cl_sample$CC[k, ]==grp)
  train_sub = train[, col_grp]
  a = rep(mean(train_sub, na.rm = T), ncol(train_sub))
  b = 1 - a
  m0 = 1
  if(cl_sample$NN[k, grp] >= 15){
    m0 = m0_selections[which(m0_selections$V2 == grp), 3]
  }
  block[[grp]] = dpmm_beta(train_sub, a, b, m0= m0, burn_in = 1000, sample_size=sample_size, thinning = 15)
}
## only store the cluster assignments to save space.
for (grp in 1:length(block)) {
  block[[grp]] = block[[grp]]$CC
}
save(list = 'block', file = paste('block_res_', as.character(k), '.RData', sep = ""))
     