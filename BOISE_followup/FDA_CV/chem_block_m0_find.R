#load('/ua/pyu58/results/col_clust_res_5.RData')

source("/ua/pyu58/codes/dpmm_beta.R")
source("/ua/pyu58/codes/Initial_beta.R")
source("/ua/pyu58/codes/Update_beta.R")
set.seed(817)
args <- commandArgs()
id = as.numeric(args[3])
grp = as.numeric(args[4])
m0 = as.numeric(args[5])

test_ids = read.table('/ua/pyu58/data/Test_IDS.txt')
test_id = test_ids[id,1]
## helper function to find m0
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

load('/ua/pyu58/data/fda_data.RData')
test = dat[rownames(dat)==test_id,]
complete_idx = which(!is.na(test))
train = dat[!row.names(dat)==test_id, complete_idx]
test = test[complete_idx]
cids = colnames(train)

sample_size = 100
cl = read.csv('/ua/pyu58/data/chemical_clustering_res.csv', header = T)
cl = cl[cl$cid %in% cids,]
cl$C = NA
tab = table(cl$clust)
ord = order(tab, decreasing = T)
for (k in 1:length(tab)) {
  cl$C[which(cl$clust == names(tab)[ord[k]])] = k
}
cl$clust = NULL
tab = table(cl$C)
for (k in 1:length(tab)){
  if(tab[k] < 10)
    break
}
cl$C[which(cl$C > k)] = k

col_grp = cl$cid[which(cl$C == grp)]
train = train[, colnames(train)%in%col_grp]
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
n = nrow(train)
prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
cl_sample = dpmm_beta(train, a, b, m0= m0, burn_in = 1500, sample_size=sample_size, thinning = 15)
diff = abs(mean(cl_sample$KK) - prior_m0_expect)

result = data.frame("grp"=grp, "m0" = m0, "Differ" = diff)

## Write results
#save(list = "cl_sample", file = paste('testid_', as.character(id), '_block_', as.character(grp), '.RData', sep = ""))
if(m0 < 10){
  write.table(result,file = paste("testID_", as.character(id), '_', as.character(grp), "_0", as.character(m0), ".txt", sep = ""),
              col.names = F,row.names = F)
} else{
  write.table(result,file = paste("testID_", as.character(id), '_', as.character(grp), "_", as.character(m0), ".txt", sep = ""),
              col.names = F,row.names = F)
}

