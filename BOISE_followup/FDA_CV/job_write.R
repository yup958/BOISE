### for m0 find
test_ids = read.table('Test_IDS.txt')
load('fda_data.RData')
total_k = 0
grps = c()
ids = c()
for (id in 31:60) {
  test_id = test_ids[id,1]
  test = dat[rownames(dat)==test_id,]
  complete_idx = which(!is.na(test))
  train = dat[!row.names(dat)==test_id, complete_idx]
  test = test[complete_idx]
  cids = colnames(train)
  cl = read.csv('chemical_clustering_res.csv', header = T)
  cl = cl[cl$cid %in% cids,]
  cl$C = NA
  tab = table(cl$clust)
  ord = order(tab, decreasing = T)
  for (k in 1:length(tab)) {
    cl$C[which(cl$clust == names(tab)[ord[k]])] = k
  }
  cl$clust = NULL
  tab = table(cl$C)
  for (K in 1:length(tab)){
    if(tab[K] < 10)
      break
  }
  total_k = total_k + K
  for (k in 1:K) {
    grps = c(grps, rep(k,25))
  }
  ids = c(ids, rep(id,25*K))
}

A = data.frame(id = ids, grp = grps)
A$m0 = rep(1:25, total_k)
write.table(A,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)

## rearrangement blocks
m0_selections = read.table('~/CHTC_Downloads/FDA_cv/prior_mass_chem.txt', header = F)
for (id in 31:60) {
  test_id = test_ids[id,1]
  test = dat[rownames(dat)==test_id,]
  complete_idx = which(!is.na(test))
  train = dat[!row.names(dat)==test_id, complete_idx]
  test = test[complete_idx]
  cids = colnames(train)
  cl = read.csv('chemical_clustering_res.csv', header = T)
  cl = cl[cl$cid %in% cids,]
  cl$C = NA
  tab = table(cl$clust)
  ord = order(tab, decreasing = T)
  for (k in 1:length(tab)) {
    cl$C[which(cl$clust == names(tab)[ord[k]])] = k
  }
  cl$clust = NULL
  tab = table(cl$C)
  for (K in 1:length(tab)){
    if(tab[K] < 10)
      break
  }
  cl$C[which(cl$C > K)] = K # total number of groups of compounds
  ## load m0s
  m0_selection = m0_selections[which(m0_selections$V1==id),]
  m0s = m0_selection$V3
  
  ## load blocks
  block = list()
  for (k in 1:K) {
    load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block_',
               as.character(k), '.RData', sep = ''))
    block[[k]] = cl_sample$CC
  }
  save(list = c('block', 'cl', 'train', 'test', 'm0s', 'test_id'), file = paste('testid_', as.character(id),
                                                                         '_block.RData', sep=''))
}

### for original boise m0 find
ids = c()
m0 = c()
for (id in 31:60) {
  ids = c(ids, rep(id, 30))
  m0 = c(m0, 1:30)
}
A = data.frame(id = ids, m0 = m0)
write.table(A,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)

### for informer search
test_ids = read.table('Test_IDS.txt')
informs = read.table('~/CHTC_Downloads/FDA_cv/orig_new_informer_11.txt')
load('fda_data.RData')
candidates = c()
ids = c()
count = 0
for  (id in 21:40) {
  test_id = test_ids[id,1]
  test = dat[rownames(dat)==test_id,]
  complete_idx = which(!is.na(test))
  count = count + length(complete_idx)
  train = dat[!row.names(dat)==test_id, complete_idx]
  test = test[complete_idx]
  cids = colnames(train)
  pre_inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
  #pre_inform=c()
  cand_num = ncol(train) 
  ids = c(ids, rep(id,cand_num - length(pre_inform)))
  candidates= c(candidates, setdiff(1:cand_num, pre_inform))
}
A = data.frame(id = ids, candidate = candidates)
write.table(A,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)

## rearrangement original clustering
m0_selections = read.table('~/CHTC_Downloads/FDA_cv/prior_mass_orig.txt', header = F)
for (id in 31:60) {
  test_id = test_ids[id,1]
  test = dat[rownames(dat)==test_id,]
  complete_idx = which(!is.na(test))
  train = dat[!row.names(dat)==test_id, complete_idx]
  test = test[complete_idx]
  cids = colnames(train)
  ## load m0
  m0_selection = m0_selections[which(m0_selections$V1==id),]
  m0 = m0_selection$V2
  
  ## load cl_sample
  load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
  save(list = c('cl_sample', 'train', 'test', 'm0', 'test_id'), 
       file = paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
}

### PCBA informer search
set.seed(817)
test_ids = sample(1:102)
id = test_ids[15]
print(id)
load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
print(ncol(train))
candidates = 30001:ncol(train)
ids = rep(id, length(candidates))
A = data.frame(id = ids, candidate = candidates)
write.table(A,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)

