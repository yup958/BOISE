setwd("~/RAwork/BOISE/BOISE_followup/FDA_prospective/")
set.seed(817)
load('fda_data.RData')
complete_indicator = read.csv('all_FDA_assays_cln_smiles_pccids.csv', header = T)
rownames(complete_indicator) = complete_indicator$PUBCHEM_CID
complete_indicator = complete_indicator[5:ncol(complete_indicator)]
### targets with different sets of compounds
# columns 1,2,3,4,8,9,10,11,12,13,14,15,16,18,19,20,21,23,24,25,26,27. Size 688 x 937
# columns 5, 6. Size 688 x 458
# columns 7, 22. Size 688 x 936
# column 17. Size 688 x 308
valid_cols = c(1, 5, 7, 17) 
valid_col = valid_cols[id]
complete_cpd = rownames(complete_indicator)[which(complete_indicator[,valid_col] == 1)]
complete_cpd = sapply(complete_cpd, function(x){return(paste('CID', as.character(x), sep = '_'))})
train = dat[, which(colnames(dat) %in% complete_cpd)]
dim(train)

## rearrange informers for original BOISE
nA = 20
informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_L2_informer_', as.character(nA), '.txt', sep=''), header = F)
chr_inform = unlist(lapply(informs$V2, function(x){return(paste(x, collapse = ','))}))
informs$V2 = chr_inform
informs$V3 = NULL
write.table(informs, file = paste('~/CHTC_Downloads/FDA_prospective/pros_orig_L2_informer_', as.character(nA), '.txt', sep=''),
            row.names = F, col.names = F)

pre_informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_L2_informer_', as.character(nA-1), '.txt', sep=''), header = F)
informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_L2_informer_', as.character(nA), '.txt', sep=''), header = F)
for (id in 1:4) {
  pre_inform = pre_informs$V2[id]
  A = informs$V2[id]
  inform = paste(pre_inform, A, collapse = ',')
  informs$V2[id] = inform
}
write.table(informs, file = paste('~/CHTC_Downloads/FDA_prospective/pros_orig_L2_informer_', as.character(nA), '.txt', sep=''),
            row.names = F, col.names = F)

## rearrange informers for block BOISE
nA = 30
informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_', as.character(nA), '.txt', sep=''), header = F)
chr_inform = unlist(lapply(informs$V2, function(x){return(paste(x, collapse = ','))}))
informs$V2 = chr_inform
informs$V3 = NULL
write.table(informs, file = paste('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_', as.character(nA), '.txt', sep=''),
            row.names = F, col.names = F)

pre_informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_', as.character(nA-1), '.txt', sep=''), header = F)
informs = read.table(paste('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_', as.character(nA), '.txt', sep=''), header = F)
for (id in 1:4) {
  pre_inform = pre_informs$V2[id]
  A = informs$V2[id]
  inform = paste(pre_inform, A, collapse = ',')
  informs$V2[id] = inform
}
write.table(informs, file = paste('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_', as.character(nA), '.txt', sep=''),
            row.names = F, col.names = F)

## fast BOISE with original clustering
library(entropy)
sample_size = 250
nA=30
informs = data.frame(V1 = 1:4, V2 = rep('',4))
system.time({
  for (id in 1:4) {
    load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(id), '.RData', sep = ''))
    inform_scores = rep(0, ncol(train))
    for (i in 1:sample_size) {
      for (j in 1:ncol(train)) {
        x = train[, j]
        complete_idx = which(!is.na(x))
        prev_entropy = entropy(cl_sample$CC[i, complete_idx])
        active_idx = which(x==1)
        inactive_idx = which(x==0)
        active_rate = mean(x, na.rm = T)
        post_entropy = 0
        if(length(active_idx) > 0){
          post_entropy = post_entropy + active_rate * entropy(cl_sample$CC[i, active_idx])
        }
        if(length(inactive_idx) > 0){
          post_entropy = post_entropy + (1-active_rate) * entropy(cl_sample$CC[i, inactive_idx])
        }
        mutual_info = prev_entropy - post_entropy
        inform_scores[j] = inform_scores[j] + mutual_info
      }
    }
    inform = order(inform_scores, decreasing = T)[1:nA]
    chr_inform = paste(inform, collapse = ' ')
    informs$V2[id] = chr_inform
  }
})
write.table(informs, file = paste('~/CHTC_Downloads/FDA_prospective/pros_fast_informer_', as.character(nA), '.txt', sep=''),
            row.names = F, col.names = F)


## Write a matrix for informers x test targets
informs = read.table('~/CHTC_Downloads/FDA_prospective/pros_fast_informer_30.txt', header = F)
inform_cids = c()
for (id in 1:4) {
  load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(id), '.RData', sep = ''))
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
  inform_cid = colnames(train)[inform]
  inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
  inform_cids = union(inform_cids, inform_cid)
}
interm_mat = complete_indicator[which(rownames(complete_indicator) %in% inform_cids), ]
grps = rep(1, 27)
grps[c(5, 6)] = 2
grps[c(7,22)] = 3
grps[17] = 4
for (col_id in 1:27) {
  interm_mat[,col_id] = NA
  grp_id = grps[col_id]
  load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(grp_id), '.RData', sep = ''))
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[grp_id]), split = ' ')))
  inform_cid = colnames(train)[inform]
  inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
  interm_mat[which(rownames(interm_mat)%in% inform_cid), col_id] = 1
}
interm_mat$PUBCHEM_CID = rownames(interm_mat)
write.csv(interm_mat, file = 'fast_Boise_interm_mat.csv', row.names = F)

## ranking
informs = read.table('~/CHTC_Downloads/FDA_prospective/pros_fast_informer_30.txt', header = F)
interm_mat = read.csv('fast_Boise_interm_mat.csv', header = T)
sample_size = 250
final_mat = complete_indicator
for (col_id in 1:27) {
  grp_id = grps[col_id]
  load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(grp_id), '.RData', sep = ''))
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[grp_id]), split = ' ')))
  inform_cid = colnames(train)[inform]
  inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
  xA = c()
  for (cid in inform_cid) {
    tmp_res = interm_mat[which(interm_mat$PUBCHEM_CID == cid), col_id+1]
    print(tmp_res)
    xA = c(xA, tmp_res)
  }
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  ## inner function for Evaluate
  P = clust_sum(cl_sample,train,sample_size, a, b)
  Score = rep(0, ncol(train))
  nA=length(inform)
  m = ncol(train)
  log_post_probs = matrix(0, 1, sample_size)
  post_thetas = matrix(0, sample_size, m)
  for (k in 1:sample_size){
    postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, a, b, m0)
    log_post_probs[k] = postls$log_post_prob
    post_thetas[k, ] = postls$post_theta
  }
  log_post_probs = log_post_probs - max(log_post_probs)
  post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
  Score = post_probs %*% post_thetas
  ## fill Score into final mat
  final_mat[,col_id] = NA
  for (row in 1:nrow(final_mat)) {
    cpd = paste('CID', rownames(final_mat)[row], sep = '_')
    if(cpd %in% colnames(train)){
      final_mat[row, col_id] = Score[which(colnames(train) == cpd)]
    }
  }
}

