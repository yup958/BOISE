setwd('~/RAwork/BOISE/BOISE_followup/FDA_CV/')
set.seed(817)
source('chem_evaluate.R')
sample_size = 100
baseline_roc = roc_results
baseline_nef = nef_results
chem_scores = read.table('~/CHTC_Downloads/FDA_cv/chem_scores.txt', header = F)
orig_scores = read.table('~/CHTC_Downloads/FDA_cv/orig_scores.txt', header = F)
colnames(chem_scores) = c('testid', 'chemid', 'score')
colnames(orig_scores) = c('testid', 'chemid', 'score')
## for original CRP clustering
for (id in 1:60) {
  load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
  # inform_scores = rep(0, ncol(train))
  # for (i in 1:sample_size) {
  #   for (j in 1:ncol(train)) {
  #     x = train[, j]
  #     complete_idx = which(!is.na(x))
  #     prev_entropy = entropy(cl_sample$CC[i, complete_idx])
  #     active_idx = which(x==1)
  #     inactive_idx = which(x==0)
  #     active_rate = mean(x, na.rm = T)
  #     post_entropy = 0
  #     if(length(active_idx) > 0){
  #       post_entropy = post_entropy + active_rate * entropy(cl_sample$CC[i, active_idx])
  #     }
  #     if(length(inactive_idx) > 0){
  #       post_entropy = post_entropy + (1-active_rate) * entropy(cl_sample$CC[i, inactive_idx])
  #     }
  #     mutual_info = prev_entropy - post_entropy
  #     inform_scores[j] = inform_scores[j] + mutual_info
  #   }
  # }
  inform_scores = orig_scores[which(orig_scores$testid == id),'score']
  nT = as.integer(ncol(train) * 0.1)
  for (nA in 1:20) {
    #inform = order(inform_scores, decreasing = T)[1:nA]
    inform = order(inform_scores)[1:nA]
    roc_name = paste('roc_', as.character(nA), sep = '')
    nef_name = paste('nef_', as.character(nA), sep = '')
    a = rep(mean(train, na.rm = T), ncol(train))
    b = 1 - a
    baseline_roc[id, roc_name] = Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
    baseline_nef[id, nef_name] = Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
  }
}
write.table(baseline_roc, file = 'fast_info_1_roc_results.txt',row.names = F)
write.table(baseline_nef, file = 'fast_nef_results.txt',row.names = F)

## for block clustering
max_block_size = 5
for (id in 1:60) {
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  # inform_scores = rep(0, ncol(train))
  # for (cpd in 1:ncol(train)) {
  #   grp = cl$C[cpd]
  #   sub_cols = which(cl$C == grp) 
  #   sub_train = train[ , sub_cols]
  #   for (i in 1:sample_size) {
  #     x = train[, cpd]
  #     complete_idx = which(!is.na(x))
  #     prev_entropy = entropy(block[[grp]][i, complete_idx])
  #     active_idx = which(x==1)
  #     inactive_idx = which(x==0)
  #     active_rate = mean(x, na.rm = T)
  #     post_entropy = 0
  #     if(length(active_idx) > 0){
  #       post_entropy = post_entropy + active_rate * entropy(block[[grp]][i, active_idx])
  #     }
  #     if(length(inactive_idx) > 0){
  #       post_entropy = post_entropy + (1-active_rate) * entropy(block[[grp]][i, inactive_idx])
  #     }
  #     mutual_info = prev_entropy - post_entropy
  #     inform_scores[cpd] = inform_scores[cpd] + mutual_info
  #   }
  # }
  inform_scores = chem_scores[which(chem_scores$testid == id),'score']
  # Evaluate
  nT = as.integer(ncol(train) * 0.1)
  inform = c()
  sub_grps = c()
  idx = 1
  for (nA in 1:30) {
    #inform = order(inform_scores, decreasing = T)[1:nA]
    # inform = order(inform_scores)[1:nA]
    tmp_inform = order(inform_scores)[idx]
    tmp_grp = cl$C[tmp_inform]
    while(length(which(sub_grps == tmp_grp)) >= max_block_size){
      idx = idx + 1
      tmp_inform = order(inform_scores)[idx]
      tmp_grp = cl$C[tmp_inform]
    }
    idx = idx + 1
    inform = c(inform, tmp_inform)
    sub_grps = c(sub_grps, tmp_grp)
    roc_name = paste('roc_', as.character(nA), sep = '')
    nef_name = paste('nef_', as.character(nA), sep = '')
    Scores = evaluate_interm(cl, inform, train, test, m0s, block, sample_size)
    xA = test[inform]
    Scores[inform[which(xA==1)]] = rep(1, sum(xA))
    Scores[inform[which(xA==0)]] = rep(0, sum(1-xA))
    ## ROCAUC
    Response = as.vector(test)
    Scores = as.vector(Scores)
    rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
    rocauc = rocobj$auc
    baseline_roc[id, roc_name] = rocauc
    ## NEF
    nT = as.integer(ncol(train) * 0.1)
    top = order(Scores,decreasing = T)[1:nT]
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = min(hit,nT)
    nef10 = ((pred_hit/nT - hit/length(test)) / (maxhit/nT - hit/length(test)) + 1)/2
    baseline_nef[id, nef_name] = nef10  
    }
}
write.table(baseline_roc, file = './results/chem_fast_info_1_roc_results.txt',row.names = F)
write.table(baseline_nef, file = './results/chem_fast_info_1_nef_results.txt',row.names = F)

