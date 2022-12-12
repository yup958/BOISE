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
informs = read.table('~/CHTC_Downloads/FDA_prospective/pros_chem_L1_informer_30.txt', header = F)
sample_size = 250
final_mat = read.csv('all_FDA_assays_cln_smiles_pccids_binarylabels_peng_format.csv', header = T)
rownames(final_mat) = final_mat$PUBCHEM_CID
final_mat = final_mat[5:ncol(final_mat)]
grps = rep(1, 27)
grps[c(5, 6)] = 2
grps[c(7,22)] = 3
grps[17] = 4
pros_nef_results = read.table('pros_nef_results.txt', header = T)
pros_roc_results = read.table('pros_roc_results.txt', header = T)
recovery_counts = read.table('hits_recovered_results.txt', header = T)
# pros_nef_results = data.frame(assays = colnames(final_mat), orig_Boise_L1 = rep(NA, 27), 
#                               orig_Boise_L2 = rep(NA, 27), block_Boise = rep(NA, 27),
#                               fast_Boise = rep(NA, 27))
# pros_roc_results = data.frame(assays = colnames(final_mat), orig_Boise_L1 = rep(NA, 27), 
#                               orig_Boise_L2 = rep(NA, 27), block_Boise = rep(NA, 27),
#                               fast_Boise = rep(NA, 27))
# recovery_counts = data.frame(assays = colnames(final_mat), total_hit = rep(NA, 27),
#                              orig_Boise_L1 = rep(NA, 27), orig_Boise_L2 = rep(NA, 27),
#                              block_Boise = rep(NA, 27), fast_Boise = rep(NA, 27),
#                              rand_Boise = rep(NA, 27))
## original BOISE evaluation
for (col_id in 1:27) {
  if (sum(final_mat[, col_id], na.rm = T) == 0){
    print(paste("Column ", as.character(col_id), " is NA", sep = '')) ## not valid assay
    next
  }
  grp_id = grps[col_id]
  load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(grp_id), '.RData', sep = ''))
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[grp_id]), split = ' ')))
  inform_cid = colnames(train)[inform]
  inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
  xA = c()
  for (cid in inform_cid) {
    tmp_res = final_mat[cid, col_id]
    #print(tmp_res)
    xA = c(xA, tmp_res)
  }
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  ## fill Score into final mat
  train_cids = colnames(train)
  train_cids = sapply(train_cids, function(s){return(unlist(strsplit(s,'_'))[2])})
  test = final_mat[train_cids, col_id]
  nT = round(0.1 * ncol(train))
  pros_nef_results[col_id, 'fast_Boise'] = Evaluate(cl_sample, inform, 'nef', 0.1,
                                                       test,train,nT,sample_size, a,b,m0)
  pros_roc_results[col_id, 'fast_Boise'] = Evaluate(cl_sample, inform, 'rocauc', 0.1,
                                                       test,train,nT,sample_size, a,b,m0)
  # nTop = 100
  # recovery_counts$total_hit[col_id] = sum(test)
  # ## inner function for Evaluate
  # P = clust_sum(cl_sample,train,sample_size, a, b)
  # Score = rep(0, ncol(train))
  # nA=length(inform)
  # m = ncol(train)
  # log_post_probs = matrix(0, 1, sample_size)
  # post_thetas = matrix(0, sample_size, m)
  # for (k in 1:sample_size){
  #   postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, a, b, m0)
  #   log_post_probs[k] = postls$log_post_prob
  #   post_thetas[k, ] = postls$post_theta
  # }
  # log_post_probs = log_post_probs - max(log_post_probs)
  # post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
  # Score = post_probs %*% post_thetas
  # top = order(Score,decreasing = T)[1:nTop]
  # pred_hit = sum(test[top])
  # recovery_counts[col_id, 'fast_Boise'] = pred_hit
}

source('~/RAwork/BOISE/BOISE_followup/FDA_CV/chem_evaluate.R')
row_sample_size = 250
for (col_id in 1:27) {
  if (sum(final_mat[, col_id], na.rm = T) == 0){
    print(paste("Column ", as.character(col_id), " is NA", sep = '')) ## not valid assay
    next
  }
  grp_id = grps[col_id]
  load(paste('~/CHTC_Downloads/FDA_prospective/testid_', as.character(grp_id), '_block.RData', sep = ''))
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[grp_id]), split = ' ')))
  inform_cid = colnames(train)[inform]
  inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
  ## fill Score into final mat
  train_cids = colnames(train)
  train_cids = sapply(train_cids, function(s){return(unlist(strsplit(s,'_'))[2])})
  test = final_mat[train_cids, col_id]
  Scores = evaluate_interm(cl, inform, train, test, m0s, block, row_sample_size)
  # xA = test[inform]
  # ## ROCAUC
  # Response = as.vector(test)
  # Scores = as.vector(Scores)
  # rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
  # rocauc = rocobj$auc
  # pros_roc_results[col_id, 'block_Boise'] = rocauc
  # ## NEF
  # nT = as.integer(ncol(train) * 0.1)
  # top = order(Scores,decreasing = T)[1:nT]
  # pred_hit = sum(test[top])
  # hit = sum(test)
  # maxhit = min(hit,nT)
  # nef10 = ((pred_hit/nT - hit/length(test)) / (maxhit/nT - hit/length(test)) + 1)/2
  # pros_nef_results[col_id, 'block_Boise'] = nef10
  nTop = 100
  recovery_counts$total_hit[col_id] = sum(test)
  top = order(Scores,decreasing = T)[1:nTop]
  pred_hit = sum(test[top])
  recovery_counts[col_id, 'block_Boise'] = pred_hit
}
# write.table(pros_nef_results, file = 'pros_nef_results.txt',
#             row.names = F)
# write.table(pros_roc_results, file = 'pros_roc_results.txt',
#             row.names = F)
write.table(recovery_counts, file = 'hits_recovered_results.txt',
            row.names = F)

## rand BOISE
set.seed(817)
pros_nef_results$rand_Boise = rep(NA, 27)
pros_roc_results$rand_Boise = rep(NA, 27)
iter = 25
nA=30
nTop = 100
for (col_id in 1:27) {
  if (sum(final_mat[, col_id], na.rm = T) == 0){
    print(paste("Column ", as.character(col_id), " is NA", sep = '')) ## not valid assay
    next
  }
  grp_id = grps[col_id]
  load(paste('~/CHTC_Downloads/FDA_prospective/pros_orig_clust_', as.character(grp_id), '.RData', sep = ''))
  # roc = 0
  # nef = 0
  hits = 0
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  ## fill Score into final mat
  train_cids = colnames(train)
  train_cids = sapply(train_cids, function(s){return(unlist(strsplit(s,'_'))[2])})
  test = final_mat[train_cids, col_id]
  nT = round(0.1 * ncol(train))
  for (i in 1:iter) {
    inform = sample(1:ncol(train), nA)
    inform_cid = colnames(train)[inform]
    inform_cid = sapply(inform_cid, function(s){return(unlist(strsplit(s,'_'))[2])})
    xA = c()
    for (cid in inform_cid) {
      tmp_res = final_mat[cid, col_id]
      #print(tmp_res)
      xA = c(xA, tmp_res)
    }
    # roc = roc + Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
    # nef = nef + Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
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
    top = order(Score,decreasing = T)[1:nTop]
    pred_hit = sum(test[top])
    hits = hits + pred_hit
  }
  # roc = roc / iter
  # nef = nef / iter
  # pros_nef_results[col_id, 'rand_Boise'] = nef
  # pros_roc_results[col_id, 'rand_Boise'] = roc
  hits = hits / iter
  recovery_counts$rand_Boise[col_id] = hits
}

### Plot
pros_nef_results = read.table('pros_nef_results.txt', header = T)
pros_roc_results = read.table('pros_roc_results.txt', header = T)
pros_nef_results = pros_nef_results[which(!is.na(pros_nef_results$orig_Boise_L1)), ]
pros_roc_results = pros_roc_results[which(!is.na(pros_roc_results$orig_Boise_L1)), ]
m=18
temp = data.frame("NEF10" = rep(0, 5 * m), "IBR" = c(rep("rand_Boise",m),rep("fast_Boise",m),
                                                     rep("Boise_L2",m),rep("block_Boise",m), rep("Boise_L1",m)))
temp$NEF10 = c(pros_nef_results$rand_Boise, pros_nef_results$fast_Boise,pros_nef_results$orig_Boise_L2,
               pros_nef_results$block_Boise, pros_nef_results$orig_Boise_L1)
temp$IBR = as.factor(temp$IBR)
temp$IBR <- factor(temp$IBR, levels =c('rand_Boise', 'fast_Boise', 'Boise_L2',
                                       'block_Boise', 'Boise_L1'))

legend_title_size = 13
legend_text_size = 13
axis_title_size = 14
axis_text_size = 13
p1 <- ggplot(temp, aes(x = IBR, y = NEF10, fill = IBR)) +  
  geom_violin( trim = F, adjust = 0.5 )+
  geom_boxplot(width = 0.05,aes(fill = "white")) + 
  scale_y_continuous(limits = c(0.55,1))+
  scale_fill_manual(values = c( "blue","green","lightblue","purple", "orange", "#000000"))+
  stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  xlab('BOISE variants') +
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

temp = data.frame("ROCAUC" = rep(0, 5 * m), "IBR" = c(rep("rand_Boise",m),rep("fast_Boise",m),
                                                     rep("Boise_L2",m),rep("block_Boise",m), rep("Boise_L1",m)))
temp$ROCAUC = c(pros_roc_results$rand_Boise, pros_roc_results$fast_Boise,pros_roc_results$orig_Boise_L2,
               pros_roc_results$block_Boise, pros_roc_results$orig_Boise_L1)
temp$IBR = as.factor(temp$IBR)
temp$IBR <- factor(temp$IBR, levels =c('rand_Boise', 'fast_Boise', 'Boise_L2',
                                       'block_Boise', 'Boise_L1'))
p2 <- ggplot(temp, aes(x = IBR, y = ROCAUC, fill = IBR)) +  
  geom_violin( trim = F, adjust = 0.5 )+
  geom_boxplot(width = 0.05,aes(fill = "white")) + 
  scale_y_continuous(limits = c(0.55,1))+
  scale_fill_manual(values = c( "blue","green","lightblue","purple", "orange", "#000000"))+
  stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  xlab('BOISE variants') +
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
myGrobs <- list(p1,p2)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 1,ncol = 2)


