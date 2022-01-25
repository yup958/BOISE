setwd('~/RAwork/BOISE/BOISE_followup/FDA_analysis/')
source('evaluate_interm.R')
set.seed(2)
load('col_clust_res_5.RData')
m0_selections = read.table('~/CHTC_Downloads/Prior_mass_FDA_separate.txt')

## train-test separation
load('fda_data_rearranged.RData')
block_res$row_block = block_res$row_block[1:5,]
block_res$col_block = block_res$col_block[1:5,]
active_ind = apply(block_res$col_block, 2, sum)
active_ind = which(active_ind > 0)
cpds = colnames(block_res$col_block[, active_ind])
dat = dat[, cpds]

count = apply(dat,1, function(x){return(sum(is.na(x)))})
test_AID = names(sort(count)[1:100])
train_AID = setdiff(rownames(dat), test_AID)
train = dat[train_AID, ]
test = dat[test_AID,]
rm(dat)

col_sample_size = 100
row_sample_size = 100

inform = c(91, 10, 747, 85, 137, 525, 62, 189, 297, 470, 638, 82, 46, 713, 59)
inform_CID = colnames(train)[inform]
nT = as.integer(ncol(train) * 0.1)

## evaluate function of the top set
Evaluate <-function(cl_sample, inform, measure, percent,
           test, train, nT, col_sample_size, row_sample_size, m0_selections){
  Scores = rep(0, ncol(train))
  complete_idx = which(!is.na(test))
  xA = test[inform]
  nA = length(inform)
  for (k in 1:col_sample_size) {
    # load block RData file
    if(k == 8){
      next
    }
    load(paste('~/CHTC_Downloads/block_res_', as.character(k), '.RData', sep = ''))
    # specific m0 priors
    m0_selection = m0_selections[which(m0_selections$V1==k),]
    ## cluster assignment. Merge all column clusters with size <= 2
    cl = cl_sample$CC[k,]
    for (grp in 1:cl_sample$KK[k]){
      if(cl_sample$NN[k, grp] <= 2)
        break
    }
    cl[which(cl > grp)] = grp
    Scores = Scores + evaluate_interm(cl, inform, train, test,
                                      m0_selection, block, row_sample_size)
  }
  Scores = Scores / (col_sample_size-1)
  
  Scores[inform[which(xA==1)]] = rep(max(Scores) + 1, sum(xA, na.rm = T))
  Scores[inform[which(xA==0)]] = rep(min(Scores) - 1, sum(1-xA, na.rm = T))
  return(Scores)
  # test = as.vector(test[complete_idx])
  # Scores = as.vector(Scores[complete_idx])
  # if(measure == "nef"){
  #   if(sum(test, na.rm = T) == 0){
  #     print('All observed responses are 0!')
  #     result = NA
  #   } else{
  #     nTop = round(ncol(train) * percent, 0)
  #     top = order(Scores,decreasing = T)[1:nTop]
  #     pred_hit = sum(test[top])
  #     hit = sum(test)
  #     maxhit = min(hit,nTop)
  #     result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  #   }
  # } else if(measure == 'selectivity_nef'){
  #   nTop = round(ncol(train) * percent, 0)
  #   top = order(Scores,decreasing = T)[1:nTop]
  #   selectivity_score = colMeans(train, na.rm = T)
  #   test = test - selectivity_score
  #   pred_hit = sum(test[top])
  #   hit = sum(test)
  #   maxhit = sum(sort(test, decreasing = T)[1:nTop])
  #   result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  # } else if(measure == "rocauc"){
  #   if(sum(test, na.rm = T) == 0){
  #     print('All observed responses are 0!')
  #     result = NA
  #   } else{
  #     rocobj = pROC::roc(response = test, predictor = Scores, quiet = TRUE)
  #     result = rocobj$auc
  #   }
  # } else{
  #   print("Criteria is not supported.")
  #   result = 0
  # }
  # return(result)
}

roc_results = read.table('block_roc_results.txt', sep = ' ', header=T)
nA_results = read.table('block_nA_results.txt', sep = ' ', header=T)
#nef_results = read.table('nef_results.txt', sep = ' ', header=T)
cache = list()

for (k in 15:15) {
  nA_name = paste('nA_', as.character(k), sep = '')
  roc_name = paste('rocauc_', as.character(k), sep = '')
  #nef_name = paste('nef_', as.character(k), sep = '')
  roc_results[, roc_name] = rep(NA, 100)
  #nef_results[, nef_name] = rep(NA, 100)
  nA_results[, nA_name] = rep(0, 100)
  tmp_inform = inform[1:k]
  for (i in 1:100) {
    print(i)
    if(sum(test[i, ], na.rm = T) == 0){
      print('All observed responses are 0!')
      next
    }
    xa = paste(test[i, tmp_inform], collapse = '')
    complete_idx = which(!is.na(test[i,]))
    valid_inform = intersect(complete_idx, tmp_inform)
    nA_results[i, nA_name] = length(valid_inform)
    if (length(valid_inform) >= 1){
      if(xa %in% names(cache)){
        Scores = cache[[xa]]
      } else{
        Scores = Evaluate(cl_sample, tmp_inform, measure='rocauc', percent=0.1,
                          test=test[i,], train, nT, col_sample_size, row_sample_size, m0_selections)
        cache[[xa]] = Scores
      }
      Response = as.vector(test[i, complete_idx])
      Scores = as.vector(Scores[complete_idx])
      rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
      roc_results[i, roc_name] = rocobj$auc
    }
  }
  print(summary(roc_results[, roc_name]))
}

write.table(nA_results, file = 'block_nA_results.txt',row.names = F)
write.table(roc_results, file = 'block_roc_results.txt',row.names = F)
write.table(nef_results, file = 'block_nef_results.txt',row.names = F)

### Plots
library(ggplot2)
results = data.frame('Informer_size' = 3:15)
roc_results = read.table('block_roc_results.txt', sep = ' ', header=T)
results[,'Boise_block'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[3:15]
roc_results = read.table('roc_results.txt', sep = ' ', header=T)
results[,'Boise_original'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[4:16]
results[, 'Baseline'] = rep(mean(baseline_roc, na.rm=T), 13)
legend_title_size = 13
legend_text_size = 13
axis_title_size = 16
axis_text_size = 15
p <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = Boise_block))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_block, color = 'block_Boise'))+
  geom_line(mapping = aes(x = Informer_size, y = Baseline, color = 'Baseline'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_original))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_original, color = 'original_Boise'))+
  scale_color_manual(name = "Methods", values = c("block_Boise" = "darkblue", 
                                                  'original_Boise' = 'green',
                                                  "Baseline" = "red"))+
  scale_x_continuous('Informer sizes', breaks = seq(3,15, by=1))+
  scale_y_continuous('ROCAUC mean', limits = c(0.82, 0.86))+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.position = c(0.1,0.8))
p
