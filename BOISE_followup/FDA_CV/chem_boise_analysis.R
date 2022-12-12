setwd('~/RAwork/BOISE/BOISE_followup/FDA_CV/')
source('chem_evaluate.R')
set.seed(817)

nA = 30
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA), '.txt', sep=''), header = F)
chr_inform = unlist(lapply(informs$V2, function(x){return(paste(x, collapse = ','))}))
informs$V2 = chr_inform
informs$V3 = NULL
write.table(informs, file = paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA), '.txt',sep=''),
            row.names = F, col.names = F)

pre_informs = read.table(paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA-1),'.txt', sep=''), header = F)
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA), '.txt',sep=''), header = F)
for (id in 1:60) {
  pre_inform = pre_informs$V2[id]
  A = informs$V2[id]
  inform = paste(pre_inform, A, collapse = ',')
  informs$V2[id] = inform
}
write.table(informs, file = paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA), '.txt',sep=''),
            row.names = F, col.names = F)


# roc_results = data.frame(ids = 1:60)
# nef_results = data.frame(ids = 1:60)
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/Informer_', as.character(nA), '.txt',sep=''), header = F)
roc_results = read.table('./results/chem_roc_results.txt', sep = ' ', header=T)
nef_results = read.table('./results/chem_nef_results.txt', sep = ' ', header=T)

k=nA
row_sample_size = 100
roc_name = paste('roc_', as.character(k), sep = '')
nef_name = paste('nef_', as.character(k), sep = '')
roc_results[, roc_name] = rep(NA, 60)
nef_results[, nef_name] = rep(NA, 60)

for (id in 1:60) {
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  Scores = evaluate_interm(cl, inform, train, test, m0s, block, row_sample_size)
  xA = test[inform]
  Scores[inform[which(xA==1)]] = rep(1, sum(xA))
  Scores[inform[which(xA==0)]] = rep(0, sum(1-xA))
  ## ROCAUC
  Response = as.vector(test)
  Scores = as.vector(Scores)
  rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
  rocauc = rocobj$auc
  roc_results[id, roc_name] = rocauc
  ## NEF
  nT = as.integer(ncol(train) * 0.1)
  top = order(Scores,decreasing = T)[1:nT]
  pred_hit = sum(test[top])
  hit = sum(test)
  maxhit = min(hit,nT)
  nef10 = ((pred_hit/nT - hit/length(test)) / (maxhit/nT - hit/length(test)) + 1)/2
  nef_results[id, nef_name] = nef10
}

write.table(roc_results, file = './results/chem_roc_results.txt',row.names = F)
write.table(nef_results, file = './results/chem_nef_results.txt',row.names = F)
  
baseline_roc = read.table('./results/chem_rdinfo_roc_results.txt', header = T)
baseline_nef = read.table('./results/chem_rdinfo_nef_results.txt', header = T)
iter = 25
row_sample_size = 100
for (id in 1:60) {
  print(id)
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  ncpd = ncol(train)
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  nT = as.integer(ncol(train) * 0.1)
  for (nA in (11:20)*10) {
    roc_name = paste('roc_', as.character(nA), sep = '')
    nef_name = paste('nef_', as.character(nA), sep = '')
    roc = 0
    nef = 0
    for (i in 1:iter) {
      inform = sample(1:ncpd, nA)
      Scores = evaluate_interm(cl, inform, train, test, m0s, block, row_sample_size)
      xA = test[inform]
      Scores[inform[which(xA==1)]] = rep(1, sum(xA))
      Scores[inform[which(xA==0)]] = rep(0, sum(1-xA))
      ## ROCAUC
      Response = as.vector(test)
      Scores = as.vector(Scores)
      rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
      tmp_roc = rocobj$auc
      ## NEF10
      nT = as.integer(ncol(train) * 0.1)
      top = order(Scores,decreasing = T)[1:nT]
      pred_hit = sum(test[top])
      hit = sum(test)
      maxhit = min(hit,nT)
      tmp_nef = ((pred_hit/nT - hit/length(test)) / (maxhit/nT - hit/length(test)) + 1)/2
      roc = roc + tmp_roc
      nef = nef + tmp_nef
    }
    roc = roc / iter
    nef = nef / iter
    baseline_roc[id, roc_name] = roc
    baseline_nef[id, nef_name] = nef
  }
}
write.table(baseline_roc, file = './results/chem_rdinfo_roc_results.txt',row.names = F)
write.table(baseline_nef, file = './results/chem_rdinfo_nef_results.txt',row.names = F)


baseline_roc = rep(NA,60)
baseline_nef = rep(NA,60)
# for naive baseline
for (id in 1:60) {
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  Scores = apply(train, 2, function(x){return(mean(x, na.rm=T))})
  ## ROCAUC
  Response = as.vector(test)
  Scores = as.vector(Scores)
  rocobj = pROC::roc(response = Response, predictor = Scores, quiet = TRUE)
  rocauc = rocobj$auc
  baseline_roc[id] = rocauc
  ## NEF
  nT = as.integer(ncol(train) * 0.1)
  top = order(Scores,decreasing = T)[1:nT]
  pred_hit = sum(test[top])
  hit = sum(test)
  maxhit = min(hit,nT)
  nef10 = ((pred_hit/nT - hit/length(test)) / (maxhit/nT - hit/length(test)) + 1)/2
  baseline_nef[id] = nef10
}
baseline_roc = baseline_roc[-58]
baseline_nef = baseline_nef[-58]

### Plots
library(ggplot2)
results = data.frame('Informer_size' = (1:30))
cols = c(11, 21, 31:48)
roc_results = read.table('./results/chem_roc_results.txt', sep = ' ', header=T)
results[,'Boise_block'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
roc_results = read.table('./results/chem_fast_roc_results.txt', sep = ' ', header=T)
results[,'Boise_fast'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
roc_results = read.table('./results/chem_fast_info_1_roc_results.txt', sep = ' ', header=T)
results[,'fast_pel1'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[cols]
roc_results = read.table('./results/chem_rdinfo_roc_results.txt', sep = ' ', header=T)
results[,'Boise_rand'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
results[, 'Baseline'] = rep(mean(baseline_roc, na.rm=T), 30)

nef_results = read.table('./results/chem_nef_results.txt', sep = ' ', header=T)
results[,'Boise_block'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
nef_results = read.table('./results/chem_fast_nef_results.txt', sep = ' ', header=T)
results[,'Boise_fast'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
nef_results = read.table('./results/chem_fast_info_1_nef_results.txt', sep = ' ', header=T)
results[,'fast_pel1'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[cols]
nef_results = read.table('./results/chem_rdinfo_nef_results.txt', sep = ' ', header=T)
results[,'Boise_rand'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:31]
results[, 'Baseline'] = rep(mean(baseline_nef, na.rm=T), 30)

legend_title_size = 13
legend_text_size = 13
axis_title_size = 14
axis_text_size = 13
p4 <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = Boise_block))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_block, color = 'block_Boise'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_fast))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_fast, color = 'fast_Boise'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_rand))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_rand, color = 'rand_Boise'))+
  geom_line(mapping = aes(x = Informer_size, y = Baseline, color = 'Baseline'))+
  scale_color_manual(name = "Methods", values = c("block_Boise" = "darkblue",
                                                  'fast_Boise' = 'purple',
                                                  'rand_Boise' = 'orange',
                                                  'orig_Boise' = 'green',
                                                  "Baseline" = "red"))+
  scale_x_continuous('Informer sizes', breaks = seq(1,30, by=2))+
  #scale_y_continuous('ROCAUC mean', limits = c(0.789, 0.849))+
  scale_y_continuous('NEF10 mean', limits = c(0.673, 0.785))+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        #legend.position = c(0.85,0.2),
        legend.position = 'none')

p3 <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = Boise_rand))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_rand, color = 'rand_Boise'))+
  geom_point(mapping = aes(x = Informer_size, y = fast_entropy))+
  geom_line(mapping = aes(x = Informer_size, y = fast_entropy, color = 'fast_entropy'))+
  geom_point(mapping = aes(x = Informer_size, y = fast_pel1))+
  geom_line(mapping = aes(x = Informer_size, y = fast_pel1, color = 'fast_pel1'))+
  scale_color_manual(name = "Methods", values = c("rand_Boise" = "red",
                                                  'fast_entropy' = 'darkblue',
                                                  'fast_pel1' = 'purple'))+
  scale_x_continuous('Informer sizes', breaks = seq(10,200, by=10))+
  scale_y_continuous('ROCAUC mean', limits = c(0.806, 0.950))+
  #scale_y_continuous('NEF mean', limits = c(0.692, 0.918))+
  labs(title = "Clustering on separate matrices")+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        #legend.position = c(0.85,0.2),
        legend.position = c(0.85,0.2))
