setwd('~/RAwork/BOISE/BOISE_followup/FDA_CV/')
set.seed(817)

nA = 15
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA), '.txt', sep=''), header = F)
chr_inform = unlist(lapply(informs$V2, function(x){return(paste(x, collapse = ','))}))
informs$V2 = chr_inform
informs$V3 = NULL
write.table(informs, file = paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA), '.txt',sep=''),
            row.names = F, col.names = F)

pre_informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA-1),'.txt', sep=''), header = F)
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA), '.txt',sep=''), header = F)
for (id in 1:60) {
  pre_inform = pre_informs$V2[id]
  A = informs$V2[id]
  inform = paste(pre_inform, A, collapse = ',')
  informs$V2[id] = inform
}
write.table(informs, file = paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA), '.txt',sep=''),
            row.names = F, col.names = F)

# roc_results = data.frame(ids = 1:60)
# nef_results = data.frame(ids = 1:60)
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(nA), '.txt',sep=''), header = F)
roc_results = read.table('./results/orig_L0_roc_results.txt', sep = ' ', header=T)
nef_results = read.table('./results/orig_L0_nef_results.txt', sep = ' ', header=T)
chem_roc_results = read.table('./results/chem_roc_results.txt', sep = ' ', header=T)
chem_nef_results = read.table('./results/chem_nef_results.txt', sep = ' ', header=T)
baseline_roc = read.table('./results/orig_L1_roc_results.txt', sep = ' ', header=T)
baseline_nef = read.table('./results/orig_L1_nef_results.txt', sep = ' ', header=T)

k=nA
sample_size = 100
roc_name = paste('roc_', as.character(k), sep = '')
nef_name = paste('nef_', as.character(k), sep = '')
roc_results[, roc_name] = rep(NA, 60)
nef_results[, nef_name] = rep(NA, 60)
#pr_results = rep(NA,60)

for (id in 1:60) {
  print(id)
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
  load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  nT = as.integer(ncol(train) * 0.1)
  rocauc = Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
  roc_results[id, roc_name] = rocauc
  nef10 = Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
  nef_results[id, nef_name] = nef10
  # pr_results[id] = Evaluate(cl_sample, inform, 'prauc', 0.1, test, train, nT,sample_size,a,b,m0)
}
rocauc
nef10
roc_results$roc_17[id]
nef_results$nef_17[id]
baseline_roc$roc_18[id]
baseline_nef$nef_18[id]

write.table(roc_results, file = './results/orig_L0_roc_results.txt',row.names = F)
write.table(nef_results, file = './results/orig_L0_nef_results.txt',row.names = F)

## other nef thresholds
sample_size = 100
# old_nef20_results = data.frame(ids = 1:60)
# old_nef30_results = data.frame(ids = 1:60)
# new_nef20_results = data.frame(ids = 1:60)
# new_nef30_results = data.frame(ids = 1:60)
load('./results/nef_2030_l012.RData')

for (nA in 1:20) {
  informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_informer_', as.character(nA), '.txt',sep=''), header = F)
  nef_name = paste('nef_', as.character(nA), sep = '')
  old_nef20_results[, nef_name] = rep(NA, 60)
  old_nef30_results[, nef_name] = rep(NA, 60)
  for (id in 1:60) {
    inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
    load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
    a = rep(mean(train, na.rm = T), ncol(train))
    b = 1 - a
    nT = as.integer(ncol(train) * 0.2)
    nef20 = Evaluate(cl_sample, inform, 'nef', 0.2, test, train, nT,sample_size,a,b,m0)
    old_nef20_results[id, nef_name] = nef20
    nT = as.integer(ncol(train) * 0.3)
    nef30 = Evaluate(cl_sample, inform, 'nef', 0.3, test, train, nT,sample_size,a,b,m0)
    old_nef30_results[id, nef_name] = nef30
  }
}

informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_new_informer_', as.character(19), '.txt',sep=''), header = F)
for (nA in 19:19) {
  print(nA)
  nef_name = paste('nef_', as.character(nA), sep = '')
  l2_nef20_results[, nef_name] = rep(NA, 60)
  l2_nef30_results[, nef_name] = rep(NA, 60)
  for (id in 1:60) {
    inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
    inform = inform[1:nA]
    load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
    a = rep(mean(train, na.rm = T), ncol(train))
    b = 1 - a
    nT = as.integer(ncol(train) * 0.2)
    nef20 = Evaluate(cl_sample, inform, 'nef', 0.2, test, train, nT,sample_size,a,b,m0)
    l2_nef20_results[id, nef_name] = nef20
    nT = as.integer(ncol(train) * 0.3)
    nef30 = Evaluate(cl_sample, inform, 'nef', 0.3, test, train, nT,sample_size,a,b,m0)
    l2_nef30_results[id, nef_name] = nef30
  }
}
l0_nef20_results = data.frame(ids = 1:60)
l0_nef30_results = data.frame(ids = 1:60)
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_L0_informer_', as.character(14), '.txt',sep=''), header = F)
for (nA in 15:15) {
  print(nA)
  nef_name = paste('nef_', as.character(nA), sep = '')
  l0_nef20_results[, nef_name] = rep(NA, 60)
  l0_nef30_results[, nef_name] = rep(NA, 60)
  for (id in 1:60) {
    inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
    inform = inform[1:nA]
    load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
    a = rep(mean(train, na.rm = T), ncol(train))
    b = 1 - a
    nT = as.integer(ncol(train) * 0.2)
    nef20 = Evaluate(cl_sample, inform, 'nef', 0.2, test, train, nT,sample_size,a,b,m0)
    l0_nef20_results[id, nef_name] = nef20
    nT = as.integer(ncol(train) * 0.3)
    nef30 = Evaluate(cl_sample, inform, 'nef', 0.3, test, train, nT,sample_size,a,b,m0)
    l0_nef30_results[id, nef_name] = nef30
  }
}
save(list = c('l0_nef20_results','l0_nef30_results','l1_nef20_results',
              'l1_nef30_results','l2_nef20_results','l2_nef30_results'), 
     file = './results/nef_2030_l012.RData')

### informers overlapping
nA = 20
chem_informs = read.table(paste('~/CHTC_Downloads/Informer_', as.character(nA), '.txt', sep=''), header = F)
orig_informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_informer_', as.character(nA), '.txt', sep=''), header = F)
load('fda_data.RData')
counts = data.frame(aid = 1:60,count = rep(0, 60))
for (id in 1:60) {
  chem_inform = as.numeric(unlist(strsplit(as.character(chem_informs$V2[id]), split = ' ')))
  orig_inform = as.numeric(unlist(strsplit(as.character(orig_informs$V2[id]), split = ' ')))
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  chem_inform_cid = colnames(train)[chem_inform]
  orig_inform_cid = colnames(train)[orig_inform]
  counts[id, 'count'] = length(intersect(chem_inform_cid, orig_inform_cid))
}
summary(counts$count) # mean 6.3, median 7
overlap_mat = matrix(0, 60, 60)
informs = chem_informs
#informs = orig_informs
for (id1 in 1:60) {
  inform = as.numeric(unlist(strsplit(as.character(informs$V2[id1]), split = ' ')))
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id1), '_block.RData',sep=''))
  inform_cid1 = colnames(train)[inform]
  for (id2 in 1:30) {
    inform = as.numeric(unlist(strsplit(as.character(informs$V2[id2]), split = ' ')))
    load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id2), '_block.RData',sep=''))
    inform_cid2 = colnames(train)[inform]
    overlap_mat[id1, id2] = length(intersect(inform_cid1, inform_cid2))
  }
}
library(plot.matrix)
plot(overlap_mat)

baseline_roc = rep(NA,60)
baseline_nef = rep(NA,60)
# for naive baseline
for (id in 1:60) {
  load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
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

# for BF_w
baseline_roc = roc_results
baseline_nef = nef_results
sim = as.data.frame(data.table::fread('Similarity_mat_FDA.csv.gz'))
rownames(sim) = sim$V1
sim$V1 = NULL
for (id in 1:30) {
  load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
  complete_cpd = colnames(train)
  tmp_sim = sim[complete_cpd,complete_cpd]
  for (nA in 1:20) {
    roc_name = paste('roc_', as.character(nA), sep = '')
    nef_name = paste('nef_', as.character(nA), sep = '')
    freq = apply(train, 2, function(x){return(mean(x, na.rm = T))})
    inform = order(freq, decreasing = T)[1:nA]
    inform_cid = colnames(train)[inform]
    xA = as.matrix(test[inform])
    w = as.matrix(tmp_sim[,inform_cid], nrow = nrow(tmp_sim), ncol = nA)
    Scores = w %*% xA
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
write.table(baseline_roc, file = 'bfw_roc_results.txt',row.names = F)
write.table(baseline_nef, file = 'bfw_nef_results.txt',row.names = F)

# for random informer set
baseline_roc = read.table('./results/rdinfo_roc_results.txt', header = T)
baseline_nef = read.table('./results/rdinfo_nef_results.txt', header = T)
iter = 25
sample_size = 100
for (id in 1:60) {
  print(id)
  load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
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
      roc = roc + Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
      nef = nef + Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
    }
    roc = roc / iter
    nef = nef / iter
    baseline_roc[id, roc_name] = roc
    baseline_nef[id, nef_name] = nef
  }
}
write.table(baseline_roc, file = './results/rdinfo_roc_results.txt',row.names = F)
write.table(baseline_nef, file = './results/rdinfo_nef_results.txt',row.names = F)


### Plots
library(ggplot2)
results = data.frame('Informer_size' = (1:20))
cols = c(11, 21:39)
roc_results = read.table('./results/orig_L1_roc_results.txt', sep = ' ', header=T)
results[1:20,'orig_Boise'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
roc_results = read.table('fast_roc_results.txt', sep = ' ', header=T)
results[,'Boise_fast'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
# roc_results = read.table('fast_info_1_roc_results.txt', sep = ' ', header=T)[-58,]
# results[,'fast_pel1'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[cols]
roc_results = read.table('./results/rdinfo_roc_results.txt', sep = ' ', header=T)
results[,'Boise_rand'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
roc_results = read.table('./results/orig_L2_roc_results.txt', sep = ' ', header = T)
results[1:15, 'Boise_L2'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
roc_results = read.table('./results/orig_L0_roc_results.txt', sep = ' ', header = T)
results[1:15, 'Boise_L0'] = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[, 'Baseline'] = rep(mean(baseline_roc, na.rm=T), 20)

nef_results = read.table('./results/orig_L1_nef_results.txt', sep = ' ', header=T)
results[1:20,'orig_Boise'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
nef_results = read.table('fast_nef_results.txt', sep = ' ', header=T)
results[,'Boise_fast'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
# nef_results = read.table('fast_info_1_nef_results.txt', sep = ' ', header=T)[-58,]
# results[,'fast_pel1'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[cols]
nef_results = read.table('./results/rdinfo_nef_results.txt', sep = ' ', header=T)
results[,'Boise_rand'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:21]
nef_results = read.table('./results/orig_L2_nef_results.txt', sep = ' ', header=T)
results[1:15,'Boise_L2'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
nef_results = read.table('./results/orig_L0_nef_results.txt', sep = ' ', header=T)
results[1:15,'Boise_L0'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[, 'Baseline'] = rep(mean(baseline_nef, na.rm=T), 20)

revnef_results = read.table('./Res/block_revnef_results.txt', sep = ' ', header=T)
results[,'Boise_block'] = apply(revnef_results, 2,function(x){return(mean(x, na.rm = T))})[2:14]
revnef_results = read.table('./Res/revnef_results.txt', sep = ' ', header=T)
results[,'Boise_original'] = apply(revnef_results, 2,function(x){return(mean(x, na.rm = T))})[4:16]

legend_title_size = 13
legend_text_size = 13
axis_title_size = 14
axis_text_size = 13
p1 <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = orig_Boise))+
  geom_line(mapping = aes(x = Informer_size, y = orig_Boise, color = 'orig_Boise'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_fast))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_fast, color = 'fast_Boise'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_rand))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_rand, color = 'rand_Boise'))+
  geom_line(mapping = aes(x = Informer_size, y = Baseline, color = 'Baseline'))+
  # geom_point(mapping = aes(x = Informer_size, y = Boise_L0))+
  # geom_line(mapping = aes(x = Informer_size, y = Boise_L0, color = 'Boise_L0'))+
  # geom_point(mapping = aes(x = Informer_size, y = Boise_L1))+
  # geom_line(mapping = aes(x = Informer_size, y = Boise_L1, color = 'Boise_L1'))+
  # geom_point(mapping = aes(x = Informer_size, y = Boise_L2))+
  # geom_line(mapping = aes(x = Informer_size, y = Boise_L2, color = 'Boise_L2'))+
  scale_color_manual(name = "Methods", values = c('fast_Boise' = 'purple',
                                                  'rand_Boise' = 'orange',
                                                  'Boise_L0' = 'darksalmon',
                                                  'orig_Boise' = 'green',
                                                  'Boise_L2' = 'blue',
                                                  "Baseline" = "red"))+
  scale_x_continuous('Informer sizes', breaks = seq(1,20, by=1))+
  scale_y_continuous('ROCAUC mean', limits = c(0.789, 0.849))+
  #scale_y_continuous('NEF10 mean', limits = c(0.673, 0.785))+
  #scale_y_continuous('revNEF mean', limits = c(0.66, 0.73))+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        #legend.position = c(0.85,0.25),
        legend.position = c(0.85,0.3))
myGrobs <- list(p1,p2)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 1)
myGrobs <- list(p1,p3,p2,p4)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 2)

p1 <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = fast_pel1))+
  geom_line(mapping = aes(x = Informer_size, y = fast_pel1, color = 'fast_pel1'))+
  geom_point(mapping = aes(x = Informer_size, y = fast_entropy))+
  geom_line(mapping = aes(x = Informer_size, y = fast_entropy, color = 'fast_entropy'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_rand))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_rand, color = 'rand_Boise'))+
  scale_color_manual(name = "Methods", values = c("rand_Boise" = "red",
                                                  'fast_entropy' = 'darkblue',
                                                  'fast_pel1' = 'purple'))+
  scale_x_continuous('Informer sizes', breaks = seq(10,200, by=10))+
  scale_y_continuous('ROCAUC mean', limits = c(0.806, 0.950))+
  #scale_y_continuous('NEF mean', limits = c(0.692, 0.918))+
  labs(title = "Clustering on whole matrix")+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        #legend.position = c(0.85,0.25),
        legend.position = c(0.85, 0.25))

## plot of different nef thresholds
results = data.frame('Informer_size' = (1:15))
nef_results = read.table('./results/orig_L0_nef_results.txt', sep = ' ', header=T)
results[1:15,'Boise_L0'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
nef_results = read.table('./results/orig_L1_nef_results.txt', sep = ' ', header=T)
results[1:15,'Boise_L1'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
nef_results = read.table('./results/orig_L2_nef_results.txt', sep = ' ', header=T)
results[1:15,'Boise_L2'] = apply(nef_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]

results[1:15,'Boise_L0'] = apply(l0_nef20_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[1:15,'Boise_L1'] = apply(l1_nef20_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[1:15,'Boise_L2'] = apply(l2_nef20_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]

results[1:15,'Boise_L0'] = apply(l0_nef30_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[1:15,'Boise_L1'] = apply(l1_nef30_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]
results[1:15,'Boise_L2'] = apply(l2_nef30_results, 2,function(x){return(mean(x, na.rm = T))})[2:16]

p1 <- ggplot(results)+
  geom_point(mapping = aes(x = Informer_size, y = Boise_L0))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_L0, color = 'Boise_L0'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_L1))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_L1, color = 'Boise_L1'))+
  geom_point(mapping = aes(x = Informer_size, y = Boise_L2))+
  geom_line(mapping = aes(x = Informer_size, y = Boise_L2, color = 'Boise_L2'))+
  scale_color_manual(name = "Methods", values = c('Boise_L0' = 'darksalmon',
                                                  'Boise_L1' = 'green',
                                                  'Boise_L2' = 'blue',
                                                  "Baseline" = "red"))+
  scale_x_continuous('Informer sizes', breaks = seq(1,20, by=1))+
  ylab('NEF10 mean')+
  #scale_y_continuous('NEF mean', limits = c(0.673, 0.789))+
  theme(axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.position = c(0.8,0.3))
myGrobs <- list(p1,p2,p3)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 3,ncol = 1)


### informer set analysis
id=1
nA = 20
informs = read.table(paste('~/CHTC_Downloads/FDA_cv/orig_informer_', as.character(nA), '.txt', sep=''), header = F)
inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
# cl_sample summary
sample_size = 100
adj_mat = matrix(0, nrow(train), nrow(train))
for (i in 1:nrow(train)) {
  for (j in 1:nrow(train)) {
    adj_mat[i,j] = sum(cl_sample$CC[,i]==cl_sample$CC[,j])
  }
}
adj_mat = adj_mat / sample_size

candidates = c()
for (i in 1:nrow(adj_mat)) {
  tar_grp = which(adj_mat[i,] == 1)
  if(length(tar_grp) > 1)
    candidates = c(candidates, i)
}
tar_grps = list()
k=1
while(length(candidates)>0){
  c = candidates[1]
  tar_grps[[k]] = which(adj_mat[c,]==1)
  candidates = setdiff(candidates, tar_grps[[k]])
  k = k+1
}
sub_mats = lapply(tar_grps, function(x){return(train[x,])})
active_rates = lapply(sub_mats, function(X){return(apply(X, 2, function(y){return(mean(y,na.rm = T))}))})
inform_mats = matrix(0, nrow=26, ncol = ncol(train))
for (i in 1:26) {
  inform_mats[i, ] = active_rates[[i]]
}
inform_mats = inform_mats[c(1,4,6,8,9),]
inform_mats[,inform]

image(adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 600, 50), side=2, line=0.5, at=seq(0, 600, 50) / 687, las=1, cex=1)
mtext(text=seq(0, 600, 50), side=1, line=0.5, at=seq(0, 600, 50) / 687, las=1, cex=1)
