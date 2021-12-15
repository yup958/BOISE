### Write/revise job list for Condor submission
A = 1:933
pel1_res = read.csv('~/CHTC_Downloads/Informer_16.txt', header = F, sep = ' ')
pel1_res = pel1_res[order(pel1_res$V2),]
pre_inform = c(25, 195, 264, 285, 525, 189, 102, 259, 32, 205, 112, 71, 103, 95,138)
inform = pel1_res$V1[1] ## 25, 195, 264, 285, 525, 189,102, 259, 32, 205, 112, 71,103,95,138
inform = c(pre_inform, inform)
inform
A = A[-which(A %in% inform)]

write.table(A,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)

pel1_res = data.frame(nA = 3:14, pel = rep(NA, 12))
for (i in 3:14) {
  tmp_df = read.csv(paste('~/CHTC_Downloads/Informer_', as.character(i), '.txt', sep = ''), header = F,sep =' ')
  tmp_df = tmp_df[order(tmp_df$V2), ]
  pel1_res$pel[i-2] = tmp_df$V2[1]
}
