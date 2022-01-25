### Write/revise job list for Condor submission
A = 1:933
pel1_res = read.csv('~/CHTC_Downloads/Block_Informer_15.txt', header = F, sep = ' ')
pel1_res = pel1_res[order(pel1_res$V2),]
pre_inform = c(91, 10, 747, 85, 137, 525, 62, 189, 297, 470, 638, 82, 46, 713)
inform = pel1_res$V1[1] ## 91, 10, 747, 85, 137, 525, 62, 189, 297, 470, 638, 82, 46, 713
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
