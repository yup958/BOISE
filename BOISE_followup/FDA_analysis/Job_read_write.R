### Write/revise job list for Condor submission
A = 1:933
pel1_res = read.csv('~/CHTC_Downloads/Block_Informer_15.txt', header = F, sep = ' ')
pel1_res = pel1_res[order(pel1_res$V2),]
pre_inform = c(91, 10, 747, 85, 137, 525, 62, 189, 297, 470, 638, 82, 46, 713,59)
inform = pel1_res$V1[1] ## 91, 10, 747, 85, 137, 525, 62, 189, 297, 470, 638, 82, 46, 713,59
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

### for m0 find
A = data.frame(grp = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20)))
A$m0 = rep(1:20, 5)

A = data.frame(grp = 1:19, m0 = c(4, 4, 4, 5, 5, 2, 2, 1, 2, 3, 4, 3, 4, 1, 1, 2, 1, 1, 1))
write.table(A,file = '~/Upload/Boise_followup/prior_mass_FDA_chemical_sep.txt',col.names = F,row.names = F)

