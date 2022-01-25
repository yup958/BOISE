nlines = sum(cl_sample$KK) * 10
job_list = data.frame('k'=rep(NA, nlines), 'grp' = rep(NA, nlines), 'm0' = rep(NA, nlines))
idx = 1
thres = 15
for (k in 91:100) {
  for (m0 in 11:20) {
    job_list[idx, ] = c(k, 1, m0)
    idx = idx + 1
  }
  for (m0 in 10:19) {
    job_list[idx, ] = c(k, 2, m0)
    idx = idx + 1
  }
  for (grp in 3:cl_sample$KK[k]) {
    if(cl_sample$NN[k, grp] < thres){
      break
    } else{
      for (m0 in 1:10) {
        job_list[idx, ] = c(k, grp, m0)
        idx = idx + 1
      }
    }
  }
}
job_list = job_list[1:(idx-1),]
write.table(job_list,file = '~/Upload/Boise_followup/job_list.txt',col.names = F,row.names = F)
