load('/ua/pyu58/results/col_clust_res_5.RData')
source("/ua/pyu58/codes/dpmm_beta.R")
source("/ua/pyu58/codes/Initial_beta.R")
source("/ua/pyu58/codes/Update_beta.R")
set.seed(2)
args <- commandArgs()
k = as.numeric(args[3])
grp = as.numeric(args[4])
m0 = as.numeric(args[5])

## helper function to find m0
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

train = t(train)
sample_size = 100
col_grp = which(cl_sample$CC[k, ]==grp)
train = train[, col_grp]
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
n = nrow(train)
prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
cl_sample = dpmm_beta(train, a, b, m0= m0, burn_in = 1000, sample_size=sample_size, thinning = 15)
diff = abs(mean(cl_sample$KK) - prior_m0_expect)

result = data.frame( "k" = k, "grp"=grp, "m0" = m0, "Differ" = diff)

## Write results
#save(list = c('train', 'cl_sample'), file = paste('col_clust_res_', as.character(m0), '.RData', sep = ""))
if(m0 < 10){
  write.table(result,file = paste("Prior_mass_", as.character(k), "_", as.character(grp), "_0", as.character(m0), ".txt", sep = ""),
              col.names = F,row.names = F)
} else{
  write.table(result,file = paste("Prior_mass_", as.character(k), "_", as.character(grp), "_", as.character(m0), ".txt", sep = ""),
              col.names = F,row.names = F)
}

