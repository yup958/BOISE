#load data
source("/ua/pyu58/codes/dpmm_beta.R")
source("/ua/pyu58/codes/Initial_beta.R")
source("/ua/pyu58/codes/Update_beta.R")

set.seed(817)
args <- commandArgs()
id = as.numeric(args[3])
m0 = as.numeric(args[4])

test_ids = read.table('/z/Comp/boise/data/Test_IDS.txt')
test_id = test_ids[id,1]

## helper function to find m0
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

## Train-test split
load('/z/Comp/boise/data/fda_data.RData')
test = dat[rownames(dat)==test_id,]
complete_idx = which(!is.na(test))
train = dat[!row.names(dat)==test_id, complete_idx]
test = test[complete_idx]
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
n = nrow(train)

prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
sample_size = 100
cl_sample = dpmm_beta(train, a, b, m0= m0, burn_in = 1000, sample_size=sample_size, thinning = 15)
diff = abs(mean(cl_sample$KK) - prior_m0_expect)
result = data.frame('k' = id, "m0" = m0, "Differ" = diff)

## Write results
# save(list = 'cl_sample', file = paste('orig_clust_res_', as.character(id), '.RData', sep = ""))
if(m0 < 10){
 write.table(result,file = paste("testID_", as.character(id),"_mass_0", as.character(m0), ".txt", sep = ""),col.names = F,row.names = F)
} else{
 write.table(result,file = paste("testID_", as.character(id), "_mass_", as.character(m0), ".txt", sep = ""),col.names = F,row.names = F)
}
