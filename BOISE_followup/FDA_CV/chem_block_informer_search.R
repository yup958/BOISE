set.seed(817)
source('/ua/pyu58/codes/chem_score_cpd.R')
source('/ua/pyu58/codes/chem_score_inform.R')
source("/ua/pyu58/codes/clust_sum.R")
source("/ua/pyu58/codes/pel2_beta.R")


## candidate inform
args <- commandArgs()
id = as.numeric(args[3])
A = as.numeric(args[4])

# load block RData file
load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))

# load pre informers
informs = read.table('~/CHTC_Downloads/FDA_cv/Informer_29.txt')
pre_inform = as.numeric(unlist(strsplit(as.character(informs$V2[id]), split = ' ')))
#pre_inform = c()
inform = c(pre_inform, A)
nA = length(inform)

## PEL1 Score for all the column clustering assignments
pel1 = 0
nT = as.integer(0.1 * ncol(train))
interm_size = 1000
row_sample_size = 100

## evaluate for all column clustering assignments
## compute scores and probability for all possible intermediate xA's
final_Scores = score_inform(cl, inform, train, m0s, block,
                              row_sample_size, interm_size, nT, simplified=T, thres = 2^16)
wt = unlist(lapply(final_Scores, function(x){return(unname(x$lg_wt))}))
wt = exp(wt - max(wt))
wt = wt / sum(wt)
pel2 = unlist(lapply(final_Scores, function(x){return(sum(1-unname(x$sc)))}))
pel1 = pel1 + sum(wt * pel2)

## Write result
result = data.frame('id' = id, 'A' = A,  "PEL" = pel1) 
if(A < 10){
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA), "_00", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else if(A < 100){
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA),"_0", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
} else{
  write.table(result,file = paste("testid_", as.character(id), "_Informer", as.character(nA),"_", as.character(A), ".txt", sep = ""),
              col.names = F,row.names = F)
}





