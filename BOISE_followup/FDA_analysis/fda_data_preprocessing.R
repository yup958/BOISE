## 1st round: largest complete data matrix
setwd('~/RAwork/BOISE/BOISE_followup/FDA_analysis/')
set.seed(2)
load('fda_data_short.RData')
assay = unique(data_short$PUBCHEM_AID)
cpd = unique(data_short$PUBCHEM_CID)
dat = matrix(NA, nrow = length(assay), ncol = length(cpd))
rownames(dat) = assay
colnames(dat) = cpd

for (i in 1:nrow(data_short)) {
  if(data_short$PUBCHEM_OUTCOME[i] == 'Active'){
    dat[data_short$PUBCHEM_AID[i], data_short$PUBCHEM_CID[i]] = 1
  } else{
    dat[data_short$PUBCHEM_AID[i], data_short$PUBCHEM_CID[i]] = 0
  }
}
mean(dat, na.rm = T)

par(mfrow = c(1,2))
row_counts = apply(dat, 1, function(x){return(mean(is.na(x)))})
hist(row_counts, breaks = 50,xlab = 'Missing rate', ylab = 'Counts', main = 'Hist of NA Rate in Assays')
col_counts = apply(dat, 2, function(x){return(mean(is.na(x)))})
hist(col_counts, breaks = 50,xlab = 'Missing rate', ylab = 'Counts', main = 'Hist of NA Rate in Compounds')

### Not working: greedy is not good on this problem.
# counts = apply(dat, 2, function(x){return(sum(!is.na(x)))})
# counts = sort(counts, decreasing = T)
# cid = names(counts[1:500])
# A = dat[, cid]
# counts = apply(A, 1, function(x){return(sum(!is.na(x)))})
# sort(counts, decreasing = T)[1]

### Helper function for initialization: max square matrix
largestSquare <- function(M, k = 1) {
  nr <- nrow(M); nc <- ncol(M)
  S <- 0*M; S[1, ] <- M[1, ]; S[, 1] <- M[, 1]
  for(i in 2:nr) 
    for(j in 2:nc)
      if (M[i, j] == 1) S[i, j] = min(S[i, j-1], S[i-1, j], S[i-1, j-1]) + 1
  o <- head(order(-S), k)
  d <- data.frame(row = row(M)[o], col = col(M)[o], mx = S[o])
  apply(d, 1, function(x) { 
    dn <- dimnames(M[x[1] - 1:x[3] + 1, x[2] - 1:x[3] + 1])
    out <- c(rownames(M) %in% dn[[1]], colnames(M) %in% dn[[2]]) + 0
    setNames(out, unlist(dimnames(M)))
  })
}

### function to find maximum complete submatrix
find_max_complete_submat <- function(A, penalty = prod(dim(dat))){
  A.na <- is.na(A) + 0
  Ainf <- ifelse(A.na, -penalty, 1) # used by f. Use the same penalization -prod(dim(dat)).
  nr <- nrow(A) # used by f
  f <- function(x) - c(x[seq(nr)] %*% Ainf %*% x[-seq(nr)])
  # starting values
  
  # Input is a square matrix of zeros and ones.
  # Output is a matrix with k columns such that first column defines the
  # largest square submatrix of ones, second defines next largest and so on.
  # Based on algorithm given here:
  # http://www.geeksforgeeks.org/maximum-size-sub-matrix-with-all-1s-in-a-binary-matrix/
  s <- seriation::seriate(A.na)
  p <- seriation::permute(A.na, s)
  # calcualte largest square submatrix in p of zeros rearranging to be in A's  order
  st <- largestSquare(1-p)[unlist(dimnames(A)), 1]
  # optimize relaxed problem
  res <- optim(st, f, lower = 0*st, upper = st^0, method = "L-BFGS-B")
  # retrieve results
  active_AID = res$par[1:nrow(A)]
  active_CID = res$par[(nrow(A)+1):(nrow(A)+ncol(A))]
  active_AID = names(active_AID)[which(active_AID>0)]
  active_CID = names(active_CID)[which(active_CID>0)]
  return(list(AID = active_AID, CID = active_CID))
}

largest_block = find_max_complete_submat(A = dat, penalty = 2)
#1:617, top 100 =564;top50 = 573; top10 = 595

A = dat[,largest_block$CID]
counts = apply(A, 1, function(x){return(sum(!is.na(x)))})
sort(counts, decreasing = T)[50]

save(list = c('dat', 'largest_block'), file = 'max_complete_data.RData')

### function to rearrange matrix into blocks. Stopping criteria: less than 20 cpds selected or > 5% missing
rearrange_mat_block <- function(M, cpd_thres = 20, mis_thres = 0.05){
  max_blocks = 100
  block_res = list(row_block = matrix(0, nrow = max_blocks, ncol = nrow(M)),
                   col_block = matrix(0, nrow = max_blocks, ncol = ncol(M)))
  colnames(block_res$row_block) = rownames(M)
  colnames(block_res$col_block) = colnames(M)
  k = 1 ## number of blocks
  A_cur = M ## Initial matrix
  while (TRUE) {
    print(k)
    tmp_res = find_max_complete_submat(A = A_cur)
    if(length(tmp_res$CID) < cpd_thres | mean(is.na(M[tmp_res$AID, tmp_res$CID])) > mis_thres | k > 99){
      break
    }
    block_res$row_block[k, tmp_res$AID] = 1
    block_res$col_block[k, tmp_res$CID] = 1
    A_cur = A_cur[ , setdiff(colnames(A_cur), tmp_res$CID)]
    k = k + 1
  }
  return(block_res)
}

block_res = rearrange_mat_block(M = dat, cpd_thres = 20, mis_thres = 0.05) ## 11 blocks for cpd>20; miss>0.05
block_res$row_block = block_res$row_block[1:11,]
block_res$col_block = block_res$col_block[1:11,]

### save results
save(list = c('dat', 'block_res'), file = 'fda_data_rearranged.RData')
load('fda_data_rearranged.RData')

### visualization
# lof function for ranking of block matrices
lof <- function(Z){
  binary_str = apply(Z, 2, function(z){
    return(paste(z, collapse = ""))
  })
  col_order = order(binary_str, decreasing = T)
  return(Z[, col_order])
}
# Get rearranged columns / rows
block_res$col_block = lof(block_res$col_block)
block_res$row_block = lof(block_res$row_block)
active_ind = apply(block_res$col_block, 2, sum)
active_ind = which(active_ind > 0)
cpds = colnames(block_res$col_block[, active_ind])
## rearrange row indices to make better visualization
end = 1
k = 1
assays = colnames(block_res$row_block)
while (k <= nrow(block_res$row_block)) {
  begin = end
  while (block_res$row_block[k, end] == 1) {
    end = end + 1
  }
  assays[begin:(end-1)] = rev(assays[begin:(end-1)])
  k = k+1
}
block_res$row_block = block_res$row_block[, assays]
## visualize blocked data
rotate <- function(x) t(apply(x, 2, rev))
A = dat[assays, cpds]
A = (1 + A) / 2
A[is.na(A)] = 0
# image(rotate(dat))
image(rotate(A), axes = F, xlab = 'Compounds', ylab = 'Targets')
mtext(text=seq(0, 600, 50), side=2, line=0.5, at=seq(0, 600, 50) / 688, las=1, cex=1)
mtext(text=seq(0, 900, 50), side=1, line=0.5, at=seq(0, 900, 50) / 933, las=1, cex=1)
print(apply(block_res$col_block, 1, sum))
print(apply(block_res$row_block, 1, sum))

col_grp5 = which(block_res$col_block[5,]>0)
cpds = colnames(block_res$col_block[,col_grp5])
row_grp5 = which(block_res$row_block[5,]>0)
assays = colnames(block_res$row_block[,row_grp5])
A = dat[assays, cpds]
counts = apply(A,2,function(x){return(sum(x, na.rm = T))})
sort(counts, decreasing = T)
save(list = c('block_res', 'dat', 'A', 'rotate'), file = 'fda_data_rearrange_w_blocks.RData')
