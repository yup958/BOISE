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
mean(is.na(dat))


### Continuous relaxation; Optimization for largest complete matrix
A = dat
A.na <- is.na(A) + 0

Ainf <- ifelse(A.na, -prod(dim(A)), 1) # used by f
nr <- nrow(dat) # used by f
f <- function(x) - c(x[seq(nr)] %*% Ainf %*% x[-seq(nr)])

# starting values

# Input is a square matrix of zeros and ones.
# Output is a matrix with k columns such that first column defines the
# largest square submatrix of ones, second defines next largest and so on.
# Based on algorithm given here:
# http://www.geeksforgeeks.org/maximum-size-sub-matrix-with-all-1s-in-a-binary-matrix/
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
s <- seriation::seriate(A.na)
p <- seriation::permute(A.na, s)
# calcualte largest square submatrix in p of zeros rearranging to be in A's  order
st <- largestSquare(1-p)[unlist(dimnames(A)), 1]

res <- optim(st, f, lower = 0*st, upper = st^0, method = "L-BFGS-B")
active_AID = res$par[1:nrow(dat)]
active_CID = res$par[(nrow(dat)+1):(nrow(dat)+ncol(dat))]
active_AID = names(active_AID)[which(active_AID>0)]
active_CID = names(active_CID)[which(active_CID>0)]

### save results
dat = dat[, active_CID]
train = dat[active_AID,]
test = dat[setdiff(rownames(dat),active_AID),]
save(list = c('dat', 'train', 'test'), file = 'max_complete_data.RData')

