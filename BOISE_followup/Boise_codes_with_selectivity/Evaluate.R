Evaluate <-
function(cl_sample, inform, measure, percent,
                     test, train, nT, sample_size, alpha, beta, m0){
  #source("npel2.R")
  #source("clust_sum.R")
  if (!require('pROC')) {
    install.packages("pROC")
    library(pROC)
  }
  if (!require('ROCR')) {
    install.packages("ROCR")
    library(ROCR)
  }
  P = clust_sum(cl_sample,train,sample_size, alpha, beta)
  Score = rep(0, ncol(train))
  xA = test[inform]
  nA=length(inform)
  m = ncol(train)
  post_probs = matrix(0, 1, sample_size)
  post_thetas = matrix(0, sample_size, m)
  post_theta_mean = matrix(0, sample_size, m)
  for (k in 1:sample_size){
    postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, alpha, beta, m0)
    post_probs[k] = postls$post_prob
    post_thetas[k, ] = postls$post_theta
    post_theta_mean[k, ] = post_theta_j(P[[k]], x0=train, xA, A=inform, alpha, beta, m0)
  }
  post_probs = post_probs / (sum(post_probs))
  Score = post_probs %*% (post_thetas - post_theta_mean)
  Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA))
  Score[inform[which(xA==0)]] = rep(min(Score) - 1, nA - sum(xA))
  test = as.vector(test)
  Score = as.vector(Score)
  if(measure == "nef"){
    nTop = round(ncol(train) * percent, 0)
    top = order(Score,decreasing = T)[1:nTop]
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = min(hit,nTop)
    result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  } else if(measure == 'selectivity_nef'){
    nTop = round(ncol(train) * percent, 0)
    top = order(Score,decreasing = T)[1:nTop]
    selectivity_score = colMeans(train, na.rm = T)
    test = test - selectivity_score
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = sum(sort(test, decreasing = T)[1:nTop])
    result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  } else if(measure == "rocauc"){
    rocobj = pROC::roc(test,Score)
    result = rocobj$auc
  } else if(measure %in% c("mat", "f")){
    pred.obj = ROCR::prediction(Score, test)
    perform.obj = ROCR::performance(pred.obj, measure)
    result = max(unlist(perform.obj@y.values),na.rm = T)
  } else{
    print("Criteria is not supported.")
    result = 0
  }
  return(result)
}
