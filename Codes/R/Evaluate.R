Evaluate <-
function(cl_sample, inform, measure,
                     test, train, nT, iter, a, b ,alpha){
  #source("npel2.R")
  #source("clust_sum.R")
  P = clust_sum(cl_sample,train,iter, a, b)
  Score = rep(0, ncol(train))
  xA = test[inform]
  nA=length(inform)
  m = ncol(train)
  post_probs = matrix(0, 1, iter)
  post_thetas = matrix(0, iter, m)
  for (k in 1:iter){
    postls = pel2_beta(P[[k]], x0=train, xA, nA, A=inform, nT=36, a, b, alpha=15)
    post_probs[k] = postls$post_prob
    post_thetas[k, ] = postls$post_theta
  }
  post_probs = post_probs / (sum(post_probs))
  Score = post_probs %*% post_thetas
  Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA))
  Score[inform[which(xA==0)]] = rep(min(Score) - 1, nA - sum(xA))
  test = as.vector(test)
  Score = as.vector(Score)
  if(measure == "nef"){
    top = order(Score,decreasing = T)[1:nT]
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = min(hit,nT)
    result = ((pred_hit/nT - hit/ncol(train)) / (maxhit/nT - hit/ncol(train)) + 1)/2
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
