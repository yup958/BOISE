pel1_beta <-
function(cl_sample, P, iter, size, A, nA, nT = 10, a, b, x0, alpha = 2){
  #source("npel2.R")
  ### Compute PEL1
  XA = matrix(0, iter*size, nA)
  for (i in 1:iter) {
    XA[(size*i-size+1):(size*i),] = cl_sample$XX[i, ,A]
  }
  XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
  tab = table(XA)
  YA = names(tab)
  l = length(tab)
  tmp_pel2 = sapply(1:l, function(x){
    xA = as.numeric(strsplit(YA[x], "")[[1]])
    post_probs = matrix(0, 1, iter)
    post_thetas = matrix(0, iter, ncol(x0))
    for (k in 1:iter){
      postls = pel2_beta(P[[k]], x0, xA, nA, A, nT, a, b, alpha)
      post_probs[k] = postls$post_prob
      post_thetas[k, ] = postls$post_theta
    }
    post_probs = post_probs / (sum(post_probs))
    Score = post_probs %*% post_thetas
    return(sum(sort(1-Score)[1:nT]))
  })
  pel1 = (tmp_pel2 * tab) / iter
  return(sum(pel1)/size)
}
