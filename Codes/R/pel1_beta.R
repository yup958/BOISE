pel1_beta <-
function(cl_sample, P, sample_size, interm_size, A, nT = 10, alpha, beta, x0, m0 = 2){
  #source("npel2.R")
  ### Compute PEL1
  nA = length(A)
  XA = matrix(0, sample_size * interm_size, nA)
  for (i in 1:sample_size) {
    XA[(interm_size * i - interm_size + 1):(interm_size * i),] = cl_sample$XX[i, ,A]
  }
  XA = apply(XA,1,function(x){return(paste(x, collapse = ""))})
  tab = table(XA)
  YA = names(tab)
  l = length(tab)
  tmp_pel2 = sapply(1:l, function(x){
    xA = as.numeric(strsplit(YA[x], "")[[1]])
    post_probs = matrix(0, 1, sample_size)
    post_thetas = matrix(0, sample_size, ncol(x0))
    for (k in 1:sample_size){
      postls = pel2_beta(P[[k]], x0, xA, A, nT, alpha, beta, m0)
      post_probs[k] = postls$post_prob
      post_thetas[k, ] = postls$post_theta
    }
    post_probs = post_probs / (sum(post_probs))
    Score = post_probs %*% post_thetas
    return(sum(sort(1-Score)[1:nT]))
  })
  pel1 = (tmp_pel2 * tab) / sample_size
  return(sum(pel1)/interm_size)
}
