## Intermediate function for multiple informer features
## Input xb, Nb, xA,A
## Output posterior theta

multi_post <- function(xb, Nb,xA,A){
  if(Nb == 1){          
    tmp <- list(K = length(xb), L = length(A), A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Vector.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- mean(la$alpha)
    btahat <- apply(la$bta, 2, mean)
    theta <- exp(alphahat + btahat)/(1 + exp(alphahat + btahat))
    return(theta)
  } else{
    tmp <- list(J = Nb, K = dim(xb)[2], L = length(A), A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Matrix.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- mean(la$alpha)
    btahat <- apply(la$bta, 2, mean)
    theta <- exp(alphahat + btahat)/(1 + exp(alphahat + btahat))
    return(theta)
  }
}

# ## Test
# xb <- dat[1,]
# Nb <- 1
# xA <- dat[100, c(10,100)]
# A <- c(10,100)
# theta <- multi_post(xb,Nb,xA,A)
