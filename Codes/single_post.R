## Intermediate function for single informer feature
## Input xb, Nb, xA,A
## Output posterior theta

single_post <- function(xb, Nb,xA,A){
  if(Nb == 1){          
    tmp <- list(K = length(xb), A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Vector_1.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- mean(la$alpha)
    btahat <- apply(la$bta, 2, mean)
    theta <- exp(alphahat + btahat)/(1 + exp(alphahat + btahat))
    return(theta)
  } else{
    tmp <- list(J = Nb, K = dim(xb)[2], A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Matrix_1.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 2000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- mean(la$alpha)
    btahat <- apply(la$bta, 2, mean)
    theta <- exp(alphahat + btahat)/(1 + exp(alphahat + btahat))
    return(theta)
  }
}

## Test
# xb <- dat[1:10,]
# Nb <- 10
# xA <- dat[100, 10]
# A <- 10
# theta <- single_post(xb,Nb,xA,A)
