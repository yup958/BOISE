## Intermediate function for multiple informer features
## Input xb, Nb, xA,A
## Output posterior theta

multi_post <- function(xb, Nb,xA,A){
  if(Nb == 1){          
    tmp <- list(K = length(xb), L = length(A), A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Vector.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 8000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- la$alpha[ ,Nb + 1]
    btahat <- la$bta
    theta <- 0
    for (i in 1:length(alphahat)) {
      theta <- theta + exp(alphahat[i] + btahat[i, ])/(1 + exp(alphahat[i] + btahat[i, ]))
    }
    return(theta/length(alphahat))
  } else{
    tmp <- list(J = Nb, K = dim(xb)[2], L = length(A), A = A, x = xb, y = xA)
    fit <- stan(file = "Informer_Rasch_Matrix.stan", data = tmp, verbose = FALSE, control = list(max_treedepth = 15), iter = 8000)
    la <- extract(fit, permuted = TRUE) # return a list of arrays 
    alphahat <- la$alpha[ ,Nb + 1]
    btahat <- la$bta
    theta <- 0
    for (i in 1:length(alphahat)) {
      theta <- theta + exp(alphahat[i] + btahat[i, ])/(1 + exp(alphahat[i] + btahat[i, ]))
    }
    return(theta/length(alphahat))
  }
}

## Test
# xb <- subdat[which(cl$C == 1),]
# Nb <- 43
# xA <- subdat[1, c(10,50)]
# A <- c(10,50)
# theta <- multi_post(xb,Nb,xA,A)
