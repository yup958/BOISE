### Indian buffet process prior
ibp_prior <- function(m0, N){
  ## m0: prior mass; N: number of targets / rows
  
  # preallocate Z with upper bound of N*m0 number of latent factors
  cols = N*m0
  Z = matrix(NA, nrow=N, ncol=cols) 
  
  # start with some dishes/assigments
  dishes = rpois(1, m0)      
  zeroes = cols - dishes   # fill in the rest of potential dishes
  Z[1,] = c(rep(1, dishes), rep(0, zeroes))
  
  for(i in 2:N){
    prev = i-1
    # esoteric line that gets the last dish sampled without a search for it
    last_previously_sampled_dish = sum(colSums(Z[1:prev, ,drop=F]) > 0)    
    
    # initialize 
    dishes_previously_sampled = matrix(0, nrow=1, ncol=last_previously_sampled_dish)
    
    # calculate probability of sampling from previous dishes
    dish_prob = colSums(Z[1:prev, 1:last_previously_sampled_dish, drop=F]) / i
    dishes_previously_sampled[1,] = rbinom(n=last_previously_sampled_dish, size=1, prob=dish_prob)
    
    # sample new dish and assign based on results
    new_dishes = rpois(1, m0/i)
    zeroes = ncol(Z) - (last_previously_sampled_dish + new_dishes)
    Z[i,] = c(dishes_previously_sampled, rep(1,new_dishes), rep(0, zeroes))
  }
  
  # return only the dimensions sampled
  last_sampled_dish = sum(colSums(Z[1:prev, ]) > 0) 
  return(Z[, 1:last_sampled_dish])
}

### test
# m0 = 1
# N = 100
# Z = ibp_prior(m0, N)
# while(sum(rowSums(Z) == 0) > 0){
#   Z = ibp_prior(m0, N)
# }
# rowSums(Z)
