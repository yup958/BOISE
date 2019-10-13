## This is a test part
set.seed(2154)
#setwd("~/Drug_Discovery/BOISE")
library(rstan)
#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load("BOISE_test2.RData")
rm(foo)
rm(u)
#load("clustering.RData")
#u <- foo$scaled.x
#dat <- 1*(u > .5 )
#source("dpmm.R")
source("pel1.R")
value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value)
A = i+1
nA = 1
nT = 10
pel_value = pel1(A, nA, nT, cl, dat)
if(nchar(value) == 1){
  write.table(pel_value,file = paste("Test_result00", value, ".txt", sep = ""),col.names = F)
} else if(nchar(value) == 2){
  write.table(pel_value,file = paste("Test_result0", value, ".txt", sep = ""),col.names = F)
} else{
  write.table(pel_value,file = paste("Test_result", value, ".txt", sep = ""),col.names = F)
}

#high_hit = which(apply(dat,2,sum)>5)
#clust1 = sample(which(cl$C == 1), 40)
#clust2 = sample(which(cl$C == 2), 20)
#clust3 = sample(which(cl$C == 3), 20)
#subdat = dat[c(clust1,clust2,clust3), high_hit]
#cl = dpmm(subdat,a = 4, aux = 4, iter = 100)
#cl$K = 3
#cl$N = c(40,20,20)
#cl$C = c(rep(1,40),rep(2,20),rep(3,20))

# pel_value[1] = 0.731215
# load("pkis1.rda")
# col_name = colnames(pkis1)
# huikun_AS <- c("5482344", "6539081", "6539107","6539108","10173796", "10308522", "11163861",
#                "11626927", "23646938", "25218600", "25218601", "25218614", "44536036", "53239967",
#                "56604034", "57391096")
# huikun <- rep(0,16)
# for (i in 1:length(huikun_AS)) {
#   huikun[i] <- which(col_name == huikun_AS[i])
# }


# pel_value[huikun]
# test <- apply(dat,2,sum)
# pel_test <- pel_value[huikun]
# test[huikun[order(pel_test)]]
