### Analysis based on Beta model
source("pel1_beta.R")
source("dpmm_beta.R")
#load data
load("clustering.RData")
u <- foo$scaled.x
dat <- 1*(u > .5 )
rm(u)
rm(foo)

#choose prior
a = apply(dat,2,sum)/100
b = 2 - a

#choose iterations, warm up step, step length, divergence alpha
warm = 200
iter = 20
step = 5
alpha = 2
#sample clustering assignments
cl_sample = dpmm_beta(a,b,x0 = dat, warm, iter, step, alpha)

#compute pel1 for nA = 1
pel1 = rep(0,366)
for (i in 1:366) {
  A = i
  nA = 1
  nT = 5
  pel1[i] = pel1_beta(cl_sample, iter, A, nA, nT, a, b, x0 = dat, alpha)
}

# Plotting
library(ggplot2)
testresult <- data.frame(pel1_value = pel1, label = 1:366, hitrate = apply(dat,2,sum))
p <- ggplot(data = testresult, mapping = aes(x = hitrate, y = pel1_value))+
  geom_point(color = "blue")
print(p)
