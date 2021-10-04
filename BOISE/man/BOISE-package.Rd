\name{BOISE-package}
\alias{BOISE-package}
\alias{BOISE}
\docType{package}
\title{
\packageTitle{BOISE}
}
\description{
\packageDescription{BOISE}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BOISE}
\packageIndices{BOISE}

}
\author{
\packageAuthor{BOISE}

Maintainer: \packageMaintainer{BOISE}
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}

\keyword{ package }
\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
library(BOISE)
#load data
data(pkis1)
## Use 2sd criteria to create binary matrix
x0 <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
#choose iterations, warm up step, step length, prior mass
burn_in = 500
sample_size = 10
thinning = 10
m0 = 15
# One experiment in leave-one-out cross validation for BOISE framework
i = 1
test = x0[i, ]
train = x0[-i, ]
alpha = rep(mean(train),ncol(train))
beta = 1 - alpha
cl_sample = dpmm_beta(train, alpha, beta, m0, burn_in, sample_size, thinning)
nA = 2
nT = 36
interm_size = 1000
inform = Boise(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0 = train, m0, mcParallel = F)
## Add 1 more informer
inform = Boise_Aug(cl_sample, sample_size, interm_size, nT, alpha, beta, train, m0,
                    inform, nAdd = 1, mcParallel = F)
## To maintain consistency, we would keep the initial clustering samples 
## and use the same assignments in all the computations.
nef.result = Evaluate(cl_sample,inform, measure = "nef", percent = 0.1, test,train,
                      nT,sample_size, alpha, beta, m0)
auc.result = Evaluate(cl_sample, inform, measure = "rocauc", percent = 0.1, test, train,
                      nT,sample_size, alpha, beta, m0)
}
