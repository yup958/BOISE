\name{dpmm_beta}
\alias{dpmm_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to do MCMC sampling for conjugate prior (Beta-binomial) clustering assignments
}
\description{
Output: A list cl_sample, including clustering assignments of size "sample_size".
}
\usage{
dpmm_beta(x0, alpha, beta, m0, burn_in = 200, sample_size = 10, thinning = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
  Original interaction data x0
}
  \item{alpha}{
  Prior hyperparameters alpha0
}
  \item{beta}{
  Prior hypereparameters beta0
}
  \item{m0}{
  Chinese restaurant prior mass m0.
}
  \item{burn_in}{
  Burn-in parameter. How many samples are dropped in the beginning.
}
  \item{sample_size}{
  Sample size of DPMM
}
  \item{thinning}{
  Thinning parameters. Gap between consecutive samples.
}
}
\details{
Typically, we choose hyperparameters alpha and beta with the mean of actives/inactives in the data set. We choose m0 with several tests that aim to minimize the gap of cluster numbers between Chinese restaurant process prior and posterior. 
}
\value{
cl_sample, is a list of clustering assignment samples. cl_sample = {KK, NN, CC}, where KK is an array  of size "iter" with each element as the number of clusters; NN is a "iter" * "row(x0)" matrix with each row as the cluster assignment of the targets; CC is summary of every cluster in each assignment. 
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
It may take hours for sampling, please try small sample size and burn-in/thinning step first. Only works for binary data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
data(pkis1)
## Use 2sd criteria to create binary matrix
x0 <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
alpha = rep(mean(x0,na.rm = T),ncol(x0))
beta = 1-alpha
cl_sample = dpmm_beta(x0, alpha, beta, m0 = 15, burn_in = 500, sample_size=10, thinning=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
