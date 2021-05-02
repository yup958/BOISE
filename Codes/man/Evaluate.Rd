\name{Evaluate}
\alias{Evaluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to evaluate selected informer set given the DPMM sampling.
}
\description{
Input: DPMM cluster sampling "cl_sample", informer set "inform", measurement criteria "measure", test set, train set, nA, nT, MCMC sample size iter, priors a and b, prior mass alpha
Output: Evaluated value under given criteria
}
\usage{
Evaluate(cl_sample, inform, measure, percent, test, train, nT, sample_size, alpha, beta, m0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
A list of clustering assignments. cl_sample = {KK, NN, CC}.
}
  \item{inform}{
Selected informer set
}
  \item{measure}{
Measurement methods. mesure = c("nef", "rocauc", "mat", "f")
}
  \item{percent}{
Percentage of Top set used in NEF evaluation. Between 0 and 1. Only used for "NEF".
}
  \item{test}{
True interaction on new target i*. Should be a vector/array
}
  \item{train}{
Equivalent to original interaction data x0.
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{sample_size}{
Sample size of Dirichlet process mixture model clustering
}
  \item{alpha}{
Prior hyperparameters alpha0
}
  \item{beta}{
Prior hyperparameters beta0
}
  \item{m0}{
Chinese restaurant prior mass m0
}
}
\details{
"nef": NEF10 metric; "rocauc": ROCAUC; "mat":MCC; "f": F1 score.
}
\value{
Evaluated ranking performance under given criteria for given informer set. The ranking is based on posterior expectation, as described in BOISE paper.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
Currently nef is only for NEF10 metric.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006813}
}
\examples{
#suppose you have the dat and cl_sample as in dpmm_beta
data(pkis1)
x0 <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
alpha = rep(mean(x0,na.rm = T),ncol(x0))
beta = 1-alpha
test = rbinom(ncol(x0), 1, 0.5)
nef = Evaluate(cl_sample, inform = c(1,2), measure = "nef", percent = 0.1, test, train = x0, nT = 10,
sample_size = 10, alpha, beta, m0 = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
