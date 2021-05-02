\name{Boise}
\alias{Boise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to identify a Bayes optimal informer set.
}
\description{
Input a binary matrix of bioactivity data and corresponding Dirichlet-process-mixture-model samples  of
clustering over targets (with dpmm_beta.R to get cl_sample object). This function 
computes and minimizes posterior expected loss and returns an informer set A of size nA compounds.
}
\usage{
Boise(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
Posterior sampling of clusterings: A list of clustering assignments. cl_sample = {KK, NN, CC}.
}
  \item{sample_size}{
Sample size of clustering assignments sampled from Dirichlet process mixture model clustering. 
}
  \item{interm_size}{
Sample size of possible intermediate outcomes xA on the new target.
}
  \item{nA}{
Length of informer set
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{alpha}{
Prior hyperparameters alpha0
}
  \item{beta}{
Prior hyperparameters beta0
}
  \item{x0}{
Initial bioactivity data.
}
  \item{m0}{
Chinese restaurant prior mass m0
}
}
\details{
x0 is the original data used throughout the process. alpha, beta, m0 and sample_size are the same parameters for \code{dpmm_beta}. interm_size is a new parameter for sample size of possible outcomes on the new target.  
}
\value{
Bayes optimal informer set A with size nA.
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michael A. Newton
}
\note{
  Please test the function for small nA. When nA >= 8 we recommend using Boise_Aug to augment informer set step by step, otherwise it may take days to complete the selection. The algorithm used in this function corresponds to Algorithm 2 and 4 in the reference paper.
}


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
cl_sample = dpmm_beta(x0, alpha, beta, m0=15, burn_in=500, sample_size=10, thinning=10)
interm_size = 1000
A = Boise(cl_sample, sample_size = 10, interm_size = 1000, nA = 2, nT = 10, alpha, beta, x0, m0 = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
