\name{Boise_Aug}
\alias{Boise_Aug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to augment an existing informer set, usually used in parallel computing.
}
\description{
Prerequisite is to sample DPMM samples with dpmm_beta.R to get cl_sample object. This function will take the
input of an old informer set and other corresponding parameters to give a 
larger BOISE informer set with size augmented by "nAdd".
}
\usage{
Boise_Aug(cl_sample, sample_size, interm_size, nT, alpha, beta, x0, m0, inform, nAdd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
A list of {K, N, C} comes from dpmm_beta.R
}
  \item{sample_size}{
Sample size of clustering assignments from Dirichlet process mixture model.
}
  \item{interm_size}{
Sample size of possible intermediate outcomes xA
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
  \item{inform}{
An old informer set with smaller size.
}
  \item{nAdd}{
Number of new informers to add into original informer set.
}
}
\details{
Parameters should be the same as in \code{Boise}. nAdd is the number of new informers added in this process.
}
\value{
An augmented informer set with size increased by "nAdd"
}
\references{
Yu, P., Ericksen, S.S., Gitter, A. and Newton, M.A. (2020) \emph{Bayes Optimal Informer Sets for Drug Discovery}. University of Wisconsin, Technical Report, 2020.
}
\author{
Peng Yu, Michale A. Newton
}
\note{
When original informer set size is large, say >=10, we recommend to choose nAdd = 1 or 2 at each time, otherwise may take really long time. We highly recommeend use parallel computing to select large informer sets. 
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
cl_sample = dpmm_beta(x0, alpha, beta, m0=15, burn_in=500, sample_size=10, thinning=10)
interm_size = 1000
A = Boise(cl_sample, sample_size = 10, interm_size = 1000, nA = 2, nT = 10, alpha, beta, x0, m0 = 15)
A_new = Boise_Aug(cl_sample,sample_size = 10, interm_size = 2000,nT = 10,alpha, beta, x0, m0 = 15, 
                    A, nAdd = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
