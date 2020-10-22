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
Boise_Aug(cl_sample, iter, size, nT, a, b, x0, alpha, inform, nAdd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
A list of {K, N, C} comes from dpmm_beta.R
}
  \item{iter}{
Sample size of DPMM clustering
}
  \item{size}{
Sample size of possible intermediate outcomes xA
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{a}{
Prior hyperparameters alpha0
}
  \item{b}{
Prior hyperparameters beta0
}
  \item{x0}{
Original interaction data.
}
  \item{alpha}{
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An augmented informer set with size increased by "nAdd"
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
When original informer set size is large, say >=10, we recommend to choose nAdd = 1 or 2 at each time, otherwise may take really long time. We highly recommeend use parallel computing to select large informer sets. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#suppose you have the dat and cl_sample as in dpmm_beta
data(pkis1)
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
a = rep(mean(dat,na.rm = T),ncol(dat))
b = 1-a
cl_sample = dpmm_beta(a,b,x0 = dat, warm=500, iter=10, step=10, alpha=15)
size = 1000
A = Boise(cl_sample, iter = 10, size = 1000, nA = 2, nT = 10, a, b, x0 = dat, alpha = 15)
A_new = Boise_Aug(cl_sample,iter = 10, size = 2000,nT = 10,a,b,x0 = dat,alpha = 15, 
                    A, nAdd = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
