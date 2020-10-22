\name{Boise}
\alias{Boise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to BOISE informer set selection, corresponding to Algorithm 2 in the paper
}
\description{
Prerequisite is to sample DPMM samples with dpmm_beta.R to get cl_sample object. This function will take the
input of cl_sample and other corresponding parameters to give a 
Bayes optimal informer set A with size nA.
}
\usage{
Boise(cl_sample, iter, size, nA, nT, a, b, x0, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl_sample}{
A list of clustering assignments. cl_sample = {KK, NN, CC}.
}
  \item{iter}{
Sample size of DPMM clustering
}
  \item{size}{
Sample size of possible intermediate outcomes xA
}
  \item{nA}{
Length of informer set
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Bayes optimal informer set A with size nA.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
  Please test the function for small nA. When nA >= 8 we recommend using Boise_Aug to augment informer set step by step, otherwise it may take days to complete the selection. 
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
