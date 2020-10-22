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
Evaluate(cl_sample, inform, measure, test, train, nT, iter, a, b, alpha)
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
  \item{test}{
True interaction on new target i*. Should be a vector/array
}
  \item{train}{
Equivalent to original interaction data x0.
}
  \item{nT}{
Length of top set, used to define loss function in BOISE.
}
  \item{iter}{
Sample size of DPMM clustering
}
  \item{a}{
Prior hyperparameters alpha0
}
  \item{b}{
Prior hyperparameters beta0
}
  \item{alpha}{
Chinese restaurant prior mass m0
}
}
\details{
"nef": NEF10 metric; "rocauc": ROCAUC; "mat":MCC; "f": F1 score.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
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
test = rbinom(ncol(dat), 1, 0.5)
nef = Evaluate(cl_sample, inform = c(1,2), measure = "nef", test, train = dat, nT = 10, iter = 10,
      a, b, alpha = 15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
